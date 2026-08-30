"""Shared scaffolding for the canonical experiments of this project.

Provides: environment capture, deterministic BLAS threading, sha256 hashing,
an append-only experiment_log.json writer, and the canonical figure style.

Deliberately dependency-light: numpy + matplotlib only.  Copied verbatim into
every experiment folder so that each folder reproduces on its own.
"""
from __future__ import annotations

import hashlib
import json
import os
import platform
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

# --- Rule 11: BLAS thread count is controlled and recorded -------------------
# Must be set before numpy is imported by anything downstream.
BLAS_THREADS = os.environ.setdefault("OMP_NUM_THREADS", "1")
for _v in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, BLAS_THREADS)

import numpy as np  # noqa: E402


def utcnow() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def sha256_file(path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def sha256_array(arr) -> str:
    return hashlib.sha256(np.ascontiguousarray(arr).tobytes()).hexdigest()


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode()).hexdigest()


# --- Rule 8: system metadata -------------------------------------------------
def capture_env(env_dir, project_root=None):
    env_dir = Path(env_dir)
    env_dir.mkdir(parents=True, exist_ok=True)

    def _sh(cmd, cwd=None):
        try:
            return subprocess.run(cmd, shell=True, capture_output=True, text=True,
                                  timeout=30, cwd=cwd).stdout.strip()
        except Exception:
            return ""

    cpu_model = ""
    for line in _sh("lscpu").splitlines():
        if line.startswith("Model name:"):
            cpu_model = line.split(":", 1)[1].strip()
            break
    try:
        import psutil  # noqa
        ram = round(psutil.virtual_memory().total / 1e9, 2)
    except Exception:
        ram = round(os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES") / 1e9, 2)

    versions = {}
    for mod in ("numpy", "scipy", "pyscf", "openfermion", "openfermionpyscf", "matplotlib"):
        try:
            versions[mod] = __import__(mod).__version__
        except Exception:
            versions[mod] = None

    blas = ""
    try:
        cfg = np.show_config(mode="dicts")
        blas = str(cfg.get("Build Dependencies", {}).get("blas", {}).get("name", ""))
    except Exception:
        blas = "unknown"

    info = {
        "os": platform.system(),
        "kernel": platform.release(),
        "distribution": _sh("lsb_release -ds") or platform.platform(),
        "machine": platform.machine(),
        "cpu_model": cpu_model,
        "cpu_cores_physical": os.cpu_count() and int(_sh("lscpu -p=Core,Socket | grep -v '^#' | sort -u | wc -l") or 0),
        "cpu_cores_logical": os.cpu_count(),
        "ram_total_gb": ram,
        "gpu": None,
        "python_version": sys.version.split()[0],
        "python_executable": sys.executable,
        "blas_implementation": blas,
        "blas_threads": BLAS_THREADS,
        "package_versions": versions,
        "captured_at": utcnow(),
    }
    (env_dir / "system_info.json").write_text(json.dumps(info, indent=2) + "\n")

    # Rule 6: exact dependency versions
    freeze = _sh(f"{sys.executable} -m pip freeze")
    (env_dir / "requirements.txt").write_text(freeze + "\n")

    # Rule 7: git state of the project
    if project_root:
        (env_dir / "git_commit.txt").write_text(_sh("git rev-parse HEAD", cwd=project_root) + "\n")
        status = _sh("git status --porcelain", cwd=project_root)
        (env_dir / "git_status.txt").write_text(status + "\n")
    return info


# --- append-only experiment log ---------------------------------------------
class ExperimentLog:
    """Append-only writer for experiment_log.json.

    Runs are never overwritten; each call to `record` appends one entry and
    rewrites the file so an interrupted session still leaves a valid document.
    """

    def __init__(self, path, **header):
        self.path = Path(path)
        if self.path.exists():
            self.doc = json.loads(self.path.read_text())
            self.doc.update({k: v for k, v in header.items() if k not in ("runs", "summary")})
        else:
            self.doc = dict(header)
            self.doc.setdefault("runs", [])
            self.doc.setdefault("summary", {})
        self.doc.setdefault("runs", [])

    def next_id(self) -> int:
        return len(self.doc["runs"]) + 1

    def record(self, *, parameters, status, seed=None, exit_code=0,
               t_start=None, t_end=None, wallclock=None, inputs=None,
               outputs=None, notes="", peak_memory_mb=None):
        entry = {
            "run_id": self.next_id(),
            "timestamp_start": t_start or utcnow(),
            "timestamp_end": t_end or utcnow(),
            "wallclock_seconds": round(wallclock, 3) if wallclock is not None else None,
            "status": status,
            "seed": seed,
            "parameters": parameters,
            "exit_code": exit_code,
            "peak_memory_mb": peak_memory_mb,
            "inputs": inputs or {"data_files": [], "data_hashes": {}},
            "outputs": outputs or {"result_files": [], "scalar_results": {}},
            "stdout_path": None,
            "stderr_path": None,
            "notes": notes,
        }
        self.doc["runs"].append(entry)
        self.flush()
        return entry

    def summarize(self, *, aggregate_metrics, key_findings, matches_hypothesis, interpretation):
        runs = self.doc["runs"]
        self.doc["summary"] = {
            "n_runs_total": len(runs),
            "n_runs_success": sum(r["status"] == "success" for r in runs),
            "n_runs_failure": sum(r["status"] == "failure" for r in runs),
            "n_runs_aborted": sum(r["status"] == "aborted" for r in runs),
            "aggregate_metrics": aggregate_metrics,
            "key_findings": key_findings,
            "matches_hypothesis": matches_hypothesis,
            "interpretation": interpretation,
        }
        self.flush()

    def flush(self):
        self.path.write_text(json.dumps(self.doc, indent=2, default=_jsonable) + "\n")


def _jsonable(o):
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    return str(o)


def peak_rss_mb():
    import resource
    return round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0, 1)


# --- figures -----------------------------------------------------------------
# Okabe-Ito, canonical order.  Index 0 (black) is reserved for the reference
# curve and is never one of several competing methods.
OKABE_ITO = ["#000000", "#E69F00", "#56B4E9", "#009E73",
             "#F0E442", "#0072B2", "#D55E00", "#CC79A7"]
MARKERS = ["o", "s", "^", "D", "v", "P", "*", "X"]

# springer-1col family: the manuscript is built with svjour3 [smallextended],
# whose \textwidth is 338.0pt = 4.677in (measured, not assumed).  Figures are
# authored at the width they PRINT so LaTeX applies no scale factor: matplotlib
# point sizes are absolute, so any \includegraphics down-scale would shrink
# every font below the 7pt floor.  WIDE is \linewidth; SINGLE is 0.75\linewidth.
# Keep these in step with the \includegraphics widths in manuscript/main.tex.
FIGSIZE_SINGLE = (3.51, 2.63)
FIGSIZE_WIDE = (4.677, 3.505)


def use_style(figsize=FIGSIZE_SINGLE):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from cycler import cycler
    try:
        import scienceplots  # noqa: F401
        plt.style.use(["science", "ieee", "no-latex"])
    except Exception:
        pass
    matplotlib.rcParams.update({
        "figure.figsize": figsize,
        "savefig.dpi": 300,
        # NOT bbox_inches="tight": trimming changes the PDF page size, and the figure
        # validator checks the page against the venue column width.  constrained_layout
        # fits the content inside the declared figsize instead.
        "savefig.bbox": None,
        "savefig.pad_inches": 0.0,
        "figure.constrained_layout.use": True,
        "figure.constrained_layout.h_pad": 0.02,
        "figure.constrained_layout.w_pad": 0.02,
        "savefig.format": "pdf",
        "axes.prop_cycle": cycler(color=OKABE_ITO),
        "font.size": 8,
        "axes.labelsize": 8,
        "legend.fontsize": 7,
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
    })
    return plt


def write_csv(path, rows, fieldnames):
    import csv
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in fieldnames})
    return path


class Timer:
    def __enter__(self):
        self.t0 = time.time()
        self.start = utcnow()
        return self

    def __exit__(self, *a):
        self.dt = time.time() - self.t0
        self.end = utcnow()
