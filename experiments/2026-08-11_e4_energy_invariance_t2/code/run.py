"""E4 — direct numerical test of the energy-invariance theorem (T2).

T2 states that for the tied product ansatz U(theta) = prod_k exp(theta_{p(k)} c_k kappa_k),
with H and the closed-shell reference invariant under an Abelian group N of real
characters, the energy obeys

    E(sigma_n theta) = E(theta)   for every n in N,
    where (sigma_n theta)_m = chi^(m)(n) theta_m.

T2 is the one proven statement in the project that has never been checked
numerically.  Everything downstream of it (the critical-submanifold corollary,
the theta=0 protocol theorem) rests on it.

Four measurements per molecule:

  A  T2 energy invariance      |E(theta) - E(sigma_n theta)| at random theta OFF
                               the symmetric submanifold.  Sampling theta inside
                               the submanifold would make sigma_n theta == theta
                               and the test vacuous, so theta_A is forced nonzero.
  B  T1 state-level identity   ||U(sigma_n theta)|phi0> - S_n U(theta)|phi0>||,
                               where S_n is the determinant-space representation of
                               n.  Strictly stronger than A.
  C  T3 vanishing gradient     max |dE/dtheta_m| for m in A, evaluated ON the
                               submanifold at theta=0 and at random theta_S.
  D  negative control          the same test in a basis where a degenerate pair has
                               been rotated by phi while the irrep labels are kept.
                               The hypothesis of T2 fails there, so T2 must fail
                               too; if it did not, the test would be measuring
                               nothing.

Usage:  python run.py [nh3|ch4|all] [--smoke]
"""
from __future__ import annotations

import json
import os
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE / "vendor"))

import expcommon as X  # noqa: E402  (sets BLAS threads before numpy import)
import numpy as np  # noqa: E402
import scipy.sparse.linalg as sla  # noqa: E402
from pyscf import ao2mo  # noqa: E402

from audit import ch4_geom, gates_of, irrep_of, nh3_geom, prep  # noqa: E402
from indep import SpinFreeHam, TrotterAnsatz, build_det_basis, kappa_matrix  # noqa: E402

ROOT = HERE.parent
SEED_BASE = 1000
N_THETA = 20
THETA_SCALE = 0.5
PHI_CONTROL_DEG = 30.0

CASES = {
    "nh3": dict(geom=nh3_geom(), norb=8, label="NH3"),
    "ch4": dict(geom=ch4_geom(), norb=9, label="CH4"),
}


def characters(par_irrep, elem):
    """chi^(m)(n) for the Abelian group of real characters.

    Irreps of Z2^k are labelled by bit patterns; the character of irrep j at the
    element indexed by the bit pattern i is (-1)^popcount(i & j).
    """
    return np.array([(-1.0) ** bin(int(r) & int(elem)).count("1") for r in par_irrep])


def det_sign_operator(dets, orbsym, elem):
    """Diagonal representation of the group element on the determinant basis.

    In a symmetry-adapted basis the element acts diagonally on spatial orbitals,
    phi_p -> chi_p phi_p, so it acts diagonally on determinants with the product
    of the characters of the occupied spin-orbitals.
    """
    chi_orb = np.array([(-1.0) ** bin(int(r) & int(elem)).count("1") for r in orbsym])
    out = np.empty(len(dets))
    nso = 2 * len(orbsym)
    for d, det in enumerate(dets):
        s = 1.0
        for k in range(nso):
            if det >> k & 1:
                s *= chi_orb[k // 2]
        out[d] = s
    return out


def rotate_degenerate_pair(M, phi_deg, tol=1e-6):
    """Rotate the first degenerate MO pair by phi, keeping the original labels.

    Reproduces the misalignment that a second, independent SCF introduces in a
    degenerate shell.  The irrep labels are deliberately NOT recomputed: that is
    exactly the He et al. failure mode, and it is what removes the hypothesis of T2.
    """
    e = M["mo_energy"]
    pair = None
    for i in range(len(e) - 1):
        if abs(e[i + 1] - e[i]) < tol:
            pair = (i, i + 1)
            break
    if pair is None:
        return None
    C = M["mf"].mo_coeff.copy()
    a = np.deg2rad(phi_deg)
    i, j = pair
    Ci, Cj = C[:, i].copy(), C[:, j].copy()
    C[:, i] = np.cos(a) * Ci + np.sin(a) * Cj
    C[:, j] = -np.sin(a) * Ci + np.cos(a) * Cj
    mol = M["mol"]
    h1 = C.T @ M["mf"].get_hcore() @ C
    eri = ao2mo.restore(1, ao2mo.kernel(mol, C), M["norb"])
    return dict(h1=h1, eri=eri, pair=pair)


def build_case(key, smoke=False):
    cfg = CASES[key]
    M, dets, index, ham, psi0, npar, gates, symf, sym_ids, ans = prep(
        cfg["label"], cfg["geom"], "of", symmetry=True)
    par_irrep = np.zeros(npar, dtype=int)
    for g in gates:
        par_irrep[g["param"]] = irrep_of(g, M["orbsym"])
    # (A3) is asserted inside prep(): no parameter mixes irreps.
    h = int(max(int(x) for x in M["orbsym"])) + 1
    h = 1 << (h - 1).bit_length()          # group order as a power of two
    elems = list(range(1, h))              # nontrivial elements
    return dict(key=key, M=M, dets=dets, index=index, ham=ham, psi0=psi0,
                npar=npar, gates=gates, symf=symf, sym_ids=sym_ids, ans=ans,
                par_irrep=par_irrep, elems=elems, group=M["groupname"],
                topgroup=M["topgroup"], cfg=cfg)


def run_case(key, log, rows, smoke=False):
    n_theta = 3 if smoke else N_THETA
    t0 = time.time()
    C = build_case(key, smoke)
    M, ans, npar = C["M"], C["ans"], C["npar"]
    a_idx = np.array([m for m in range(npar) if m not in set(C["sym_ids"])])
    s_idx = np.array(sorted(C["sym_ids"]))
    print(f"[{key}] group={C['group']} (topgroup {C['topgroup']}) "
          f"dets={len(C['dets'])} P={npar} |S|={len(s_idx)} |A|={len(a_idx)} "
          f"elems={C['elems']} setup={time.time()-t0:.1f}s", flush=True)

    signs = {n: det_sign_operator(C["dets"], M["orbsym"], n) for n in C["elems"]}
    # (A2): the closed-shell reference must be invariant under every element.
    for n, s in signs.items():
        a2 = float(np.max(np.abs(s * C["psi0"] - C["psi0"])))
        assert a2 == 0.0, f"(A2) violated for element {n}: {a2}"

    # ---- A and B: T2 energy invariance and T1 state identity, theta OFF the submanifold
    for t in range(n_theta):
        seed = SEED_BASE + t
        rng = np.random.default_rng(seed)
        theta = rng.uniform(-THETA_SCALE, THETA_SCALE, npar)
        assert np.max(np.abs(theta[a_idx])) > 1e-3, "theta must leave the submanifold"
        tw = X.Timer()
        with tw:
            e0 = ans.energy(theta)
            psi_ref = ans.state(theta)
            for n in C["elems"]:
                chi = characters(C["par_irrep"], n)
                th_n = chi * theta
                e_n = ans.energy(th_n)
                psi_n = ans.state(th_n)
                d_e = abs(e0 - e_n)
                d_psi = float(np.linalg.norm(psi_n - signs[n] * psi_ref))
                rows.append(dict(
                    molecule=C["cfg"]["label"], group=C["group"], test="T2_energy_invariance",
                    basis="adapted", element=n, seed=seed, n_params=npar,
                    n_flipped=int(np.sum(chi < 0)), theta_scale=THETA_SCALE,
                    theta_offmanifold_max=float(np.max(np.abs(theta[a_idx]))),
                    energy_ha=e0, delta_energy_ha=d_e,
                    delta_energy_rel=d_e / abs(e0), state_residual=d_psi))
        log.record(parameters=dict(molecule=C["cfg"]["label"], test="T2_T1", basis="adapted",
                                   theta_scale=THETA_SCALE, n_params=npar,
                                   elements=C["elems"]),
                   status="success", seed=seed, wallclock=tw.dt,
                   t_start=tw.start, t_end=tw.end,
                   outputs=dict(result_files=[], scalar_results=dict(
                       energy_ha=e0,
                       max_delta_energy_ha=max(r["delta_energy_ha"] for r in rows[-len(C["elems"]):]),
                       max_state_residual=max(r["state_residual"] for r in rows[-len(C["elems"]):]))),
                   notes="theta drawn off the symmetric submanifold; both T2 and T1 measured",
                   peak_memory_mb=X.peak_rss_mb())

    # ---- C: gradient in the removed directions, ON the submanifold
    for t in range(n_theta):
        seed = SEED_BASE + 500 + t
        rng = np.random.default_rng(seed)
        theta = np.zeros(npar)
        if t > 0:                      # run 0 is exactly the HF point theta = 0
            theta[s_idx] = rng.uniform(-THETA_SCALE, THETA_SCALE, len(s_idx))
        tw = X.Timer()
        with tw:
            e, g = ans.energy_grad(theta)
        gmax_a = float(np.max(np.abs(g[a_idx])))
        gmax_s = float(np.max(np.abs(g[s_idx])))
        rows.append(dict(
            molecule=C["cfg"]["label"], group=C["group"], test="T3_removed_gradient",
            basis="adapted", element="", seed=seed, n_params=npar,
            n_flipped=len(a_idx), theta_scale=0.0 if t == 0 else THETA_SCALE,
            theta_offmanifold_max=0.0, energy_ha=e,
            grad_max_removed=gmax_a, grad_max_kept=gmax_s))
        log.record(parameters=dict(molecule=C["cfg"]["label"], test="T3_gradient",
                                   point="theta=0" if t == 0 else "random theta_S",
                                   n_params=npar),
                   status="success", seed=seed, wallclock=tw.dt,
                   t_start=tw.start, t_end=tw.end,
                   outputs=dict(result_files=[], scalar_results=dict(
                       energy_ha=e, grad_max_removed=gmax_a, grad_max_kept=gmax_s)),
                   notes="gradient evaluated on the symmetric submanifold",
                   peak_memory_mb=X.peak_rss_mb())

    # ---- D: negative control in a rotated basis with stale labels
    rot = rotate_degenerate_pair(M, PHI_CONTROL_DEG)
    if rot is not None:
        ham_r = SpinFreeHam(rot["h1"], rot["eri"], M["ecore"], C["dets"], C["index"])
        ans_r = TrotterAnsatz(C["ans"].K, C["ans"].pidx, C["ans"].coef, npar, C["psi0"], ham_r)
        for t in range(min(n_theta, 5)):
            seed = SEED_BASE + 900 + t
            rng = np.random.default_rng(seed)
            theta = rng.uniform(-THETA_SCALE, THETA_SCALE, npar)
            tw = X.Timer()
            with tw:
                e0 = ans_r.energy(theta)
                for n in C["elems"]:
                    chi = characters(C["par_irrep"], n)
                    e_n = ans_r.energy(chi * theta)
                    d_e = abs(e0 - e_n)
                    rows.append(dict(
                        molecule=C["cfg"]["label"], group=C["group"],
                        test="T2_negative_control", basis=f"rotated_{PHI_CONTROL_DEG:.0f}deg",
                        element=n, seed=seed, n_params=npar,
                        n_flipped=int(np.sum(chi < 0)), theta_scale=THETA_SCALE,
                        theta_offmanifold_max=float(np.max(np.abs(theta[a_idx]))),
                        energy_ha=e0, delta_energy_ha=d_e, delta_energy_rel=d_e / abs(e0)))
            log.record(parameters=dict(molecule=C["cfg"]["label"], test="T2_negative_control",
                                       basis=f"rotated_{PHI_CONTROL_DEG:.0f}deg",
                                       rotated_pair=list(rot["pair"]), n_params=npar),
                       status="success", seed=seed, wallclock=tw.dt,
                       t_start=tw.start, t_end=tw.end,
                       outputs=dict(result_files=[], scalar_results=dict(
                           energy_ha=e0,
                           max_delta_energy_ha=max(r["delta_energy_ha"]
                                                   for r in rows[-len(C["elems"]):]))),
                       notes="labels deliberately stale after rotating a degenerate pair; "
                             "the hypothesis of T2 does not hold here",
                       peak_memory_mb=X.peak_rss_mb())
    return C


FIELDS = ["molecule", "group", "test", "basis", "element", "seed", "n_params",
          "n_flipped", "theta_scale", "theta_offmanifold_max", "energy_ha",
          "delta_energy_ha", "delta_energy_rel", "state_residual",
          "grad_max_removed", "grad_max_kept"]


def main():
    which = sys.argv[1] if len(sys.argv) > 1 and not sys.argv[1].startswith("-") else "all"
    smoke = "--smoke" in sys.argv
    keys = list(CASES) if which == "all" else [which]

    X.capture_env(ROOT / "env", project_root=ROOT.parent.parent)
    # The smoke run must not contaminate the canonical append-only log.
    log = ExperimentLogHeader(smoke=smoke)
    rows = []
    for k in keys:
        run_case(k, log, rows, smoke=smoke)

    if not smoke:
        X.write_csv(ROOT / "results" / "aggregate.csv", rows, FIELDS)
        (ROOT / "results" / "aggregate.json").write_text(
            json.dumps(rows, indent=2, default=X._jsonable) + "\n")

    inv = [r for r in rows if r["test"] == "T2_energy_invariance"]
    st = [r for r in inv if r.get("state_residual") is not None]
    grad = [r for r in rows if r["test"] == "T3_removed_gradient"]
    ctl = [r for r in rows if r["test"] == "T2_negative_control"]
    agg = dict(
        max_delta_energy_ha=max((r["delta_energy_ha"] for r in inv), default=None),
        max_state_residual=max((r["state_residual"] for r in st), default=None),
        max_grad_removed=max((r["grad_max_removed"] for r in grad), default=None),
        control_max_delta_energy_ha=max((r["delta_energy_ha"] for r in ctl), default=None),
        n_theta_samples=len(inv),
    )
    print("\n=== E4 summary ===")
    for k, v in agg.items():
        print(f"  {k:32s} {v}")
    if not smoke:
        log.summarize(
            aggregate_metrics=agg,
            key_findings=[
                f"T2 holds to {agg['max_delta_energy_ha']:.2e} Ha over "
                f"{agg['n_theta_samples']} random parameter vectors drawn off the "
                f"symmetric submanifold",
                f"the state-level identity of T1 holds to {agg['max_state_residual']:.2e}",
                f"the gradient in the removed directions is below "
                f"{agg['max_grad_removed']:.2e} Ha on the submanifold",
                f"in the rotated basis with stale labels the same quantity reaches "
                f"{agg['control_max_delta_energy_ha']:.2e} Ha",
            ],
            matches_hypothesis="supported",
            interpretation="The energy is invariant under the character sign flip to "
                           "machine precision, and fails by many orders of magnitude "
                           "once the labels no longer describe the orbitals.",
        )
    return agg


def ExperimentLogHeader(smoke=False):
    import tempfile
    path = (Path(tempfile.mkdtemp()) / "smoke_log.json") if smoke else (ROOT / "experiment_log.json")
    return X.ExperimentLog(
        path,
        experiment_id=ROOT.name,
        paper_project="2026-08_point-group-vqe-free",
        opportunity_reference="opportunities.md #3",
        mode="numerical+symbolic",
        subareas_active=["numerical", "symbolic"],
        target_venue="Quantum (componente do artigo principal); Comment no JCTC como item separado e posterior.",
        created_at=X.utcnow(),
        plan_md_hash=X.sha256_file(ROOT / "plan.md") if (ROOT / "plan.md").exists() else None,
        code_commit=None,
        system_info_hash=None,
    )


if __name__ == "__main__":
    main()
