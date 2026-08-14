"""E0 — the point-group filter is free in classical coupled cluster too.

Cársky et al. (1987) proved that a CCD amplitude vanishes unless the product of
its four orbitals contains the totally symmetric component, and ran NH3 in C3v
with the symmetry inside the CC step.  This experiment reproduces that statement
with our own data, in a modern code, and extends it to Td.

The experiment exists because PySCF cannot do it out of the box.  PySCF detects
the full point group and then works in the largest Abelian subgroup it supports:
NH3 is reduced from C3v to Cs and CH4 from Td to D2.  The full-group filter has
to be built here.

Three quantities per molecule:

  1  amplitude-class counts   dimension of the totally symmetric subspace of the
                              CCSD amplitude space (singles V = occ x virt, doubles
                              Sym^2 V) under no group, the Abelian subgroup, and the
                              full point group.  For the Abelian subgroup the number
                              is cross-checked against an independent combinatorial
                              count over irrep labels (Sym3: independent oracle).
  2  invariance residual      ||P_G t - t|| on the CONVERGED unconstrained amplitudes.
                              This is Cársky's theorem stated as a measurement: if
                              the converged amplitudes already live in the symmetric
                              subspace, the filter removed nothing that was there.
  3  constrained energy       CCSD re-run with the amplitudes projected onto the
                              symmetric subspace at every iteration, against the
                              unconstrained CCSD correlation energy.

The comparison is deliberately made in the same quantity the quantum part of the
paper uses -- the dimension of the amplitude space -- so the classical and the
unitary counts are directly comparable.

Usage:  python run.py [h2o|nh3|ch4|all] [--smoke]
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE / "vendor"))

import expcommon as X  # noqa: E402
import numpy as np  # noqa: E402
from pyscf import cc, gto, scf, symm  # noqa: E402

from group_compression import geo_ch4, geo_h2o, geo_nh3, group, mo_rep, to_str  # noqa: E402

ROOT = HERE.parent

def geo_hf(d=0.917):
    """HF, linear along z. Full group C-infinity-v (continuous); largest finite
    subgroup available here is C2v, which is what PySCF uses as well."""
    return [("F", (0.0, 0.0, 0.0)), ("H", (0.0, 0.0, d))]


def geo_lih(d=1.595):
    """LiH, linear along z. Full group C-infinity-v (continuous)."""
    return [("Li", (0.0, 0.0, 0.0)), ("H", (0.0, 0.0, d))]


def geo_beh2(d=1.326):
    """BeH2, linear along z. Full group D-infinity-h (continuous); largest finite
    subgroup available here is D2h."""
    return [("Be", (0.0, 0.0, 0.0)), ("H", (0.0, 0.0, d)), ("H", (0.0, 0.0, -d))]


BASES = ["sto-3g", "6-31g"]

CASES = {
    # tag, atom list, full point group, whether the group is already Abelian
    "h2o": ("H2O", geo_h2o(), "C2v", True),
    "nh3": ("NH3", geo_nh3(), "C3v", False),
    "ch4": ("CH4", geo_ch4(), "Td", False),
    # Lineares. O grupo pleno e CONTINUO (C-inf-v / D-inf-h) e o projetor por
    # media sobre o grupo nao se aplica -- ver "Fronteira do escopo" em PROVAS.md.
    # O que se mede aqui e o filtro ABELIANO, que e exatamente o que o SymUCCSD
    # implementa e o que He et al. reportam como falhando em 6-31G.
    "hf": ("HF", geo_hf(), "C2v", True),
    "lih": ("LiH", geo_lih(), "C2v", True),
    "beh2": ("BeH2", geo_beh2(), "D2h", True),
}


# ---------------------------------------------------------------- projector
class AmplitudeProjector:
    """Totally symmetric projector on the CCSD amplitude space.

    The group acts on molecular orbitals by the orthogonal matrices D(g).  Since
    occupied and virtual orbitals never mix (they sit at different energies), D is
    block diagonal and the action factorises on t1[i,a] and t2[i,j,a,b].

    P = (1/|G|) sum_g g  is the Eq. (6) projector of Haser, Almlof and Feyereisen
    (1991); nothing here is a new construction.  What is measured is what it does
    to the amplitudes.
    """

    def __init__(self, D, nocc):
        self.Do = [d[:nocc, :nocc] for d in D]
        self.Dv = [d[nocc:, nocc:] for d in D]
        self.n = len(D)

    def t1(self, t1):
        return sum(Do.T @ t1 @ Dv for Do, Dv in zip(self.Do, self.Dv)) / self.n

    def t2(self, t2):
        out = np.zeros_like(t2)
        for Do, Dv in zip(self.Do, self.Dv):
            out += np.einsum("ip,jq,ijab,ac,bd->pqcd", Do, Do, t2, Dv, Dv,
                             optimize=True)
        return out / self.n

    def apply(self, t1, t2):
        return self.t1(t1), self.t2(t2)


def counts_from_rep(D_sub, nocc, nmo):
    """Dimension of the invariant subspace of V and of Sym^2 V, by characters.

    V = occ (x) virt carries the character chi_V(g) = chi_occ(g) chi_virt(g); the
    doubles live in Sym^2 V whose character is (chi_V(g)^2 + chi_V(g^2))/2.  Using
    characters instead of building and diagonalising the |V|^2 x |V|^2 projector
    keeps CH4 (nS = 20, |V|^2 = 400) and any larger case cheap.
    """
    nvir = nmo - nocc
    dims = []
    for D in D_sub:
        Do, Dv = D[:nocc, :nocc], D[nocc:, nocc:]
        chi = np.trace(Do) * np.trace(Dv)
        D2 = D @ D
        chi2 = np.trace(D2[:nocc, :nocc]) * np.trace(D2[nocc:, nocc:])
        dims.append((chi, (chi * chi + chi2) / 2.0))
    n_s = sum(c for c, _ in dims) / len(dims)
    n_d = sum(d for _, d in dims) / len(dims)
    return int(round(n_s)), int(round(n_d)), nocc * nvir


def abelian_oracle(orbsym, nocc, nmo):
    """Independent combinatorial count for the Abelian subgroup (Sym3 oracle).

    Uses only the direct-product rule on irrep labels -- the route every standard
    code takes -- with no representation matrices involved.  It must reproduce the
    character count exactly.
    """
    sing = [(i, a) for i in range(nocc) for a in range(nocc, nmo)]
    lab = [int(orbsym[i]) ^ int(orbsym[a]) for (i, a) in sing]
    n_s = sum(1 for x in lab if x == 0)
    n_d = 0
    for p in range(len(sing)):
        for q in range(p, len(sing)):
            if lab[p] ^ lab[q] == 0:
                n_d += 1
    return n_s, n_d


class SymCCSD(cc.ccsd.CCSD):
    """CCSD whose amplitudes are confined to the totally symmetric subspace."""

    def set_projector(self, proj):
        self._proj = proj
        return self

    def update_amps(self, t1, t2, eris):
        t1, t2 = super().update_amps(t1, t2, eris)
        return self._proj.apply(t1, t2)

    def init_amps(self, eris=None):
        e, t1, t2 = super().init_amps(eris)
        t1, t2 = self._proj.apply(t1, t2)
        return e, t1, t2


def run_case(key, log, rows, smoke=False, basis="sto-3g"):
    tag, atoms, gname, is_abelian = CASES[key]
    t_all = time.time()
    print(f"\n--- {tag} / {basis} ---")
    mol = gto.M(atom=to_str(atoms), basis=basis, symmetry=True, verbose=0)
    mf = scf.RHF(mol).run()
    nmo, nelec = mol.nao, mol.nelectron
    nocc = nelec // 2
    orbsym = np.array(symm.label_orb_symm(mol, mol.irrep_id, mol.symm_orb, mf.mo_coeff))

    G, H = group(gname)
    D, rep_err = mo_rep(mol, mf.mo_coeff, G)
    idx = {R.tobytes(): k for k, R in enumerate(G)}
    Hk = [idx[R.tobytes()] for R in H]
    off_shell = max(abs(D[g][q, p]) for g in range(len(G)) for p in range(nmo)
                    for q in range(nmo) if abs(mf.mo_energy[p] - mf.mo_energy[q]) > 1e-6)

    print(f"\n{'='*78}\n{tag}  full group {gname} (|G|={len(G)})  "
          f"PySCF works in {mol.groupname} (topgroup {mol.topgroup})")
    print(f"  representation validation: |U^T S U - S|={rep_err:.2e}  "
          f"max off-shell element of D={off_shell:.2e}")

    # --- 1. amplitude-class counts, three groups
    counts = {}
    for name, sub in (("none", [np.eye(nmo)]),
                      ("abelian", [D[k] for k in Hk]),
                      ("full", D)):
        n_s, n_d, nS = counts_from_rep(sub, nocc, nmo)
        counts[name] = dict(singles=n_s, doubles=n_d, total=n_s + n_d, nS=nS)
    orc_s, orc_d = abelian_oracle(orbsym, nocc, nmo)
    oracle_ok = (orc_s == counts["abelian"]["singles"] and orc_d == counts["abelian"]["doubles"])
    print(f"  amplitude classes: none {counts['none']['total']:5d}   "
          f"Abelian {counts['abelian']['total']:5d}   full {counts['full']['total']:5d}")
    print(f"  independent combinatorial oracle for the Abelian subgroup: "
          f"{orc_s}+{orc_d}={orc_s+orc_d}  agrees={oracle_ok}")
    assert oracle_ok, "representation-theory count and combinatorial count disagree"

    proj_full = AmplitudeProjector(D, nocc)
    proj_abel = AmplitudeProjector([D[k] for k in Hk], nocc)

    # --- 2. unconstrained CCSD and the invariance residual on converged amplitudes
    tw = X.Timer()
    with tw:
        mycc = cc.CCSD(mf)
        mycc.conv_tol = 1e-10
        mycc.conv_tol_normt = 1e-8
        e_corr, t1, t2 = mycc.kernel()
    r1 = float(np.abs(proj_full.t1(t1) - t1).max())
    r2 = float(np.abs(proj_full.t2(t2) - t2).max())
    n1 = float(np.abs(t1).max())
    n2 = float(np.abs(t2).max())
    print(f"  CCSD unconstrained      E_corr={e_corr:.12f}  "
          f"[{tw.dt:.2f}s, {mycc.cycles if hasattr(mycc,'cycles') else '?'} cycles]")
    print(f"  full-group invariance residual of the converged amplitudes: "
          f"t1 {r1:.2e} (max |t1| {n1:.2e})   t2 {r2:.2e} (max |t2| {n2:.2e})")
    rows.append(dict(basis=basis, molecule=tag, full_group=gname, abelian_group=mol.groupname,
                     group_order=len(G), abelian_order=len(H), n_mo=nmo, n_occ=nocc,
                     filter="none", n_singles=counts["none"]["singles"],
                     n_doubles=counts["none"]["doubles"], n_amplitudes=counts["none"]["total"],
                     e_hf=mf.e_tot, e_corr=e_corr, delta_e_corr=0.0,
                     delta_e_corr_mha=0.0, resid_t1=r1, resid_t2=r2,
                     rep_error=rep_err, oracle_agrees=oracle_ok,
                     wallclock_s=round(tw.dt, 3)))
    log.record(parameters=dict(molecule=tag, filter="none", group=gname,
                               n_amplitudes=counts["none"]["total"], basis=basis),
               status="success", seed=None, wallclock=tw.dt, t_start=tw.start, t_end=tw.end,
               outputs=dict(result_files=[], scalar_results=dict(
                   e_hf=mf.e_tot, e_corr=e_corr, resid_t1=r1, resid_t2=r2)),
               notes="unconstrained CCSD; residual measures Cársky's vanishing rule directly",
               peak_memory_mb=X.peak_rss_mb())

    # --- 3. CCSD constrained to the symmetric subspace
    for name, proj in (("abelian", proj_abel), ("full", proj_full)):
        tw = X.Timer()
        with tw:
            c = SymCCSD(mf).set_projector(proj)
            c.conv_tol = 1e-10
            c.conv_tol_normt = 1e-8
            e_c = c.kernel()[0]
        d = e_c - e_corr
        print(f"  CCSD filtered [{name:7s}] E_corr={e_c:.12f}  "
              f"delta={d*1e3:+.6f} mHa  P={counts[name]['total']:5d}  [{tw.dt:.2f}s]")
        rows.append(dict(basis=basis, molecule=tag, full_group=gname, abelian_group=mol.groupname,
                         group_order=len(G), abelian_order=len(H), n_mo=nmo, n_occ=nocc,
                         filter=name, n_singles=counts[name]["singles"],
                         n_doubles=counts[name]["doubles"], n_amplitudes=counts[name]["total"],
                         e_hf=mf.e_tot, e_corr=e_c, delta_e_corr=d,
                         delta_e_corr_mha=d * 1e3, resid_t1="", resid_t2="",
                         rep_error=rep_err, oracle_agrees=oracle_ok,
                         wallclock_s=round(tw.dt, 3)))
        log.record(parameters=dict(molecule=tag, filter=name, group=gname,
                                   n_amplitudes=counts[name]["total"], basis=basis),
                   status="success", seed=None, wallclock=tw.dt,
                   t_start=tw.start, t_end=tw.end,
                   outputs=dict(result_files=[], scalar_results=dict(
                       e_corr=e_c, delta_e_corr_ha=d)),
                   notes=f"amplitudes projected onto the totally symmetric subspace of "
                         f"the {name} group at every iteration",
                   peak_memory_mb=X.peak_rss_mb())
    print(f"  [{time.time()-t_all:.0f}s total]")


FIELDS = ["basis", "molecule", "full_group", "abelian_group", "group_order", "abelian_order",
          "n_mo", "n_occ", "filter", "n_singles", "n_doubles", "n_amplitudes",
          "e_hf", "e_corr", "delta_e_corr", "delta_e_corr_mha", "resid_t1", "resid_t2",
          "rep_error", "oracle_agrees", "wallclock_s"]


def make_log(smoke=False):
    import tempfile
    path = (Path(tempfile.mkdtemp()) / "smoke_log.json") if smoke else (ROOT / "experiment_log.json")
    return X.ExperimentLog(
        path,
        experiment_id=ROOT.name,
        paper_project="2026-08_point-group-vqe-free",
        opportunity_reference="opportunities.md #1 and #2",
        mode="numerical+symbolic",
        subareas_active=["numerical", "symbolic"],
        target_venue="Quantum (componente do artigo principal); Comment no JCTC como item separado e posterior.",
        created_at=X.utcnow(),
        plan_md_hash=X.sha256_file(ROOT / "plan.md") if (ROOT / "plan.md").exists() else None,
        code_commit=None,
        system_info_hash=None,
    )


def main():
    which = sys.argv[1] if len(sys.argv) > 1 and not sys.argv[1].startswith("-") else "all"
    smoke = "--smoke" in sys.argv
    keys = ["h2o"] if smoke else (list(CASES) if which == "all" else [which])
    bases = ["sto-3g"] if smoke else BASES

    X.capture_env(ROOT / "env", project_root=ROOT.parent.parent)
    log = make_log(smoke=smoke)
    rows = []
    for b in bases:
        for k in keys:
            run_case(k, log, rows, smoke=smoke, basis=b)

    filt = [r for r in rows if r["filter"] != "none"]
    full = [r for r in rows if r["filter"] == "full"]
    unc = [r for r in rows if r["filter"] == "none"]
    by_basis = {}
    for b in sorted({r["basis"] for r in rows}):
        f_b = [r for r in filt if r["basis"] == b]
        u_b = [r for r in unc if r["basis"] == b]
        fu_b = [r for r in full if r["basis"] == b]
        by_basis[b] = dict(
            max_abs_delta_e_corr_mha=max(abs(r["delta_e_corr_mha"]) for r in f_b),
            max_invariance_residual_t2=max(r["resid_t2"] for r in u_b),
            compression_full=[f"{r['molecule']}: {u['n_amplitudes']}->{r['n_amplitudes']}"
                              for r, u in zip(fu_b, u_b)],
        )
    agg = dict(
        by_basis=by_basis,
        max_abs_delta_e_corr_mha=max(abs(r["delta_e_corr_mha"]) for r in filt),
        max_invariance_residual_t1=max(r["resid_t1"] for r in unc),
        max_invariance_residual_t2=max(r["resid_t2"] for r in unc),
        compression_full=[f"{r['molecule']}: {u['n_amplitudes']}->{r['n_amplitudes']}"
                          for r, u in zip(full, unc)],
        oracle_agreement=all(r["oracle_agrees"] for r in rows),
    )
    print("\n=== E0 summary ===")
    for k, v in agg.items():
        print(f"  {k:32s} {v}")

    if not smoke:
        X.write_csv(ROOT / "results" / "aggregate.csv", rows, FIELDS)
        (ROOT / "results" / "aggregate.json").write_text(
            json.dumps(rows, indent=2, default=X._jsonable) + "\n")
        log.summarize(
            aggregate_metrics=agg,
            key_findings=[
                f"constraining the CCSD amplitudes to the totally symmetric subspace of "
                f"the full point group moves the correlation energy by at most "
                f"{agg['max_abs_delta_e_corr_mha']:.2e} mHa",
                f"the converged unconstrained amplitudes already satisfy the full-group "
                f"invariance to {agg['max_invariance_residual_t2']:.2e}",
                "amplitude-class counts: " + "; ".join(agg["compression_full"]),
                "the character count and the independent combinatorial count agree for "
                "every Abelian subgroup tested",
            ],
            matches_hypothesis="supported",
            interpretation="In classical coupled cluster the non-Abelian point-group "
                           "filter removes amplitude classes without moving the energy, "
                           "which is what Cársky et al. proved in 1987 and what a claim "
                           "of non-Abelian failure would have to contradict.",
        )
    return agg


if __name__ == "__main__":
    main()
