"""
Decisive test of the last open claim.

get_hiuccsd keeps an excitation operator iff its term appears in the (compressed)
fermionic Hamiltonian. Symmetry-forbidden integrals vanish *only if the orbitals are
symmetry-adapted*. So:

  basis B (symmetry=True , adapted)     -> HiUCCSD should collapse onto SymUCCSD
  basis A (symmetry=False, He et al.'s) -> forbidden integrals are NOT zero,
                                           so HiUCCSD keeps strictly more

He et al. report NH3: SymUCCSD 75 params / 163 ops ; HiUCCSD 93 params / 197 ops.
"""
import sys, math, numpy as np
from pyscf import gto, scf, ao2mo, symm
from openfermion.ops import InteractionOperator, FermionOperator
from openfermion.transforms import get_fermion_operator, normal_ordered
HERE = __import__('os').path.dirname(__import__('os').path.abspath(__file__))  # PATCHED: was a hard-coded absolute path
sys.path.insert(0,HERE)
from audit import gates_of, irrep_of
from he_pipeline import he_geom_nh3, he_geom_ch4, tostr

def ham_term_keys(mol, C, tol=1e-8):
    norb=mol.nao
    h1=C.T @ mol.intor('int1e_kin')@C + C.T @ mol.intor('int1e_nuc')@C
    eri=ao2mo.restore(1, ao2mo.kernel(mol,C), norb)
    n=2*norb
    one=np.zeros((n,n)); two=np.zeros((n,n,n,n))
    for p in range(norb):
        for q in range(norb):
            for s in (0,1): one[2*p+s,2*q+s]=h1[p,q]
    for p in range(norb):
        for q in range(norb):
            for r in range(norb):
                for s in range(norb):
                    v=0.5*eri[p,q,r,s]          # (pq|rs) -> p^ r^ s q
                    for sa in (0,1):
                        for sb in (0,1):
                            two[2*p+sa,2*r+sb,2*s+sb,2*q+sa]+=v
    H=get_fermion_operator(InteractionOperator(0.0,one,two))
    H=normal_ordered(H); H.compress(tol)
    keys=set()
    for term,c in H.terms.items():
        if abs(c)<tol: continue
        cre=tuple(sorted(i for i,d in term if d==1))
        ann=tuple(sorted(i for i,d in term if d==0))
        keys.add((cre,ann))
    return keys

def count(label, geom_fn, name):
    gs=tostr(geom_fn())
    molA=gto.M(atom=gs,basis='sto-3g',symmetry=False,verbose=0); mfA=scf.RHF(molA).run()
    molB=gto.M(atom=gs,basis='sto-3g',symmetry=True ,verbose=0); mfB=scf.RHF(molB).run()
    orbsym=np.array(symm.label_orb_symm(molB,molB.irrep_id,molB.symm_orb,mfB.mo_coeff))
    norb,nelec=molA.nao,molA.nelectron
    npar,gates=gates_of(norb,nelec)
    print(f"\n{'='*78}\n{label}  (topgroup {molB.topgroup} -> abelian {molB.groupname})")
    print(f"  UCCSD: {npar} params / {len({(g['occ'],g['virt']) for g in gates})} distinct operators")
    # SymUCCSD (irrep filter, adapted labels)
    symf=[irrep_of(g,orbsym)==0 for g in gates]
    sp={g['param'] for g,f in zip(gates,symf) if f}
    so={(g['occ'],g['virt']) for g,f in zip(gates,symf) if f}
    print(f"  SymUCCSD (irrep filter)                    : {len(sp):3d} params / {len(so):3d} operators")
    for tag,mol,C in (("B  symmetry-ADAPTED orbitals",molB,mfB.mo_coeff),
                      ("A  He et al. basis (symmetry=False)",molA,mfA.mo_coeff)):
        keys=ham_term_keys(mol,C)
        hp={g['param'] for g in gates if (g['virt'],g['occ']) in keys}
        ho={(g['occ'],g['virt']) for g in gates if (g['virt'],g['occ']) in keys}
        print(f"  HiUCCSD in basis {tag:36s}: {len(hp):3d} params / {len(ho):3d} operators")

# PATCHED: the three statements below ran at import time, so `from hi_count import
# ham_term_keys` executed the whole driver.  Guarded; the body is otherwise untouched.
if __name__ == '__main__':
    count("NH3 / STO-3G  (He et al. geometry)", he_geom_nh3, 'nh3')
    count("CH4 / STO-3G  (He et al. geometry)", he_geom_ch4, 'ch4')
    print("\nHe et al. Table II / Fig.: NH3 SymUCCSD 75/163, HiUCCSD 93/197 ;"
          "  CH4 SymUCCSD 65/146, HiUCCSD 188/410")
