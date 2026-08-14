"""Hipotese: a corrupcao exige DUAS camadas degeneradas.
Com uma so, trocar seus dois membros e' apenas renomear x<->y globalmente -- inocuo.
Com duas ou mais, uma troca independente em uma delas nao e' mais renomeacao global."""
import numpy as np, math
from pyscf import gto, scf
def nh3(d=1.09,a=107.3):
    t=math.radians(a); r=d*math.sin(t/2); z=-d*math.cos(t/2)
    return [['N',[0,0,0]]]+[['H',[r*math.cos(2*math.pi*k/3),r*math.sin(2*math.pi*k/3),z]] for k in range(3)]
def ch4(d=1.09):
    v=[[1,1,1],[1,-1,-1],[-1,1,-1],[-1,-1,1]]
    return [['C',[0,0,0]]]+[['H',list(d*np.array(x)/np.linalg.norm(x))] for x in v]
def h2o(d=1.02):
    a=np.radians(104.5); return [['O',[0,0,0]],['H',[d,0,0]],['H',[d*np.cos(a),d*np.sin(a),0]]]
def c2h4(d1=1.34,d2=1.08,a=120.0):
    r=a/180*np.pi; x,y=d2*np.cos(r),d2*np.sin(r); hc=d1/2
    return [['C',[-hc,0,0]],['C',[hc,0,0]],['H',[-hc-x,y,0]],['H',[-hc-x,-y,0]],['H',[hc+x,y,0]],['H',[hc+x,-y,0]]]
M=[('HF',[['H',[0,0,0]],['F',[0,0,1.0]]]),('LiH',[['H',[0,0,0]],['Li',[0,0,1.55]]]),
   ('H2O',h2o()),('BeH2',[['Be',[0,0,0]],['H',[0,0,1.32]],['H',[0,0,-1.32]]]),
   ('NH3',nh3()),('CH4',ch4()),('N2',[['N',[0,0,0]],['N',[0,0,1.2]]]),
   ('CO',[['C',[0,0,0]],['O',[0,0,1.2]]]),('NaH',[['H',[0,0,0]],['Na',[0,0,1.64]]]),('C2H4',c2h4())]
gs=lambda g:"; ".join(f"{a} {x[0]:.10f} {x[1]:.10f} {x[2]:.10f}" for a,x in g)
FAIL={'NH3','CH4','N2','CO','NaH'}
for basis in ('sto-3g','6-31g'):
    print(f"\n===== base {basis} =====")
    print(f"{'mol':6s} {'grupo':6s} {'camadas degeneradas':>20s} {'>=2?':>5s} {'falha (He)':>11s} {'bate':>5s}")
    print("-"*62)
    for name,geo in M:
        try:
            m=gto.M(atom=gs(geo),basis=basis,symmetry=True,verbose=0); mf=scf.RHF(m).run()
        except Exception as e:
            print(f"{name:6s} erro: {str(e)[:40]}"); continue
        e_=mf.mo_energy; sh=[]; k=0
        while k < len(e_):
            j=k
            while j+1 < len(e_) and abs(e_[j+1]-e_[k])<1e-6: j+=1
            if j>k: sh.append((k,j))
            k=j+1
        n=len(sh); pred=n>=2
        obs = name in FAIL if basis=='sto-3g' else None
        ok = ('YES' if pred==obs else 'no') if obs is not None else '-'
        print(f"{name:6s} {m.groupname:6s} {n:20d} {str(pred):>5s} {str(obs):>11s} {ok:>5s}")
