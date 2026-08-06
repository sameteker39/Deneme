"""Literature-based iterative EFA item purification."""
import numpy as np, pandas as pd, sys, warnings
warnings.filterwarnings('ignore')
sys.path.insert(0,'/home/user/Deneme/efa')
from efa_core import smooth, efa, fit_indices, communalities
from factor_analyzer.factor_analyzer import calculate_kmo

CRIT = dict(h2=.30, load=.40, cross_diff=.20, cross_sec=.30, msa=.60, redund=.85, minitems=3)

def diagnose(R, cols, k, N, rot='oblimin'):
    Rs = smooth(R, 0.02)
    L, Phi = efa(Rs, k, rot, 'minres')
    h2 = communalities(L, Phi)
    A = np.abs(L)
    order = np.argsort(-A, axis=1)
    prim = A[np.arange(len(cols)), order[:,0]]
    sec  = A[np.arange(len(cols)), order[:,1]] if k>1 else np.zeros(len(cols))
    fac  = order[:,0]+1
    try: msa = calculate_kmo(pd.DataFrame(np.linalg.inv(np.linalg.inv(Rs))))[0]
    except Exception: msa = np.full(len(cols), np.nan)
    d = pd.DataFrame({'Madde':cols,'Faktor':fac,'Yuk':prim.round(3),'2.Yuk':sec.round(3),
                      'Fark':(prim-sec).round(3),'h2':h2.round(3)})
    fi = fit_indices(Rs, L, k, N, Phi)
    return d, L, Phi, fi

def msa_vec(R):
    Ri = np.linalg.inv(smooth(R,0.02))
    D = np.sqrt(np.diag(Ri)); P = -Ri/np.outer(D,D); np.fill_diagonal(P,1)
    r2 = (R**2).sum(0)-1; p2 = (P**2).sum(0)-1
    return r2/(r2+p2)

def purify(R0, cols0, k, N, verbose=True, crit=CRIT, protect=()):
    R = R0.copy(); cols = list(cols0); log = []
    while True:
        d, L, Phi, fi = diagnose(R, cols, k, N)
        m = msa_vec(R); d['MSA'] = m.round(3)
        # ihlaller
        cand = []
        for i,row in d.iterrows():
            c = row['Madde']
            if c in protect: continue
            if row['h2'] < crit['h2']:          cand.append((c,'dusuk ortak varyans (h2=%.2f)'%row['h2'], 1, row['h2']))
            elif row['Yuk'] < crit['load']:     cand.append((c,'dusuk faktor yuku (%.2f)'%row['Yuk'], 2, row['Yuk']))
            elif row['2.Yuk']>=crit['cross_sec'] and row['Fark']<crit['cross_diff']:
                                                cand.append((c,'binisik yuk (%.2f/%.2f, fark=%.2f)'%(row['Yuk'],row['2.Yuk'],row['Fark']), 3, row['Fark']))
            elif row['MSA'] < crit['msa']:      cand.append((c,'dusuk MSA (%.2f)'%row['MSA'], 4, row['MSA']))
        if not cand:
            vc = d['Faktor'].value_counts()
            d.attrs['zayif_faktor'] = [int(f) for f in vc.index if vc[f] < crit['minitems']]
            return R, cols, d, L, Phi, fi, log
        cand.sort(key=lambda t:(t[2], t[3]))
        drop, reason = cand[0][0], cand[0][1]
        j = cols.index(drop)
        keep = [i for i in range(len(cols)) if i!=j]
        R = R[np.ix_(keep,keep)]; cols = [cols[i] for i in keep]
        log.append((drop, reason, len(cols)))
        if verbose: print("   - %-4s cikarildi: %-42s (kalan %d)"%(drop, reason, len(cols)))
