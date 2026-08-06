import numpy as np, pandas as pd
from factor_analyzer import FactorAnalyzer

def smooth(R, floor=0.01):
    w,V = np.linalg.eigh(R)
    if w.min() >= floor: return R
    w = np.clip(w, floor, None)
    A = V@np.diag(w)@V.T
    d = np.sqrt(np.diag(A))
    return A/np.outer(d,d)

def fit_indices(R, L, k, N, Phi=None):
    """psych::fa uyumlu uyum iyiligi indeksleri. Egik donduruide S = L Phi L' + Psi."""
    p = R.shape[0]
    if Phi is None: Phi = np.eye(L.shape[1])
    Cm = L@Phi@L.T
    psi = np.clip(1 - np.diag(Cm), 1e-4, None)
    S = Cm + np.diag(psi)
    Si = np.linalg.inv(S)
    F = np.log(np.linalg.det(S)) - np.log(np.linalg.det(R)) + np.trace(R@Si) - p
    F = max(F, 0.0)
    df = ((p-k)**2 - p - k)/2
    mult = (N - 1 - (2*p+5)/6 - 2*k/3)
    chi2 = mult*F
    Fn = -np.log(np.linalg.det(R)); dfn = p*(p-1)/2
    chi2n = (N - 1 - (2*p+5)/6)*Fn
    d_ = max(chi2-df, 0); dn_ = max(chi2n-dfn, 0)
    rmsea = np.sqrt(d_/(df*(N-1))) if df>0 else np.nan
    cfi = 1 - d_/max(dn_, d_, 1e-12)
    tli = ((chi2n/dfn)-(chi2/df))/((chi2n/dfn)-1) if df>0 else np.nan
    res = R - S
    iu = np.tril_indices(p, -1)
    srmr = np.sqrt((res[iu]**2).mean())
    from scipy.stats import chi2 as _c2
    pv = 1-_c2.cdf(chi2, df) if df>0 else np.nan
    # RMSEA %90 GA
    def ncp(target, lo, hi):
        from scipy.optimize import brentq
        f=lambda nc: _c2.cdf(chi2, df, nc)-target
        try:
            if f(0)<0: return 0.0
            return brentq(f,0,max(chi2*4,100))
        except Exception: return np.nan
    lo=ncp(.95,0,0); hi=ncp(.05,0,0)
    rl = np.sqrt(lo/(df*(N-1))) if df>0 and lo==lo else np.nan
    rh = np.sqrt(hi/(df*(N-1))) if df>0 and hi==hi else np.nan
    return dict(chi2=chi2, df=df, p=pv, chi2_df=chi2/df if df>0 else np.nan,
                RMSEA=rmsea, RMSEA_lo=rl, RMSEA_hi=rh, CFI=cfi, TLI=tli, SRMR=srmr,
                TLI_=tli)

def communalities(L, Phi=None):
    if Phi is None: Phi = np.eye(L.shape[1])
    return np.diag(L@Phi@L.T)

def efa(R, k, rot='oblimin', method='minres'):
    """MINRES/ML cikarim + GPA egik dondurme (Phi desenle tam tutarli)."""
    from rotate import gpa_oblique, gpa_orth
    fa = FactorAnalyzer(n_factors=k, rotation=None, method=method, is_corr_matrix=True)
    fa.fit(R)
    A = fa.loadings_
    if rot in ('varimax',):
        L, Phi = gpa_orth(A, 'varimax')
    elif rot is None or rot == 'none':
        L, Phi = A, np.eye(k)
    else:
        L, Phi = gpa_oblique(A, rot)
    # isaret normalizasyonu: her faktorun agirlikli yuk toplami pozitif olsun
    for c in range(k):
        if L[:, c].sum() < 0:
            L[:, c] *= -1
            Phi[c, :] *= -1; Phi[:, c] *= -1
    np.fill_diagonal(Phi, 1.0)
    # varyansa gore siralama
    ss = (L**2).sum(0); o = np.argsort(-ss)
    L = L[:, o]; Phi = Phi[np.ix_(o, o)]
    return L, Phi

def map_test(R, kmax=12):
    """Velicer MAP (orijinal ve 4. kuvvet)."""
    p=R.shape[0]; w,V=np.linalg.eigh(R); o=np.argsort(-w); w=w[o]; V=V[:,o]
    out=[]
    d=np.diag(R).copy()
    f2=[np.sum(R[np.triu_indices(p,1)]**2)/(p*(p-1)/2)]
    f4=[np.sum(R[np.triu_indices(p,1)]**4)/(p*(p-1)/2)]
    for m in range(1,kmax+1):
        A=V[:,:m]*np.sqrt(w[:m])
        Rp=R-A@A.T
        dg=np.sqrt(np.abs(np.diag(Rp))); dg[dg<1e-8]=1e-8
        Rs=Rp/np.outer(dg,dg)
        iu=np.triu_indices(p,1)
        f2.append(np.sum(Rs[iu]**2)/(p*(p-1)/2))
        f4.append(np.sum(Rs[iu]**4)/(p*(p-1)/2))
    return np.array(f2), np.array(f4)
