"""Polikorik korelasyon matrisi (iki asamali ML; Olsson, 1979).
Bivariate normal CDF: Gauss-Legendre ile r uzerinde integral."""
import numpy as np
from scipy.stats import norm
from scipy.optimize import minimize_scalar

_GLX, _GLW = np.polynomial.legendre.leggauss(80)

def bvn_cdf(h, k, r):
    """P(Z1<=h, Z2<=k) for standard bivariate normal with corr r. h,k arrays."""
    h = np.asarray(h, float); k = np.asarray(k, float)
    out = norm.cdf(h) * norm.cdf(k)
    if abs(r) < 1e-12:
        return out
    t = 0.5 * r * (_GLX + 1.0)          # 0..r
    w = 0.5 * r * _GLW
    H = h[..., None]; K = k[..., None]
    one_t2 = 1.0 - t**2
    expo = -(H**2 - 2*t*H*K + K**2) / (2*one_t2)
    with np.errstate(over='ignore', invalid='ignore'):
        integ = np.exp(expo) / np.sqrt(one_t2)
    integ = np.nan_to_num(integ, nan=0.0, posinf=0.0)
    return out + (integ * w).sum(-1) / (2*np.pi)

def _thresholds(x, K):
    n = len(x)
    cum = np.cumsum(np.bincount(x, minlength=K)[:K]) / n
    tau = norm.ppf(np.clip(cum[:-1], 1e-6, 1-1e-6))
    return np.concatenate(([-np.inf], tau, [np.inf]))

def polychoric_pair(x, y, ta, tb):
    """x,y: 0-indexed integer codes (pairwise complete)."""
    Ka, Kb = len(ta)-1, len(tb)-1
    obs = np.zeros((Ka, Kb))
    np.add.at(obs, (x, y), 1)
    A = np.repeat(ta[:, None], Kb+1, 1)   # (Ka+1, Kb+1)
    B = np.repeat(tb[None, :], Ka+1, 0)
    def negll(r):
        C = bvn_cdf(np.clip(A, -8, 8), np.clip(B, -8, 8), r)
        C[np.isinf(A) & (A < 0)] = 0.0; C[np.isinf(B) & (B < 0)] = 0.0
        C[np.isinf(A) & (A > 0)] = norm.cdf(np.clip(B[0], -8, 8))[None, :][0] if False else C[np.isinf(A)&(A>0)]
        # cell probs by inclusion-exclusion
        P = C[1:, 1:] - C[:-1, 1:] - C[1:, :-1] + C[:-1, :-1]
        P = np.clip(P, 1e-12, 1.0)
        return -(obs * np.log(P)).sum()
    res = minimize_scalar(negll, bounds=(-0.999, 0.999), method='bounded',
                          options={'xatol': 1e-5})
    return res.x

def polychoric_matrix(df, verbose=False):
    cols = list(df.columns)
    p = len(cols)
    codes = {}
    ths = {}
    for c in cols:
        v = df[c].values
        u = np.array(sorted(np.unique(v[~np.isnan(v)])))
        codes[c] = (v, u)
        ok = ~np.isnan(v)
        idx = np.searchsorted(u, v[ok]).astype(int)
        ths[c] = _thresholds(idx, len(u))
    R = np.eye(p)
    for i in range(p):
        for j in range(i+1, p):
            ci, cj = cols[i], cols[j]
            vi, ui = codes[ci]; vj, uj = codes[cj]
            ok = (~np.isnan(vi)) & (~np.isnan(vj))
            xi = np.searchsorted(ui, vi[ok]).astype(int)
            xj = np.searchsorted(uj, vj[ok]).astype(int)
            R[i, j] = R[j, i] = polychoric_pair(xi, xj, ths[ci], ths[cj])
        if verbose: print('.', end='', flush=True)
    return R, cols

def nearest_pd(R):
    w, V = np.linalg.eigh(R)
    if w.min() > 1e-8: return R
    w = np.clip(w, 1e-6, None)
    A = V @ np.diag(w) @ V.T
    d = np.sqrt(np.diag(A))
    return A / np.outer(d, d)
