"""Gradient Projection Algorithm ile egik/dik dondurme (Jennrich 2002;
Bernaards & Jennrich 2005). Phi desen matrisiyle tam tutarli uretilir."""
import numpy as np

def _quartimin(L, gam=0.0):
    k = L.shape[1]
    N = np.ones((k,k)) - np.eye(k)
    L2 = L**2
    X = L2 @ N
    if gam != 0:
        p = L.shape[0]
        C = np.ones((p,p))/p
        X = (np.eye(p) - gam*C) @ X
    return (L2*X).sum()/4.0, L*X

def _geomin(L, eps=0.01):
    p,k = L.shape
    L2 = L**2 + eps
    pro = np.exp(np.log(L2).sum(1)/k)
    f = pro.sum()
    Gq = (2.0/k)*(L/L2)*pro[:,None]
    return f, Gq

def _varimax_crit(L):
    p,k = L.shape
    L2 = L**2
    C = L2 - L2.mean(0)
    f = -(C**2).sum()/4.0
    return f, -L*C

CRITS = {'oblimin': lambda L: _quartimin(L,0.0), 'quartimin': lambda L: _quartimin(L,0.0),
         'geomin': _geomin, 'varimax': _varimax_crit}

def gpa_oblique(A, crit='oblimin', normalize=True, maxit=1000, eps=1e-6):
    A = np.asarray(A, float); p,k = A.shape
    if k == 1: return A.copy(), np.ones((1,1))
    hw = np.sqrt((A**2).sum(1)) if normalize else np.ones(p)
    hw[hw < 1e-8] = 1e-8
    An = A/hw[:,None] if normalize else A
    f_crit = CRITS[crit]
    T = np.eye(k); al = 1.0
    Ti = np.linalg.inv(T); L = An @ Ti.T
    f, Gq = f_crit(L)
    G = -(L.T @ Gq @ Ti).T
    for it in range(maxit):
        Gp = G - T @ np.diag((T*G).sum(0))
        s = np.sqrt((Gp**2).sum())
        if s < eps: break
        al *= 2
        for _ in range(20):
            X = T - al*Gp
            v = 1/np.sqrt((X**2).sum(0))
            Tt = X*v
            try: Ti = np.linalg.inv(Tt)
            except np.linalg.LinAlgError: al /= 2; continue
            L = An @ Ti.T
            ft, Gqt = f_crit(L)
            if ft < f - 0.5*s*s*al: break
            al /= 2
        T = Tt; f = ft
        G = -(L.T @ Gqt @ Ti).T
    Phi = T.T @ T
    if normalize: L = L*hw[:,None]
    return L, Phi

def gpa_orth(A, crit='varimax', normalize=True, maxit=1000, eps=1e-6):
    A = np.asarray(A,float); p,k = A.shape
    if k == 1: return A.copy(), np.eye(1)
    hw = np.sqrt((A**2).sum(1)) if normalize else np.ones(p)
    hw[hw<1e-8]=1e-8
    An = A/hw[:,None] if normalize else A
    f_crit = CRITS[crit]; T=np.eye(k); al=1.0
    L = An@T; f,Gq = f_crit(L); G = An.T@Gq
    for it in range(maxit):
        M = T.T@G; S=(M+M.T)/2; Gp = G - T@S
        s = np.sqrt((Gp**2).sum())
        if s<eps: break
        al*=2
        for _ in range(20):
            X = T-al*Gp
            U,_,Vt = np.linalg.svd(X, full_matrices=False)
            Tt = U@Vt
            L = An@Tt; ft,Gqt = f_crit(L)
            if ft < f-0.5*s*s*al: break
            al/=2
        T=Tt; f=ft; G = An.T@Gqt
    if normalize: L = L*hw[:,None]
    return L, np.eye(k)
