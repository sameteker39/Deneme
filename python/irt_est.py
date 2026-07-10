"""
IRTree vs EM-IRT vs Standart 2PL — Python MML-EM kestirimi (mirt olmadan).
Eksik veri (TIMSS kitapcik tasarimi) desteklenir.
"""
import numpy as np

def _p2pl(a, b, th):
    return 1.0 / (1.0 + np.exp(-a * (th - b)))

def fit_2pl_mml(X, n_grid=41, lo=-4.0, hi=4.0, max_iter=300, tol=1e-4, seed=0):
    """
    Unidim 2PL, marjinal ML (Bock-Aitkin EM). X: N x J {0,1,nan}.
    Doner: a (J,), b (J,), theta_EAP (N,), loglik.
    """
    X = np.asarray(X, float)
    N, J = X.shape
    obs = ~np.isnan(X)
    Xf = np.where(obs, X, 0.0)
    th = np.linspace(lo, hi, n_grid)
    W = np.exp(-0.5 * th**2); W /= W.sum()          # standart normal a priori
    a = np.ones(J); b = np.nanmean(X, 0); b = np.clip(-(b - 0.5) * 2, -3, 3)  # kaba baslangic
    ll_old = -np.inf
    for it in range(max_iter):
        # --- E-step ---
        P = _p2pl(a[None, :], b[None, :], th[:, None])      # (G,J)
        P = np.clip(P, 1e-9, 1 - 1e-9)
        logP, log1P = np.log(P), np.log(1 - P)
        # kisi x grid log-likelihood (yalniz gozlenen maddeler)
        LL = (Xf @ logP.T) + ((obs * (1 - Xf)) @ log1P.T)   # (N,G) ; gozlenmeyen -> katki 0
        # duzeltme: gozlenmeyen maddeler icin logP katkisi eklenmemeli -> Xf=0 & obs=0 halleder
        LL += np.log(W)[None, :]
        m = LL.max(1, keepdims=True)
        post = np.exp(LL - m); post /= post.sum(1, keepdims=True)   # (N,G)
        ll = (m.ravel() + np.log(np.exp(LL - m).sum(1))).sum()
        # --- M-step: her madde icin beklenen sayimlar ---
        # n_jq = sum_i post_iq*obs_ij ; r_jq = sum_i post_iq*x_ij*obs_ij
        nq = post.T @ obs.astype(float)      # (G,J)
        rq = post.T @ (Xf * obs)             # (G,J)
        for j in range(J):
            aj, bj = a[j], b[j]
            for _ in range(50):
                p = np.clip(_p2pl(aj, bj, th), 1e-9, 1 - 1e-9)
                w = nq[:, j] * p * (1 - p)
                res = rq[:, j] - nq[:, j] * p
                # gradyan (a,b): dP/da=(th-b)p(1-p), dP/db=-a p(1-p)
                g_a = np.sum((th - bj) * res)
                g_b = np.sum(-aj * res)
                H_aa = -np.sum(w * (th - bj)**2)
                H_bb = -np.sum(w * aj**2)
                H_ab = np.sum(w * aj * (th - bj)) - np.sum(res * (th - bj)) * 0  # yaklasik
                Hab = -np.sum(w * (-aj) * (th - bj))
                H = np.array([[H_aa, Hab], [Hab, H_bb]])
                grad = np.array([g_a, g_b])
                try:
                    step = np.linalg.solve(H, grad)
                except np.linalg.LinAlgError:
                    break
                aj_new = aj - step[0]; bj_new = bj - step[1]
                aj_new = np.clip(aj_new, 0.05, 5.0); bj_new = np.clip(bj_new, -6, 6)
                if abs(aj_new - aj) + abs(bj_new - bj) < 1e-6:
                    aj, bj = aj_new, bj_new; break
                aj, bj = aj_new, bj_new
            a[j], b[j] = aj, bj
        if abs(ll - ll_old) < tol:
            break
        ll_old = ll
    theta_eap = post @ th
    return dict(a=a, b=b, theta=theta_eap, loglik=ll, n_iter=it + 1)


def fit_irtree_2d(X1, X2, n_grid=25, lo=-4.0, hi=4.0, max_iter=200, tol=1e-3):
    """
    Iki korelasyonlu boyutlu IRTree:
      X1 (N,J): RG dugumu (hizli tahmin 0/1) -> boyut 1
      X2 (N,J): ACC dugumu (dogruluk 0/1, effortful) -> boyut 2
    Basit yapi: her madde tek boyuta yuklenir; boyutlar korelasyonlu N2(0,Sigma).
    Doner: a1,b1 (RG), a2,b2 (ACC), theta1,theta2 (EAP), rho.
    """
    X1 = np.asarray(X1, float); X2 = np.asarray(X2, float)
    N, J = X1.shape
    o1 = ~np.isnan(X1); o2 = ~np.isnan(X2)
    X1f = np.where(o1, X1, 0.0); X2f = np.where(o2, X2, 0.0)
    g = np.linspace(lo, hi, n_grid)
    T1, T2 = np.meshgrid(g, g, indexing='ij')
    t1 = T1.ravel(); t2 = T2.ravel(); Q = t1.size
    a1 = np.ones(J); b1 = np.clip(-(np.nanmean(X1, 0) - .5) * 2, -3, 3)
    a2 = np.ones(J); b2 = np.clip(-(np.nanmean(X2, 0) - .5) * 2, -3, 3)
    rho = 0.0
    for it in range(max_iter):
        # a priori (korelasyonlu normal)
        det = 1 - rho**2
        quad = (t1**2 - 2 * rho * t1 * t2 + t2**2) / det
        logW = -0.5 * quad - 0.5 * np.log(det)
        # olabilirlik
        P1 = np.clip(_p2pl(a1[None, :], b1[None, :], t1[:, None]), 1e-9, 1 - 1e-9)  # (Q,J)
        P2 = np.clip(_p2pl(a2[None, :], b2[None, :], t2[:, None]), 1e-9, 1 - 1e-9)
        LL = (X1f @ np.log(P1).T + (o1 * (1 - X1f)) @ np.log(1 - P1).T)
        LL += (X2f @ np.log(P2).T + (o2 * (1 - X2f)) @ np.log(1 - P2).T)
        LL += logW[None, :]
        m = LL.max(1, keepdims=True)
        post = np.exp(LL - m); post /= post.sum(1, keepdims=True)     # (N,Q)
        # M-step maddeler
        def upd(a, b, Xf, o, tg):
            nq = post.T @ o.astype(float); rq = post.T @ (Xf * o)
            for j in range(J):
                aj, bj = a[j], b[j]
                for _ in range(30):
                    p = np.clip(_p2pl(aj, bj, tg), 1e-9, 1 - 1e-9)
                    w = nq[:, j] * p * (1 - p); res = rq[:, j] - nq[:, j] * p
                    g_a = np.sum((tg - bj) * res); g_b = np.sum(-aj * res)
                    H_aa = -np.sum(w * (tg - bj)**2); H_bb = -np.sum(w * aj**2)
                    Hab = np.sum(w * aj * (tg - bj))
                    H = np.array([[H_aa, Hab], [Hab, H_bb]]); grad = np.array([g_a, g_b])
                    try: step = np.linalg.solve(H, grad)
                    except np.linalg.LinAlgError: break
                    aj = np.clip(aj - step[0], 0.05, 5.0); bj = np.clip(bj - step[1], -6, 6)
                a[j], b[j] = aj, bj
            return a, b
        a1, b1 = upd(a1, b1, X1f, o1, t1)
        a2, b2 = upd(a2, b2, X2f, o2, t2)
        # M-step rho (varyanslar 1 sabit)
        E_t1t2 = np.sum(post * (t1 * t2)[None, :]) / N
        rho_new = np.clip(E_t1t2, -0.95, 0.95)
        conv = abs(rho_new - rho)
        rho = rho_new
        if conv < tol and it > 5:
            break
    th1 = post @ t1; th2 = post @ t2
    return dict(a1=a1, b1=b1, a2=a2, b2=b2, theta1=th1, theta2=th2, rho=rho, n_iter=it + 1)
