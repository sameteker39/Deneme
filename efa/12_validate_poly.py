import numpy as np, pandas as pd, sys, warnings; warnings.filterwarnings('ignore')
sys.path.insert(0,'/home/user/Deneme/efa')
from polychor import polychoric_matrix
rng=np.random.default_rng(11)
print("Polikorik tahminci dogrulamasi: surekli veriden kesilerek uretilen sirali veri\n")
for rho in [.3,.5,.7]:
    for cuts,lbl in [([-1.5,-0.5,0.5,1.5],'simetrik'), ([-2.2,-1.6,-1.0,-0.2],'carpik(tavan)'),
                     ([-2.8,-2.3,-1.9,-1.2],'asiri carpik')]:
        S=np.array([[1,rho],[rho,1]])
        X=rng.multivariate_normal([0,0],S,3000)
        O=np.digitize(X[:,0],cuts)+1; P=np.digitize(X[:,1],cuts)+1
        df=pd.DataFrame({'a':O.astype(float),'b':P.astype(float)})
        R,_=polychoric_matrix(df)
        pear=np.corrcoef(O,P)[0,1]
        mn=min(pd.Series(O).value_counts().min(),pd.Series(P).value_counts().min())
        print("gercek r=%.2f %-14s -> polikorik=%.3f (hata %+.3f) | pearson=%.3f | min hucre n=%d"%(
            rho,lbl,R[0,1],R[0,1]-rho,pear,mn))
print("\n--- Kucuk orneklem + seyrek hucre (n=310, asiri carpik) ---")
for rep in range(5):
    S=np.array([[1,.5],[.5,1]]); X=rng.multivariate_normal([0,0],S,310)
    cuts=[-2.8,-2.4,-2.0,-1.3]
    O=np.digitize(X[:,0],cuts)+1; P=np.digitize(X[:,1],cuts)+1
    df=pd.DataFrame({'a':O.astype(float),'b':P.astype(float)})
    R,_=polychoric_matrix(df)
    mn=min(pd.Series(O).value_counts().min(),pd.Series(P).value_counts().min())
    print("  gercek .50 -> polikorik=%.3f  pearson=%.3f  min hucre=%d"%(R[0,1],np.corrcoef(O,P)[0,1],mn))
