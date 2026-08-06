import pickle, numpy as np, pandas as pd, sys, warnings; warnings.filterwarnings('ignore')
sys.path.insert(0,'.')
from items import lab, TXT, ITEM_NO
from efa_core import smooth, efa, communalities, map_test
pd.set_option('display.width',250)
d=pickle.load(open('D.pkl','rb')); R=smooth(d['R'],0.02); cols=d['cols']
w=np.linalg.eigvalsh(R)[::-1]
print("Ozdegerler:",np.round(w[:9],2)," Kaiser>1:",(w>1).sum())
f2,f4=map_test(R,10); print("MAP -> k=%d / %d"%(np.argmin(f2),np.argmin(f4)))
for k in [4,5]:
    L,Phi=efa(R,k,'oblimin','minres'); h2=communalities(L,Phi)
    print("\n"+"="*118); print("k=%d (31 madde)"%k); print("="*118)
    A=np.abs(L); pri=A.argmax(1)
    for f in range(k):
        idx=[i for i in range(31) if pri[i]==f]
        idx.sort(key=lambda i:-A[i,f])
        print("\n--- F%d ---"%(f+1))
        for i in idx:
            oth="".join(" F%d=%+.2f"%(j+1,L[i,j]) for j in range(k) if j!=f and abs(L[i,j])>=.30)
            print("   %+.2f  %-58s h2=%.2f%s"%(L[i,f],lab(cols[i],56),h2[i],oth))
    print("\nPhi:\n",pd.DataFrame(Phi.round(2),index=['F%d'%(i+1) for i in range(k)],
                                  columns=['F%d'%(i+1) for i in range(k)]).to_string())
