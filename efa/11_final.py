import pickle, numpy as np, pandas as pd, sys, warnings; warnings.filterwarnings('ignore')
sys.path.insert(0,'/home/user/Deneme/efa')
from final_engine import *
from efa_core import smooth, efa, fit_indices, communalities, map_test
d=pickle.load(open('Rc.pkl','rb')); R=d['Rraw']; cols=d['cols']; N=306
pd.set_option('display.width',200)
store={}
for k0 in [3,4,5]:
    print("\n"+"#"*76); print("# Baslangic k = %d"%k0); print("#"*76)
    Rc, cc, k = R, cols, k0
    for it in range(5):
        print("  [tur %d] k=%d, %d madde"%(it+1,k,len(cc)))
        r=purify_at_k(Rc,cc,k,N)
        Rc,cc=r['R'],r['cols']
        # zayif faktor kontrolu -> k dusur
        vc=r['d']['Faktor'].value_counts()
        if (vc<CRIT['minitems']).any() and k>2:
            print("     ! faktorde <3 madde -> k=%d'e dusuruluyor"%(k-1)); k-=1; continue
        # boyutlulugu yeniden degerlendir
        f2,f4=map_test(smooth(Rc,0.02),min(10,len(cc)//3))
        kmap=int(np.argmin(f2))
        if kmap!=k and 2<=kmap<=6 and it<3:
            print("     ! arindirilmis sette MAP -> k=%d, yeniden arindiriliyor"%kmap); k=kmap; continue
        break
    fi=r['fi']
    print("  SONUC: k=%d, %d madde | RMSEA=%.3f CFI=%.3f TLI=%.3f SRMR=%.3f"%(
        k,len(cc),fi['RMSEA'],fi['CFI'],fi['TLI'],fi['SRMR']))
    print("  atilan:", ", ".join(x[0] for x in r['log']))
    store[k0]=dict(r=r,k=k,cols=cc,R=Rc)
pickle.dump(store,open('final.pkl','wb'))
