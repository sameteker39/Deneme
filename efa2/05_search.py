import pickle, numpy as np, pandas as pd, sys, warnings; warnings.filterwarnings('ignore')
sys.path.insert(0,'.')
from polychor import polychoric_matrix
from efa_core import smooth, efa, communalities
from final_engine import purify_at_k, solution, CRIT
pd.set_option('display.width',250)
d=pickle.load(open('D.pkl','rb')); R=d['R']; cols=d['cols']; N=202
out={}
for k in [2,3,4,5]:
    print("\n"+"#"*80); print("# k=%d icin arindirma"%k); print("#"*80)
    r=purify_at_k(R,cols,k,N)
    dd=r['d']; vc=dd['Faktor'].value_counts()
    out[k]=r
    print("  -> %d madde | faktor buyuklukleri %s | yakinsadi=%s"%(
        len(r['cols']), vc.sort_index().to_dict(), r['ok']))
    print("  atilan:", ", ".join(x[0] for x in r['log']) if r['log'] else "yok")
pickle.dump(out,open('purified.pkl','wb'))
