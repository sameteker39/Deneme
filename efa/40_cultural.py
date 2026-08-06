import pickle, numpy as np, pandas as pd, sys, warnings; warnings.filterwarnings('ignore')
sys.path.insert(0,'.')
from items import TXT, ITEM_NO, REVERSE
from polychor import polychoric_matrix
from efa_core import smooth
d=pickle.load(open('D.pkl','rb')); df=d['df']; dfc=d['dfc']   # M9 zaten ters kodlanmis
# ORIJINAL 32'lik numaralandirma -> veri sutunu
def col_of(orig):
    if orig==23: return None                  # kontrol maddesi, veride yok
    return 'M%d'%(orig if orig<=22 else orig-1)
CUT_ORIG=[5,6,9,13,18,19,20,23,24,25,26,32]
cut=[col_of(o) for o in CUT_ORIG]
print("Kulturel uygunluk gerekcesiyle cikarilan maddeler:")
for o in CUT_ORIG:
    c=col_of(o)
    if c is None: print("  orij.%-2d -> (kontrol maddesi, veride yok)"%o)
    else: print("  orij.%-2d -> %-4s (yeni no m%d)  %s"%(o,c,ITEM_NO[c],TXT[ITEM_NO[c]]))
cut=[c for c in cut if c]
keep=[c for c in df.columns if c not in cut]
print("\nKalan %d madde:"%len(keep))
for c in keep: print("  %-4s (orij.%2d)  %s"%(c, ITEM_NO[c] if ITEM_NO[c]<=22 else ITEM_NO[c]+1, TXT[ITEM_NO[c]]))
sub=df[keep]; tot=sub.sum(axis=1)
print("\n=== Kalan ters madde denetimi ===")
bad=[]
for c in keep:
    n=ITEM_NO[c]; r=sub[c].corr(tot-sub[c])
    if n in REVERSE: print("  %-4s (ters madde) MTK=%+.3f -> %s"%(c,r,"zaten ters kodlu, islem gerekmiyor" if r>0 else "TERS KODLANIYOR"))
    if r<0: bad.append(c)
print("Negatif MTK'li kalan madde:", bad if bad else "yok -> tum maddeler ayni yonde")
R,_=polychoric_matrix(dfc[keep])
print("\nPolikorik: min ozdeger=%.4f  max|r|=%.3f"%(np.linalg.eigvalsh(R).min(),
      np.abs(R[np.triu_indices(len(keep),1)]).max()))
pickle.dump({'keep':keep,'R':R,'sub':sub,'dfc':dfc[keep]},open('K20.pkl','wb'))
