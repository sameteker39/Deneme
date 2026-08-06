import pickle, numpy as np, pandas as pd, sys, warnings; warnings.filterwarnings('ignore')
sys.path.insert(0,'/home/user/Deneme/efa')
d=pickle.load(open('Rc.pkl','rb')); df=d['df']; cols=d['cols']
b=pickle.load(open('boot.pkl','rb')); sd=b['sd']; iu=b['iu']
# madde basina ortalama bootstrap SE
SE=np.zeros((31,31)); SE[iu]=sd; SE=SE+SE.T
mse=SE.sum(0)/30
t=pd.DataFrame({
 'Ort':df.mean().round(2), 'SS':df.std().round(2),
 'Tepe_kat_%':(100*df.apply(lambda s: s.value_counts().max()/s.notna().sum())).round(1),
 'Carpiklik':df.skew().round(2), 'Basiklik':df.kurt().round(2),
 'Min_hucre':df.apply(lambda s: int(s.value_counts().reindex([1,2,3,4,5]).fillna(0).min())),
 'Boot_SE':mse.round(3)}, index=cols)
t['Asiri_tavan']=(t['Tepe_kat_%']>=85)&(t['Carpiklik'].abs()>3)
print(t.sort_values('Tepe_kat_%',ascending=False).to_string())
print("\nOlcut: tepe kategoride >=%85 VE |carpiklik|>3")
print("Isaretlenen maddeler:", t.index[t['Asiri_tavan']].tolist())
print("\nBu maddelerin ort. bootstrap SE'si: %.3f | digerleri: %.3f"%(
    mse[t['Asiri_tavan'].values].mean(), mse[~t['Asiri_tavan'].values].mean()))
