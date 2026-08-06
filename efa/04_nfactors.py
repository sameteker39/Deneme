import pickle, numpy as np, pandas as pd, sys
sys.path.insert(0,'/home/user/Deneme/efa')
from efa_core import smooth, map_test, efa, fit_indices
d=pickle.load(open('R.pkl','rb')); R=smooth(d['R']); cols=d['cols']; df=d['df']; N=int(round(d['n_mean']))
p=R.shape[0]
print("N (ort. ikili) =",N)
w=np.linalg.eigvalsh(R)[::-1]
print("\n--- Ozdegerler (polikorik) ---")
print(pd.DataFrame({'Faktor':range(1,16),'Ozdeger':w[:15].round(3),
      'Aciklanan %':(100*w[:15]/p).round(2),'Kumulatif %':(100*np.cumsum(w[:15])/p).round(2)}).to_string(index=False))
print("Kaiser (>1):",(w>1).sum())

f2,f4=map_test(R,12)
print("\n--- Velicer MAP ---")
print(pd.DataFrame({'k':range(0,13),'MAP(kare)':f2.round(5),'MAP(4.kuvvet)':f4.round(5)}).to_string(index=False))
print("MAP min (kare) -> k =",int(np.argmin(f2)), "| 4.kuvvet -> k =",int(np.argmin(f4)))

# Ampirik Kaiser Olcutu (Braeken & van Assen, 2017)
J,Nn=p,N
ek=[]; ref=[]
lam=w.copy(); l_prev=1.0
for j in range(1,13):
    ref_j=max(((1+np.sqrt(J/Nn))**2)*(J-np.sum(lam[:j-1]))/(J-j+1), 1)
    ref.append(ref_j); ek.append(lam[j-1]>ref_j)
print("\n--- Ampirik Kaiser Olcutu ---")
print(pd.DataFrame({'k':range(1,13),'Ozdeger':lam[:12].round(3),'Referans':np.round(ref,3),'Tut':ek}).to_string(index=False))
print("EKC -> k =", int(np.argmin(ek)) if False in ek else 12)
