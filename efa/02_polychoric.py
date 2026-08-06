import pandas as pd, numpy as np, sys, time
sys.path.insert(0,'/home/user/Deneme/efa')
from polychor import polychoric_matrix, nearest_pd
pd.set_option('display.width',250); pd.set_option('display.max_columns',60)

SRC='/root/.claude/uploads/e45f5564-4a80-51fe-a001-9fa87c7cfa26/6a6ea3ae-FAfA_veri_31madde.xlsx'
df = pd.read_excel(SRC)

# --- M9 incelemesi ---
P = df.corr()
print("M9'un diger maddelerle Pearson korelasyonlari:")
s = P['M9'].drop('M9').sort_values()
print(s.round(2).to_string())
print("\nNegatif korelasyon sayisi: %d / 30 | ortalama r = %.3f" % ((s<0).sum(), s.mean()))
