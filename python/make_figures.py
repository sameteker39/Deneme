"""
APA uyumlu figurler (Ingilizce) — gercek mirt ciktilarindan.
Girdi:  outputs/T1_madde_parametreleri.csv, outputs/T5_duyarlilik.csv
Cikti:  figs/Figure1_item_agreement.png, Figure2_prevalence_validity.png,
        Figure3_sensitivity.png  (600 dpi)
Tasarim: gri-tonlamaya uygun, renk-korligune guvenli palet (dogrulanmis),
kimlik cizgisi, minimal eksen, panel harfleri.
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import rcParams

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
IN = os.path.join(ROOT, "outputs")
OUT = os.path.join(ROOT, "figs")
os.makedirs(OUT, exist_ok=True)

ip = pd.read_csv(os.path.join(IN, "T1_madde_parametreleri.csv"))
t5 = pd.read_csv(os.path.join(IN, "T5_duyarlilik.csv"))

rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 9.5,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.edgecolor": "#444444", "axes.linewidth": 0.8,
    "xtick.color": "#444444", "ytick.color": "#444444",
    "axes.labelcolor": "#222222", "text.color": "#222222",
    "figure.facecolor": "white", "axes.facecolor": "white",
})
POINT, IDENT = "#2a78d6", "#9a9a9a"
C_RG, C_SOL, C_RHO = "#eda100", "#1baf7a", "#4a3aa7"
labels = ["NT10\n(wins.)", "NT10\n(unwins.)", "NT5", "Fixed 5 s", "MRTQ"]
xr = np.arange(len(labels))


def panel(ax, l):
    ax.text(-0.16, 1.06, l, transform=ax.transAxes, fontsize=11,
            fontweight="bold", va="top", ha="left")


def scatter(ax, x, y, xl, yl):
    lo, hi = min(x.min(), y.min()), max(x.max(), y.max())
    pad = (hi - lo) * 0.06; lo -= pad; hi += pad
    ax.plot([lo, hi], [lo, hi], ls=(0, (5, 4)), color=IDENT, lw=1, zorder=1)
    ax.scatter(x, y, s=26, facecolor=POINT, edgecolor="white",
               linewidth=0.5, alpha=0.85, zorder=3)
    ax.text(0.05, 0.93, f"$r$ = {np.corrcoef(x, y)[0,1]:.3f}",
            transform=ax.transAxes, fontsize=9.5, va="top")
    ax.set_xlim(lo, hi); ax.set_ylim(lo, hi); ax.set_aspect("equal", "box")
    ax.set_xlabel(xl); ax.set_ylabel(yl); ax.tick_params(length=3)


# ===== FIGURE 1 =====
fig, ax = plt.subplots(2, 2, figsize=(6.8, 6.8))
scatter(ax[0, 0], ip.b_std, ip.b_em, "Difficulty (Standard 2PL)", "Difficulty (EM-IRT)")
scatter(ax[0, 1], ip.b_em, ip.b_tree, "Difficulty (EM-IRT)", "Difficulty (IRTree)")
scatter(ax[1, 0], ip.a_std, ip.a_em, "Discrimination (Standard 2PL)", "Discrimination (EM-IRT)")
scatter(ax[1, 1], ip.a_em, ip.a_tree, "Discrimination (EM-IRT)", "Discrimination (IRTree)")
for a, l in zip(ax.ravel(), ["(a)", "(b)", "(c)", "(d)"]):
    panel(a, l)
fig.tight_layout(w_pad=2.5, h_pad=2.2)
fig.savefig(os.path.join(OUT, "Figure1_item_agreement.png"), dpi=600, bbox_inches="tight")

# ===== FIGURE 2 =====
fig, ax = plt.subplots(1, 2, figsize=(7.6, 3.5))
ax[0].bar(xr, t5["RG_orani"], width=0.6, color=POINT, edgecolor="white", linewidth=0.8)
for i, v in enumerate(t5["RG_orani"]):
    ax[0].text(i, v + 0.004, f"{v:.3f}", ha="center", fontsize=8)
ax[0].set_xticks(xr); ax[0].set_xticklabels(labels, fontsize=8)
ax[0].set_ylabel("Rapid-guessing rate"); ax[0].set_ylim(0, 0.17)
ax[0].tick_params(length=3); panel(ax[0], "(a)")
w = 0.38
ax[1].bar(xr - w / 2, t5["RG_dogruluk"], w, label="Rapid guess", color=C_RG, edgecolor="white", linewidth=0.8)
ax[1].bar(xr + w / 2, t5["caba_dogruluk"], w, label="Solution", color=C_SOL, edgecolor="white", linewidth=0.8)
ax[1].axhline(0.25, ls=(0, (5, 4)), color="#d1495b", lw=1)
ax[1].text(len(labels) - 1.5, 0.27, "chance = .25", color="#d1495b", fontsize=8)
ax[1].set_xticks(xr); ax[1].set_xticklabels(labels, fontsize=8)
ax[1].set_ylabel("Proportion correct"); ax[1].set_ylim(0, 0.65)
ax[1].tick_params(length=3); panel(ax[1], "(b)")
ax[1].legend(frameon=False, fontsize=8, loc="upper right", handlelength=1.1)
fig.tight_layout(w_pad=2.5)
fig.savefig(os.path.join(OUT, "Figure2_prevalence_validity.png"), dpi=600, bbox_inches="tight")

# ===== FIGURE 3 =====
fig, ax = plt.subplots(1, 2, figsize=(7.6, 3.5))
mk = {"a_kor_ANA": ("Discrimination", "o", "#2a78d6"),
      "b_kor_ANA": ("Difficulty", "s", "#1baf7a"),
      "theta_kor_ANA": ("Ability", "^", "#eda100")}
for col, (lab, m, c) in mk.items():
    ax[0].plot(xr, t5[col], marker=m, color=c, lw=1.4, ms=6, label=lab,
               markeredgecolor="white", markeredgewidth=0.6)
ax[0].set_xticks(xr); ax[0].set_xticklabels(labels, fontsize=8)
ax[0].set_ylabel("Correlation with primary NT10"); ax[0].set_ylim(0.90, 1.005)
ax[0].tick_params(length=3); panel(ax[0], "(a)")
ax[0].legend(frameon=False, fontsize=8, loc="lower left", handlelength=1.4)
ax[1].bar(xr, t5["rho_RG_ACC"], width=0.6, color=C_RHO, edgecolor="white", linewidth=0.8)
for i, v in enumerate(t5["rho_RG_ACC"]):
    ax[1].text(i, v - 0.015, f"{v:.2f}", ha="center", va="top", fontsize=8, color="#222222")
ax[1].axhline(0, color="#444", lw=0.8)
ax[1].set_xticks(xr); ax[1].set_xticklabels(labels, fontsize=8)
ax[1].set_ylabel("Disengagement–ability correlation"); ax[1].set_ylim(-0.52, 0.03)
ax[1].tick_params(length=3); panel(ax[1], "(b)")
fig.tight_layout(w_pad=2.5)
fig.savefig(os.path.join(OUT, "Figure3_sensitivity.png"), dpi=600, bbox_inches="tight")
print("Figurler yazildi:", OUT)
