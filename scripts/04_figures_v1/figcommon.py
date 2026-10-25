"""Shared style + data paths for all manuscript figures (publication, light mode, print-safe)."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, pandas as pd

U = "/mnt/user-data/uploads/Protein and RNA/"
ANA = "/home/claude/ana/"
CONF = "/home/claude/mirror/results/confirmatory/"
OUT = "/home/claude/figs/out/"

# validated categorical order (light mode); gray = null / reference
BLUE, ORANGE, AQUA, YELLOW, MAGENTA, GREEN, VIOLET, RED = (
    "#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948")
GRAY, LIGHTGRAY, INK, INK2 = "#9a9893", "#d9d7d2", "#0b0b0b", "#52514e"
SET_COLORS = {"Response": BLUE, "RNA": ORANGE, "Protein_MS": AQUA, "RPPA": YELLOW,
              "Response+RNA": MAGENTA, "Response+Protein_MS": GREEN}
SET_LABELS = {"Response": "Tahoe response", "RNA": "Baseline RNA", "Protein_MS": "Baseline protein (MS)",
              "RPPA": "Baseline protein (RPPA)", "Response+RNA": "Response + RNA",
              "Response+Protein_MS": "Response + protein (MS)"}
SETS = list(SET_COLORS)

def style():
    plt.rcParams.update({
        "font.family": "Liberation Sans", "font.size": 7, "axes.titlesize": 7.5, "axes.labelsize": 7,
        "xtick.labelsize": 6.5, "ytick.labelsize": 6.5, "legend.fontsize": 6.5, "axes.linewidth": 0.6,
        "xtick.major.width": 0.6, "ytick.major.width": 0.6, "xtick.major.size": 2.5, "ytick.major.size": 2.5,
        "axes.spines.top": False, "axes.spines.right": False, "axes.edgecolor": INK2,
        "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2, "text.color": INK,
        "axes.titleweight": "bold", "axes.titlelocation": "left", "legend.frameon": False,
        "savefig.dpi": 300, "figure.dpi": 150, "pdf.fonttype": 42, "ps.fonttype": 42,
        "lines.linewidth": 1.2, "patch.linewidth": 0.5})

def panel(ax, letter, title=None, x=-0.16, y=1.06):
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=10, fontweight="bold", va="bottom", ha="left")
    if title: ax.set_title(title, pad=4)

def rank(x): return pd.Series(x).rank().values
def spear(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float); ok = ~(np.isnan(a) | np.isnan(b))
    if ok.sum() < 3 or np.std(a[ok]) == 0 or np.std(b[ok]) == 0: return np.nan
    return np.corrcoef(rank(a[ok]), rank(b[ok]))[0, 1]

def save(fig, name):
    import os; os.makedirs(OUT, exist_ok=True)
    fig.savefig(OUT + name + ".png", dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(OUT + name + ".pdf", bbox_inches="tight", facecolor="white")
    fig.savefig(OUT + name + ".tif", dpi=300, bbox_inches="tight", facecolor="white", pil_kwargs={"compression": "tiff_lzw"})
    print("saved", name)
