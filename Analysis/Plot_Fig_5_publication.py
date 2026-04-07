"""
Plot_Fig_5_publication.py
─────────────────────────
Generates publication-quality Figure 5 (distribution of simulation
runtimes) from the benchmark timing CSV.

Outputs (written to the same directory as this script):
  Figure_5.eps   – vector graphic (journal preferred format)
  Figure_5.tif   – 600 dpi RGB TIFF, LZW-compressed (combination artwork)
  Figure_5.pdf   – vector PDF (transparency-safe backup)

Journal specification applied:
  - Vector graphics: EPS
  - Combination artwork (histogram + reference lines): TIFF ≥ 600 dpi
  - Fonts embedded (Type 42 / TrueType)
  - Color consistent with Figures 3 and 4 (Wong 2011 palette)

Usage:
  python Analysis/Plot_Fig_5_publication.py

Requirements: see Analysis/requirements.txt (matplotlib, pandas, numpy,
              scipy, Pillow)
"""

import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy import stats

# ── Paths ─────────────────────────────────────────────────────────────────────
SCRIPT_DIR = Path(__file__).parent.resolve()
DATA_FILE  = SCRIPT_DIR.parent / "data" / "Figure_5" / "Simulation_Runtimes_Minutes.csv"
OUT_DIR    = SCRIPT_DIR

if not DATA_FILE.exists():
    sys.exit(f"ERROR: data file not found: {DATA_FILE}\n"
             f"       Check that data/Figure_5/ exists in the repository root.")

# ── Load data ─────────────────────────────────────────────────────────────────
df       = pd.read_csv(DATA_FILE)
runtimes = pd.to_numeric(df["Runtime_minutes"]).values

n      = len(runtimes)
mean   = runtimes.mean()
sd     = runtimes.std(ddof=1)
median = np.median(runtimes)

# ── Global matplotlib settings ────────────────────────────────────────────────
plt.rcParams.update({
    "font.family":        "serif",
    "font.serif":         ["Times New Roman", "DejaVu Serif"],
    "font.size":          10,
    "axes.labelsize":     11,
    "axes.titlesize":     11,
    "xtick.labelsize":    9,
    "ytick.labelsize":    9,
    "legend.fontsize":    9,
    "legend.frameon":     True,
    "legend.framealpha":  0.9,
    "legend.edgecolor":   "0.7",
    "axes.linewidth":     0.8,
    "xtick.major.width":  0.8,
    "ytick.major.width":  0.8,
    "xtick.minor.width":  0.5,
    "ytick.minor.width":  0.5,
    "xtick.direction":    "in",
    "ytick.direction":    "in",
    "xtick.top":          True,
    "ytick.right":        True,
    "lines.linewidth":    1.5,
    "pdf.fonttype":       42,
    "ps.fonttype":        42,
    "savefig.bbox":       "tight",
    "savefig.pad_inches": 0.05,
})

# ── Figure ────────────────────────────────────────────────────────────────────
# Double-column width (174 mm = 6.85 in); height by 3:2 ratio
FIG_W, FIG_H = 6.85, 6.85 / 1.5

fig, ax = plt.subplots(figsize=(FIG_W, FIG_H))

# Histogram — 30 bins
ax.hist(
    runtimes, bins=30,
    color="#0072B2",        # Wong blue — consistent with Figures 3 & 4
    edgecolor="white",
    linewidth=0.5,
    zorder=3
)

# Mean reference line
ax.axvline(
    mean, color="#D55E00", linewidth=1.4, linestyle="--", zorder=4,
    label=f"Mean = {mean:.1f}\u202fmin"
)

# ±1 SD shaded region
ax.axvspan(
    mean - sd, mean + sd,
    alpha=0.12, color="#D55E00", linewidth=0, zorder=2,
    label=f"\u00b1\u202f1\u202fSD  ({mean - sd:.1f}\u2013{mean + sd:.1f}\u202fmin)"
)

ax.set_xlabel("Runtime (minutes)", labelpad=4)
ax.set_ylabel("Number of simulations", labelpad=4)

ax.xaxis.set_major_locator(ticker.MultipleLocator(5))
ax.xaxis.set_minor_locator(ticker.MultipleLocator(1))
ax.yaxis.set_major_locator(ticker.MaxNLocator(integer=True, nbins=6))
ax.yaxis.set_minor_locator(ticker.AutoMinorLocator(2))

ax.set_xlim(runtimes.min() - 1, runtimes.max() + 1)
ax.set_ylim(bottom=0)

leg = ax.legend(loc="upper right", borderpad=0.6, handlelength=1.8,
                labelspacing=0.4, handletextpad=0.5)
leg.get_frame().set_linewidth(0.6)

fig.tight_layout()

# ── Save ──────────────────────────────────────────────────────────────────────
def save(fig, stem):
    eps  = OUT_DIR / f"{stem}.eps"
    tiff = OUT_DIR / f"{stem}.tif"
    pdf  = OUT_DIR / f"{stem}.pdf"
    fig.savefig(eps,  format="eps")
    fig.savefig(pdf,  format="pdf")
    fig.savefig(tiff, format="tiff", dpi=600,
                pil_kwargs={"compression": "tiff_lzw"})
    # Convert to RGB — alpha channel not accepted by all journal systems
    from PIL import Image
    img = Image.open(tiff).convert("RGB")
    img.save(tiff, format="TIFF", dpi=(600, 600), compression="tiff_lzw")
    print(f"  EPS  → {eps}")
    print(f"  TIFF → {tiff}  (600 dpi, RGB, LZW)")
    print(f"  PDF  → {pdf}")

print("Saving Figure 5 …")
save(fig, "Figure_5")
plt.close(fig)
print("Done.")
