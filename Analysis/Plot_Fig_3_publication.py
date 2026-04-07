"""
Plot_Fig_3_publication.py
─────────────────────────
Generates publication-quality Figure 3 (mean species counts over time)
from the four ensemble CSV files produced by the analysis pipeline.

Outputs (written to the same directory as this script):
  Figure_3.eps   – vector graphic (journal preferred format)
  Figure_3.tif   – 600 dpi RGB TIFF, LZW-compressed (combination artwork)
  Figure_3.pdf   – vector PDF (transparency-safe backup)

Journal specification applied:
  - Vector graphics: EPS
  - Combination artwork (line + shading): TIFF ≥ 600 dpi
  - Fonts embedded (Type 42 / TrueType)
  - Colorblind-safe palette (Wong 2011)
  - Shaded bands: mean ± SEM across n = 5 independent simulations

Usage:
  python Analysis/Plot_Fig_3_publication.py

Requirements: see Analysis/requirements.txt (matplotlib, pandas, numpy, Pillow)
"""

import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# ── Paths ─────────────────────────────────────────────────────────────────────
SCRIPT_DIR = Path(__file__).parent.resolve()
DATA_DIR   = SCRIPT_DIR.parent / "data" / "Figure_3"
OUT_DIR    = SCRIPT_DIR

# Verify data files exist before proceeding
DATA_FILES = {
    "AggregateProne": DATA_DIR / "Appending_AggregateProne_Count.csv",
    "Native":         DATA_DIR / "Appending_Native_Count.csv",
    "Oligomer":       DATA_DIR / "Appending_Oligomer_Count.csv",
    "Aggregate":      DATA_DIR / "Appending_Aggregate_Count.csv",
}
for label, path in DATA_FILES.items():
    if not path.exists():
        sys.exit(f"ERROR: data file not found: {path}\n"
                 f"       Check that data/Figure_3/ exists in the repository root.")

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

# ── Colorblind-safe palette (Wong 2011) ───────────────────────────────────────
COLORS = {
    "Native":         "#0072B2",   # blue
    "AggregateProne": "#E69F00",   # orange
    "Oligomer":       "#009E73",   # green
    "Aggregate":      "#D55E00",   # vermilion
}
LINE_STYLES = {
    "Native":         "-",
    "AggregateProne": "--",
    "Oligomer":       "-.",
    "Aggregate":      ":",
}

# ── Data loading ──────────────────────────────────────────────────────────────
MAX_TIMESTEP = 300

def load(path, max_t=MAX_TIMESTEP):
    """Return (timesteps, mean, SEM, n_simulations) filtered to max_t."""
    df   = pd.read_csv(path)
    df   = df[df.iloc[:, 0] <= max_t].copy()
    t    = df.iloc[:, 0].values
    cols = df.iloc[:, 1:]
    n    = cols.shape[1]
    return t, cols.mean(axis=1).values, cols.sem(axis=1).values, n

t, m_nat, s_nat, n = load(DATA_FILES["Native"])
_, m_ap,  s_ap,  _ = load(DATA_FILES["AggregateProne"])
_, m_oli, s_oli, _ = load(DATA_FILES["Oligomer"])
_, m_agg, s_agg, _ = load(DATA_FILES["Aggregate"])

# ── Figure ────────────────────────────────────────────────────────────────────
# Double-column width (174 mm = 6.85 in); height by 3:2 ratio
FIG_W, FIG_H = 6.85, 6.85 / 1.5

fig, ax = plt.subplots(figsize=(FIG_W, FIG_H))

series = [
    ("Native",         t, m_nat, s_nat),
    ("AggregateProne", t, m_ap,  s_ap),
    ("Oligomer",       t, m_oli, s_oli),
    ("Aggregate",      t, m_agg, s_agg),
]
for label, ts, mean, sem in series:
    ax.plot(ts, mean,
            color=COLORS[label], linestyle=LINE_STYLES[label],
            linewidth=1.6, label=label, zorder=3)
    ax.fill_between(ts, mean - sem, mean + sem,
                    color=COLORS[label], alpha=0.15,
                    linewidth=0, zorder=2)

ax.set_xlim(0, MAX_TIMESTEP)
ax.set_ylim(bottom=0)
ax.xaxis.set_major_locator(ticker.MultipleLocator(50))
ax.xaxis.set_minor_locator(ticker.MultipleLocator(10))
ax.yaxis.set_major_locator(ticker.MaxNLocator(integer=True, nbins=6))
ax.yaxis.set_minor_locator(ticker.AutoMinorLocator(2))

ax.set_xlabel("Timestep", labelpad=4)
ax.set_ylabel(f"Mean monomer count (\u00b1\u202fSEM, n\u202f=\u202f{n})", labelpad=4)

leg = ax.legend(loc="upper right", borderpad=0.6, handlelength=2.0,
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

print("Saving Figure 3 …")
save(fig, "Figure_3")
plt.close(fig)
print("Done.")
