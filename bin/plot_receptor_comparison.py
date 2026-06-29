#!/usr/bin/env python3
"""
plot_receptor_comparison.py — Per-receptor overlaid comparison of MD trajectories.

For each receptor (ACR157, QCL936, XP273) generates a figure where each subplot
contains one line per treatment (GORE4 / SKTI / BEN) on the same axes.

Layout: 3 rows × 2 columns
  [Backbone RMSD]  [Ligand RMSD]
  [Contacts]       [H-bonds]
  [Radius of Gyration] [Ligand SASA]

Usage (from repo root on server):
  python3 bin/plot_receptor_comparison.py
  python3 bin/plot_receptor_comparison.py --receptors ACR157 QCL936
  python3 bin/plot_receptor_comparison.py --outdir figures/ --window-ns 5
"""

import argparse
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib as mpl

ROOT = Path(__file__).parent.parent

mpl.rcParams.update({
    'font.size':      10,
    'axes.titlesize': 11,
    'axes.labelsize': 10,
    'legend.fontsize': 8,
    'figure.dpi':    150,
})

# ── Colors and labels ──────────────────────────────────────────────────────

TREATMENT_STYLE = {
    "GORE4": {"color": "#1a6faf", "lw": 1.8},
    "SKTI":  {"color": "#2ca02c", "lw": 1.8},
    "BEN":   {"color": "#d62728", "lw": 1.8},
}

TREATMENT_FULL = {
    "GORE4": "GORE4",
    "SKTI":  "SKTI",
    "BEN":   "BEN (benzamidine)",
}

# ── Metrics (xvg filename, y-axis label, subplot title) ───────────────────

METRICS = [
    ("rmsd_backbone.xvg", "RMSD (nm)",  "Backbone RMSD"),
    ("rmsd_ligante.xvg",  "RMSD (nm)",  "Ligand RMSD"),
    ("numcont.xvg",       "N. atoms",   "Receptor–ligand contacts"),
    ("hbond.xvg",         "N. H-bonds", "Hydrogen bonds"),
    ("gyrate.xvg",        "Rg (nm)",    "Radius of Gyration"),
    ("sasa_ligante.xvg",  "SASA (nm²)", "Ligand SASA"),
]

RECEPTORS = ["ACR157", "QCL936", "XP273"]


# ── XVG reader ─────────────────────────────────────────────────────────────

def load_xvg(path: Path) -> np.ndarray | None:
    if not path or not path.exists():
        return None
    data = []
    for ln in open(path):
        if ln.startswith(("#", "@")):
            continue
        try:
            data.append([float(x) for x in ln.split()])
        except ValueError:
            pass
    arr = np.array(data) if data else None
    return arr if arr is not None and arr.ndim == 2 and arr.shape[1] >= 2 else None


# ── Rolling statistics ─────────────────────────────────────────────────────

def rolling_stats(values: np.ndarray, window: int):
    n = len(values)
    w = max(1, min(window, n))
    v = values.astype(float)
    pad = w // 2
    vp = np.pad(v, (pad, w - pad - 1), mode="edge")
    cs = np.cumsum(np.insert(vp, 0, 0.0))
    mean = (cs[w:] - cs[:-w]) / w
    cs2 = np.cumsum(np.insert(vp ** 2, 0, 0.0))
    var = (cs2[w:] - cs2[:-w]) / w - mean ** 2
    std = np.sqrt(np.maximum(var, 0.0))
    return mean[:n], std[:n]


# ── Path discovery ─────────────────────────────────────────────────────────

def find_analise(receptor: str, treatment: str) -> Path | None:
    """Locate the analise/ directory for a receptor/treatment pair."""
    if treatment == "BEN":
        candidates = [
            ROOT / "results" / f"{receptor}-BEN" / f"{receptor}-BEN" / "analise",
            ROOT / "results" / f"{receptor}-BEN" / "analise",
        ]
        for c in candidates:
            if c.is_dir():
                return c
        return None

    # GORE4 / SKTI: search in *_NEW directories
    suffix_patterns = [f"{receptor}-{treatment}_NEW", f"{receptor}-{treatment}-NEW"]
    for pat in suffix_patterns:
        for base in sorted(ROOT.glob(pat)):
            if not base.is_dir():
                continue
            for analise in sorted(base.glob("**/analise")):
                if analise.is_dir() and (analise / "rmsd_backbone.xvg").exists():
                    return analise
    return None


# ── Per-receptor figure ────────────────────────────────────────────────────

def plot_receptor(receptor: str, outdir: Path, window_ns: float = 5.0):
    print(f"\n[{receptor}]")

    # Discover data directories
    analise = {}
    for treatment in ("GORE4", "SKTI", "BEN"):
        d = find_analise(receptor, treatment)
        if d:
            print(f"  {treatment:5s} → {d.relative_to(ROOT)}")
            analise[treatment] = d
        else:
            print(f"  {treatment:5s} → not found")

    if not analise:
        print("  [SKIP] no data found")
        return

    n_rows, n_cols = 3, 2
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(13, n_rows * 3.8),
                             constrained_layout=True)
    fig.suptitle(
        f"{receptor} — MD trajectories: GORE4 vs SKTI vs BEN",
        fontsize=13, fontweight="bold",
    )

    for idx, (xvg_fname, ylabel, title) in enumerate(METRICS):
        ax = axes[idx // n_cols, idx % n_cols]
        plotted = False

        for treatment, adir in analise.items():
            data = load_xvg(adir / xvg_fname)
            if data is None:
                continue

            t = data[:, 0]
            v = data[:, 1]
            dt = float(t[1] - t[0]) if len(t) > 1 else 0.01
            win = max(3, int(window_ns / dt))
            rm, rs = rolling_stats(v, win)

            style = TREATMENT_STYLE[treatment]
            color = style["color"]
            lw    = style["lw"]
            label = f"{TREATMENT_FULL[treatment]}  (mean {v.mean():.3f})"

            ax.plot(t, v,         lw=0.3,  color=color, alpha=0.15)
            ax.fill_between(t, rm - rs, rm + rs, alpha=0.12, color=color)
            ax.plot(t, rm,        lw=lw,   color=color, label=label)
            plotted = True

        ax.set_xlabel("Time (ns)")
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.grid(alpha=0.25)
        if plotted:
            ax.legend(framealpha=0.85)

    p = outdir / f"{receptor}_comparison.png"
    fig.savefig(p, dpi=200, bbox_inches="tight")
    plt.close(fig)
    size_kb = p.stat().st_size // 1024
    print(f"  [OK] {p.name}  ({size_kb} KB)")


# ── Main ──────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--receptors", nargs="+", default=RECEPTORS,
                    help=f"Receptors to process (default: {RECEPTORS})")
    ap.add_argument("--outdir",    default="figures",
                    help="Output directory (default: figures/)")
    ap.add_argument("--window-ns", type=float, default=5.0,
                    help="Rolling mean window in ns (default: 5.0)")
    args = ap.parse_args()

    outdir = ROOT / args.outdir
    outdir.mkdir(parents=True, exist_ok=True)
    print(f"Output: {outdir}")

    for receptor in args.receptors:
        plot_receptor(receptor, outdir, window_ns=args.window_ns)

    print("\n[OK] done")


if __name__ == "__main__":
    main()
