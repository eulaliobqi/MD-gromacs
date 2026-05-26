#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Gera painel de figuras a partir dos .xvg do GROMACS e (opcionalmente)
dos resultados MM-GBSA do gmx_MMPBSA.

Saídas:
  painel_completo.png  — 6 painéis (sem MMGBSA)
  painel_final.png     — 8 painéis (com MMGBSA, via --mmgbsa-csv)
  + 1 PNG individual por métrica

Uso:
  plot_results.py --analise-dir <dir> [--titulo T] [--mmgbsa-csv F] [--output F]
"""
import argparse, os, sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.gridspec as gridspec

mpl.rcParams['font.size'] = 10
mpl.rcParams['axes.titlesize'] = 11


# ── Leitores ──────────────────────────────────────────────────────────────────

def load_xvg(path):
    if not path or not os.path.exists(path):
        return None
    data = []
    for ln in open(path):
        if ln.startswith(('#', '@')):
            continue
        try:
            data.append([float(x) for x in ln.split()])
        except ValueError:
            pass
    return np.array(data) if data else None


def load_mmgbsa_csv(path):
    """Lê o CSV do gmx_MMPBSA; retorna dict {colname: array} ou None."""
    if not path or not os.path.exists(path):
        return None
    try:
        with open(path) as f:
            lines = [l.rstrip() for l in f if l.strip() and not l.startswith('#')]
        if len(lines) < 2:
            return None
        header = [h.strip() for h in lines[0].split(',')]
        rows = []
        for ln in lines[1:]:
            try:
                rows.append([float(x) for x in ln.split(',')])
            except ValueError:
                pass
        if not rows:
            return None
        arr = np.array(rows)
        return {h: arr[:, i] for i, h in enumerate(header) if i < arr.shape[1]}
    except Exception:
        return None


# ── Helpers de plot ───────────────────────────────────────────────────────────

def plot_line(ax, data, ylabel, color, title=None, hline=None, xlabel="Tempo (ns)"):
    if data is None:
        ax.text(0.5, 0.5, "(arquivo não encontrado)", ha='center', va='center',
                transform=ax.transAxes, fontsize=9, color='gray')
        ax.set_xticks([])
        ax.set_yticks([])
        return
    ax.plot(data[:, 0], data[:, 1], lw=0.8, color=color)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if title:
        media = data[:, 1].mean()
        std   = data[:, 1].std()
        ax.set_title(f"{title}\n(média: {media:.3f} ± {std:.3f})")
    if hline is not None:
        ax.axhline(hline, ls='--', color='red', alpha=0.5)
    ax.grid(alpha=0.3)


def plot_rmsf(ax, data):
    if data is None:
        ax.text(0.5, 0.5, "(rmsf não encontrado)", ha='center', va='center',
                transform=ax.transAxes, fontsize=9, color='gray')
        return
    ax.bar(data[:, 0], data[:, 1], width=1.0, color='seagreen', alpha=0.8)
    ax.bar(data[-5:, 0], data[-5:, 1], width=1.0, color='crimson',
           label='Ligante (5 últ.)')
    ax.legend(fontsize=8)
    ax.set_xlabel("Resíduo")
    ax.set_ylabel("RMSF (nm)")
    ax.set_title("Flutuação por resíduo (RMSF)")
    ax.grid(alpha=0.3, axis='y')


def plot_mmgbsa_total(ax, mmgbsa):
    """ΔG total ao longo do tempo com média ± desvio."""
    total = None
    for key in ('TOTAL', 'DELTA total', 'delta_total', 'DeltaG'):
        if key in mmgbsa:
            total = mmgbsa[key]
            break
    if total is None:
        ax.text(0.5, 0.5, "(TOTAL não encontrado no CSV)", ha='center', va='center',
                transform=ax.transAxes, fontsize=9, color='gray')
        return
    frames = np.arange(len(total))
    media, std = total.mean(), total.std()
    ax.plot(frames, total, lw=0.8, color='darkred', alpha=0.7)
    ax.axhline(media, ls='--', color='red', lw=1.2,
               label=f'Média: {media:.2f} kcal/mol')
    ax.fill_between(frames, media - std, media + std,
                    alpha=0.15, color='red')
    ax.set_xlabel("Frame")
    ax.set_ylabel("ΔG bind (kcal/mol)")
    ax.set_title(f"Energia de Ligação MM-GBSA\n(média: {media:.2f} ± {std:.2f} kcal/mol)")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)


def plot_mmgbsa_components(ax, mmgbsa):
    """Componentes VdW, Elet, GB, SA como barras com média ± DP."""
    components = {}
    for key, label in [('VDWAALS', 'VdW'), ('EEL', 'Elet'),
                       ('EGB', 'GB solvat.'), ('ESURF', 'SA')]:
        if key in mmgbsa:
            components[label] = mmgbsa[key]
    if not components:
        ax.text(0.5, 0.5, "(componentes não encontrados)", ha='center', va='center',
                transform=ax.transAxes, fontsize=9, color='gray')
        return
    labels = list(components.keys())
    means  = [v.mean() for v in components.values()]
    stds   = [v.std()  for v in components.values()]
    colors = ['steelblue', 'coral', 'mediumpurple', 'gold'][:len(labels)]
    bars = ax.bar(labels, means, yerr=stds, color=colors, alpha=0.85,
                  capsize=4, edgecolor='black', linewidth=0.5)
    ax.axhline(0, color='black', lw=0.8)
    ax.set_ylabel("Energia (kcal/mol)")
    ax.set_title("Componentes MM-GBSA\n(média ± DP)")
    for bar, m in zip(bars, means):
        ax.text(bar.get_x() + bar.get_width() / 2,
                bar.get_height() + (2 if m >= 0 else -6),
                f"{m:.1f}", ha='center', fontsize=8)
    ax.grid(alpha=0.3, axis='y')


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--analise-dir", required=True)
    ap.add_argument("--titulo",      default="Dinâmica Molecular")
    ap.add_argument("--mmgbsa-csv",  default=None)
    ap.add_argument("--output",      default="painel_completo.png")
    args = ap.parse_args()

    D = args.analise_dir
    xvg = {
        "rmsd_bb":  load_xvg(os.path.join(D, "rmsd_backbone.xvg")),
        "rmsd_lig": load_xvg(os.path.join(D, "rmsd_ligante.xvg")),
        "rmsf":     load_xvg(os.path.join(D, "rmsf_residuos.xvg")),
        "rg":       load_xvg(os.path.join(D, "gyrate.xvg")),
        "ncont":    load_xvg(os.path.join(D, "numcont.xvg")),
        "hbond":    load_xvg(os.path.join(D, "hbond.xvg")),
    }

    mmgbsa = load_mmgbsa_csv(args.mmgbsa_csv) if args.mmgbsa_csv else None
    has_mmgbsa = mmgbsa is not None

    nrows = 4 if has_mmgbsa else 3
    fig, axes = plt.subplots(nrows, 2, figsize=(14, nrows * 4))
    fig.suptitle(args.titulo, fontsize=14, fontweight='bold')

    # Row 1
    media_bb = xvg["rmsd_bb"][:, 1].mean() if xvg["rmsd_bb"] is not None else None
    plot_line(axes[0, 0], xvg["rmsd_bb"],  "RMSD (nm)",  "navy",
              "RMSD do Backbone", hline=media_bb)
    plot_line(axes[0, 1], xvg["rmsd_lig"], "RMSD (nm)",  "darkorange",
              "RMSD do Ligante (peptídio)")

    # Row 2
    plot_rmsf(axes[1, 0], xvg["rmsf"])
    plot_line(axes[1, 1], xvg["rg"],       "Rg (nm)",    "purple",
              "Raio de Giro")

    # Row 3
    plot_line(axes[2, 0], xvg["ncont"],    "N. átomos",  "teal",
              "Contatos receptor–ligante (<0.4 nm)")
    plot_line(axes[2, 1], xvg["hbond"],    "N. H-bonds", "indianred",
              "Pontes de hidrogênio receptor–ligante")

    # Row 4 — MM-GBSA (optional)
    if has_mmgbsa:
        plot_mmgbsa_total(axes[3, 0], mmgbsa)
        plot_mmgbsa_components(axes[3, 1], mmgbsa)

    plt.tight_layout()
    out = os.path.join(D, args.output)
    plt.savefig(out, dpi=150, bbox_inches='tight')
    print(f"Salvo: {out}")
    plt.close()

    # PNGs individuais
    specs = [
        ("rmsd_bb",  "RMSD (nm)",    "navy"),
        ("rmsd_lig", "RMSD (nm)",    "darkorange"),
        ("rg",       "Rg (nm)",      "purple"),
        ("ncont",    "N. átomos",    "teal"),
        ("hbond",    "N. H-bonds",   "indianred"),
    ]
    for key, ylabel, color in specs:
        if xvg[key] is None:
            continue
        fig, ax = plt.subplots(figsize=(9, 5))
        ax.plot(xvg[key][:, 0], xvg[key][:, 1], lw=0.8, color=color)
        ax.set_xlabel("Tempo (ns)")
        ax.set_ylabel(ylabel)
        ax.set_title(key)
        ax.grid(alpha=0.3)
        plt.tight_layout()
        plt.savefig(os.path.join(D, f"{key}.png"), dpi=150, bbox_inches='tight')
        plt.close()

    if xvg["rmsf"] is not None:
        fig, ax = plt.subplots(figsize=(9, 5))
        plot_rmsf(ax, xvg["rmsf"])
        plt.tight_layout()
        plt.savefig(os.path.join(D, "rmsf.png"), dpi=150, bbox_inches='tight')
        plt.close()

    print("Pronto.")


if __name__ == "__main__":
    main()
