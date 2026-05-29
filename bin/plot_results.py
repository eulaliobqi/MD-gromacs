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
        "rmsd_bb":      load_xvg(os.path.join(D, "rmsd_backbone.xvg")),
        "rmsd_lig":     load_xvg(os.path.join(D, "rmsd_ligante.xvg")),
        "rmsf":         load_xvg(os.path.join(D, "rmsf_residuos.xvg")),
        "rg":           load_xvg(os.path.join(D, "gyrate.xvg")),
        "ncont":        load_xvg(os.path.join(D, "numcont.xvg")),
        "hbond":        load_xvg(os.path.join(D, "hbond.xvg")),
        "sasa_prot":    load_xvg(os.path.join(D, "sasa_protein.xvg")),
        "sasa_lig":     load_xvg(os.path.join(D, "sasa_ligante.xvg")),
        "dist_ser":     load_xvg(os.path.join(D, "dist_ser.xvg")),
        "dist_his":     load_xvg(os.path.join(D, "dist_his.xvg")),
        "dist_asp":     load_xvg(os.path.join(D, "dist_asp.xvg")),
    }

    mmgbsa = load_mmgbsa_csv(args.mmgbsa_csv) if args.mmgbsa_csv else None
    has_mmgbsa = mmgbsa is not None

    has_sasa  = xvg["sasa_prot"] is not None or xvg["sasa_lig"] is not None
    has_triad = any(xvg[k] is not None for k in ("dist_ser", "dist_his", "dist_asp"))

    nrows = 3 + has_sasa + has_triad + (2 if has_mmgbsa else 0)
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

    next_row = 3

    # Row SASA
    if has_sasa:
        ax_sasa = axes[next_row, 0]
        if xvg["sasa_prot"] is not None:
            ax_sasa.plot(xvg["sasa_prot"][:, 0], xvg["sasa_prot"][:, 1],
                         lw=0.8, color="steelblue", label="Receptor")
        if xvg["sasa_lig"] is not None:
            ax2 = ax_sasa.twinx()
            ax2.plot(xvg["sasa_lig"][:, 0], xvg["sasa_lig"][:, 1],
                     lw=0.8, color="tomato", label="Ligante")
            ax2.set_ylabel("SASA ligante (nm²)", color="tomato")
            ax2.tick_params(axis='y', labelcolor='tomato')
            ax2.legend(loc='upper right', fontsize=8)
        ax_sasa.set_xlabel("Tempo (ns)")
        ax_sasa.set_ylabel("SASA receptor (nm²)", color="steelblue")
        ax_sasa.tick_params(axis='y', labelcolor='steelblue')
        ax_sasa.set_title("Área Acessível ao Solvente (SASA)")
        ax_sasa.legend(loc='upper left', fontsize=8)
        ax_sasa.grid(alpha=0.3)

        # SASA ligante individual no col 1
        plot_line(axes[next_row, 1], xvg["sasa_lig"], "SASA (nm²)", "tomato",
                  "SASA do Ligante (enterramento)")
        next_row += 1

    # Row Tríade
    if has_triad:
        ax_tri = axes[next_row, 0]
        cores_triad = {"dist_ser": ("Ser195", "forestgreen"),
                       "dist_his": ("His57",  "royalblue"),
                       "dist_asp": ("Asp102", "crimson")}
        for key, (label, cor) in cores_triad.items():
            if xvg[key] is not None:
                ax_tri.plot(xvg[key][:, 0], xvg[key][:, 1],
                            lw=0.8, color=cor, label=label, alpha=0.85)
        ax_tri.axhline(0.5, ls='--', color='gray', alpha=0.5, label='0.5 nm')
        ax_tri.set_xlabel("Tempo (ns)")
        ax_tri.set_ylabel("Distância mínima (nm)")
        ax_tri.set_title("Distâncias à Tríade Catalítica")
        ax_tri.legend(fontsize=8)
        ax_tri.grid(alpha=0.3)

        # média por resíduo da tríade como barras no col 1
        ax_bar = axes[next_row, 1]
        labels_bar, means_bar, stds_bar, colors_bar = [], [], [], []
        for key, (label, cor) in cores_triad.items():
            if xvg[key] is not None:
                labels_bar.append(label)
                means_bar.append(xvg[key][:, 1].mean())
                stds_bar.append(xvg[key][:, 1].std())
                colors_bar.append(cor)
        if labels_bar:
            bars = ax_bar.bar(labels_bar, means_bar, yerr=stds_bar,
                              color=colors_bar, alpha=0.8, capsize=5,
                              edgecolor='black', linewidth=0.5)
            ax_bar.axhline(0.5, ls='--', color='gray', alpha=0.5, label='0.5 nm')
            for bar, m in zip(bars, means_bar):
                ax_bar.text(bar.get_x() + bar.get_width() / 2,
                            bar.get_height() + 0.01,
                            f"{m:.3f} nm", ha='center', fontsize=9)
            ax_bar.set_ylabel("Distância média (nm)")
            ax_bar.set_title("Ocupação Média da Tríade Catalítica\n(média ± DP)")
            ax_bar.legend(fontsize=8)
            ax_bar.grid(alpha=0.3, axis='y')
        next_row += 1

    # Row MM-GBSA (optional)
    if has_mmgbsa:
        plot_mmgbsa_total(axes[next_row, 0], mmgbsa)
        plot_mmgbsa_components(axes[next_row, 1], mmgbsa)

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

    # SASA individual
    if xvg["sasa_lig"] is not None:
        fig, ax = plt.subplots(figsize=(9, 5))
        plot_line(ax, xvg["sasa_lig"], "SASA (nm²)", "tomato",
                  "SASA do Ligante no Complexo")
        plt.tight_layout()
        plt.savefig(os.path.join(D, "sasa_ligante.png"), dpi=150, bbox_inches='tight')
        plt.close()

    if xvg["sasa_prot"] is not None:
        fig, ax = plt.subplots(figsize=(9, 5))
        plot_line(ax, xvg["sasa_prot"], "SASA (nm²)", "steelblue",
                  "SASA do Receptor")
        plt.tight_layout()
        plt.savefig(os.path.join(D, "sasa_protein.png"), dpi=150, bbox_inches='tight')
        plt.close()

    # Tríade individual
    triad_keys = [("dist_ser", "Ser195", "forestgreen"),
                  ("dist_his", "His57",  "royalblue"),
                  ("dist_asp", "Asp102", "crimson")]
    any_triad = any(xvg[k] is not None for k, *_ in triad_keys)
    if any_triad:
        fig, ax = plt.subplots(figsize=(9, 5))
        for key, label, cor in triad_keys:
            if xvg[key] is not None:
                ax.plot(xvg[key][:, 0], xvg[key][:, 1],
                        lw=0.8, color=cor, label=label)
        ax.axhline(0.5, ls='--', color='gray', alpha=0.5, label='0.5 nm')
        ax.set_xlabel("Tempo (ns)")
        ax.set_ylabel("Distância mínima (nm)")
        ax.set_title("Distâncias à Tríade Catalítica")
        ax.legend()
        ax.grid(alpha=0.3)
        plt.tight_layout()
        plt.savefig(os.path.join(D, "triad_distances.png"), dpi=150, bbox_inches='tight')
        plt.close()

    print("Pronto.")


if __name__ == "__main__":
    main()
