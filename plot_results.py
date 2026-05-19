#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Gera figuras principais a partir dos .xvg produzidos por run_analyses.sh.

Saidas:
  painel_completo.png  - 6 paineis (RMSD bb, RMSD lig, RMSF, Rg, contatos, H-bonds)
  + 1 PNG individual por metrica

Uso:
  python3 plot_results.py --analise-dir <dir> [--titulo TITLE]
"""
import argparse, os, sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams['font.size'] = 10
mpl.rcParams['axes.titlesize'] = 11


def load_xvg(path):
    if not os.path.exists(path): return None
    data = []
    for ln in open(path):
        if ln.startswith(('#', '@')): continue
        try:
            data.append([float(x) for x in ln.split()])
        except ValueError:
            pass
    return np.array(data) if data else None


def plot_line(ax, data, ylabel, color, title=None, hline=None):
    if data is None:
        ax.text(0.5, 0.5, "(arquivo nao encontrado)", ha='center', va='center',
                transform=ax.transAxes, fontsize=10, color='gray')
        ax.set_xticks([]); ax.set_yticks([])
        return
    ax.plot(data[:,0], data[:,1], lw=0.8, color=color)
    ax.set_xlabel("Tempo (ns)"); ax.set_ylabel(ylabel)
    if title:
        media = data[:,1].mean(); std = data[:,1].std()
        ax.set_title(f"{title} (media: {media:.3f} ± {std:.3f})")
    if hline is not None:
        ax.axhline(hline, ls='--', color='red', alpha=0.4)
    ax.grid(alpha=0.3)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--analise-dir", required=True)
    ap.add_argument("--titulo", default="Dinamica Molecular")
    args = ap.parse_args()

    D = args.analise_dir
    files = {
        "rmsd_bb":   os.path.join(D, "rmsd_backbone.xvg"),
        "rmsd_lig":  os.path.join(D, "rmsd_ligante.xvg"),
        "rmsf":      os.path.join(D, "rmsf_residuos.xvg"),
        "rg":        os.path.join(D, "gyrate.xvg"),
        "ncont":     os.path.join(D, "numcont.xvg"),
        "hbond":     os.path.join(D, "hbond.xvg"),
    }
    data = {k: load_xvg(p) for k, p in files.items()}

    # ---------- Painel 6-em-1 ----------
    fig, axes = plt.subplots(3, 2, figsize=(14, 12))
    fig.suptitle(args.titulo, fontsize=14, fontweight='bold')

    # 1) RMSD backbone
    if data["rmsd_bb"] is not None:
        media = data["rmsd_bb"][:,1].mean()
        plot_line(axes[0,0], data["rmsd_bb"], "RMSD (nm)", "navy",
                  "RMSD do Backbone", hline=media)

    # 2) RMSD ligante
    plot_line(axes[0,1], data["rmsd_lig"], "RMSD (nm)", "darkorange",
              "RMSD do Ligante (peptidio)")

    # 3) RMSF
    ax = axes[1,0]
    if data["rmsf"] is not None:
        d = data["rmsf"]
        ax.bar(d[:,0], d[:,1], width=1.0, color="seagreen", alpha=0.8)
        # destaca 5 ultimos residuos (geralmente o ligante)
        n = len(d)
        ax.bar(d[-5:,0], d[-5:,1], width=1.0, color="crimson",
               label="Ligante (5 ultimos)")
        ax.legend()
        ax.set_xlabel("Residuo"); ax.set_ylabel("RMSF (nm)")
        ax.set_title("Flutuacao por residuo (RMSF)")
        ax.grid(alpha=0.3, axis='y')
    else:
        ax.text(0.5,0.5,"(rmsf nao encontrado)", ha='center', transform=ax.transAxes)

    # 4) Raio de giro
    plot_line(axes[1,1], data["rg"], "Rg (nm)", "purple", "Raio de Giro")

    # 5) Contatos (numero de atomos em contato < 0.4 nm)
    plot_line(axes[2,0], data["ncont"], "N. atomos", "teal",
              "Contatos receptor-ligante (<0.4 nm)")

    # 6) Pontes de hidrogenio
    plot_line(axes[2,1], data["hbond"], "N. H-bonds", "indianred",
              "Pontes de hidrogenio receptor-ligante")

    plt.tight_layout()
    out = os.path.join(D, "painel_completo.png")
    plt.savefig(out, dpi=150, bbox_inches='tight')
    print(f"Salvo: {out}")
    plt.close()

    # PNGs individuais
    for k, arr in data.items():
        if arr is None: continue
        fig, ax = plt.subplots(figsize=(9, 5))
        if k == "rmsf":
            ax.bar(arr[:,0], arr[:,1], width=1.0, color="seagreen", alpha=0.8)
            ax.bar(arr[-5:,0], arr[-5:,1], width=1.0, color="crimson")
            ax.set_xlabel("Residuo"); ax.set_ylabel("RMSF (nm)")
        else:
            ax.plot(arr[:,0], arr[:,1], lw=0.8)
            ax.set_xlabel("Tempo (ns)")
            ax.set_ylabel({"rmsd_bb":"RMSD (nm)","rmsd_lig":"RMSD (nm)",
                           "rg":"Rg (nm)","ncont":"N. atomos","hbond":"N. H-bonds"}[k])
        ax.set_title(k)
        ax.grid(alpha=0.3)
        out = os.path.join(D, f"{k}.png")
        plt.savefig(out, dpi=150, bbox_inches='tight')
        plt.close()

    print("Pronto.")


if __name__ == "__main__":
    main()
