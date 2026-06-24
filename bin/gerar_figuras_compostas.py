#!/usr/bin/env python3
"""
gerar_figuras_compostas.py — Monta figuras compostas para o artigo.

Descobre automaticamente as PNGs individuais de cada sistema e cria:
  Figura 1 — Painéis DM série GORE4 (RMSD/contacts/H-bonds/SASA)
  Figura 2 — Distâncias tríade série GORE4
  Figura 3 — Painéis DM série SKTI
  Figura 4 — Distâncias tríade série SKTI
  Figura 5 — Evento de dissociação série BEN
  Figura 6 — Mapas de contato resíduo×resíduo
  Figura 7 — ProLIF fingerprints de persistência

Uso (na raiz do repositório no servidor):
  python3 bin/gerar_figuras_compostas.py [--outdir figures/]

Saídas: figures/Figure_1.png ... figures/Figure_7.png  (300 dpi, 180 mm wide)
"""

import argparse
import sys
from pathlib import Path
import warnings
warnings.filterwarnings("ignore")

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.image as mpimg
    from matplotlib.gridspec import GridSpec
    import numpy as np
except ImportError:
    print("ERRO: instalar matplotlib  (mamba install -n md-gromacs matplotlib)")
    sys.exit(1)

ROOT = Path(__file__).parent.parent

# ── Mapeamento de sistemas ─────────────────────────────────────────────────
# Cada entrada: (label_artigo, paths_candidatos_para_painel, paths_candidatos_triad)

GORE4_SYSTEMS = [
    ("QCL936-GORE4\n(estável)",
     ["results-MD/GORE4/QCL936/painel_completo.png",
      "ACR157-GORE4_NEW/MD/qcl936-gore4-c3/painel_completo.png"],
     ["results-MD/GORE4/QCL936/triad_distances.png",
      "ACR157-GORE4_NEW/MD/qcl936-gore4-c3/triad_distances.png"]),

    ("ACR157-GORE4\n(estável)",
     ["results-MD/GORE4/ACR157/painel_completo.png",
      "ACR157-GORE4_NEW/MD/acr157-gore4-c1/painel_completo.png"],
     ["results-MD/GORE4/ACR157/triad_distances.png",
      "ACR157-GORE4_NEW/MD/acr157-gore4-c1/triad_distances.png"]),

    ("XP273-GORE4\n(estável)",
     ["results-MD/GORE4/XP273/painel_completo.png",
      "XP273-GORE4_NEW/MD/xp273-gore4-c1/painel_completo.png"],
     ["results-MD/GORE4/XP273/triad_distances.png",
      "XP273-GORE4_NEW/MD/xp273-gore4-c1/triad_distances.png"]),

    ("XP352-GORE4\n(dissociação)",
     ["results-MD/GORE4/XP352/painel_completo.png",
      "XP352-GORE4_NEW/MD/xp352-gore4-c4/painel_completo.png"],
     ["results-MD/GORE4/XP352/triad_distances.png",
      "XP352-GORE4_NEW/MD/xp352-gore4-c4/triad_distances.png"]),
]

SKTI_SYSTEMS = [
    ("ACR157-SKTI\n(canônico 4/4)",
     ["results-MD/SKTI/ACR157/painel_completo.png",
      "ACR157-SKTI_NEW/MD/acr157-skti-c2/painel_completo.png"],
     ["results-MD/SKTI/ACR157/triad_distances.png",
      "ACR157-SKTI_NEW/MD/acr157-skti-c2/triad_distances.png"]),

    ("QCL936-SKTI\n(His+Ser+S1)",
     ["results-MD/SKTI/QCL936/painel_completo.png",
      "QCL936-SKTI_NEW/MD/qcl936-skti-c2/painel_completo.png"],
     ["results-MD/SKTI/QCL936/triad_distances.png",
      "QCL936-SKTI_NEW/MD/qcl936-skti-c2/triad_distances.png"]),

    ("XP273-SKTI\n(Tyr+Asp+S1)",
     ["results-MD/SKTI/XP273/painel_completo.png",
      "XP273-SKTI_NEW/MD/xp273-skti-c2/painel_completo.png"],
     ["results-MD/SKTI/XP273/triad_distances.png",
      "XP273-SKTI_NEW/MD/xp273-skti-c2/triad_distances.png"]),
]

BEN_SYSTEMS = [
    ("XP273-BEN\n(~80 ns, S1=Ile229)",
     ["results/XP273-BEN/XP273-BEN/painel_completo.png",
      "results-MD/BEN/XP273/painel_completo.png"]),

    ("ACR157-BEN\n(~95 ns, S1=Ile205)",
     ["results/ACR157-BEN/ACR157-BEN/painel_completo.png",
      "results-MD/BEN/ACR157/painel_completo.png"]),

    ("XP352-BEN\n(~125 ns, S1=Asp262)",
     ["results/XP352-BEN/XP352-BEN/painel_completo.png",
      "results-MD/BEN/XP352/painel_completo.png"]),

    ("QCL936-BEN\n(~150 ns, S1=Asp241)",
     ["results/QCL936-BEN/QCL936-BEN/painel_completo.png",
      "results-MD/BEN/QCL936/painel_completo.png"]),
]

CONTACT_MAP_SYSTEMS = [
    # (label, path_candidates)
    ("QCL936-GORE4",  ["results-MD/GORE4/QCL936/contact_map/contact_map.png",
                       "ACR157-GORE4_NEW/MD/qcl936-gore4-c3/contact_map/contact_map.png"]),
    ("ACR157-GORE4",  ["results-MD/GORE4/ACR157/contact_map/contact_map.png",
                       "ACR157-GORE4_NEW/MD/acr157-gore4-c1/contact_map/contact_map.png"]),
    ("XP273-GORE4",   ["results-MD/GORE4/XP273/contact_map/contact_map.png",
                       "XP273-GORE4_NEW/MD/xp273-gore4-c1/contact_map/contact_map.png"]),
    ("ACR157-SKTI",   ["results-MD/SKTI/ACR157/contact_map/contact_map.png",
                       "ACR157-SKTI_NEW/MD/acr157-skti-c2/contact_map/contact_map.png"]),
    ("QCL936-SKTI",   ["results-MD/SKTI/QCL936/contact_map/contact_map.png",
                       "QCL936-SKTI_NEW/MD/qcl936-skti-c2/contact_map/contact_map.png"]),
    ("XP273-SKTI",    ["results-MD/SKTI/XP273/contact_map/contact_map.png",
                       "XP273-SKTI_NEW/MD/xp273-skti-c2/contact_map/contact_map.png"]),
]

PROLIF_SYSTEMS = [
    ("QCL936-GORE4",  ["results-MD/GORE4/QCL936/prolif/prolif_persistence.png",
                       "ACR157-GORE4_NEW/MD/qcl936-gore4-c3/prolif/prolif_persistence.png"]),
    ("ACR157-GORE4",  ["results-MD/GORE4/ACR157/prolif/prolif_persistence.png",
                       "ACR157-GORE4_NEW/MD/acr157-gore4-c1/prolif/prolif_persistence.png"]),
    ("XP273-GORE4",   ["results-MD/GORE4/XP273/prolif/prolif_persistence.png",
                       "XP273-GORE4_NEW/MD/xp273-gore4-c1/prolif/prolif_persistence.png"]),
    ("ACR157-SKTI",   ["results-MD/SKTI/ACR157/prolif/prolif_persistence.png",
                       "ACR157-SKTI_NEW/MD/acr157-skti-c2/prolif/prolif_persistence.png"]),
    ("QCL936-SKTI",   ["results-MD/SKTI/QCL936/prolif/prolif_persistence.png",
                       "QCL936-SKTI_NEW/MD/qcl936-skti-c2/prolif/prolif_persistence.png"]),
    ("XP273-SKTI",    ["results-MD/SKTI/XP273/prolif/prolif_persistence.png",
                       "XP273-SKTI_NEW/MD/xp273-skti-c2/prolif/prolif_persistence.png"]),
]


# ── Helpers ────────────────────────────────────────────────────────────────

def resolve(candidates):
    """Retorna primeiro caminho existente ou None."""
    for c in candidates:
        p = ROOT / c
        if p.exists():
            return p
    return None


def load_img(path):
    if path is None:
        return None
    try:
        return mpimg.imread(str(path))
    except Exception:
        return None


def placeholder(ax, label, msg="arquivo não encontrado\n(rodar no servidor)"):
    ax.set_facecolor("#f0f0f0")
    ax.text(0.5, 0.5, f"{label}\n\n{msg}", ha="center", va="center",
            transform=ax.transAxes, fontsize=9, color="#555555",
            bbox=dict(boxstyle="round,pad=0.4", fc="white", ec="#aaaaaa"))
    ax.axis("off")


def panel_letter(i):
    return chr(ord("A") + i)


def add_panel_label(ax, letter, fontsize=13):
    ax.text(-0.04, 1.05, f"({letter})", transform=ax.transAxes,
            fontsize=fontsize, fontweight="bold", va="top", ha="right")


def make_composite(system_list, key_fn, n_cols, figsize, fig_title,
                   letter_labels=True):
    """Cria figura composta com as imagens de cada sistema."""
    n = len(system_list)
    n_rows = (n + n_cols - 1) // n_cols
    fig, axes = plt.subplots(n_rows, n_cols, figsize=figsize,
                              constrained_layout=True)
    if n_rows == 1 and n_cols == 1:
        axes = np.array([[axes]])
    elif n_rows == 1:
        axes = axes[np.newaxis, :]
    elif n_cols == 1:
        axes = axes[:, np.newaxis]

    for idx, entry in enumerate(system_list):
        r, c = divmod(idx, n_cols)
        ax = axes[r, c]
        label = entry[0]
        candidates = key_fn(entry)
        path = resolve(candidates)
        img = load_img(path)
        if img is not None:
            ax.imshow(img, aspect="auto")
            ax.axis("off")
            ax.set_title(label, fontsize=10, pad=4, fontweight="bold")
        else:
            placeholder(ax, label)
        if letter_labels:
            add_panel_label(ax, panel_letter(idx))

    # Ocultar eixos excedentes
    for idx in range(n, n_rows * n_cols):
        r, c = divmod(idx, n_cols)
        axes[r, c].axis("off")

    fig.suptitle(fig_title, fontsize=12, fontweight="bold", y=1.01)
    return fig


# ── Figura 1 — Painéis DM GORE4 ───────────────────────────────────────────
def fig1(outdir):
    fig = make_composite(
        GORE4_SYSTEMS,
        key_fn=lambda e: e[1],
        n_cols=2,
        figsize=(14, 10),
        fig_title="Figura 1 — Parâmetros dinâmicos dos complexos série GORE4 (100 ns)"
    )
    out = outdir / "Figure_1.png"
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  [OK] Figura 1: {out}")


# ── Figura 2 — Tríade GORE4 ────────────────────────────────────────────────
def fig2(outdir):
    fig = make_composite(
        GORE4_SYSTEMS,
        key_fn=lambda e: e[2],
        n_cols=2,
        figsize=(14, 10),
        fig_title="Figura 2 — Distâncias mínimas ao sítio catalítico — série GORE4"
    )
    out = outdir / "Figure_2.png"
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  [OK] Figura 2: {out}")


# ── Figura 3 — Painéis DM SKTI ────────────────────────────────────────────
def fig3(outdir):
    fig = make_composite(
        SKTI_SYSTEMS,
        key_fn=lambda e: e[1],
        n_cols=3,
        figsize=(18, 6),
        fig_title="Figura 3 — Parâmetros dinâmicos dos complexos série SKTI (100 ns)"
    )
    out = outdir / "Figure_3.png"
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  [OK] Figura 3: {out}")


# ── Figura 4 — Tríade SKTI ────────────────────────────────────────────────
def fig4(outdir):
    fig = make_composite(
        SKTI_SYSTEMS,
        key_fn=lambda e: e[2],
        n_cols=3,
        figsize=(18, 6),
        fig_title="Figura 4 — Distâncias mínimas ao sítio catalítico — série SKTI"
    )
    out = outdir / "Figure_4.png"
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  [OK] Figura 4: {out}")


# ── Figura 5 — Dissociação BEN ────────────────────────────────────────────
def fig5(outdir):
    fig = make_composite(
        BEN_SYSTEMS,
        key_fn=lambda e: e[1],
        n_cols=2,
        figsize=(14, 10),
        fig_title="Figura 5 — Dissociação da benzamidina (BEN) em quatro isoformas (200 ns)"
    )
    out = outdir / "Figure_5.png"
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  [OK] Figura 5: {out}")


# ── Figura 6 — Mapas de contato ───────────────────────────────────────────
def fig6(outdir):
    fig = make_composite(
        CONTACT_MAP_SYSTEMS,
        key_fn=lambda e: e[1],
        n_cols=3,
        figsize=(18, 12),
        fig_title="Figura 6 — Mapas de frequência de contato resíduo×resíduo (dist < 0,4 nm)"
    )
    out = outdir / "Figure_6.png"
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  [OK] Figura 6: {out}")


# ── Figura 7 — ProLIF persistence ─────────────────────────────────────────
def fig7(outdir):
    fig = make_composite(
        PROLIF_SYSTEMS,
        key_fn=lambda e: e[1],
        n_cols=3,
        figsize=(18, 12),
        fig_title="Figura 7 — Fingerprints ProLIF: persistência temporal de interações (%)"
    )
    out = outdir / "Figure_7.png"
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  [OK] Figura 7: {out}")


# ── Main ──────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--outdir", default="figures",
                        help="Diretório de saída das figuras (default: figures/)")
    parser.add_argument("--only", type=int, nargs="+",
                        help="Gerar apenas as figuras indicadas, ex: --only 1 2 6")
    args = parser.parse_args()

    outdir = ROOT / args.outdir
    outdir.mkdir(parents=True, exist_ok=True)
    print(f"Gerando figuras em: {outdir}")
    print()

    figs = {1: fig1, 2: fig2, 3: fig3, 4: fig4, 5: fig5, 6: fig6, 7: fig7}
    to_run = args.only if args.only else list(figs.keys())

    for n in to_run:
        if n in figs:
            figs[n](outdir)
        else:
            print(f"  [SKIP] Figura {n}: não definida")

    print()
    print("=== Resumo ===")
    for p in sorted(outdir.glob("Figure_*.png")):
        size_kb = p.stat().st_size // 1024
        print(f"  {p.name:20s}  {size_kb:5d} KB")
    print()
    print("[OK] gerar_figuras_compostas.py concluído")
    print("Para embedder no Word: python3 bin/embed_figuras_docx.py")


if __name__ == "__main__":
    main()
