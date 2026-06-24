#!/usr/bin/env python3
"""
gerar_figuras_compostas.py — Monta figuras compostas para o artigo.

Auto-descobre os PNGs individuais de cada sistema nas subpastas *_NEW/MD/,
results-MD/, results/ e gera:

  Figura 1 — Painéis DM série GORE4 (RMSD/contacts/H-bonds/SASA)
  Figura 2 — Distâncias tríade série GORE4
  Figura 3 — Painéis DM série SKTI
  Figura 4 — Distâncias tríade série SKTI
  Figura 5 — Evento de dissociação série BEN
  Figura 6 — Mapas de contato resíduo×resíduo
  Figura 7 — ProLIF fingerprints de persistência

Uso (na raiz do repositório no servidor):
  python3 bin/gerar_figuras_compostas.py [--outdir figures/] [--only 1 2 6] [--scan]

Saídas: figures/Figure_N.png  (300 dpi)
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
    import numpy as np
except ImportError:
    print("ERRO: instalar matplotlib  (mamba install -n md-gromacs matplotlib)")
    sys.exit(1)

ROOT = Path(__file__).parent.parent


# ── Auto-descoberta de arquivos ────────────────────────────────────────────

_SKIP_DIRS = {"figures", "local_viz", "bin", "__pycache__", "work"}

def _discover(filename):
    """
    Retorna dict {sid_normalizado: Path} para todos os `filename` encontrados.
    Ignora Nextflow work/ dirs (hashes ilegíveis) e dirs utilitários.
    """
    result = {}
    for p in sorted(ROOT.glob(f"**/{filename}")):
        parts = p.parts
        if any(x in _SKIP_DIRS for x in parts):
            continue
        parent = p.parent
        if parent.name in ("analise", "contact_map", "prolif"):
            sid_dir = parent.parent
        else:
            sid_dir = parent
        sid = sid_dir.name.lower().replace("-", "").replace("_", "")
        if sid not in result:
            result[sid] = p
    return result


def _discover_individual_plots():
    """
    Para cada SID com analise/, descobre plots individuais de DM.
    Retorna dict {sid: {fname: Path}}.
    """
    INDIVIDUAL = ["rmsd_bb.png", "rmsd_lig.png", "ncont.png", "hbond.png",
                  "sasa_ligante.png", "rg.png"]
    per_sid = {}
    for fname in INDIVIDUAL:
        for p in sorted(ROOT.glob(f"**/analise/{fname}")):
            if any(x in _SKIP_DIRS for x in p.parts):
                continue
            sid = p.parent.parent.name.lower().replace("-", "").replace("_", "")
            per_sid.setdefault(sid, {})[fname] = p
    return per_sid


def _compose_painel_from_individuals(plots_for_sid):
    """
    Monta painel 2×2 com plots individuais de analise/ quando painel_completo.png ausente.
    Retorna fig matplotlib ou None se nenhum arquivo existe.
    """
    ORDER  = ["rmsd_bb.png",  "ncont.png", "hbond.png", "rmsd_lig.png"]
    LABELS = ["RMSD backbone", "Contatos (# átomos)", "H-bonds", "RMSD ligante"]
    imgs = [(lbl, plots_for_sid.get(fname)) for lbl, fname in zip(LABELS, ORDER)]
    if not any(p is not None and p.exists() for _, p in imgs):
        return None
    fig, axes = plt.subplots(2, 2, figsize=(9, 7), constrained_layout=True)
    for ax, (lbl, p) in zip(axes.flat, imgs):
        if p and p.exists():
            try:
                ax.imshow(mpimg.imread(str(p)), aspect="auto")
                ax.axis("off")
                ax.set_title(lbl, fontsize=9, pad=3)
            except Exception:
                placeholder(ax, lbl)
        else:
            placeholder(ax, lbl, "arquivo ausente")
    return fig


def _norm(name):
    return name.lower().replace("-", "").replace("_", "").replace("\n", "")


# ── Catálogo de sistemas (label, padrão sid) ───────────────────────────────

GORE4_ENTRIES = [
    ("QCL936-GORE4\n(estável)",    ["qcl936gore4c3", "qcl936gore4c1", "qcl936gore4"]),
    ("ACR157-GORE4\n(estável)",    ["acr157gore4c1", "acr157gore4"]),
    ("XP273-GORE4\n(estável)",     ["xp273gore4c1", "xp273gore4"]),
    ("XP352-GORE4\n(dissociação)", ["xp352gore4c4", "xp352gore4c1", "xp352gore4"]),
]

SKTI_ENTRIES = [
    ("ACR157-SKTI\n(canônico 4/4)", ["acr157sktic2", "acr157sktic1", "acr157skti"]),
    ("QCL936-SKTI\n(His+Ser+S1)",   ["qcl936sktic2", "qcl936sktic1", "qcl936skti"]),
    ("XP273-SKTI\n(Tyr+Asp+S1)",    ["xp273sktic2", "xp273sktic1", "xp273skti"]),
]

BEN_ENTRIES = [
    ("XP273-BEN\n(~80 ns, S1=Ile229)",  ["xp273ben", "xp273benxp273ben"]),
    ("ACR157-BEN\n(~95 ns, S1=Ile205)", ["acr157ben", "acr157benacr157ben"]),
    ("XP352-BEN\n(~125 ns, S1=Asp262)", ["xp352ben", "xp352benxp352ben"]),
    ("QCL936-BEN\n(~150 ns, S1=Asp241)",["qcl936ben", "qcl936benqcl936ben"]),
]

STABLE_ENTRIES = GORE4_ENTRIES[:3] + SKTI_ENTRIES  # 6 sistemas estáveis para Fig 6-7


def _find(index, patterns):
    """Busca no index por qualquer um dos patterns."""
    for pat in patterns:
        if pat in index:
            return index[pat]
    # Tentativa parcial: prefixo
    for key, path in index.items():
        for pat in patterns:
            if pat in key or key in pat:
                return path
    return None


# ── Funções de renderização ────────────────────────────────────────────────

def placeholder(ax, label, msg="arquivo não encontrado\n(rodar no servidor)"):
    ax.set_facecolor("#f0f0f0")
    ax.text(0.5, 0.5, f"{label}\n\n{msg}", ha="center", va="center",
            transform=ax.transAxes, fontsize=9, color="#555555",
            bbox=dict(boxstyle="round,pad=0.4", fc="white", ec="#aaaaaa"))
    ax.axis("off")


def add_panel_label(ax, letter, fontsize=13):
    ax.text(-0.04, 1.05, f"({letter})", transform=ax.transAxes,
            fontsize=fontsize, fontweight="bold", va="top", ha="right")


def render_panel(ax, img_path, label, letter=None):
    if img_path and img_path.exists():
        try:
            img = mpimg.imread(str(img_path))
            ax.imshow(img, aspect="auto")
            ax.axis("off")
            ax.set_title(label, fontsize=10, pad=4, fontweight="bold")
        except Exception as e:
            placeholder(ax, label, f"erro ao carregar:\n{e}")
    else:
        placeholder(ax, label)
    if letter:
        add_panel_label(ax, letter)


def make_composite(entries, index, n_cols, figsize, fig_title, idx_individual=None):
    """
    Monta figura composta. Para cada painel, tenta:
    1. painel_completo.png do index
    2. Fallback: compõe 2×2 a partir de plots individuais (rmsd_bb, ncont, hbond, rmsd_lig)
    3. Placeholder cinza se nada encontrado
    """
    n = len(entries)
    n_rows = (n + n_cols - 1) // n_cols
    # Cada célula da grade principal recebe uma sub-imagem do painel
    fig, axes = plt.subplots(n_rows, n_cols, figsize=figsize, constrained_layout=True)
    if n_rows == 1 and n_cols == 1:
        axes = np.array([[axes]])
    elif n_rows == 1:
        axes = axes[np.newaxis, :]
    elif n_cols == 1:
        axes = axes[:, np.newaxis]

    for idx, (label, patterns) in enumerate(entries):
        r, c = divmod(idx, n_cols)
        ax = axes[r, c]
        path = _find(index, patterns)

        # Fallback para plots individuais
        if (path is None or not path.exists()) and idx_individual:
            # Encontra o SID do primeiro pattern que tenha plots individuais
            for pat in patterns:
                if pat in idx_individual:
                    sub_fig = _compose_painel_from_individuals(idx_individual[pat])
                    if sub_fig is not None:
                        # Salva fig temporária em memória e carrega como imagem
                        import io
                        buf = io.BytesIO()
                        sub_fig.savefig(buf, format="png", dpi=100, bbox_inches="tight")
                        plt.close(sub_fig)
                        buf.seek(0)
                        img = mpimg.imread(buf)
                        ax.imshow(img, aspect="auto")
                        ax.axis("off")
                        ax.set_title(label, fontsize=10, pad=4, fontweight="bold")
                        ax.text(-0.04, 1.05, f"({chr(ord('A') + idx)})",
                                transform=ax.transAxes, fontsize=13, fontweight="bold",
                                va="top", ha="right")
                        break
            else:
                render_panel(ax, path, label, chr(ord("A") + idx))
        else:
            render_panel(ax, path, label, chr(ord("A") + idx))

    for idx in range(n, n_rows * n_cols):
        r, c = divmod(idx, n_cols)
        axes[r, c].axis("off")

    fig.suptitle(fig_title, fontsize=12, fontweight="bold", y=1.01)
    return fig


# ── Funções por figura ─────────────────────────────────────────────────────

def fig1(outdir, idx_painel, idx_ind):
    fig = make_composite(GORE4_ENTRIES, idx_painel, n_cols=2, figsize=(14, 10),
                         fig_title="Figura 1 — Parâmetros dinâmicos dos complexos série GORE4 (100 ns)",
                         idx_individual=idx_ind)
    p = outdir / "Figure_1.png"
    fig.savefig(p, dpi=300, bbox_inches="tight"); plt.close(fig)
    print(f"  [OK] {p.name}  ({p.stat().st_size//1024} KB)")


def fig2(outdir, idx_triad):
    fig = make_composite(GORE4_ENTRIES, idx_triad, n_cols=2, figsize=(14, 10),
                         fig_title="Figura 2 — Distâncias mínimas ao sítio catalítico — série GORE4")
    p = outdir / "Figure_2.png"
    fig.savefig(p, dpi=300, bbox_inches="tight"); plt.close(fig)
    print(f"  [OK] {p.name}  ({p.stat().st_size//1024} KB)")


def fig3(outdir, idx_painel, idx_ind):
    fig = make_composite(SKTI_ENTRIES, idx_painel, n_cols=3, figsize=(18, 6),
                         fig_title="Figura 3 — Parâmetros dinâmicos dos complexos série SKTI (100 ns)",
                         idx_individual=idx_ind)
    p = outdir / "Figure_3.png"
    fig.savefig(p, dpi=300, bbox_inches="tight"); plt.close(fig)
    print(f"  [OK] {p.name}  ({p.stat().st_size//1024} KB)")


def fig4(outdir, idx_triad):
    fig = make_composite(SKTI_ENTRIES, idx_triad, n_cols=3, figsize=(18, 6),
                         fig_title="Figura 4 — Distâncias mínimas ao sítio catalítico — série SKTI")
    p = outdir / "Figure_4.png"
    fig.savefig(p, dpi=300, bbox_inches="tight"); plt.close(fig)
    print(f"  [OK] {p.name}  ({p.stat().st_size//1024} KB)")


def fig5(outdir, idx_painel, idx_ind):
    fig = make_composite(BEN_ENTRIES, idx_painel, n_cols=2, figsize=(14, 10),
                         fig_title="Figura 5 — Dissociação da benzamidina em quatro isoformas (200 ns)",
                         idx_individual=idx_ind)
    p = outdir / "Figure_5.png"
    fig.savefig(p, dpi=300, bbox_inches="tight"); plt.close(fig)
    print(f"  [OK] {p.name}  ({p.stat().st_size//1024} KB)")


def fig6(outdir, idx_cmap):
    fig = make_composite(STABLE_ENTRIES, idx_cmap, n_cols=3, figsize=(18, 12),
                         fig_title="Figura 6 — Mapas de frequência de contato resíduo×resíduo (dist < 0,4 nm)")
    p = outdir / "Figure_6.png"
    fig.savefig(p, dpi=300, bbox_inches="tight"); plt.close(fig)
    print(f"  [OK] {p.name}  ({p.stat().st_size//1024} KB)")


def fig7(outdir, idx_prolif):
    fig = make_composite(STABLE_ENTRIES, idx_prolif, n_cols=3, figsize=(18, 12),
                         fig_title="Figura 7 — Fingerprints ProLIF: persistência temporal de interações (%)")
    p = outdir / "Figure_7.png"
    fig.savefig(p, dpi=300, bbox_inches="tight"); plt.close(fig)
    print(f"  [OK] {p.name}  ({p.stat().st_size//1024} KB)")


# ── Scan para debug ────────────────────────────────────────────────────────

def cmd_scan():
    print("=== Descoberta automática de arquivos ===")
    for fname, desc in [
        ("painel_completo.png", "Painéis DM"),
        ("triad_distances.png", "Distâncias tríade"),
        ("contact_map.png",     "Contact maps"),
        ("prolif_persistence.png", "ProLIF"),
    ]:
        idx = _discover(fname)
        print(f"\n{desc} ({fname}):  {len(idx)} encontrados")
        for sid, path in sorted(idx.items()):
            print(f"  {sid:40s}  {path.relative_to(ROOT)}")

    idx_ind = _discover_individual_plots()
    print(f"\nPlots individuais (analise/rmsd_bb.png etc.):  {len(idx_ind)} SIDs")
    for sid, plots in sorted(idx_ind.items()):
        print(f"  {sid:40s}  {sorted(plots.keys())}")


# ── Main ──────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--outdir", default="figures")
    parser.add_argument("--only", type=int, nargs="+",
                        help="Gerar apenas figuras indicadas, ex: --only 1 2")
    parser.add_argument("--scan", action="store_true",
                        help="Listar todos os PNGs encontrados (diagnóstico)")
    args = parser.parse_args()

    if args.scan:
        cmd_scan()
        return

    outdir = ROOT / args.outdir
    outdir.mkdir(parents=True, exist_ok=True)

    # Constrói índices uma vez
    print("Descobrindo arquivos no repositório...")
    idx_painel = _discover("painel_completo.png")
    idx_triad  = _discover("triad_distances.png")
    idx_cmap   = _discover("contact_map.png")
    idx_prolif = _discover("prolif_persistence.png")
    idx_ind    = _discover_individual_plots()

    counts = {
        "painel_completo": len(idx_painel),
        "triad_distances": len(idx_triad),
        "contact_map": len(idx_cmap),
        "prolif": len(idx_prolif),
        "individual_plots (sids)": len(idx_ind),
    }
    print(f"  Encontrados: {counts}")
    print(f"  Saídas em  : {outdir}")
    print()

    figs = {
        1: lambda: fig1(outdir, idx_painel, idx_ind),
        2: lambda: fig2(outdir, idx_triad),
        3: lambda: fig3(outdir, idx_painel, idx_ind),
        4: lambda: fig4(outdir, idx_triad),
        5: lambda: fig5(outdir, idx_painel, idx_ind),
        6: lambda: fig6(outdir, idx_cmap),
        7: lambda: fig7(outdir, idx_prolif),
    }
    to_run = args.only if args.only else list(figs.keys())
    for n in to_run:
        if n in figs:
            figs[n]()

    print()
    print("[OK] gerar_figuras_compostas.py concluído")
    print("Para embedder no Word: python3 bin/embed_figuras_docx.py")


if __name__ == "__main__":
    main()
