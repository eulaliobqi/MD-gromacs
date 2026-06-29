#!/usr/bin/env python3
"""
plot_pharmacophore_profile.py — Perfil farmacofórico por complexo a partir dos dados ProLIF.

Para cada complexo lê prolif_summary.csv (gerado por run_prolif.py) e gera:
  - {sample_id}_pharmacophore.png  — perfil individual (barras + donut)

Além disso, cria uma figura comparativa:
  - pharmacophore_comparison.png   — heatmap de persistência por tipo de interação

Uso:
  python3 bin/plot_pharmacophore_profile.py [--outdir figures/] [--min-pct 5]

Dependência: pandas, matplotlib (já no md-gromacs env)
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib as mpl

ROOT = Path(__file__).parent.parent

mpl.rcParams.update({
    "font.size":       11,
    "axes.titlesize":  13,
    "axes.labelsize":  11,
    "legend.fontsize": 9,
})

# ── Farmacóforos ───────────────────────────────────────────────────────────

INTERACTION_TYPES = [
    "HBDonor", "HBAcceptor",
    "Hydrophobic",
    "Ionic", "Cationic", "Anionic",
    "PiStacking", "PiCation", "CationPi",
    "VdWContact",
]

PHARMA_CLASS = {
    "HBDonor":    "H-Bond",
    "HBAcceptor": "H-Bond",
    "Hydrophobic": "Hydrophobic",
    "Ionic":       "Electrostatic",
    "Cationic":    "Electrostatic",
    "Anionic":     "Electrostatic",
    "PiStacking":  "Aromatic",
    "PiCation":    "Aromatic",
    "CationPi":    "Aromatic",
    "VdWContact":  "vdW Contact",
}

CLASS_COLOR = {
    "H-Bond":        "#1f77b4",
    "Hydrophobic":   "#ff7f0e",
    "Electrostatic": "#d62728",
    "Aromatic":      "#9467bd",
    "vdW Contact":   "#7f7f7f",
    "Other":         "#2ca02c",
}

# Ordem das classes no donut
CLASS_ORDER = ["H-Bond", "Hydrophobic", "Electrostatic", "Aromatic", "vdW Contact", "Other"]


# ── Utilitários ────────────────────────────────────────────────────────────

def extract_itype(col: str) -> str:
    """Extrai o tipo de interação de um nome de coluna ProLIF."""
    for itype in INTERACTION_TYPES:
        if itype in col:
            return itype
    return "Other"


def extract_resname(col: str) -> str:
    """Retorna rótulo legível sem o tipo de interação."""
    col = str(col)
    for itype in INTERACTION_TYPES:
        col = col.replace(itype, "")
    return col.strip(".-_").replace("--", "-").strip("-")


def load_summary(csv_path: Path) -> pd.DataFrame | None:
    try:
        df = pd.read_csv(csv_path)
        if "interaction" not in df.columns or "persistence_pct" not in df.columns:
            return None
        df["itype"] = df["interaction"].apply(extract_itype)
        df["class"] = df["itype"].apply(lambda x: PHARMA_CLASS.get(x, "Other"))
        df["residue"] = df["interaction"].apply(extract_resname)
        df["color"]  = df["class"].apply(lambda x: CLASS_COLOR.get(x, CLASS_COLOR["Other"]))
        return df.sort_values("persistence_pct", ascending=False).reset_index(drop=True)
    except Exception as e:
        print(f"  [WARN] erro lendo {csv_path}: {e}")
        return None


def find_summaries() -> dict[str, Path]:
    """Descobre todos os prolif_summary.csv no repositório."""
    results = {}
    for csv in sorted(ROOT.glob("**/prolif/prolif_summary.csv")):
        if any(x in csv.parts for x in {"work", "__pycache__", "figures"}):
            continue
        # SID = nome do diretório de nível superior acima de prolif/
        # Estrutura típica: results/{ID}/{ID}/prolif/ ou {ID}_NEW/MD/{id}/prolif/
        parts = csv.parts
        prolif_idx = parts.index("prolif")
        # Tenta achar o ID do sistema subindo até 3 níveis
        for depth in (1, 2, 3):
            if prolif_idx - depth >= 0:
                candidate = parts[prolif_idx - depth]
                # Prioriza nomes que pareçam IDs de sistema
                if any(x in candidate.upper() for x in ("GORE4", "SKTI", "BEN", "ACR", "QCL", "XP")):
                    results[candidate] = csv
                    break
        else:
            results[parts[prolif_idx - 1]] = csv
    return results


# ── Figura individual ─────────────────────────────────────────────────────

def plot_single(sample_id: str, df: pd.DataFrame, outdir: Path, min_pct: float = 5.0):
    df_filt = df[df["persistence_pct"] >= min_pct].copy()
    top = df_filt.head(20)

    # Agregação por classe para donut
    class_agg = (
        df_filt.groupby("class")["persistence_pct"]
        .sum()
        .reindex(CLASS_ORDER)
        .dropna()
    )

    fig, (ax_bar, ax_donut) = plt.subplots(
        1, 2,
        figsize=(16, max(6, len(top) * 0.45 + 2)),
        gridspec_kw={"width_ratios": [3, 1.2]},
        constrained_layout=True,
    )
    fig.suptitle(
        f"{sample_id} — Pharmacophore Profile (MD 100 ns)",
        fontsize=14, fontweight="bold",
    )

    # ── Barras horizontais ─────────────────────────────────────────────────
    if top.empty:
        ax_bar.text(0.5, 0.5, "No interactions ≥ min threshold",
                    ha="center", va="center", transform=ax_bar.transAxes, fontsize=12)
        ax_bar.axis("off")
    else:
        y = np.arange(len(top))
        bars = ax_bar.barh(y, top["persistence_pct"], color=top["color"],
                           edgecolor="white", linewidth=0.5, height=0.7)
        ax_bar.set_yticks(y)
        ax_bar.set_yticklabels(top["residue"], fontsize=9)
        ax_bar.invert_yaxis()
        ax_bar.set_xlabel("Persistence (%)", fontsize=11)
        ax_bar.set_title("Residue–ligand interactions", fontsize=12, pad=6)
        ax_bar.axvline(20, color="gray", lw=0.8, ls="--", alpha=0.7)
        ax_bar.axvline(50, color="red",  lw=0.8, ls="--", alpha=0.7)
        ax_bar.set_xlim(0, 105)
        ax_bar.grid(axis="x", alpha=0.3)

        # Valor na barra
        for bar, val in zip(bars, top["persistence_pct"]):
            ax_bar.text(val + 0.8, bar.get_y() + bar.get_height() / 2,
                        f"{val:.0f}%", va="center", fontsize=8)

        # Legenda de classes
        handles = [
            mpatches.Patch(color=CLASS_COLOR[cls], label=cls)
            for cls in CLASS_ORDER if cls in top["class"].values
        ]
        ax_bar.legend(handles=handles, title="Pharmacophore class",
                      fontsize=9, title_fontsize=9,
                      loc="lower right", framealpha=0.9)

    # ── Donut ─────────────────────────────────────────────────────────────
    if class_agg.empty:
        ax_donut.axis("off")
    else:
        colors = [CLASS_COLOR.get(c, CLASS_COLOR["Other"]) for c in class_agg.index]
        wedges, texts, autotexts = ax_donut.pie(
            class_agg.values,
            labels=class_agg.index,
            colors=colors,
            autopct="%1.0f%%",
            startangle=90,
            pctdistance=0.75,
            wedgeprops={"width": 0.55, "edgecolor": "white", "linewidth": 1.5},
        )
        for t in texts:
            t.set_fontsize(9)
        for at in autotexts:
            at.set_fontsize(8)
            at.set_color("white")
            at.set_fontweight("bold")
        ax_donut.set_title("By class (Σ persistence)", fontsize=11, pad=8)

    p = outdir / f"{sample_id}_pharmacophore.png"
    fig.savefig(p, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  [OK] {p.name}  ({p.stat().st_size // 1024} KB)")


# ── Figura comparativa (heatmap) ──────────────────────────────────────────

def plot_comparison(all_data: dict[str, pd.DataFrame], outdir: Path, min_pct: float = 5.0):
    if len(all_data) < 2:
        return

    # Pivot: linhas = sistemas, colunas = classes farmacofóricas
    rows = []
    for sid, df in all_data.items():
        class_sum = (
            df[df["persistence_pct"] >= min_pct]
            .groupby("class")["persistence_pct"]
            .sum()
            .reindex(CLASS_ORDER)
            .fillna(0)
        )
        class_sum.name = sid
        rows.append(class_sum)

    hm = pd.DataFrame(rows).fillna(0)
    # Filtra colunas com algum valor
    hm = hm.loc[:, (hm > 0).any()]

    fig, ax = plt.subplots(
        figsize=(max(8, len(hm.columns) * 1.5), max(5, len(hm) * 0.55 + 2)),
        constrained_layout=True,
    )
    fig.suptitle(
        "Pharmacophore Profile — All complexes (Σ persistence by class, %)",
        fontsize=13, fontweight="bold",
    )

    im = ax.imshow(hm.values, aspect="auto", cmap="YlOrRd", vmin=0)
    plt.colorbar(im, ax=ax, label="Σ persistence (%)", shrink=0.6)

    ax.set_xticks(range(len(hm.columns)))
    ax.set_xticklabels(hm.columns, rotation=30, ha="right", fontsize=11)
    ax.set_yticks(range(len(hm.index)))
    ax.set_yticklabels(hm.index, fontsize=10)

    for i in range(len(hm.index)):
        for j in range(len(hm.columns)):
            val = hm.values[i, j]
            if val > 0:
                ax.text(j, i, f"{val:.0f}", ha="center", va="center",
                        fontsize=9, color="black" if val < hm.values.max() * 0.6 else "white",
                        fontweight="bold")

    p = outdir / "pharmacophore_comparison.png"
    fig.savefig(p, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  [OK] {p.name}  ({p.stat().st_size // 1024} KB)")


# ── Main ──────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--outdir",  default="figures",
                    help="Diretório de saída (default: figures/)")
    ap.add_argument("--min-pct", type=float, default=5.0,
                    help="Threshold mínimo de persistência %% para exibir (default: 5)")
    args = ap.parse_args()

    outdir = ROOT / args.outdir
    outdir.mkdir(parents=True, exist_ok=True)
    print(f"Output: {outdir}\n")

    summaries = find_summaries()
    if not summaries:
        print("[WARN] Nenhum prolif_summary.csv encontrado.")
        print("  Execute run_prolif.py para cada complexo antes deste script.")
        return

    print(f"Complexos encontrados: {len(summaries)}")
    for sid, path in sorted(summaries.items()):
        print(f"  {sid:30s}  {path.relative_to(ROOT)}")
    print()

    all_data = {}
    for sid, csv_path in sorted(summaries.items()):
        df = load_summary(csv_path)
        if df is None or df.empty:
            print(f"  [SKIP] {sid} — dados vazios ou inválidos")
            continue
        print(f"[{sid}]  {len(df)} interações")
        plot_single(sid, df, outdir, min_pct=args.min_pct)
        all_data[sid] = df

    if len(all_data) >= 2:
        print("\nGerando figura comparativa...")
        plot_comparison(all_data, outdir, min_pct=args.min_pct)

    print("\n[OK] done")


if __name__ == "__main__":
    main()
