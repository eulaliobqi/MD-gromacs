#!/usr/bin/env python3
"""
Analisa resultados de docking Vina de BEN contra 4 tripsinas Spodoptera.
- Extrai scores de todos os modos de cada receptor
- Calcula distâncias da melhor pose aos resíduos catalíticos
- Gera: ben_scores.csv, ben_distances.csv, ben_docking_panel.png

Uso:
    python3 bin/analyze_ben_docking.py --ben-dir spodoptera-ben
"""

import argparse
import os
import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd


RECEPTORS = ["ACR157", "QCL936", "XP273", "XP352"]

# Resíduos catalíticos por receptor (número PDB, label)
TRIAD_MAP = {
    "ACR157": [(69, "His69"), (114, "Asp114"), (211, "Ser211"), (205, "Ile205/S1")],
    "QCL936": [(92, "His92"),  (142, "Asp142"), (247, "Ser247"), (241, "Asp241/S1")],
    "XP273":  [(83, "Tyr83"),  (132, "Asp132"), (234, "Ser234"), (229, "Ile229/S1")],
    "XP352":  [(112,"Arg112"),(166, "Asp166"), (268, "Ser268"), (262, "Asp262/S1")],
}


def parse_vina_scores(pdbqt_path):
    """Extrai (modo, afinidade, rmsd_lb, rmsd_ub) do arquivo _out.pdbqt.

    Suporta dois formatos:
    - Vina >= 1.2: linha 'REMARK VINA RESULT: aff rmsd_lb rmsd_ub'
    - Vina customizado (design-inibidores): sem REMARK, usa .log associado
    """
    scores = []

    # Formato 1: REMARK VINA RESULT no PDBQT (Vina 1.2+)
    with open(pdbqt_path) as f:
        for line in f:
            m = re.match(r"REMARK VINA RESULT:\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)", line)
            if m:
                mode = len(scores) + 1
                scores.append((mode, float(m.group(1)), float(m.group(2)), float(m.group(3))))
    if scores:
        return scores

    # Formato 2: log capturado via stdout (Vina customizado, sem REMARK)
    log_path = str(pdbqt_path).replace("_out.pdbqt", "_vina.log")
    if os.path.exists(log_path):
        with open(log_path) as f:
            for line in f:
                m = re.match(r"^\s+(\d+)\s+([-\d.]+)\s+[-\d.]+\s+[-\d.]+", line)
                if m:
                    mode = int(m.group(1))
                    scores.append((mode, float(m.group(2)), 0.0, 0.0))
    return scores


def parse_vina_models(pdbqt_path):
    """Lê todos os modelos (poses) do arquivo _out.pdbqt."""
    models = []
    current = []
    with open(pdbqt_path) as f:
        for line in f:
            if line.startswith("MODEL"):
                current = []
            elif line.startswith("ENDMDL"):
                if current:
                    models.append(current)
            else:
                current.append(line)
    return models


def extract_coords_from_model(lines):
    """Extrai coordenadas (x,y,z) de todas as linhas ATOM/HETATM de um modelo."""
    coords = []
    for line in lines:
        if line.startswith("ATOM") or line.startswith("HETATM"):
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                coords.append(np.array([x, y, z]))
            except ValueError:
                continue
    return coords


def parse_receptor_residue_coords(pdb_path, residue_nums):
    """Extrai todos os átomos pesados de cada resíduo catalítico no receptor."""
    res_atoms = {r: [] for r in residue_nums}
    with open(pdb_path) as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            atom_name = line[12:16].strip()
            if atom_name.startswith("H"):
                continue
            try:
                resnum = int(line[22:26].strip())
            except ValueError:
                continue
            if resnum in res_atoms:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                res_atoms[resnum].append(np.array([x, y, z]))
    return res_atoms


def min_distance(ligand_coords, receptor_atom_list):
    """Distância mínima entre qualquer átomo do ligante e qualquer átomo do resíduo."""
    if not ligand_coords or not receptor_atom_list:
        return float("nan")
    lig = np.array(ligand_coords)
    rec = np.array(receptor_atom_list)
    dists = np.linalg.norm(lig[:, None, :] - rec[None, :, :], axis=2)
    return float(np.min(dists))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--ben-dir", default="spodoptera-ben", help="Diretório base")
    args = parser.parse_args()

    ben_dir = Path(args.ben_dir)
    res_dir = ben_dir / "results"

    # ── 1. Extrair scores ──────────────────────────────────────────────────────
    all_scores = []
    for rec in RECEPTORS:
        out_pdbqt = res_dir / f"{rec}_ben_out.pdbqt"
        if not out_pdbqt.exists():
            print(f"AVISO: {out_pdbqt} não encontrado — pulando {rec}", file=sys.stderr)
            continue
        scores = parse_vina_scores(out_pdbqt)
        for mode, aff, rlb, rub in scores:
            all_scores.append({"receptor": rec, "mode": mode,
                                "affinity_kcal_mol": aff,
                                "rmsd_lb": rlb, "rmsd_ub": rub})

    df_scores = pd.DataFrame(all_scores)
    scores_csv = res_dir / "ben_scores.csv"
    df_scores.to_csv(scores_csv, index=False)
    print(f"Scores salvos: {scores_csv}")

    # ── 2. Calcular distâncias (melhor pose modo 1) ────────────────────────────
    dist_rows = []
    best_scores = {}

    for rec in RECEPTORS:
        out_pdbqt = res_dir / f"{rec}_ben_out.pdbqt"
        receptor_pdb = ben_dir / f"{rec}-final.pdb"
        if not out_pdbqt.exists() or not receptor_pdb.exists():
            continue

        models = parse_vina_models(out_pdbqt)
        if not models:
            print(f"AVISO: sem modelos em {out_pdbqt}", file=sys.stderr)
            continue

        best_model_lines = models[0]
        lig_coords = extract_coords_from_model(best_model_lines)

        triad = TRIAD_MAP[rec]
        res_nums = [r for r, _ in triad]
        rec_atoms = parse_receptor_residue_coords(receptor_pdb, res_nums)

        row = {"receptor": rec}
        scores_this = df_scores[df_scores["receptor"] == rec]
        if not scores_this.empty:
            best_aff = scores_this[scores_this["mode"] == 1]["affinity_kcal_mol"].values
            row["best_affinity"] = float(best_aff[0]) if len(best_aff) else float("nan")
            best_scores[rec] = row["best_affinity"]

        for resnum, label in triad:
            dist = min_distance(lig_coords, rec_atoms.get(resnum, []))
            row[label] = round(dist / 10.0, 3)  # Å → nm

        dist_rows.append(row)

    df_dist = pd.DataFrame(dist_rows)
    dist_csv = res_dir / "ben_distances.csv"
    df_dist.to_csv(dist_csv, index=False)
    print(f"Distâncias salvas: {dist_csv}")

    # ── 3. Painel ──────────────────────────────────────────────────────────────
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    fig.suptitle("BEN Docking — AutoDock Vina | 4 Tripsinas Spodoptera", fontsize=13, weight="bold")

    # Subplot 1: scores (modo 1 + distribuição dos 9 modos)
    ax1 = axes[0]
    recs_with_data = [r for r in RECEPTORS if r in best_scores]
    colors = ["#2196F3", "#4CAF50", "#FF9800", "#9C27B0"][:len(recs_with_data)]

    if not df_scores.empty:
        for i, rec in enumerate(recs_with_data):
            sub = df_scores[df_scores["receptor"] == rec]["affinity_kcal_mol"].values
            x_pos = np.ones(len(sub)) * i + np.random.uniform(-0.15, 0.15, len(sub))
            ax1.scatter(x_pos, sub, color=colors[i], alpha=0.5, s=30, zorder=2)
            ax1.bar(i, best_scores[rec], color=colors[i], alpha=0.8, width=0.5,
                    label=f"{rec}: {best_scores[rec]:.1f}", zorder=1)

    ax1.set_xticks(range(len(recs_with_data)))
    ax1.set_xticklabels(recs_with_data, fontsize=11)
    ax1.set_ylabel("Afinidade (kcal/mol)", fontsize=11)
    ax1.set_title("Melhor score (barra) + 9 modos (pontos)", fontsize=10)
    ax1.legend(fontsize=9, loc="lower right")
    ax1.axhline(0, color="gray", linewidth=0.5, linestyle="--")
    ax1.invert_yaxis()  # mais negativo = melhor

    # Subplot 2: heatmap distâncias BEN→tríade
    ax2 = axes[1]
    if not df_dist.empty:
        triad_cols = [c for c in df_dist.columns if c not in ("receptor", "best_affinity")]
        mat = df_dist.set_index("receptor")[triad_cols].reindex(RECEPTORS)
        mat_vals = mat.values.astype(float)

        cmap = plt.cm.RdYlGn_r
        norm = mcolors.Normalize(vmin=0.2, vmax=1.5)
        im = ax2.imshow(mat_vals, cmap=cmap, norm=norm, aspect="auto")
        plt.colorbar(im, ax=ax2, label="Distância mínima (nm)")

        ax2.set_xticks(range(len(triad_cols)))
        ax2.set_xticklabels(triad_cols, rotation=35, ha="right", fontsize=9)
        ax2.set_yticks(range(len(RECEPTORS)))
        ax2.set_yticklabels(RECEPTORS, fontsize=10)
        ax2.set_title("Distância BEN → sítio ativo (nm)\nVerde < 0.35 nm = contato", fontsize=10)

        # Anotações
        for i in range(mat_vals.shape[0]):
            for j in range(mat_vals.shape[1]):
                val = mat_vals[i, j]
                if not np.isnan(val):
                    color = "white" if val > 0.9 or val < 0.3 else "black"
                    ax2.text(j, i, f"{val:.2f}", ha="center", va="center",
                             fontsize=8, color=color, weight="bold")

        # Linha de corte 0.35 nm
        ax2.axvline(-0.5, color="none")

    plt.tight_layout()
    panel_path = res_dir / "ben_docking_panel.png"
    plt.savefig(panel_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Painel salvo: {panel_path}")

    # ── 4. Resumo terminal ─────────────────────────────────────────────────────
    print("\n=== RESUMO BEN DOCKING ===")
    print(f"{'Receptor':<10} {'Best (kcal/mol)':<18} {'Contatos <0.35nm'}")
    print("-" * 50)
    for _, row in df_dist.iterrows():
        rec = row["receptor"]
        aff = row.get("best_affinity", float("nan"))
        triad_cols = [c for c in row.index if c not in ("receptor", "best_affinity")]
        n_contacts = sum(1 for c in triad_cols if not np.isnan(row[c]) and row[c] < 0.35)
        print(f"{rec:<10} {aff:<18.2f} {n_contacts}/4")


if __name__ == "__main__":
    main()
