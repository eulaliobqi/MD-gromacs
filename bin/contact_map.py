#!/usr/bin/env python3
"""
Mapa de contato residuo(receptor) × residuo(ligante) ao longo da trajetória MD.

Calcula a frequência (fração de frames) em que cada par de resíduos mantém
ao menos um par de átomos a distância < cutoff. Saídas:
  contact_map.csv  — matriz de frequência (0-1)
  contact_map.png  — heatmap (receptor_residues × ligand_residues)

Uso:
    python3 bin/contact_map.py \
        --tpr  results/ACR157-GORE4/ACR157-GORE4/prod/md.tpr \
        --xtc  results/ACR157-GORE4/ACR157-GORE4/prod/md_fit.xtc \
        --ndx  results/ACR157-GORE4/ACR157-GORE4/analise/lig.ndx \
        --outdir results/ACR157-GORE4/ACR157-GORE4/contact_map \
        --sample-id ACR157-GORE4 \
        [--cutoff 0.4] [--start-ns 0] [--end-ns 100] [--stride 10] [--min-freq 0.05]
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


# ── Dependências ───────────────────────────────────────────────────────────────

def _check_deps():
    try:
        import MDAnalysis
        from MDAnalysis.lib import distances as mda_dist
    except ImportError:
        sys.exit(
            "ERROR: MDAnalysis não instalado.\n"
            "  mamba install -n md-gromacs mdanalysis seaborn\n"
        )
    try:
        import seaborn
    except ImportError:
        sys.exit(
            "ERROR: seaborn não instalado.\n"
            "  mamba install -n md-gromacs seaborn\n"
        )

_check_deps()
import MDAnalysis as mda
from MDAnalysis.lib import distances as mda_dist
import seaborn as sns


# ── Parser do arquivo de índice GROMACS (.ndx) ────────────────────────────────

def parse_ndx(ndx_path: str) -> dict:
    """Retorna dict {nome_grupo: [índices 0-based]}."""
    groups, current = {}, None
    with open(ndx_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith('['):
                current = line.strip('[]').strip()
                groups[current] = []
            elif current and line:
                groups[current].extend(int(x) - 1 for x in line.split())
    return groups


def _generate_gro_ref(tpr: str, xtc: str, gro_out: str):
    """Extrai frame 0 como GRO usando gmx_mpi trjconv (fallback quando TPR muito novo)."""
    import subprocess
    for gmx in ('mpirun -np 1 gmx_mpi', 'gmx_mpi', 'gmx'):
        cmd = f'echo "0" | {gmx} trjconv -s "{tpr}" -f "{xtc}" -dump 0 -o "{gro_out}" -nobackup 2>/dev/null'
        ret = subprocess.run(cmd, shell=True)
        if ret.returncode == 0 and Path(gro_out).exists():
            print(f"[contact_map] Topologia GRO gerada: {gro_out}")
            return
    raise RuntimeError(
        "Não foi possível gerar GRO de referência.\n"
        "Tente manualmente: echo '0' | mpirun -np 1 gmx_mpi trjconv "
        f"-s {tpr} -f {xtc} -dump 0 -o frame0_ref.gro"
    )


def load_universe(tpr: str, xtc: str) -> mda.Universe:
    """Carrega Universe do MDAnalysis. Se TPR versão não suportada, gera GRO como fallback."""
    try:
        return mda.Universe(tpr, xtc)
    except (ValueError, NotImplementedError) as e:
        err = str(e)
        if 'tpx version' not in err.lower() and 'not support' not in err.lower():
            raise
        # TPR muito novo para esta versão do MDAnalysis — usa GRO como topologia
        print(f"[contact_map] AVISO: {err}")
        print("[contact_map] Gerando topologia GRO alternativa via gmx_mpi trjconv...")
        gro_ref = str(Path(tpr).parent / 'frame0_ref.gro')
        if not Path(gro_ref).exists():
            _generate_gro_ref(tpr, xtc, gro_ref)
        print(f"[contact_map] Carregando com GRO: {gro_ref}")
        return mda.Universe(gro_ref, xtc)


# ── Mapa de contato vetorizado ────────────────────────────────────────────────

def compute_contact_map(
    tpr: str, xtc: str, ndx: str,
    cutoff_nm: float = 0.4,
    start_ns: float = 0.0,
    end_ns: float | None = None,
    stride: int = 10,
) -> tuple[np.ndarray, list, list, int]:
    """
    Retorna:
        freq_matrix  — shape (n_rec_res, n_lig_res), valores 0-1
        rec_labels   — lista de labels dos resíduos do receptor
        lig_labels   — lista de labels dos resíduos do ligante
        n_frames     — número de frames analisados
    """
    cutoff_ang = cutoff_nm * 10.0  # MDAnalysis usa Ångströms

    u = load_universe(tpr, xtc)
    ndx_groups = parse_ndx(ndx)

    for grp in ('Receptor', 'Ligante'):
        if grp not in ndx_groups:
            sys.exit(f"ERROR: grupo '{grp}' não encontrado em {ndx}. "
                     f"Grupos disponíveis: {list(ndx_groups.keys())}")

    rec_sel = u.atoms[ndx_groups['Receptor']]
    lig_sel = u.atoms[ndx_groups['Ligante']]

    # Ordena resíduos únicos (por resindex interno do MDAnalysis — estável)
    rec_res_sorted = sorted(set(zip(rec_sel.resindices, rec_sel.resids,
                                     rec_sel.resnames)),
                             key=lambda x: x[0])
    lig_res_sorted = sorted(set(zip(lig_sel.resindices, lig_sel.resids,
                                     lig_sel.resnames)),
                             key=lambda x: x[0])

    rec_residx_list = [x[0] for x in rec_res_sorted]
    lig_residx_list = [x[0] for x in lig_res_sorted]
    rec_labels = [f"{x[2]}{x[1]}" for x in rec_res_sorted]
    lig_labels = [f"{x[2]}{x[1]}" for x in lig_res_sorted]

    n_rec_res = len(rec_residx_list)
    n_lig_res = len(lig_residx_list)

    # Mapeamento átomo → índice de posição na matriz
    residx2row_rec = {ri: i for i, ri in enumerate(rec_residx_list)}
    residx2row_lig = {ri: i for i, ri in enumerate(lig_residx_list)}
    rec_atom_row = np.array([residx2row_rec[r] for r in rec_sel.resindices])
    lig_atom_row = np.array([residx2row_lig[r] for r in lig_sel.resindices])

    n_lig_atoms = len(lig_sel)

    print(f"[contact_map] Receptor: {n_rec_res} resíduos ({len(rec_sel)} átomos)")
    print(f"[contact_map] Ligante:  {n_lig_res} resíduos ({len(lig_sel)} átomos)")
    print(f"[contact_map] Cutoff: {cutoff_nm} nm | Stride: {stride} frames")
    if end_ns:
        print(f"[contact_map] Janela temporal: {start_ns}–{end_ns} ns")

    freq_sum = np.zeros((n_rec_res, n_lig_res), dtype=np.float32)
    n_frames = 0

    for ts in u.trajectory[::stride]:
        t_ns = ts.time / 1000.0
        if t_ns < start_ns:
            continue
        if end_ns is not None and t_ns > end_ns:
            break

        # Matriz de distâncias átomo × átomo (Å)
        d = mda_dist.distance_array(rec_sel.positions, lig_sel.positions)
        contacts = (d < cutoff_ang).astype(np.int8)  # (n_rec_atoms, n_lig_atoms)

        # Redução 1: agrupa por resíduo do receptor
        # tmp[i_rec_res, j_lig_atom] = max dos contatos de todos os átomos de rec_res_i
        tmp = np.zeros((n_rec_res, n_lig_atoms), dtype=np.int8)
        np.maximum.at(tmp, rec_atom_row, contacts)

        # Redução 2: agrupa por resíduo do ligante
        frame_contact = np.zeros((n_lig_res, n_rec_res), dtype=np.int8)
        np.maximum.at(frame_contact, lig_atom_row, tmp.T)
        freq_sum += frame_contact.T.astype(np.float32)

        n_frames += 1
        if n_frames % 100 == 0:
            print(f"[contact_map]   frame {n_frames} | t = {t_ns:.1f} ns")

    if n_frames == 0:
        sys.exit("ERROR: nenhum frame processado. Verifique --start-ns e --end-ns.")

    freq_matrix = freq_sum / n_frames
    print(f"[contact_map] Total de frames: {n_frames}")
    return freq_matrix, rec_labels, lig_labels, n_frames


# ── Visualização ──────────────────────────────────────────────────────────────

def plot_contact_map(
    freq_matrix: np.ndarray,
    rec_labels: list,
    lig_labels: list,
    sample_id: str,
    n_frames: int,
    cutoff_nm: float,
    min_freq: float,
    outdir: Path,
):
    # Filtra resíduos com pelo menos min_freq de contato
    mask_rec = freq_matrix.max(axis=1) >= min_freq
    mask_lig = freq_matrix.max(axis=0) >= min_freq

    if mask_rec.sum() == 0:
        mask_rec[:] = True
    if mask_lig.sum() == 0:
        mask_lig[:] = True

    df_full = pd.DataFrame(freq_matrix, index=rec_labels, columns=lig_labels)
    df_plot = df_full.loc[np.array(rec_labels)[mask_rec],
                           np.array(lig_labels)[mask_lig]]

    fig_h = max(5, df_plot.shape[0] * 0.25)
    fig_w = max(4, df_plot.shape[1] * 0.55)

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    sns.heatmap(df_plot,
                ax=ax,
                cmap='YlOrRd',
                vmin=0, vmax=1,
                cbar_kws={'label': 'Frequência de contato'},
                linewidths=0.2,
                linecolor='lightgray',
                annot=(df_plot.shape[0] * df_plot.shape[1] <= 200))

    ax.set_title(
        f'{sample_id} — Mapa de Contato (cutoff {cutoff_nm} nm)\n'
        f'{n_frames} frames | resíduos com freq ≥ {min_freq:.0%}',
        fontsize=12,
    )
    ax.set_xlabel('Ligante', fontsize=10)
    ax.set_ylabel('Receptor', fontsize=10)
    ax.tick_params(axis='x', rotation=45, labelsize=8)
    ax.tick_params(axis='y', rotation=0, labelsize=8)
    plt.tight_layout()

    png = outdir / 'contact_map.png'
    fig.savefig(png, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"[contact_map] Salvo: {png}")
    return df_full


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(description='Contact map receptor × ligante')
    ap.add_argument('--tpr',       required=True)
    ap.add_argument('--xtc',       required=True)
    ap.add_argument('--ndx',       required=True)
    ap.add_argument('--outdir',    required=True)
    ap.add_argument('--sample-id', required=True)
    ap.add_argument('--cutoff',    type=float, default=0.4,
                    help='Distância de cutoff em nm (padrão: 0.4)')
    ap.add_argument('--start-ns',  type=float, default=0.0)
    ap.add_argument('--end-ns',    type=float, default=None)
    ap.add_argument('--stride',    type=int,   default=10,
                    help='Usar 1 a cada N frames (padrão: 10)')
    ap.add_argument('--min-freq',  type=float, default=0.05,
                    help='Frequência mínima para exibir no heatmap (padrão: 0.05)')
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    freq_matrix, rec_labels, lig_labels, n_frames = compute_contact_map(
        tpr=args.tpr, xtc=args.xtc, ndx=args.ndx,
        cutoff_nm=args.cutoff,
        start_ns=args.start_ns,
        end_ns=args.end_ns,
        stride=args.stride,
    )

    # Salva CSV completo
    df_full = plot_contact_map(
        freq_matrix, rec_labels, lig_labels,
        sample_id=args.sample_id,
        n_frames=n_frames,
        cutoff_nm=args.cutoff,
        min_freq=args.min_freq,
        outdir=outdir,
    )
    csv_path = outdir / 'contact_map.csv'
    df_full.to_csv(csv_path, float_format='%.4f')
    print(f"[contact_map] Salvo: {csv_path}")

    # Top contatos
    flat = [
        (rec_labels[i], lig_labels[j], freq_matrix[i, j])
        for i in range(len(rec_labels))
        for j in range(len(lig_labels))
        if freq_matrix[i, j] >= args.min_freq
    ]
    flat.sort(key=lambda x: -x[2])
    print(f"\n[contact_map] Top 20 contatos (receptor × ligante × freq):")
    for r, l, f in flat[:20]:
        print(f"  {r:12s} × {l:12s}  {f:.3f} ({f*100:.1f}%)")

    # Salva top contacts CSV
    df_top = pd.DataFrame(flat, columns=['receptor_residue', 'ligand_residue', 'frequency'])
    df_top.to_csv(outdir / 'top_contacts.csv', index=False, float_format='%.4f')
    print(f"\n[contact_map] Concluído: {n_frames} frames analisados")


if __name__ == '__main__':
    main()
