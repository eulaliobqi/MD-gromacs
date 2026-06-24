#!/usr/bin/env python3
"""
ProLIF interaction fingerprints para complexos receptor-ligante da DM.

Compatível com:
  - GORE4  (peptídeo 5 aa, chain B)
  - BEN    (benzamidina, HETATM resname BNZ/B)
  - SKTI   (proteína 177 aa — usa contact-map residue mode)

Saídas (em --outdir):
  prolif_barcode.png       — heatmap por frame (branco=sem interação, colorido=presente)
  prolif_persistence.png   — barplot de persistência (% frames) por par resíduo:tipo
  prolif_fingerprint.csv   — matriz completa frame × interação
  prolif_summary.csv       — persistência por interação (top 30)

Instalação:
  mamba install -n md-gromacs prolif mdanalysis

Uso:
    python3 bin/run_prolif.py \
        --tpr  results/ACR157-GORE4/ACR157-GORE4/prod/md.tpr \
        --xtc  results/ACR157-GORE4/ACR157-GORE4/prod/md_fit.xtc \
        --ndx  results/ACR157-GORE4/ACR157-GORE4/analise/lig.ndx \
        --outdir results/ACR157-GORE4/ACR157-GORE4/prolif \
        --sample-id ACR157-GORE4 \
        [--start-ns 0] [--end-ns 100] [--stride 10]
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
    missing = []
    try:
        import MDAnalysis
    except ImportError:
        missing.append("mdanalysis")
    try:
        import prolif
    except ImportError:
        missing.append("prolif")
    if missing:
        sys.exit(
            f"ERROR: pacotes ausentes: {', '.join(missing)}\n"
            f"  mamba install -n md-gromacs {' '.join(missing)}\n"
        )

_check_deps()
import MDAnalysis as mda
import prolif as plf


# ── Utilitários ───────────────────────────────────────────────────────────────

def parse_ndx(ndx_path: str) -> dict:
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


def _find_existing_topology(tpr: str) -> str | None:
    """Busca GRO/PDB existente no diretório do TPR (frame0_ref.gro, *.gro, *.pdb)."""
    d = Path(tpr).parent
    candidates = (
        [d / 'frame0_ref.gro'] +
        sorted(d.glob('*.gro')) +
        [d.parent / 'complexo' / 'complexo.pdb'] +
        sorted(d.glob('*.pdb'))
    )
    for c in candidates:
        if c.exists() and c.stat().st_size > 1000:
            return str(c)
    return None


def _generate_gro_ref(tpr: str, xtc: str, gro_out: str):
    """Tenta gerar GRO de referência via subprocess (gmx editconf/trjconv)."""
    import subprocess

    def try_cmd(cmd: str) -> bool:
        r = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        return r.returncode == 0 and Path(gro_out).exists()

    for gmx in ('mpirun -np 1 gmx_mpi', 'gmx_mpi', 'gmx'):
        if try_cmd(f'{gmx} editconf -f "{tpr}" -o "{gro_out}" -nobackup 2>&1'):
            print(f"[prolif] Topologia GRO gerada via editconf: {gro_out}")
            return

    for gmx in ('mpirun -np 1 gmx_mpi', 'gmx_mpi', 'gmx'):
        if try_cmd(
            f'echo "System" | {gmx} trjconv '
            f'-s "{tpr}" -f "{xtc}" -b 0 -e 0 -o "{gro_out}" -nobackup 2>&1'
        ):
            print(f"[prolif] Topologia GRO gerada via trjconv: {gro_out}")
            return

    raise RuntimeError(f"Não foi possível gerar GRO de referência para {tpr}")


def load_universe(tpr: str, xtc: str) -> mda.Universe:
    """Carrega Universe. Se TPR v138 não suportado: busca GRO/PDB existente, depois tenta gerar."""
    try:
        return mda.Universe(tpr, xtc)
    except (ValueError, NotImplementedError) as e:
        err = str(e)
        if 'tpx version' not in err.lower() and 'not support' not in err.lower():
            raise

        print(f"[prolif] TPR versão não suportada: {err}")

        # 1. Busca topologia alternativa já existente (md.gro, *.gro, etc.)
        alt = _find_existing_topology(tpr)
        if alt:
            print(f"[prolif] Usando topologia alternativa: {alt}")
            try:
                return mda.Universe(alt, xtc)
            except Exception as e2:
                print(f"[prolif]   falhou com {alt}: {e2}")

        # 2. Tenta gerar GRO via editconf/trjconv
        gro_ref = str(Path(tpr).parent / 'frame0_ref.gro')
        try:
            _generate_gro_ref(tpr, xtc, gro_ref)
            return mda.Universe(gro_ref, xtc)
        except RuntimeError:
            pass

        raise RuntimeError(
            f"MDAnalysis não suporta TPR v138. Corrija com:\n"
            f"  mamba update -n md-gromacs mdanalysis   # ou:\n"
            f"  mpirun -np 1 gmx_mpi editconf -f {tpr} -o {Path(tpr).parent}/frame0_ref.gro"
        )


def trajectory_slice(u, start_ns: float, end_ns: float | None, stride: int):
    """Retorna iterável de ts selecionados por tempo."""
    frames = []
    for ts in u.trajectory[::stride]:
        t_ns = ts.time / 1000.0
        if t_ns < start_ns:
            continue
        if end_ns is not None and t_ns > end_ns:
            break
        frames.append(ts.frame)
    return frames


# ── ProLIF fingerprint ─────────────────────────────────────────────────────────

def run_prolif(
    tpr: str, xtc: str, ndx: str,
    start_ns: float, end_ns: float | None, stride: int,
    sample_id: str, outdir: Path,
):
    u = load_universe(tpr, xtc)

    # ProLIF usa RDKitConverter que exige atributo 'elements'.
    # GRO não contém elementos — deduz a partir dos nomes dos átomos.
    try:
        elems = u.atoms.elements
        if len(elems) == 0 or not any(elems):
            raise AttributeError("elements vazios")
    except AttributeError:
        import warnings
        from MDAnalysis.topology import guessers
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            guessed = guessers.guess_types(u.atoms.names)
        u.add_TopologyAttr('elements', guessed)
        print("[prolif] 'elements' deduzidos de atom names (GRO não contém elementos)")

    ndx_groups = parse_ndx(ndx)

    for grp in ('Receptor', 'Ligante'):
        if grp not in ndx_groups:
            sys.exit(f"ERROR: grupo '{grp}' não encontrado em {ndx}.")

    rec_ag = u.atoms[ndx_groups['Receptor']]
    lig_ag = u.atoms[ndx_groups['Ligante']]

    # Verifica presença de H explícito no ligante (necessário para HBDonor/HBAcceptor)
    try:
        n_H = sum(1 for a in lig_ag.atoms
                  if a.name.startswith('H') or getattr(a, 'element', '') == 'H')
        if n_H == 0:
            print("[prolif] WARN: ligante sem H explícito — HBDonor/HBAcceptor serão omitidos")
    except Exception:
        n_H = -1  # não foi possível checar

    print(f"[prolif] {sample_id}")
    print(f"[prolif] Receptor: {len(rec_ag)} átomos, Ligante: {len(lig_ag)} átomos")
    print(f"[prolif] Stride: {stride} frames | Janela: {start_ns}–{end_ns or 'fim'} ns")

    # Converte ns → frame index (10 frames/ns = 100 ps/frame padrão do pipeline)
    frames_per_ns = 10
    start_frame = max(0, int(start_ns * frames_per_ns))
    end_frame   = int(end_ns * frames_per_ns) if end_ns is not None else None
    n_frames_est = (((end_frame or len(u.trajectory)) - start_frame) // stride) + 1
    print(f"[prolif] ~{n_frames_est} frames para análise (frames {start_frame}→{end_frame or 'fim'})")

    if n_frames_est <= 0:
        sys.exit("ERROR: nenhum frame na janela temporal especificada.")

    # ProLIF Fingerprint — usa interações disponíveis na versão instalada
    # Detecta dinamicamente para evitar NameError em versões diferentes
    PREFERRED = [
        "HBDonor", "HBAcceptor",
        "Hydrophobic",
        "Cationic", "Anionic",      # ProLIF 2.x (renomeados de "Ionic")
        "Ionic",                     # ProLIF 1.x
        "PiStacking", "CationPi", "PiCation",
        "VdWContact",
    ]
    try:
        # Testa cada interação individualmente para filtrar as não disponíveis
        avail = []
        rejected = []
        for name in PREFERRED:
            try:
                plf.Fingerprint(interactions=[name])
                avail.append(name)
            except (NameError, ValueError):
                rejected.append(name)
        if rejected:
            print(f"[prolif] WARN: interações indisponíveis nesta versão ProLIF: {rejected}")
        fp = plf.Fingerprint(interactions=avail) if avail else plf.Fingerprint()
        print(f"[prolif] Interações ativas: {avail or 'default (ProLIF fallback)'}")
    except Exception:
        fp = plf.Fingerprint()
        print("[prolif] Usando set de interações default do ProLIF")

    # ProLIF 2.x+ não aceita start/stop/step como kwargs — fatia a trajetória antes
    traj_slice = u.trajectory[start_frame:end_frame:stride]

    def _run_fp(force_h=False):
        """Executa fp.run(); se force_h, bypassa checagem de H via monkey-patch."""
        if force_h:
            # Ligante sem H explícito (ex: lig.ndx só com átomos pesados).
            # plf.Molecule.from_mda(ag, force=True) pula a checagem mas mantém
            # a conversão — HBDonor/HBAcceptor não serão detectados.
            _orig = plf.Molecule.from_mda
            plf.Molecule.from_mda = lambda ag, **kw: _orig(ag, force=True, **kw)
            try:
                fp.run(traj_slice, lig_ag, rec_ag, progress=True)
            finally:
                plf.Molecule.from_mda = _orig
        else:
            fp.run(traj_slice, lig_ag, rec_ag, progress=True)

    try:
        _run_fp(force_h=False)
    except TypeError:
        # ProLIF 1.x fallback
        fp.run(u.trajectory, lig_ag, rec_ag,
               start=start_frame, stop=end_frame, step=stride, verbose=True)
    except AttributeError as e:
        if "hydrogen" not in str(e).lower():
            raise
        print(f"[prolif] Ligante sem H — aplicando force=True (BEN/molécula pequena)")
        try:
            _run_fp(force_h=True)
        except Exception as e2:
            print(f"[prolif] AVISO: ProLIF não pôde processar {sample_id}: {e2}")
            return

    print("[prolif] Fingerprint calculado")

    # DataFrame de resultados
    df = fp.to_dataframe()
    if df.empty:
        print("[prolif] AVISO: nenhuma interação detectada. "
              "Verificar se o ligante tem atributos de elemento corretos no .tpr.")
        return

    # ProLIF 2.x usa MultiIndex nas colunas (lig_res, prot_res, interaction_type).
    # Achata para strings "ProtRes-LigRes-Tipo" antes de qualquer operação.
    if isinstance(df.columns, pd.MultiIndex):
        df.columns = [
            '-'.join(str(lvl) for lvl in col if str(lvl) != 'nan').strip('-')
            for col in df.columns
        ]

    df.to_csv(outdir / 'prolif_fingerprint.csv', float_format='%.0f')
    print(f"[prolif] Salvo: prolif_fingerprint.csv ({df.shape[0]} frames × {df.shape[1]} interações)")

    # Persistência (% de frames com cada interação)
    persistence = df.mean() * 100
    persistence = persistence[persistence > 0].sort_values(ascending=False)
    df_persist = pd.DataFrame({
        'interaction': persistence.index,
        'persistence_pct': persistence.values,
    })
    df_persist.to_csv(outdir / 'prolif_summary.csv', index=False, float_format='%.2f')

    # ── Plot 1: Barcode (heatmap frame × interação) ─────────────────────────
    df_bool = df.copy().astype(float)
    top_inter = persistence.head(40).index if len(persistence) >= 40 else persistence.index
    df_bar = df_bool[top_inter]

    if not df_bar.empty:
        n_frames = len(df)
        fig_w = max(10, len(top_inter) * 0.35)
        fig_h = max(4, n_frames * 0.015)
        fig, ax = plt.subplots(figsize=(fig_w, min(fig_h, 12)))
        ax.imshow(df_bar.values.T, aspect='auto', cmap='Blues',
                  vmin=0, vmax=1, interpolation='nearest')
        ax.set_yticks(range(len(top_inter)))
        ax.set_yticklabels([str(x).replace(', ', '\n') for x in top_inter], fontsize=6)
        ax.set_xlabel('Frame', fontsize=9)
        ax.set_title(f'{sample_id} — ProLIF Barcode (top {len(top_inter)} interações)', fontsize=11)
        plt.tight_layout()
        fig.savefig(outdir / 'prolif_barcode.png', dpi=150, bbox_inches='tight')
        plt.close(fig)
        print(f"[prolif] Salvo: prolif_barcode.png")

    # ── Plot 2: Persistência (barplot) ─────────────────────────────────────
    top30 = persistence.head(30)
    if not top30.empty:
        fig, ax = plt.subplots(figsize=(max(6, len(top30) * 0.45), 5))
        colors_map = {
            'HBDonor':    '#1f77b4',
            'HBAcceptor': '#aec7e8',
            'Hydrophobic': '#ff7f0e',
            'Ionic':       '#d62728',
            'PiStacking':  '#9467bd',
            'PiCation':    '#8c564b',
            'VdWContact':  '#7f7f7f',
        }
        bar_colors = [
            colors_map.get(str(inter).split('.')[-1], '#2ca02c')
            for inter in top30.index
        ]
        ax.bar(range(len(top30)), top30.values, color=bar_colors)
        ax.set_xticks(range(len(top30)))
        ax.set_xticklabels(
            [str(x).replace(', ', '\n') for x in top30.index],
            rotation=60, ha='right', fontsize=7,
        )
        ax.set_ylabel('Persistência (%)', fontsize=10)
        ax.set_title(
            f'{sample_id} — Persistência de Interações ProLIF\n'
            f'({len(df)} frames | {start_ns}–{end_ns or "fim"} ns)',
            fontsize=11,
        )
        ax.axhline(20, color='gray', lw=0.8, ls='--', label='20%')
        ax.axhline(50, color='red',  lw=0.8, ls='--', label='50%')
        ax.legend(fontsize=8)
        plt.tight_layout()
        fig.savefig(outdir / 'prolif_persistence.png', dpi=150, bbox_inches='tight')
        plt.close(fig)
        print(f"[prolif] Salvo: prolif_persistence.png")

    # Resumo no terminal
    print(f"\n[prolif] Top 20 interações mais persistentes:")
    for inter, pct in top30.head(20).items():
        print(f"  {str(inter):50s}  {pct:.1f}%")

    print(f"\n[prolif] Concluído para {sample_id}")


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(description='ProLIF fingerprints receptor-ligante')
    ap.add_argument('--tpr',       required=True)
    ap.add_argument('--xtc',       required=True)
    ap.add_argument('--ndx',       required=True)
    ap.add_argument('--outdir',    required=True)
    ap.add_argument('--sample-id', required=True)
    ap.add_argument('--start-ns',  type=float, default=0.0)
    ap.add_argument('--end-ns',    type=float, default=None)
    ap.add_argument('--stride',    type=int,   default=10)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    run_prolif(
        tpr=args.tpr, xtc=args.xtc, ndx=args.ndx,
        start_ns=args.start_ns,
        end_ns=args.end_ns,
        stride=args.stride,
        sample_id=args.sample_id,
        outdir=outdir,
    )


if __name__ == '__main__':
    main()
