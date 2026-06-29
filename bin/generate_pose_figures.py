#!/usr/bin/env python3
"""
generate_pose_figures.py — Extrai frame representativo da DM e renderiza com PyMOL.

Para cada complexo com trajetória disponível:
  1. Extrai o frame de menor RMSD backbone à estrutura média da janela 70–100 ns
  2. Salva receptor.pdb + ligand.pdb separados
  3. Gera script PyMOL (.pml) com estilo de publicação
  4. Tenta executar PyMOL em modo headless (pymol -cq)
     Se PyMOL não estiver disponível, deixa os .pml para execução manual.

Saídas (em --outdir, padrão: poses/):
  {complex}/receptor.pdb
  {complex}/ligand.pdb
  {complex}/pose.pml
  {complex}/pose.png        ← se PyMOL disponível

Uso:
  mamba activate md-gromacs
  python3 bin/generate_pose_figures.py [--outdir poses] [--start-ns 70] [--target-ns 90]
"""

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).parent.parent

_SKIP = {"work", "__pycache__", "figures", "poses", "bin", "local_viz"}


# ── Utilitários NDX ────────────────────────────────────────────────────────

def parse_ndx(ndx_path: Path) -> dict:
    groups, current = {}, None
    with open(ndx_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith("["):
                current = line.strip("[]").strip()
                groups[current] = []
            elif current and line:
                groups[current].extend(int(x) - 1 for x in line.split())
    return groups


# ── Descoberta de trajetórias ──────────────────────────────────────────────

def find_trajectories() -> list[dict]:
    """
    Descobre pares (tpr/gro, xtc, ndx, prolif_summary) no repositório.
    Retorna lista de dicts com as chaves: sample_id, tpr, xtc, ndx, prolif_csv.
    """
    hits = []
    seen = set()

    for xtc in sorted(ROOT.glob("**/prod/md_fit.xtc")):
        if any(x in xtc.parts for x in _SKIP):
            continue
        prod_dir = xtc.parent

        # Topologia: tpr preferido; fallback GRO
        tpr = prod_dir / "md.tpr"
        if not tpr.exists():
            gros = list(prod_dir.glob("*.gro"))
            tpr = gros[0] if gros else None
        if not tpr or not tpr.exists():
            continue

        # NDX em analise/
        analise = prod_dir.parent / "analise"
        ndx = analise / "lig.ndx"
        if not ndx.exists():
            continue

        # ProLIF summary (opcional) — busca por glob para cobrir variações de path
        prolif_csv = None
        for candidate in sorted(prod_dir.parent.glob("**/prolif_summary.csv")):
            prolif_csv = candidate
            break
        if prolif_csv is None:
            # Fallback: sobe um nível (estrutura BEN: results/{ID}/{ID}/prod)
            for candidate in sorted(prod_dir.parent.parent.glob("**/prolif_summary.csv")):
                prolif_csv = candidate
                break

        # Sample ID = diretório acima de prod/
        sid = prod_dir.parent.name
        if sid in seen:
            continue
        seen.add(sid)

        hits.append({
            "sample_id": sid,
            "tpr":       tpr,
            "xtc":       xtc,
            "ndx":       ndx,
            "prolif_csv": prolif_csv,
        })

    return hits


# ── Extração de frame representativo ──────────────────────────────────────

def extract_representative_frame(
    tpr: Path, xtc: Path, ndx: Path,
    start_ns: float, target_ns: float,
    outdir: Path,
) -> tuple[Path, Path, float]:
    """
    Carrega trajetória, acha frame de menor RMSD backbone à média
    dentro da janela [start_ns, target_ns].
    Salva receptor.pdb e ligand.pdb em outdir.
    Retorna (receptor_pdb, ligand_pdb, tempo_ns).
    """
    import MDAnalysis as mda
    from MDAnalysis.analysis import rms, align

    ndx_groups = parse_ndx(ndx)
    for grp in ("Receptor", "Ligante"):
        if grp not in ndx_groups:
            raise RuntimeError(f"Grupo '{grp}' não encontrado em {ndx}")

    # Carrega Universe; tenta fallback para GRO se TPR não suportado
    try:
        u = mda.Universe(str(tpr), str(xtc))
    except (ValueError, NotImplementedError):
        gros = list(tpr.parent.glob("*.gro"))
        if not gros:
            raise
        u = mda.Universe(str(gros[0]), str(xtc))

    # Deduz elementos se ausentes (GRO não contém)
    try:
        if not any(u.atoms.elements):
            raise AttributeError
    except AttributeError:
        import warnings
        from MDAnalysis.topology import guessers
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            guessed = guessers.guess_types(u.atoms.names)
        u.add_TopologyAttr("elements", guessed)

    rec_ag = u.atoms[ndx_groups["Receptor"]]
    lig_ag = u.atoms[ndx_groups["Ligante"]]
    bb = rec_ag.select_atoms("backbone")

    # Coleta frames na janela [start_ns, target_ns]
    t0 = start_ns * 1000
    t1 = target_ns * 1000
    frame_indices = []
    frame_times   = []
    positions_bb  = []

    for ts in u.trajectory:
        if ts.time < t0:
            continue
        if ts.time > t1:
            break
        frame_indices.append(ts.frame)
        frame_times.append(ts.time)
        positions_bb.append(bb.positions.copy())

    if not frame_indices:
        # Fallback: último frame disponível
        u.trajectory[-1]
        frame_indices = [u.trajectory.ts.frame]
        frame_times   = [u.trajectory.ts.time]
        positions_bb  = [bb.positions.copy()]

    # Frame com menor RMSD backbone à média do período
    pos_arr = np.array(positions_bb)
    avg_pos = pos_arr.mean(axis=0)
    rmsd_vals = np.sqrt(((pos_arr - avg_pos) ** 2).mean(axis=(1, 2)))
    best_idx = int(np.argmin(rmsd_vals))
    best_frame = frame_indices[best_idx]
    actual_ns  = frame_times[best_idx] / 1000.0

    u.trajectory[best_frame]
    print(f"  Frame representativo: {best_frame}  (t = {actual_ns:.1f} ns, "
          f"RMSD à média = {rmsd_vals[best_idx]:.3f} Å)")

    # Salva PDBs
    rec_pdb = outdir / "receptor.pdb"
    lig_pdb = outdir / "ligand.pdb"
    rec_ag.write(str(rec_pdb))
    lig_ag.write(str(lig_pdb))

    return rec_pdb, lig_pdb, actual_ns


# ── Resíduos-chave do ProLIF ───────────────────────────────────────────────

INTERACTION_TYPES = [
    "HBDonor", "HBAcceptor", "Hydrophobic",
    "Ionic", "Cationic", "Anionic",
    "PiStacking", "PiCation", "CationPi",
]

def load_key_residues(prolif_csv: Path | None, top_n: int = 8) -> list[str]:
    """Retorna lista de números de resíduos (strings) das top interações ProLIF."""
    if not prolif_csv:
        return []
    import re
    try:
        import pandas as pd
        df = pd.read_csv(prolif_csv)
        df = df.sort_values("persistence_pct", ascending=False)
        residues = []
        for inter in df["interaction"].head(top_n * 3):
            # Extrai número de resíduo (ex: "ASP189" → "189", "SER195.A" → "195")
            nums = re.findall(r"(\d{2,})", str(inter))
            for n in nums:
                if n not in residues and int(n) < 9999:
                    residues.append(n)
                    if len(residues) >= top_n:
                        break
            if len(residues) >= top_n:
                break
        return residues
    except Exception:
        return []


# ── Script PyMOL ───────────────────────────────────────────────────────────

# Template para pequenos ligantes / peptídeos curtos (< 30 resíduos)
PYMOL_TEMPLATE_SMALL = """\
# PyMOL publication figure — {sample_id}  (t = {time_ns:.1f} ns)
load {rec_pdb}, receptor
load {lig_pdb}, ligand
bg_color white

hide everything, receptor
show cartoon, receptor
color gray75, receptor
set cartoon_fancy_helices, 1
set cartoon_rect_length, 1.2

hide everything, ligand
show sticks, ligand
color cyan, ligand
set stick_radius, 0.18
show spheres, ligand and name CA
set sphere_scale, 0.25, ligand and name CA

{key_res_block}
distance hbonds, receptor, ligand, 3.5, mode=2
set dash_color, red
set dash_width, 2.5
set dash_gap, 0.35
set dash_radius, 0.06
hide labels, hbonds

create pocket, byres receptor within 5 of ligand
show surface, pocket
color tv_blue, pocket
set transparency, 0.60, pocket

orient ligand
zoom ligand, 10
rotate y, 15
rotate x, -10

set ray_shadows, 1
set ray_opaque_background, 1
set antialias, 2
set ambient, 0.35
set direct, 0.65
set specular, 0.20

ray 2400, 1800
png {png_out}, dpi=300, ray=1
quit
"""

# Template para inibidores proteicos grandes (≥ 30 resíduos) como SKTI
PYMOL_TEMPLATE_PROTEIN = """\
# PyMOL publication figure — {sample_id}  (t = {time_ns:.1f} ns)
load {rec_pdb}, receptor
load {lig_pdb}, ligand
bg_color white

hide everything, receptor
show cartoon, receptor
color gray75, receptor
set cartoon_fancy_helices, 1

hide everything, ligand
show cartoon, ligand
color cyan, ligand
set cartoon_transparency, 0.0, ligand

{key_res_block}
distance hbonds, receptor, ligand, 3.5, mode=2
set dash_color, red
set dash_width, 2.5
set dash_gap, 0.35
set dash_radius, 0.06
hide labels, hbonds

orient ligand
zoom ligand, 15
rotate y, 15
rotate x, -10

set ray_shadows, 1
set ray_opaque_background, 1
set antialias, 2
set ambient, 0.35
set direct, 0.65
set specular, 0.20

ray 2400, 1800
png {png_out}, dpi=300, ray=1
quit
"""

PYMOL_TEMPLATE = PYMOL_TEMPLATE_SMALL  # default, redefinido em write_pymol_script

KEY_RES_WITH = """\
select key_res, receptor and resi {res_sel}
show sticks, key_res
color tv_yellow, key_res and name C*
set stick_radius, 0.15, key_res
"""

KEY_RES_NONE = "# (nenhum resíduo-chave disponível)"


def _count_residues_pdb(pdb_path: Path) -> int:
    """Conta resíduos únicos no PDB pelo campo ATOM."""
    residues = set()
    try:
        for line in open(pdb_path):
            if line.startswith(("ATOM", "HETATM")):
                residues.add((line[21], line[22:26].strip()))
    except Exception:
        pass
    return len(residues)


def write_pymol_script(
    sample_id: str,
    rec_pdb: Path,
    lig_pdb: Path,
    key_residues: list[str],
    outdir: Path,
    time_ns: float,
) -> Path:
    png_out = outdir / "pose.png"
    res_sel = "+".join(key_residues) if key_residues else ""
    key_block = KEY_RES_WITH.format(res_sel=res_sel) if res_sel else KEY_RES_NONE

    # Escolhe template conforme tamanho do ligante
    n_res = _count_residues_pdb(lig_pdb)
    template = PYMOL_TEMPLATE_PROTEIN if n_res >= 30 else PYMOL_TEMPLATE_SMALL
    lig_type = "protein inhibitor" if n_res >= 30 else "small ligand/peptide"
    print(f"  Ligante: {n_res} resíduos → {lig_type}")

    script = template.format(
        sample_id=sample_id,
        time_ns=time_ns,
        rec_pdb=str(rec_pdb.resolve()),
        lig_pdb=str(lig_pdb.resolve()),
        key_res_block=key_block,
        top_n=len(key_residues),
        png_out=str(png_out.resolve()),
    )
    pml_path = outdir / "pose.pml"
    pml_path.write_text(script)
    return pml_path


# ── Execução PyMOL ─────────────────────────────────────────────────────────

def run_pymol(pml_path: Path) -> bool:
    for cmd in ("pymol", "pymol3", "/usr/bin/pymol", "/opt/pymol/bin/pymol"):
        if shutil.which(cmd):
            try:
                r = subprocess.run(
                    [cmd, "-cq", str(pml_path)],
                    check=True, timeout=180,
                    capture_output=True, text=True,
                )
                return True
            except subprocess.CalledProcessError as e:
                # PyMOL escreve erros em stdout, não em stderr
                out = ((e.stdout or "") + (e.stderr or "")).strip()
                if out:
                    tail = "\n".join(out.splitlines()[-8:])
                    print(f"  [PyMOL output]\n{tail}")
            except subprocess.TimeoutExpired:
                print("  [PyMOL] timeout após 180 s")
    return False


# ── Main ───────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--outdir",    default="poses",
                    help="Diretório raiz de saída (default: poses/)")
    ap.add_argument("--start-ns",  type=float, default=70.0,
                    help="Início da janela para frame representativo (default: 70 ns)")
    ap.add_argument("--target-ns", type=float, default=100.0,
                    help="Fim da janela (default: 100 ns)")
    ap.add_argument("--top-res",   type=int,   default=8,
                    help="N.º de resíduos-chave a destacar (default: 8)")
    args = ap.parse_args()

    outdir_root = ROOT / args.outdir
    outdir_root.mkdir(parents=True, exist_ok=True)
    print(f"Output: {outdir_root}\n")

    trajectories = find_trajectories()
    if not trajectories:
        print("[WARN] Nenhuma trajetória encontrada (prod/md_fit.xtc + analise/lig.ndx).")
        return

    print(f"Complexos encontrados: {len(trajectories)}")
    for t in trajectories:
        print(f"  {t['sample_id']}")
    print()

    has_pymol = any(shutil.which(c) for c in ("pymol", "pymol3"))
    if not has_pymol:
        print("[INFO] PyMOL não encontrado — scripts .pml serão gerados mas não executados.")
        print("  Instale: mamba install -c schrodinger pymol-bundle\n")

    for entry in trajectories:
        sid = entry["sample_id"]
        print(f"[{sid}]")
        out = outdir_root / sid
        out.mkdir(exist_ok=True)

        try:
            rec_pdb, lig_pdb, time_ns = extract_representative_frame(
                entry["tpr"], entry["xtc"], entry["ndx"],
                start_ns=args.start_ns,
                target_ns=args.target_ns,
                outdir=out,
            )
        except Exception as e:
            print(f"  [ERRO] extração de frame: {e}")
            continue

        key_res = load_key_residues(entry["prolif_csv"], top_n=args.top_res)
        print(f"  Resíduos-chave: {key_res or '(nenhum prolif_summary.csv)'}")

        pml = write_pymol_script(sid, rec_pdb, lig_pdb, key_res, out, time_ns)
        print(f"  Script PyMOL: {pml.relative_to(ROOT)}")

        if has_pymol:
            ok = run_pymol(pml)
            if ok:
                png = out / "pose.png"
                print(f"  [OK] {png.relative_to(ROOT)}  ({png.stat().st_size // 1024} KB)")
            else:
                print("  [WARN] PyMOL retornou erro — verifique o .pml manualmente")
        else:
            print("  [SKIP] render — execute manualmente: pymol -cq " + str(pml.relative_to(ROOT)))

        print()

    print("[OK] done")
    if not has_pymol:
        print("\nPara instalar PyMOL no servidor:")
        print("  mamba install -n md-gromacs -c schrodinger pymol-bundle")
        print("  # ou:")
        print("  mamba install -n md-gromacs -c conda-forge pymol-open-source")


if __name__ == "__main__":
    main()
