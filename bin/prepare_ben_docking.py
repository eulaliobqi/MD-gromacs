#!/usr/bin/env python3
"""
Prepara receptor PDBQT e arquivo de config Vina para docking de BEN.
Usa os 4 resíduos catalíticos para calcular o centro da caixa de busca.

Uso:
    python3 bin/prepare_ben_docking.py \
        --pdb spodoptera-ben/ACR157-final.pdb \
        --name ACR157 \
        --outdir spodoptera-ben \
        --triad 69 114 211 205
"""

import argparse
import os
import subprocess
import sys
from pathlib import Path


def parse_pdb_ca(pdb_path, residue_numbers):
    """Extrai coordenadas do átomo Cα (ou N para glicina) dos resíduos especificados."""
    coords = {}
    with open(pdb_path) as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            atom_name = line[12:16].strip()
            resnum = int(line[22:26].strip())
            if resnum in residue_numbers and atom_name == "CA":
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                coords[resnum] = (x, y, z)
    missing = set(residue_numbers) - set(coords.keys())
    if missing:
        # fallback: tentar Cβ
        with open(pdb_path) as f:
            for line in f:
                if not line.startswith("ATOM"):
                    continue
                atom_name = line[12:16].strip()
                resnum = int(line[22:26].strip())
                if resnum in missing and atom_name == "CB":
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    coords[resnum] = (x, y, z)
    still_missing = set(residue_numbers) - set(coords.keys())
    if still_missing:
        print(f"AVISO: resíduos sem Cα/Cβ: {still_missing}", file=sys.stderr)
    return coords


def centroid(coords_dict):
    vals = list(coords_dict.values())
    cx = sum(v[0] for v in vals) / len(vals)
    cy = sum(v[1] for v in vals) / len(vals)
    cz = sum(v[2] for v in vals) / len(vals)
    return cx, cy, cz


def convert_receptor(pdb_path, pdbqt_path):
    """Converte receptor PDB → PDBQT com obabel (modo receptor, adiciona H polares)."""
    cmd = [
        "obabel",
        "-ipdb", str(pdb_path),
        "-opdbqt",
        "-O", str(pdbqt_path),
        "-h",         # adiciona hidrogênios polares (necessário para Vina)
        "-xr",        # modo receptor: torsões fixas
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print("ERRO obabel receptor:", result.stderr, file=sys.stderr)
        sys.exit(1)
    print(f"  Receptor PDBQT: {pdbqt_path}")


def write_vina_conf(conf_path, receptor_pdbqt, ligand_pdbqt, center, box_size=22.0,
                    exhaustiveness=32, num_modes=9, cpu=8):
    cx, cy, cz = center
    conf = f"""receptor = {receptor_pdbqt}
ligand = {ligand_pdbqt}

center_x = {cx:.3f}
center_y = {cy:.3f}
center_z = {cz:.3f}

size_x = {box_size:.1f}
size_y = {box_size:.1f}
size_z = {box_size:.1f}

exhaustiveness = {exhaustiveness}
num_modes = {num_modes}
cpu = {cpu}
energy_range = 3
"""
    with open(conf_path, "w") as f:
        f.write(conf)
    print(f"  Config Vina:    {conf_path}")
    print(f"  Box center:     ({cx:.2f}, {cy:.2f}, {cz:.2f})  size={box_size} Å")


def main():
    parser = argparse.ArgumentParser(description="Prepara receptor e config Vina para BEN docking")
    parser.add_argument("--pdb",    required=True, help="PDB do receptor (sem ligante)")
    parser.add_argument("--name",   required=True, help="ID do receptor (ex: ACR157)")
    parser.add_argument("--outdir", required=True, help="Diretório base (spodoptera-ben/)")
    parser.add_argument("--triad",  required=True, nargs=4, type=int,
                        help="4 números de resíduo do sítio ativo")
    parser.add_argument("--box",    type=float, default=22.0, help="Tamanho da caixa em Å (default 22)")
    parser.add_argument("--cpu",    type=int, default=8,   help="Threads para Vina (default 8)")
    parser.add_argument("--exhaustiveness", type=int, default=32, help="Exhaustiveness Vina (default 32)")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    pdbqt_dir = outdir / "receptors_pdbqt"
    conf_dir  = outdir / "configs"
    res_dir   = outdir / "results"
    for d in [pdbqt_dir, conf_dir, res_dir]:
        d.mkdir(parents=True, exist_ok=True)

    print(f"\n=== {args.name} ===")

    # 1. Extrair coordenadas dos resíduos catalíticos
    coords = parse_pdb_ca(args.pdb, args.triad)
    if len(coords) < 2:
        print("ERRO: menos de 2 resíduos encontrados no PDB.", file=sys.stderr)
        sys.exit(1)

    for rnum, xyz in sorted(coords.items()):
        print(f"  Res {rnum}: Cα = ({xyz[0]:.2f}, {xyz[1]:.2f}, {xyz[2]:.2f})")

    # 2. Centroide
    center = centroid(coords)

    # 3. Converter receptor para PDBQT
    pdbqt_path = pdbqt_dir / f"{args.name}.pdbqt"
    convert_receptor(args.pdb, pdbqt_path)

    # 4. Escrever config Vina
    ligand_pdbqt = str(outdir / "ligands" / "ben.pdbqt")
    out_pdbqt    = str(res_dir / f"{args.name}_ben_out.pdbqt")
    conf_path    = conf_dir / f"vina_{args.name}.conf"
    write_vina_conf(
        conf_path,
        receptor_pdbqt=str(pdbqt_path),
        ligand_pdbqt=ligand_pdbqt,
        center=center,
        box_size=args.box,
        exhaustiveness=args.exhaustiveness,
        cpu=args.cpu,
    )
    # adicionar out ao conf (sem log = — versão Vina customizada não suporta --log)
    with open(conf_path, "a") as f:
        f.write(f"\nout = {out_pdbqt}\n")

    print(f"  Pronto. Rode: vina --config {conf_path}")


if __name__ == "__main__":
    main()
