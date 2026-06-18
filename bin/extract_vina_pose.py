#!/usr/bin/env python3
"""
Extrai modo N de arquivo PDBQT de saída do AutoDock Vina e salva como PDB.

Uso:
    python3 bin/extract_vina_pose.py \
        --pdbqt spodoptera-ben/results/ACR157_ben_out.pdbqt \
        --out spodoptera-ben/ben_poses/ACR157_ben_pose1.pdb \
        --mode 1 --resname BEN
"""
import argparse
import sys
from pathlib import Path


def extract_mode(pdbqt_path, mode=1):
    atoms = []
    current_mode = 0
    in_target = False
    with open(pdbqt_path) as f:
        for line in f:
            if line.startswith("MODEL"):
                current_mode += 1
                in_target = (current_mode == mode)
            elif line.startswith("ENDMDL"):
                if in_target:
                    break
            elif in_target and (line.startswith("ATOM") or line.startswith("HETATM")):
                atoms.append(line)
    return atoms


def write_pdb(atoms, out_path, residue_name="BEN", chain="B", resnum=1):
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    lines = []
    for i, line in enumerate(atoms, 1):
        # PDBQT atom line: cols 0-53 igual PDB, cols 54+ diferem (charges, tipo AD4)
        atom_name = line[12:16]
        # Determinar elemento: último campo não-vazio ou primeiro char do nome
        elem = line[77:79].strip() if len(line) > 77 else ""
        if not elem:
            elem = line[12:16].strip().lstrip("0123456789")[:1]

        # Construir linha PDB padrão
        pdb_line = (
            f"HETATM{i:5d} {atom_name}"
            f"{residue_name:>3s} {chain}{resnum:4d}    "
            f"{line[30:38]}{line[38:46]}{line[46:54]}"
            f"  1.00  0.00          {elem:>2s}\n"
        )
        lines.append(pdb_line)
    lines.append("TER\n")
    lines.append("END\n")

    with open(out_path, "w") as f:
        f.writelines(lines)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdbqt",   required=True)
    parser.add_argument("--out",     required=True)
    parser.add_argument("--mode",    type=int, default=1)
    parser.add_argument("--resname", default="BEN")
    args = parser.parse_args()

    atoms = extract_mode(args.pdbqt, args.mode)
    if not atoms:
        print(f"ERRO: modo {args.mode} não encontrado em {args.pdbqt}", file=sys.stderr)
        return 1

    write_pdb(atoms, args.out, residue_name=args.resname)
    print(f"Pose {args.mode} extraída → {args.out} ({len(atoms)} átomos)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
