#!/usr/bin/env python3
"""
Mescla topologia de receptor (pdb2gmx) com ligante pequeno (ACPYPE/GAFF2).

Operações:
1. Mescla receptor.gro + ligand.gro → complexo.gro (com renumeração de átomos)
2. Renumera resíduo do ligante para sequência após receptor
3. Patcheia receptor.top para incluir ligand.itp e adicionar ligante em [molecules]

Uso:
    python3 bin/merge_small_molecule_topology.py \\
        --protein-gro receptor.gro \\
        --ligand-gro  BEN.acpype/BEN_GMX.gro \\
        --protein-top receptor.top \\
        --ligand-itp  BEN.acpype/BEN_GMX.itp \\
        --ligand-mol  BEN \\
        --out-gro     complexo.gro \\
        --out-top     topol.top
"""
import argparse
import re
import sys
from pathlib import Path


def read_gro(gro_path):
    """Lê arquivo GRO → (titulo, natoms, linhas_átomos, linha_box)."""
    text = Path(gro_path).read_text()
    lines = text.splitlines()
    title  = lines[0]
    natoms = int(lines[1].strip())
    atoms  = lines[2 : 2 + natoms]
    box    = lines[2 + natoms]
    return title, natoms, atoms, box


def last_resnum(atom_lines):
    """Retorna número do último resíduo de uma lista de linhas GRO."""
    last = 1
    for line in atom_lines:
        try:
            last = int(line[:5].strip())
        except ValueError:
            pass
    return last


def merge_gro(protein_gro, ligand_gro, out_gro):
    """Mescla dois GROs num único arquivo com átomos renumerados."""
    _, n_prot, prot_atoms, box = read_gro(protein_gro)
    _, n_lig, lig_atoms, _    = read_gro(ligand_gro)

    last_res = last_resnum(prot_atoms)
    lig_res  = last_res + 1

    # Atualiza número de resíduo do ligante (col 0-4) para seguir receptor
    lig_renumbered = []
    for line in lig_atoms:
        lig_renumbered.append(f"{lig_res:5d}" + line[5:])

    combined = prot_atoms + lig_renumbered

    # Renumera átomos sequencialmente (col 15-19)
    final = []
    for i, line in enumerate(combined, 1):
        final.append(line[:15] + f"{i % 100000:5d}" + line[20:])

    total = n_prot + n_lig
    with open(out_gro, "w") as f:
        f.write("Complex protein + BEN (GAFF2)\n")
        f.write(f"{total:5d}\n")
        for line in final:
            f.write(line + "\n")
        f.write(box + "\n")

    print(f"GRO mesclado: {n_prot} prot + {n_lig} BEN = {total} átomos | BEN = resíduo {lig_res}")
    return lig_res


def patch_topology(protein_top, ligand_itp, ligand_mol, out_top):
    """
    Insere #include do ITP do ligante no topol.top do receptor.
    Adiciona entrada do ligante em [molecules].
    """
    content = Path(protein_top).read_text()
    itp_name = Path(ligand_itp).name

    # Insere #include após a linha de include do forcefield
    ff_pat = r'(#include\s+"[^"]+\.ff/forcefield\.itp")'
    if re.search(ff_pat, content):
        content = re.sub(ff_pat, rf'\1\n#include "{itp_name}"', content, count=1)
    else:
        # Fallback: antes da primeira [ moleculetype ]
        content = re.sub(
            r'(\[\s*moleculetype\s*\])',
            rf'#include "{itp_name}"\n\n\1',
            content, count=1
        )

    # Adiciona ligante em [molecules] (ao final do arquivo)
    content = content.rstrip("\n") + f"\n{ligand_mol:<20} 1\n"

    Path(out_top).write_text(content)
    print(f"TOP patcheado: #include {itp_name!r} + {ligand_mol} 1 → {out_top}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--protein-gro", required=True)
    ap.add_argument("--ligand-gro",  required=True)
    ap.add_argument("--protein-top", required=True)
    ap.add_argument("--ligand-itp",  required=True)
    ap.add_argument("--ligand-mol",  default="BEN")
    ap.add_argument("--out-gro",     required=True)
    ap.add_argument("--out-top",     required=True)
    args = ap.parse_args()

    merge_gro(args.protein_gro, args.ligand_gro, args.out_gro)
    patch_topology(args.protein_top, args.ligand_itp, args.ligand_mol, args.out_top)
    return 0


if __name__ == "__main__":
    sys.exit(main())
