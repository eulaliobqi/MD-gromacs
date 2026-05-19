#!/usr/bin/env bash
# Wrapper standalone para re-executar análises sem refazer a DM.
# Uso: run_analyses.sh <BASE_DIR>
# Onde BASE_DIR tem: prep/complexo.pdb, prod/md.tpr, prod/md_fit.xtc
set -euo pipefail

BASE="${1:-${HOME}/gromacs/dm_complex}"
[[ ! -d "$BASE/prod" ]] && { echo "ERRO: $BASE/prod nao existe"; exit 1; }
[[ ! -s "$BASE/prod/md.tpr"     ]] && { echo "ERRO: md.tpr nao encontrado";     exit 1; }
[[ ! -s "$BASE/prod/md_fit.xtc" ]] && { echo "ERRO: md_fit.xtc nao encontrado"; exit 1; }

GMX=$(command -v gmx_mpi 2>/dev/null || command -v gmx 2>/dev/null)
[[ -z "$GMX" ]] && { echo "ERRO: gmx nao encontrado"; exit 1; }

PDB_COMPLEX="$BASE/prep/complexo.pdb"
[[ ! -s "$PDB_COMPLEX" ]] && { echo "ERRO: $PDB_COMPLEX nao encontrado"; exit 1; }

mkdir -p "$BASE/analise"
cd "$BASE/analise"

ln -sf "$BASE/prod/md.tpr"     sistema.tpr
ln -sf "$BASE/prod/md_fit.xtc" trajetoria.xtc

LIG_FIRST=$(awk '/^ATOM/ && substr($0,22,1)=="B" {print substr($0,23,4)+0; exit}' "$PDB_COMPLEX")
LIG_LAST=$(awk  '/^ATOM/ && substr($0,22,1)=="B" {r=substr($0,23,4)+0} END{print r}' "$PDB_COMPLEX")

[[ -z "$LIG_FIRST" || -z "$LIG_LAST" ]] && { echo "ERRO: cadeia B nao detectada"; exit 1; }
echo "[ANALISE] Ligante: residuos $LIG_FIRST a $LIG_LAST"

N_DEFAULT=$("$GMX" make_ndx -f sistema.tpr -o _default.ndx 2>&1 \
    | grep -cE "^ *[0-9]+ +[A-Za-z]" || echo "10")
LIG_IDX=$N_DEFAULT
REC_IDX=$((N_DEFAULT + 1))
rm -f _default.ndx

"$GMX" make_ndx -f sistema.tpr -o lig.ndx << EOF >/dev/null 2>&1
r ${LIG_FIRST}-${LIG_LAST}
name ${LIG_IDX} Ligante
1 & ! ${LIG_IDX}
name ${REC_IDX} Receptor
q
EOF

[[ ! -s rmsd_backbone.xvg ]] && \
    printf 'Backbone\nBackbone\n' | "$GMX" rms \
        -s sistema.tpr -f trajetoria.xtc -n lig.ndx -o rmsd_backbone.xvg -tu ns

[[ ! -s rmsd_ligante.xvg ]] && \
    printf 'Ligante\nLigante\n' | "$GMX" rms \
        -s sistema.tpr -f trajetoria.xtc -n lig.ndx -o rmsd_ligante.xvg -tu ns

[[ ! -s rmsf_residuos.xvg ]] && \
    printf 'Backbone\n' | "$GMX" rmsf \
        -s sistema.tpr -f trajetoria.xtc -n lig.ndx -o rmsf_residuos.xvg -res -fit

[[ ! -s gyrate.xvg ]] && \
    printf 'Protein\n' | "$GMX" gyrate \
        -s sistema.tpr -f trajetoria.xtc -n lig.ndx -o gyrate.xvg -tu ns

rm -f numcont.xvg mindist.xvg
printf 'Receptor\nLigante\n' | "$GMX" mindist \
    -s sistema.tpr -f trajetoria.xtc -n lig.ndx \
    -od mindist.xvg -on numcont.xvg -d 0.4 -tu ns

rm -f hbond.xvg
printf 'Receptor\nLigante\n' | "$GMX" hbond \
    -s sistema.tpr -f trajetoria.xtc -n lig.ndx -num hbond.xvg -tu ns

echo "=== Analises concluidas: $BASE/analise ==="
ls -la *.xvg
