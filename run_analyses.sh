#!/usr/bin/env bash
# Pos-processamento (PBC + fit) + analises completas
# Uso: bash run_analyses.sh <BASE_DIR>
set -euo pipefail

BASE="${1:-${HOME}/gromacs/dm_complex}"
[[ ! -d "$BASE/prod" ]] && { echo "ERRO: $BASE/prod nao existe"; exit 1; }
[[ ! -s "$BASE/prod/md.tpr" ]] && { echo "ERRO: $BASE/prod/md.tpr nao existe"; exit 1; }
[[ ! -s "$BASE/prod/md.xtc" ]] && { echo "ERRO: $BASE/prod/md.xtc nao existe"; exit 1; }

GMX=$(command -v gmx_mpi 2>/dev/null || command -v gmx 2>/dev/null)
[[ -z "$GMX" ]] && { echo "ERRO: gmx nao encontrado"; exit 1; }

# ============================================================
# ETAPA 8 - Pos-processamento
# ============================================================
cd "$BASE/prod"
if [[ ! -s md_fit.xtc ]]; then
    [[ ! -s md_nojump.xtc ]] && \
        echo 0 | "$GMX" trjconv -s md.tpr -f md.xtc -o md_nojump.xtc -pbc nojump
    [[ ! -s md_nopbc.xtc ]] && \
        echo -e '1\n0' | "$GMX" trjconv -s md.tpr -f md_nojump.xtc \
            -o md_nopbc.xtc -pbc mol -center -ur compact
    echo -e '4\n0' | "$GMX" trjconv -s md.tpr -f md_nopbc.xtc \
        -o md_fit.xtc -fit rot+trans
fi

# ============================================================
# ETAPA 9 - Analises
# ============================================================
mkdir -p "$BASE/analise"
cd "$BASE/analise"
ln -sf "$BASE/prod/md.tpr"     sistema.tpr
ln -sf "$BASE/prod/md_fit.xtc" trajetoria.xtc

# Detectar residuos do ligante (chain B) via PDB
PDB_COMPLEX="$BASE/prep/complexo.pdb"
[[ ! -s "$PDB_COMPLEX" ]] && { echo "ERRO: $PDB_COMPLEX nao existe"; exit 1; }

LIG_FIRST=$(awk '/^ATOM/ && substr($0,22,1)=="B" {print substr($0,23,4)+0; exit}' "$PDB_COMPLEX")
LIG_LAST=$(awk '/^ATOM/ && substr($0,22,1)=="B" {r=substr($0,23,4)+0} END{print r}' "$PDB_COMPLEX")

[[ -z "$LIG_FIRST" || -z "$LIG_LAST" ]] && {
    echo "ERRO: nao detectei chain B"; exit 1; }
[[ "$LIG_FIRST" -gt "$LIG_LAST" ]] && {
    echo "ERRO: faixa invertida $LIG_FIRST > $LIG_LAST"; exit 1; }
echo "[ANALISE] Ligante: residuos $LIG_FIRST a $LIG_LAST"

# === DETECTAR N_DEFAULT (numero de grupos default no .tpr) ===
N_DEFAULT=$(echo q | "$GMX" make_ndx -f sistema.tpr -o /tmp/_default.ndx 2>&1 | \
            grep -cE "^ *[0-9]+ +[A-Za-z]")
LIG_IDX=$N_DEFAULT
REC_IDX=$((N_DEFAULT + 1))
echo "[INDEX] Grupos default: $N_DEFAULT (Ligante->idx $LIG_IDX, Receptor->idx $REC_IDX)"
rm -f /tmp/_default.ndx

# Cria lig.ndx do zero
rm -f lig.ndx
"$GMX" make_ndx -f sistema.tpr -o lig.ndx <<EOF >/dev/null 2>&1
r ${LIG_FIRST}-${LIG_LAST}
name ${LIG_IDX} Ligante
1 & ! ${LIG_IDX}
name ${REC_IDX} Receptor
q
EOF

# Verifica grupos
echo "[INDEX] Grupos finais no lig.ndx:"
echo q | "$GMX" make_ndx -f sistema.tpr -n lig.ndx 2>&1 | \
    grep -E "Ligante|Receptor" | head -3

# RMSD backbone
[[ ! -s rmsd_backbone.xvg ]] && \
    printf 'Backbone\nBackbone\n' | "$GMX" rms -s sistema.tpr -f trajetoria.xtc \
        -n lig.ndx -o rmsd_backbone.xvg -tu ns
echo "[OK] rmsd_backbone.xvg"

# RMSD do ligante
[[ ! -s rmsd_ligante.xvg ]] && \
    printf 'Ligante\nLigante\n' | "$GMX" rms -s sistema.tpr -f trajetoria.xtc \
        -n lig.ndx -o rmsd_ligante.xvg -tu ns
echo "[OK] rmsd_ligante.xvg"

# RMSF
[[ ! -s rmsf_residuos.xvg ]] && \
    printf 'Backbone\n' | "$GMX" rmsf -s sistema.tpr -f trajetoria.xtc \
        -n lig.ndx -o rmsf_residuos.xvg -res -fit
echo "[OK] rmsf_residuos.xvg"

# Raio de giro
[[ ! -s gyrate.xvg ]] && \
    printf 'Protein\n' | "$GMX" gyrate -s sistema.tpr -f trajetoria.xtc \
        -n lig.ndx -o gyrate.xvg -tu ns
echo "[OK] gyrate.xvg"

# Contatos Receptor vs Ligante <0.4 nm
rm -f numcont.xvg mindist.xvg
printf 'Receptor\nLigante\n' | "$GMX" mindist -s sistema.tpr -f trajetoria.xtc \
    -n lig.ndx -od mindist.xvg -on numcont.xvg -d 0.4 -tu ns 2>&1 | tail -3
echo "[OK] numcont.xvg"

# H-bonds
rm -f hbond.xvg
printf 'Receptor\nLigante\n' | "$GMX" hbond -s sistema.tpr -f trajetoria.xtc \
    -n lig.ndx -num hbond.xvg -tu ns 2>&1 | tail -3
echo "[OK] hbond.xvg"

echo ""
echo "=== Arquivos em $BASE/analise: ==="
ls -la *.xvg