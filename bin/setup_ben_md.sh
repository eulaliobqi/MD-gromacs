#!/bin/bash
# Prepara o ambiente para DM da benzamidina (BEN):
#   1. Extrai melhor pose Vina (modo 1) de cada receptor → spodoptera-ben/ben_poses/
#   2. Gera 4 samplesheets em assets/ com paths absolutos
#
# Executar a partir de ~/gromacs/MD-gromacs
# Requer: conda activate protein_design_env  (tem obabel)
#          python3 bin/extract_vina_pose.py
set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BEN_DIR="${REPO_DIR}/spodoptera-ben"
POSES_DIR="${BEN_DIR}/ben_poses"
ASSETS_DIR="${REPO_DIR}/assets"
RESULTS_DIR="${BEN_DIR}/results"

mkdir -p "${POSES_DIR}" "${ASSETS_DIR}"

declare -A TRIAD=(
    [ACR157]="69,114,211,205"
    [QCL936]="92,142,247,241"
    [XP273]="83,132,234,229"
    [XP352]="112,166,268,262"
)

RECEPTORS=(ACR157 QCL936 XP273 XP352)

echo "=== SETUP BEN MD PIPELINE ==="
echo "Repo: ${REPO_DIR}"
echo ""

# 1. Verificar que resultados do docking existem
for REC in "${RECEPTORS[@]}"; do
    PDBQT="${RESULTS_DIR}/${REC}_ben_out.pdbqt"
    if [ ! -f "${PDBQT}" ]; then
        echo "ERRO: ${PDBQT} não encontrado — rode bin/run_ben_vina.sh primeiro"
        exit 1
    fi
done
echo "[OK] Arquivos Vina encontrados em ${RESULTS_DIR}/"

# 2. Extrair melhor pose (modo 1) de cada receptor
echo ""
echo "--- Extraindo poses Vina (modo 1) ---"
for REC in "${RECEPTORS[@]}"; do
    PDBQT="${RESULTS_DIR}/${REC}_ben_out.pdbqt"
    OUT_PDB="${POSES_DIR}/${REC}_ben_pose1.pdb"
    python3 "${REPO_DIR}/bin/extract_vina_pose.py" \
        --pdbqt "${PDBQT}" \
        --out   "${OUT_PDB}" \
        --mode  1 \
        --resname BEN
done

# 3. Verificar receptores finais
echo ""
echo "--- Verificando receptores pós-MD ---"
for REC in "${RECEPTORS[@]}"; do
    RECPDB="${BEN_DIR}/${REC}-final.pdb"
    if [ ! -f "${RECPDB}" ]; then
        echo "ERRO: ${RECPDB} não encontrado — copie os receptores para spodoptera-ben/"
        exit 1
    fi
    NATOM=$(grep -c "^ATOM" "${RECPDB}" || echo 0)
    echo "  ${REC}-final.pdb: ${NATOM} átomos ATOM"
done

# 4. Gerar samplesheets
echo ""
echo "--- Gerando samplesheets ---"
for REC in "${RECEPTORS[@]}"; do
    IFS=',' read -r T1 T2 T3 T4 <<< "${TRIAD[$REC]}"
    SS="${ASSETS_DIR}/samplesheet_${REC,,}_ben.csv"
    cat > "${SS}" << CSVEOF
sample_id,receptor,ligand,triad_1,triad_2,triad_3,triad_4
${REC}-BEN,${BEN_DIR}/${REC}-final.pdb,${POSES_DIR}/${REC}_ben_pose1.pdb,${T1},${T2},${T3},${T4}
CSVEOF
    echo "  Criado: ${SS}"
done

# 5. Verificar ACPYPE
echo ""
echo "--- Verificando ACPYPE ---"
if command -v acpype &>/dev/null; then
    ACPYPE_VER=$(acpype --version 2>&1 | head -1 || echo "desconhecida")
    echo "[OK] ACPYPE disponível: ${ACPYPE_VER}"
elif conda run -n md-gromacs acpype --version &>/dev/null 2>&1; then
    echo "[OK] ACPYPE disponível no ambiente md-gromacs"
else
    echo "[AVISO] ACPYPE não encontrado — instalar com:"
    echo "         conda run -n md-gromacs pip install acpype"
fi

# 6. Mostrar comandos de execução
echo ""
echo "=== COMANDOS PARA EXECUTAR ==="
echo ""
echo "# Executar em screen (nunca em background direto):"
echo "screen -S ben_md"
echo ""
echo "# Entrar no ambiente e ir ao diretório:"
echo "conda activate md-gromacs"
echo "cd ${REPO_DIR}"
echo "git pull"
echo ""
for REC in "${RECEPTORS[@]}"; do
    SS="${ASSETS_DIR}/samplesheet_${REC,,}_ben.csv"
    OUTDIR="${REPO_DIR}/results/${REC}-BEN"
    echo "# ${REC}-BEN:"
    echo "nextflow run main_ben.nf --input ${SS} --outdir ${OUTDIR} -profile server"
    echo ""
done
echo "# Para rodar todos de uma vez (sequencial):"
echo "for SS in ${ASSETS_DIR}/samplesheet_*_ben.csv; do"
echo "    ID=\$(basename \$SS .csv | sed 's/samplesheet_//')"
echo "    nextflow run main_ben.nf --input \$SS --outdir ${REPO_DIR}/results/\${ID^^}-BEN -profile server"
echo "done"
echo ""
echo "=== SETUP CONCLUÍDO ==="
