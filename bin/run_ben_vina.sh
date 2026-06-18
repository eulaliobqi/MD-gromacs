#!/usr/bin/env bash
# Roda AutoDock Vina de BEN contra os 4 receptores Spodoptera.
# Pré-requisito: prepare_ben_docking.py já executado para todos os receptores.
# Uso: bash bin/run_ben_vina.sh [--outdir spodoptera-ben]

set -euo pipefail

OUTDIR="${1:-spodoptera-ben}"
CONF_DIR="${OUTDIR}/configs"
RES_DIR="${OUTDIR}/results"

RECEPTORS=(ACR157 QCL936 XP273 XP352)

# Verifica vina
if ! command -v vina &>/dev/null; then
    echo "ERRO: 'vina' não encontrado. Instale com:"
    echo "  mamba install -c conda-forge autodock-vina -n md-gromacs -y"
    exit 1
fi

echo "=== Docking BEN × Tripsinas Spodoptera ==="
echo "Vina: $(vina --version 2>&1 | head -1)"
echo ""

mkdir -p "${RES_DIR}"

for RECEPTOR in "${RECEPTORS[@]}"; do
    CONF="${CONF_DIR}/vina_${RECEPTOR}.conf"
    if [[ ! -f "${CONF}" ]]; then
        echo "AVISO: config não encontrado: ${CONF} — pulando ${RECEPTOR}"
        continue
    fi

    OUT="${RES_DIR}/${RECEPTOR}_ben_out.pdbqt"
    LOG="${RES_DIR}/${RECEPTOR}_ben_vina.log"

    echo "[$(date '+%H:%M:%S')] Iniciando ${RECEPTOR}..."
    vina --config "${CONF}" 2>&1 | tee "${LOG}"

    if [[ -f "${OUT}" ]]; then
        BEST=$(grep "REMARK VINA RESULT" "${OUT}" | head -1 | awk '{print $4}')
        echo "[$(date '+%H:%M:%S')] ${RECEPTOR} concluído — melhor score: ${BEST} kcal/mol"
    else
        echo "AVISO: arquivo de saída não gerado para ${RECEPTOR}"
    fi
    echo ""
done

echo "=== Docking concluído. Resultados em ${RES_DIR}/ ==="
echo "Próximo: python3 bin/analyze_ben_docking.py --ben-dir ${OUTDIR}"
