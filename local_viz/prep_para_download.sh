#!/usr/bin/env bash
# Roda no SERVIDOR: prepara pacote leve com trajetoria reduzida + scripts PyMOL
# Uso: bash prep_para_download.sh <BASE_DIR> [STRIDE]
#   STRIDE: pular N frames (default 50 = ~200 frames de 100 ns)
set -euo pipefail

BASE="${1:?Use: bash prep_para_download.sh <BASE_DIR> [STRIDE]}"
STRIDE="${2:-50}"
GMX=$(command -v gmx_mpi 2>/dev/null || command -v gmx 2>/dev/null)
[[ -z "$GMX" ]] && { echo "ERRO: gmx nao encontrado"; exit 1; }

PKG="$BASE/visualizacao_local"
mkdir -p "$PKG"
cd "$PKG"

echo "[1/4] Gerando trajetoria reduzida (skip $STRIDE)..."
echo 0 | "$GMX" trjconv -s "$BASE/prod/md.tpr" -f "$BASE/prod/md_fit.xtc" \
    -o trajetoria.xtc -skip "$STRIDE"

echo "[2/4] Extraindo estrutura inicial (frame 0)..."
echo 0 | "$GMX" trjconv -s "$BASE/prod/md.tpr" -f "$BASE/prod/md_fit.xtc" \
    -o estrutura.pdb -dump 0

echo "[3/4] Extraindo estrutura final..."
echo 0 | "$GMX" trjconv -s "$BASE/prod/md.tpr" -f "$BASE/prod/md_fit.xtc" \
    -o estrutura_final.pdb -dump 100000

echo "[4/4] Copiando scripts PyMOL..."
cp "$(dirname "$0")"/{visualizar.pml,gerar_filme.pml,README_LOCAL.md} . 2>/dev/null || true

echo ""
echo "==============================================="
echo "  Pacote pronto em: $PKG"
echo "==============================================="
ls -lh "$PKG"
echo ""
TOT=$(du -sh "$PKG" | awk '{print $1}')
echo "  Tamanho total: $TOT"
echo ""
echo "  Para baixar para o PC local:"
echo "    scp -r usuario@servidor:$PKG ~/Downloads/"
echo ""
echo "  Ou use rsync/SFTP/etc"
