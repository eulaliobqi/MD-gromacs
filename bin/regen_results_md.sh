#!/usr/bin/env bash
# regen_results_md.sh — Regenerates all PNGs in results-MD/ using updated plot_results.py
# Run from the repo root on the server after pulling the updated scripts.
#
# Usage:
#   bash bin/regen_results_md.sh [--results-dir <dir>] [--window-ns N]
#
# Requires: conda/mamba env md-gromacs with matplotlib installed
#   mamba activate md-gromacs

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(dirname "$SCRIPT_DIR")"
RESULTS_ROOT="${REPO_ROOT}/results"
OUTROOT="${REPO_ROOT}/results-MD"
WINDOW_NS=5
PLOT_SCRIPT="${SCRIPT_DIR}/plot_results.py"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --results-dir) RESULTS_ROOT="$2"; shift 2 ;;
        --window-ns)   WINDOW_NS="$2";    shift 2 ;;
        *) echo "Unknown argument: $1"; exit 1 ;;
    esac
done

echo "=== regen_results_md.sh ==="
echo "  results dir : $RESULTS_ROOT"
echo "  output root : $OUTROOT"
echo "  window-ns   : $WINDOW_NS"
echo

# Mapping: results/{ID}/analise/ → results-MD/{TRYPSIN}/{PEPTIDE}/
# Supported IDs: ACR157-GORE4, QCL936-GORE4, XP273-GORE4, XP352-GORE4
#                ACR157-SKTI,  QCL936-SKTI,  XP273-SKTI
#                ACR157-BEN,   QCL936-BEN,   XP273-BEN,  XP352-BEN

declare -A ID_TO_OUTDIR=(
    ["ACR157-GORE4"]="GORE4/ACR157"
    ["QCL936-GORE4"]="GORE4/QCL936"
    ["XP273-GORE4"]="GORE4/XP273"
    ["XP352-GORE4"]="GORE4/XP352"
    ["ACR157-SKTI"]="SKTI/ACR157"
    ["QCL936-SKTI"]="SKTI/QCL936"
    ["XP273-SKTI"]="SKTI/XP273"
    ["XP352-SKTI"]="SKTI/XP352"
    ["ACR157-BEN"]="BENZ/ACR157"
    ["QCL936-BEN"]="BENZ/QCL936"
    ["XP273-BEN"]="BENZ/XP273"
    ["XP352-BEN"]="BENZ/XP352"
)

ok=0
skipped=0

for sample_id in "${!ID_TO_OUTDIR[@]}"; do
    outdir_rel="${ID_TO_OUTDIR[$sample_id]}"

    # Try exact match first, then any subdirectory with analise/
    analise_dir=""
    if [[ -d "${RESULTS_ROOT}/${sample_id}/analise" ]]; then
        analise_dir="${RESULTS_ROOT}/${sample_id}/analise"
    elif [[ -d "${RESULTS_ROOT}/${sample_id}/${sample_id}/analise" ]]; then
        analise_dir="${RESULTS_ROOT}/${sample_id}/${sample_id}/analise"
    else
        # Glob for any matching nested path
        for d in "${RESULTS_ROOT}/${sample_id}"*/analise "${RESULTS_ROOT}"/*"${sample_id}"*/analise; do
            if [[ -d "$d" ]]; then
                analise_dir="$d"
                break
            fi
        done
    fi

    if [[ -z "$analise_dir" ]]; then
        echo "[SKIP] $sample_id — analise/ not found under $RESULTS_ROOT/$sample_id"
        ((skipped++)) || true
        continue
    fi

    outdir="${OUTROOT}/${outdir_rel}"
    mkdir -p "$outdir"

    titulo="${sample_id} - MD 100 ns"

    echo "[RUN ] $sample_id"
    echo "       analise: $analise_dir"
    echo "       output : $outdir"

    python3 "$PLOT_SCRIPT" \
        --analise-dir "$analise_dir" \
        --titulo "$titulo" \
        --window-ns "$WINDOW_NS" \
        --output "${outdir}/painel_completo.png"

    # Copy individual PNGs from analise_dir to outdir
    for png in rmsd_bb rmsd_lig rmsf rg ncont hbond sasa_ligante sasa_protein triad_distances; do
        src="${analise_dir}/${png}.png"
        if [[ -f "$src" ]]; then
            cp "$src" "${outdir}/${png}.png"
        fi
    done

    echo "       [OK]"
    ((ok++)) || true
done

echo
echo "=== Done: $ok regenerated, $skipped skipped ==="
echo "Commit the updated PNGs:"
echo "  git add results-MD/ && git commit -m 'chore(figures): regenerate all PNGs in English'"
