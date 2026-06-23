#!/usr/bin/env bash
# run_advanced_analyses.sh — orquestra contact maps, ProLIF e MM-GBSA
# para todos os complexos já simulados.
#
# Uso (a partir da raiz do repositório no servidor, com md-gromacs ativo):
#   conda activate md-gromacs
#   bash bin/run_advanced_analyses.sh [--only contact_map|prolif|mmgbsa] [--dry-run]
#
# Sistemas cobertos:
#   GORE4 series : QCL936-GORE4, ACR157-GORE4, XP273-GORE4
#   SKTI series  : ACR157-SKTI, QCL936-SKTI, XP273-SKTI
#   BEN series   : ACR157-BEN, QCL936-BEN, XP273-BEN, XP352-BEN
#     (BEN: apenas fase ligada — end-ns conforme dissociação observada)
#
# Requer:
#   mamba install -n md-gromacs mdanalysis seaborn prolif

set -euo pipefail

ONLY=""
DRY_RUN=0
STRIDE=10
CUTOFF=0.4
RESULTS_BASE=""   # vazio = auto-detectar

while [[ $# -gt 0 ]]; do
    case "$1" in
        --only)          ONLY="$2";         shift 2 ;;
        --dry-run)       DRY_RUN=1;         shift   ;;
        --stride)        STRIDE="$2";       shift 2 ;;
        --results-base)  RESULTS_BASE="$2"; shift 2 ;;
        *) echo "Argumento desconhecido: $1"; exit 1 ;;
    esac
done

# Auto-detecta o diretório base de resultados
if [[ -z "$RESULTS_BASE" ]]; then
    for candidate in results results-MD work/results; do
        if [[ -d "$candidate" ]]; then
            RESULTS_BASE="$candidate"
            break
        fi
    done
    if [[ -z "$RESULTS_BASE" ]]; then
        echo "ERRO: nenhum diretório de resultados encontrado (results/, results-MD/)."
        echo "Use: bash bin/run_advanced_analyses.sh --results-base <caminho>"
        exit 1
    fi
fi
echo "Diretório de resultados: ${RESULTS_BASE}/"

PROLIF_OK=0
python3 -c "import prolif" 2>/dev/null && PROLIF_OK=1

echo "=== Análises avançadas MD-GROMACS ==="
echo "Modo  : ${ONLY:-'contact_map + prolif + mmgbsa'}"
echo "ProLIF: $([ $PROLIF_OK -eq 1 ] && echo OK || echo 'NÃO INSTALADO — pulando')"
echo "Stride: ${STRIDE} | Cutoff: ${CUTOFF} nm"
[[ "$DRY_RUN" -eq 1 ]] && echo "[DRY-RUN ativo]"
echo ""

# ────────────────────────────────────────────────────────────────────────────
do_contact_map() {
    local rdir="$1" sid="$2" start_ns="$3" end_ns="$4"
    local tpr="${rdir}/${sid}/prod/md.tpr"
    [[ -f "$tpr" ]] || { echo "  [SKIP] ${sid}: md.tpr ausente"; return; }
    local outdir="${rdir}/${sid}/contact_map"
    local xtc="${rdir}/${sid}/prod/md_fit.xtc"
    local ndx="${rdir}/${sid}/analise/lig.ndx"
    local end_arg=()
    [[ -n "$end_ns" ]] && end_arg=(--end-ns "$end_ns")
    echo "  [contact_map] ${sid}  ${start_ns}–${end_ns:-fim} ns"
    if [[ "$DRY_RUN" -eq 0 ]]; then
        python3 bin/contact_map.py \
            --tpr "$tpr" --xtc "$xtc" --ndx "$ndx" \
            --outdir "$outdir" --sample-id "$sid" \
            --cutoff "$CUTOFF" --start-ns "$start_ns" \
            "${end_arg[@]}" --stride "$STRIDE"
    fi
}

do_prolif() {
    local rdir="$1" sid="$2" start_ns="$3" end_ns="$4"
    [[ $PROLIF_OK -eq 1 ]] || { echo "  [SKIP] ${sid}: ProLIF ausente"; return; }
    local tpr="${rdir}/${sid}/prod/md.tpr"
    [[ -f "$tpr" ]] || { echo "  [SKIP] ${sid}: md.tpr ausente"; return; }
    local outdir="${rdir}/${sid}/prolif"
    local xtc="${rdir}/${sid}/prod/md_fit.xtc"
    local ndx="${rdir}/${sid}/analise/lig.ndx"
    local end_arg=()
    [[ -n "$end_ns" ]] && end_arg=(--end-ns "$end_ns")
    echo "  [prolif]      ${sid}  ${start_ns}–${end_ns:-fim} ns"
    if [[ "$DRY_RUN" -eq 0 ]]; then
        python3 bin/run_prolif.py \
            --tpr "$tpr" --xtc "$xtc" --ndx "$ndx" \
            --outdir "$outdir" --sample-id "$sid" \
            --start-ns "$start_ns" \
            "${end_arg[@]}" --stride "$STRIDE"
    fi
}

do_mmgbsa() {
    local rdir="$1" sid="$2" start_ns="$3" end_ns="$4"
    local tpr="${rdir}/${sid}/prod/md.tpr"
    [[ -f "$tpr" ]] || { echo "  [SKIP] ${sid}: md.tpr ausente"; return; }
    local end_arg=()
    [[ -n "$end_ns" ]] && end_arg=(--end-ns "$end_ns")
    echo "  [mmgbsa]      ${sid}  ${start_ns}–${end_ns:-fim} ns"
    if [[ "$DRY_RUN" -eq 0 ]]; then
        bash bin/run_mmgbsa_standalone.sh \
            --results-dir "$rdir" --sample-id "$sid" \
            --start-ns "$start_ns" \
            "${end_arg[@]}"
    fi
}

should_run() { [[ -z "$ONLY" || "$ONLY" == "$1" ]]; }

# ════════════════════════════════════════════════════════════════════════════
# SÉRIE GORE4 — contact map + ProLIF: 0→100 ns | MM-GBSA: 50→100 ns
# ════════════════════════════════════════════════════════════════════════════
echo "── SÉRIE GORE4 ─────────────────────────────────────────────────────────"
for sid in "QCL936-GORE4" "ACR157-GORE4" "XP273-GORE4"; do
    rdir="${RESULTS_BASE}/${sid}"
    [[ -d "$rdir" ]] || { echo "  [SKIP] $rdir não encontrado"; continue; }
    should_run contact_map && do_contact_map "$rdir" "$sid" 0 ""
    should_run prolif       && do_prolif      "$rdir" "$sid" 0 ""
    should_run mmgbsa       && do_mmgbsa      "$rdir" "$sid" 50 ""
done

# ════════════════════════════════════════════════════════════════════════════
# SÉRIE SKTI — mesma lógica (proteína grande: ProLIF válido mas mais lento)
# ════════════════════════════════════════════════════════════════════════════
echo ""
echo "── SÉRIE SKTI ─────────────────────────────────────────────────────────"
for sid in "ACR157-SKTI" "QCL936-SKTI" "XP273-SKTI"; do
    rdir="${RESULTS_BASE}/${sid}"
    [[ -d "$rdir" ]] || { echo "  [SKIP] $rdir não encontrado"; continue; }
    should_run contact_map && do_contact_map "$rdir" "$sid" 0 ""
    should_run prolif       && do_prolif      "$rdir" "$sid" 0 ""
    should_run mmgbsa       && do_mmgbsa      "$rdir" "$sid" 50 ""
done

# ════════════════════════════════════════════════════════════════════════════
# SÉRIE BEN — fase ligada apenas (0 → dissociação observada)
# S1=Ile (neutro):  XP273-BEN ~80 ns | ACR157-BEN ~95 ns
# S1=Asp (−1):      XP352-BEN ~125 ns | QCL936-BEN ~150 ns
# ════════════════════════════════════════════════════════════════════════════
echo ""
echo "── SÉRIE BEN (fase ligada) ─────────────────────────────────────────────"
declare -A BEN_DISSOC=(
    ["XP273-BEN"]=80
    ["ACR157-BEN"]=95
    ["XP352-BEN"]=125
    ["QCL936-BEN"]=150
)

for sid in "XP273-BEN" "ACR157-BEN" "XP352-BEN" "QCL936-BEN"; do
    rdir="${RESULTS_BASE}/${sid}"
    [[ -d "$rdir" ]] || { echo "  [SKIP] $rdir não encontrado"; continue; }
    dissoc="${BEN_DISSOC[$sid]}"
    should_run contact_map && do_contact_map "$rdir" "$sid" 0 "$dissoc"
    should_run prolif       && do_prolif      "$rdir" "$sid" 0 "$dissoc"
    should_run mmgbsa       && do_mmgbsa      "$rdir" "$sid" 5 "$dissoc"
done

# ── Resumo de saídas ──────────────────────────────────────────────────────────
echo ""
echo "=== Resumo de saídas ==="
printf "  %-22s  %-12s  %-8s  %-8s\n" "Sistema" "contact_map" "prolif" "mmgbsa"
printf "  %-22s  %-12s  %-8s  %-8s\n" "-------" "-----------" "------" "------"

for sid in "QCL936-GORE4" "ACR157-GORE4" "XP273-GORE4" \
           "ACR157-SKTI" "QCL936-SKTI" "XP273-SKTI" \
           "XP273-BEN" "ACR157-BEN" "XP352-BEN" "QCL936-BEN"; do
    rdir="${RESULTS_BASE}/${sid}"
    [[ -d "$rdir" ]] || continue
    cm="$rdir/$sid/contact_map/contact_map.png"
    pr="$rdir/$sid/prolif/prolif_persistence.png"
    mg="$rdir/$sid/mmgbsa/FINAL_RESULTS_MMGBSA.dat"
    printf "  %-22s  %-12s  %-8s  %-8s\n" \
        "$sid" \
        "$([ -f "$cm" ] && echo OK || echo '---')" \
        "$([ -f "$pr" ] && echo OK || echo '---')" \
        "$([ -f "$mg" ] && echo OK || echo '---')"
done

echo ""
echo "[OK] run_advanced_analyses.sh concluído"
