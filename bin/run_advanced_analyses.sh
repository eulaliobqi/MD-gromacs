#!/usr/bin/env bash
# run_advanced_analyses.sh — contact maps, ProLIF e MM-GBSA para todos os complexos.
#
# Uso:
#   conda activate md-gromacs
#   bash bin/run_advanced_analyses.sh [--only contact_map|prolif|mmgbsa] [--dry-run]
#
# Auto-descobre TODOS os complexos concluídos buscando prod/md.tpr.
# Classifica como BEN ou PEPTÍDEO pelo nome do diretório.
# BEN: analisa só a fase ligada (0→dissoc_ns conforme tabela embutida).

set -euo pipefail

ONLY=""
DRY_RUN=0
STRIDE=10
CUTOFF=0.4

while [[ $# -gt 0 ]]; do
    case "$1" in
        --only)    ONLY="$2";   shift 2 ;;
        --dry-run) DRY_RUN=1;   shift   ;;
        --stride)  STRIDE="$2"; shift 2 ;;
        *) echo "Argumento desconhecido: $1"; exit 1 ;;
    esac
done

PROLIF_OK=0
python3 -c "import prolif" 2>/dev/null && PROLIF_OK=1

echo "=== Análises avançadas MD-GROMACS ==="
echo "Modo  : ${ONLY:-'contact_map + prolif + mmgbsa'}"
echo "ProLIF: $([ $PROLIF_OK -eq 1 ] && echo OK || echo 'NÃO INSTALADO — pulando')"
echo "Stride: ${STRIDE} | Cutoff: ${CUTOFF} nm"
[[ "$DRY_RUN" -eq 1 ]] && echo "[DRY-RUN ativo]"
echo ""

# ── Tempo de dissociação para BEN (fase ligada) ────────────────────────────
ben_dissoc() {
    # Retorna tempo de dissociação em ns dado o sample_id
    local sid="${1^^}"   # uppercase para match case-insensitive
    case "$sid" in
        *XP273*BEN*|*BEN*XP273*) echo 80  ;;
        *ACR157*BEN*|*BEN*ACR157*) echo 95  ;;
        *XP352*BEN*|*BEN*XP352*) echo 125 ;;
        *QCL936*BEN*|*BEN*QCL936*) echo 150 ;;
        *) echo ""  ;;   # não-BEN: sem limite
    esac
}

is_ben() { [[ "${1^^}" == *BEN* ]]; }

# ── Funções de análise ─────────────────────────────────────────────────────
do_contact_map() {
    local rdir="$1" sid="$2" start_ns="$3" end_ns="$4"
    local tpr="${rdir}/${sid}/prod/md.tpr"
    [[ -f "$tpr" ]] || { echo "  [SKIP] ${sid}: md.tpr ausente em ${rdir}/${sid}/prod/"; return; }
    local xtc="${rdir}/${sid}/prod/md_fit.xtc"
    local ndx="${rdir}/${sid}/analise/lig.ndx"
    local outdir="${rdir}/${sid}/contact_map"
    local end_arg=()
    [[ -n "$end_ns" ]] && end_arg=(--end-ns "$end_ns")
    echo "  [contact_map] ${sid}  ${start_ns}–${end_ns:-fim} ns  →  ${outdir}"
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
    [[ $PROLIF_OK -eq 1 ]] || return 0
    local tpr="${rdir}/${sid}/prod/md.tpr"
    [[ -f "$tpr" ]] || return 0
    local xtc="${rdir}/${sid}/prod/md_fit.xtc"
    local ndx="${rdir}/${sid}/analise/lig.ndx"
    local outdir="${rdir}/${sid}/prolif"
    local end_arg=()
    [[ -n "$end_ns" ]] && end_arg=(--end-ns "$end_ns")
    echo "  [prolif]      ${sid}  ${start_ns}–${end_ns:-fim} ns"
    if [[ "$DRY_RUN" -eq 0 ]]; then
        python3 bin/run_prolif.py \
            --tpr "$tpr" --xtc "$xtc" --ndx "$ndx" \
            --outdir "$outdir" --sample-id "$sid" \
            --start-ns "$start_ns" "${end_arg[@]}" --stride "$STRIDE"
    fi
}

do_mmgbsa() {
    local rdir="$1" sid="$2" start_ns="$3" end_ns="$4"
    local tpr="${rdir}/${sid}/prod/md.tpr"
    [[ -f "$tpr" ]] || return 0
    local end_arg=()
    [[ -n "$end_ns" ]] && end_arg=(--end-ns "$end_ns")
    echo "  [mmgbsa]      ${sid}  ${start_ns}–${end_ns:-fim} ns"
    if [[ "$DRY_RUN" -eq 0 ]]; then
        bash bin/run_mmgbsa_standalone.sh \
            --results-dir "$rdir" --sample-id "$sid" \
            --start-ns "$start_ns" "${end_arg[@]}"
    fi
}

should_run() { [[ -z "$ONLY" || "$ONLY" == "$1" ]]; }

# ── Auto-descoberta de todos os complexos concluídos ──────────────────────
echo "Buscando complexos concluídos (prod/md.tpr)..."
mapfile -t TPR_LIST < <(find . -name "md.tpr" -path "*/prod/md.tpr" 2>/dev/null | sort)

if [[ ${#TPR_LIST[@]} -eq 0 ]]; then
    echo "ERRO: nenhum prod/md.tpr encontrado no diretório atual."
    echo "Verifique se está na raiz do repositório: $(pwd)"
    exit 1
fi

echo "Encontrados: ${#TPR_LIST[@]} complexos"
echo ""

# ── Processa cada complexo descoberto ─────────────────────────────────────
declare -A SEEN_SIDS   # evita duplicatas se find retornar múltiplos tpr

for tpr_path in "${TPR_LIST[@]}"; do
    # tpr_path: ./results/XP273-BEN/XP273-BEN/prod/md.tpr
    #       ou: ./ACR157-GORE4_NEW/MD/acr157-gore4-c1/prod/md.tpr
    sim_dir=$(dirname "$(dirname "$tpr_path")")   # diretório com prod/ e analise/
    sid=$(basename "$sim_dir")                    # ex: XP273-BEN ou acr157-gore4-c1
    rdir=$(dirname "$sim_dir")                    # ex: ./results/XP273-BEN ou ./ACR157-GORE4_NEW/MD
    rdir="${rdir#./}"                             # remove ./ inicial

    # Verifica ndx
    ndx="${sim_dir#./}/analise/lig.ndx"
    ndx="${sim_dir}/analise/lig.ndx"
    if [[ ! -f "$ndx" ]]; then
        echo "  [SKIP] ${sid}: lig.ndx ausente (analise/ ainda em andamento?)"
        continue
    fi

    # Evita duplicatas
    [[ -n "${SEEN_SIDS[$sid]+x}" ]] && continue
    SEEN_SIDS["$sid"]=1

    # Determina parâmetros por série
    if is_ben "$sid"; then
        dissoc=$(ben_dissoc "$sid")
        end_ns="${dissoc}"
        mmgbsa_start=5
        serie="BEN"
    else
        end_ns=""
        mmgbsa_start=50
        serie="PEPTÍDEO"
    fi

    echo "── ${sid}  [${serie}]  rdir=${rdir}"
    should_run contact_map && do_contact_map "$rdir" "$sid" 0 "$end_ns"
    should_run prolif       && do_prolif      "$rdir" "$sid" 0 "$end_ns"
    should_run mmgbsa       && do_mmgbsa      "$rdir" "$sid" "$mmgbsa_start" "$end_ns"
    echo ""
done

# ── Resumo de saídas ──────────────────────────────────────────────────────
echo "=== Resumo de saídas ==="
printf "  %-30s  %-12s  %-8s  %-8s\n" "Sistema" "contact_map" "prolif" "mmgbsa"
printf "  %-30s  %-12s  %-8s  %-8s\n" "-------" "-----------" "------" "------"

for tpr_path in "${TPR_LIST[@]}"; do
    sim_dir=$(dirname "$(dirname "$tpr_path")")
    sid=$(basename "$sim_dir")
    cm="${sim_dir}/contact_map/contact_map.png"
    pr="${sim_dir}/prolif/prolif_persistence.png"
    mg="${sim_dir}/mmgbsa/FINAL_RESULTS_MMGBSA.dat"
    printf "  %-30s  %-12s  %-8s  %-8s\n" \
        "$sid" \
        "$([ -f "$cm" ] && echo OK || echo '---')" \
        "$([ -f "$pr" ] && echo OK || echo '---')" \
        "$([ -f "$mg" ] && echo OK || echo '---')"
done

echo ""
echo "[OK] run_advanced_analyses.sh concluído"
