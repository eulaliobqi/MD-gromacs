#!/bin/bash
# Executa gmx_MMPBSA manualmente fora do Nextflow, com fix automático do bug
# de índices SS bonds no tleap (gmx_MMPBSA 1.6.5).
#
# Uso:
#   bash scripts/run_mmgbsa.sh <work_dir> <sys_name> <time_ns> <saltcon>
#
# Exemplo:
#   bash scripts/run_mmgbsa.sh ~/gromacs/MD-gromacs/work/bc/05e6b7* \
#       ACR157-GORE4 100 0.15

set -euo pipefail

NF_WORK="${1:?Informe o diretório de trabalho Nextflow do MMGBSA}"
SYS_NAME="${2:?Informe o nome do sistema (ex: ACR157-GORE4)}"
TIME_NS="${3:-100}"
SALTCON="${4:-0.15}"

TOTAL_FRAMES=$(( TIME_NS * 100 ))
INTERVAL=$(python3 -c "print(max(1, $TOTAL_FRAMES // 200))")
OUTDIR="$(dirname "$NF_WORK")/../mmgbsa_${SYS_NAME}"

# ── Preparar diretório de saída ──────────────────────────────────────────────
mkdir -p "$OUTDIR"
cd "$OUTDIR"

echo "=== Diretório: $OUTDIR ==="
echo "=== Sistema: $SYS_NAME | ${TIME_NS} ns | saltcon=${SALTCON} ==="

# ── Vincular arquivos de entrada ─────────────────────────────────────────────
ln -sf "$(realpath "$NF_WORK"/md.tpr)" md.tpr
ln -sf "$(realpath "$NF_WORK"/md_fit.xtc)" md_fit.xtc
ln -sf "$(realpath "$NF_WORK"/lig.ndx)" lig.ndx

echo "Arquivos vinculados: md.tpr  md_fit.xtc  lig.ndx"

# ── mmgbsa.in ────────────────────────────────────────────────────────────────
cat > mmgbsa.in << MEOF
&general
sys_name="${SYS_NAME}",
startframe=1,
endframe=${TOTAL_FRAMES},
interval=${INTERVAL},
verbose=0,
/
&gb
igb=2,
saltcon=${SALTCON},
/
MEOF

# ── Wrapper tleap: corrige bug gmx_MMPBSA (índices SS bonds COM_OUT errados) ─
mkdir -p bin_patch
cat > bin_patch/tleap << 'WEOF'
#!/usr/bin/env python3
import sys, os, re, subprocess
args = sys.argv[1:]
for i, a in enumerate(args):
    if a == '-f' and i+1 < len(args):
        f = args[i+1]
        if os.path.exists(f):
            t = open(f).read()
            rec = re.findall(r'bond REC_OUT[.](\d+)[.]SG REC_OUT[.](\d+)[.]SG', t)
            com = re.findall(r'bond COM_OUT[.](\d+)[.]SG COM_OUT[.](\d+)[.]SG', t)
            if rec and com:
                for w, r in zip(com, rec):
                    t = t.replace(
                        'bond COM_OUT.{}.SG COM_OUT.{}.SG'.format(w[0], w[1]),
                        'bond COM_OUT.{}.SG COM_OUT.{}.SG'.format(r[0], r[1]), 1)
                open(f, 'w').write(t)
                print('[tleap-wrapper] SS bond indices corrigidos:', list(zip(com, rec)))
        break
for d in os.environ.get('PATH', '').split(':'):
    if 'bin_patch' in d:
        continue
    for name in ('tleap', 'teLeap'):
        c = os.path.join(d, name)
        if os.path.isfile(c) and os.access(c, os.X_OK):
            sys.exit(subprocess.run([c] + args).returncode)
sys.exit(1)
WEOF
chmod +x bin_patch/tleap

# ── Executar gmx_MMPBSA ───────────────────────────────────────────────────────
echo "=== Iniciando gmx_MMPBSA (pode demorar 30-60 min) ==="
CURRENT_DIR="$PWD"
mamba run -n mmgbsa-env bash -c "
export PATH=${CURRENT_DIR}/bin_patch:\$PATH
gmx_MMPBSA -O \
    -i mmgbsa.in \
    -cs md.tpr \
    -ct md_fit.xtc \
    -ci lig.ndx \
    -cg Receptor Ligante \
    -o FINAL_RESULTS_MMGBSA.dat \
    -eo mmgbsa_results.csv \
    -nogui
" 2>&1 | tee mmgbsa.log

# ── Verificar saída ───────────────────────────────────────────────────────────
if [ -f FINAL_RESULTS_MMGBSA.dat ]; then
    echo "=== gmx_MMPBSA concluído com sucesso ==="
    echo "Gerando painel de resultados..."
    mamba run -n md-gromacs plot_results.py \
        --analise-dir . \
        --titulo "${SYS_NAME} — DM ${TIME_NS} ns [MM-GBSA]" \
        --mmgbsa-csv mmgbsa_results.csv \
        --output painel_mmgbsa.png
    echo "=== Pronto! Resultados em: $OUTDIR ==="
else
    echo "=== ERRO: gmx_MMPBSA falhou. Veja: $OUTDIR/mmgbsa.log ==="
    tail -20 mmgbsa.log
    exit 1
fi
