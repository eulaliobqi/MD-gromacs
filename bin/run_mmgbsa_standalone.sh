#!/usr/bin/env bash
# run_mmgbsa_standalone.sh — executa MM-GBSA para um sistema já simulado,
# reproduzindo a lógica do módulo MMGBSA_ROBUST sem Nextflow.
#
# Uso (a partir da raiz do repositório no servidor):
#   bash bin/run_mmgbsa_standalone.sh \
#       --results-dir results/ACR157-GORE4 \
#       --sample-id   ACR157-GORE4 \
#       [--start-ns   50] \
#       [--end-ns     100] \
#       [--saltcon    0.15]
#
# Saídas em: {results-dir}/{sample-id}/mmgbsa/
#   FINAL_RESULTS_MMGBSA.dat
#   mmgbsa_results.csv
#   decomp_results.csv
#   mmgbsa.log
#
# Requisitos:
#   conda env: mmgbsa-env  (gmx_MMPBSA, python=3.11, tleap)
#   conda env: md-gromacs  (gmx_mpi — para gmx check)

set -euo pipefail

# ── Defaults ──────────────────────────────────────────────────────────────────
RESULTS_DIR=""
SAMPLE_ID=""
START_NS=50
END_NS=""      # vazio = usar todos os frames até o fim
SALTCON=0.15
GMX_CMD="mpirun -np 1 gmx_mpi"

# ── Parse de argumentos ───────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
    case "$1" in
        --results-dir) RESULTS_DIR="$2"; shift 2 ;;
        --sample-id)   SAMPLE_ID="$2";   shift 2 ;;
        --start-ns)    START_NS="$2";    shift 2 ;;
        --end-ns)      END_NS="$2";      shift 2 ;;
        --saltcon)     SALTCON="$2";     shift 2 ;;
        --gmx-cmd)     GMX_CMD="$2";     shift 2 ;;
        *) echo "Argumento desconhecido: $1"; exit 1 ;;
    esac
done

if [[ -z "$RESULTS_DIR" || -z "$SAMPLE_ID" ]]; then
    echo "ERRO: --results-dir e --sample-id são obrigatórios"
    exit 1
fi

# ── Caminhos dos arquivos de entrada ─────────────────────────────────────────
PROD_DIR="${RESULTS_DIR}/${SAMPLE_ID}/prod"
ANAL_DIR="${RESULTS_DIR}/${SAMPLE_ID}/analise"
OUT_DIR="${RESULTS_DIR}/${SAMPLE_ID}/mmgbsa"

TPR="${PROD_DIR}/md.tpr"
XTC="${PROD_DIR}/md_fit.xtc"
NDX="${ANAL_DIR}/lig.ndx"

echo "=== MM-GBSA standalone: ${SAMPLE_ID} ==="
echo "TPR : ${TPR}"
echo "XTC : ${XTC}"
echo "NDX : ${NDX}"

for f in "$TPR" "$XTC" "$NDX"; do
    if [[ ! -f "$f" ]]; then
        echo "ERRO: arquivo não encontrado: $f"
        exit 1
    fi
done

# ── Verifica grupos Receptor/Ligante em lig.ndx ───────────────────────────────
if ! grep -q "Receptor" "$NDX"; then
    echo "ERRO: grupo 'Receptor' não encontrado em $NDX"
    exit 1
fi
if ! grep -q "Ligante" "$NDX"; then
    echo "ERRO: grupo 'Ligante' não encontrado em $NDX"
    exit 1
fi
echo "[OK] grupos 'Receptor' e 'Ligante' confirmados em lig.ndx"

# ── Conta frames reais na trajetória ─────────────────────────────────────────
echo "Contando frames em $XTC ..."
TOTAL_FRAMES=$(mamba run -n md-gromacs ${GMX_CMD} check -f "$XTC" 2>&1 \
    | grep -E "^Last frame" | awk '{print $3}' || echo "")

if [[ -z "$TOTAL_FRAMES" ]] || ! [[ "$TOTAL_FRAMES" =~ ^[0-9]+$ ]]; then
    echo "AVISO: não foi possível contar frames via gmx check — estimando 1000"
    TOTAL_FRAMES=1000
fi
echo "[OK] Total de frames: ${TOTAL_FRAMES}"

# Calcula frames a partir dos nanosegundos (detecta dt da trajetória)
# Assume 100 ps/frame (padrão do pipeline = 10 frames/ns)
FRAMES_PER_NS=10

START_FRAME=$(python3 -c "print(max(1, int(${START_NS} * ${FRAMES_PER_NS})))")
if [[ -n "$END_NS" ]]; then
    END_FRAME=$(python3 -c "print(min(${TOTAL_FRAMES}, int(${END_NS} * ${FRAMES_PER_NS})))")
else
    END_FRAME=$TOTAL_FRAMES
fi

echo "[OK] Janela de análise: ${START_NS} ns (frame ${START_FRAME}) → frame ${END_FRAME}"
echo "[OK] ${SALTCON} M NaCl"

# ── Prepara diretório de saída ────────────────────────────────────────────────
# Resolve paths absolutos ANTES de cd (paths relativos deixam de existir após cd)
ABS_TPR=$(realpath "$TPR")
ABS_XTC=$(realpath "$XTC")
ABS_NDX=$(realpath "$NDX")

# gmx_MMPBSA 1.6.x só aceita .tpr ou .pdb para -cs.
# mmgbsa-env tem GROMACS 2025.4 (TPR v137) mas TPR foram escritos por GROMACS 2026.0 (v138).
# Solução: converter md.gro → PDB via MDAnalysis (sem dependência de versão GROMACS).
GRO="${PROD_DIR}/md.gro"
PDB_CS="${PROD_DIR}/complex_mmgbsa.pdb"
ABS_CS="$ABS_TPR"   # default (falha se v138 vs v137)

if [[ -f "$GRO" ]]; then
    ABS_GRO=$(realpath "$GRO")
    if [[ ! -f "$PDB_CS" ]]; then
        echo "Convertendo md.gro → PDB via MDAnalysis (evita incompatibilidade TPR v138 vs v137)..."
        mamba run -n md-gromacs python3 -c "
import MDAnalysis as mda, warnings
warnings.filterwarnings('ignore')
u = mda.Universe('${ABS_GRO}')
u.atoms.write('${PDB_CS}')
print('[OK] PDB escrito: ${PDB_CS}')
" 2>&1 | grep -v DeprecationWarning || true
    fi
    if [[ -f "$PDB_CS" ]]; then
        ABS_CS=$(realpath "$PDB_CS")
        echo "[OK] Usando PDB como -cs: $PDB_CS"
    else
        echo "[AVISO] Conversão GRO→PDB falhou; tentando TPR diretamente (pode falhar por versão)"
    fi
else
    echo "[AVISO] md.gro ausente — usando TPR (pode falhar por incompatibilidade v138 vs v137)"
fi

mkdir -p "$OUT_DIR"
cd "$OUT_DIR"

# ── Gera mmgbsa.in ────────────────────────────────────────────────────────────
cat > mmgbsa.in << MEOF
&general
sys_name="${SAMPLE_ID}",
startframe=${START_FRAME},
endframe=${END_FRAME},
interval=1,
verbose=2,
/
&gb
igb=2,
saltcon=${SALTCON},
/
&decomp
idecomp=2,
dec_verbose=1,
/
MEOF
echo "[OK] mmgbsa.in gerado (startframe=${START_FRAME}, endframe=${END_FRAME})"

# ── Cria wrapper tleap (corrige bug SS bonds COM_OUT do gmx_MMPBSA 1.6.x) ────
mkdir -p bin_patch

cat > bin_patch/tleap << 'WEOF'
#!/usr/bin/env python3
# Wrapper tleap — corrige SS bonds em COM_OUT/LIG_OUT para complexos
# proteína-peptídeo E proteína-proteína (ex: SKTI como ligante com SS bonds).
# Bug gmx_MMPBSA 1.6.x: COM_OUT usa índices errados (não adiciona offset do receptor
# no LIGANTE, ou usa offset errado no RECEPTOR).
import sys, os, re, subprocess

LOG = os.path.join(os.getcwd(), 'tleap_wrapper.log')

def wlog(msg):
    with open(LOG, 'a') as f:
        f.write(msg + '\n')

def count_pdb_residues(pdb_path):
    """Conta resíduos únicos no PDB (sequência de ATOM/HETATM)."""
    seen = set()
    try:
        with open(pdb_path) as f:
            for line in f:
                if line.startswith(('ATOM','HETATM')):
                    resseq = line[22:26].strip()
                    chain  = line[21]
                    seen.add((chain, resseq))
        return len(seen)
    except Exception:
        return 0

args = sys.argv[1:]
wlog(f"[tleap-wrapper] chamado com: {' '.join(args)}")

for i, a in enumerate(args):
    if a == '-f' and i + 1 < len(args):
        fpath = args[i + 1]
        if not os.path.exists(fpath):
            wlog(f"[tleap-wrapper] arquivo não encontrado: {fpath}")
            break

        content   = open(fpath).read()
        work_dir  = os.path.dirname(os.path.abspath(fpath))

        rec_bonds = re.findall(r'bond\s+REC_OUT\.(\d+)\.SG\s+REC_OUT\.(\d+)\.SG', content)
        com_bonds = re.findall(r'bond\s+COM_OUT\.(\d+)\.SG\s+COM_OUT\.(\d+)\.SG', content)
        lig_bonds = re.findall(r'bond\s+LIG_OUT\.(\d+)\.SG\s+LIG_OUT\.(\d+)\.SG', content)

        wlog(f"[tleap-wrapper] REC SS bonds ({len(rec_bonds)}): {rec_bonds}")
        wlog(f"[tleap-wrapper] COM SS bonds ({len(com_bonds)}): {com_bonds}")
        wlog(f"[tleap-wrapper] LIG SS bonds ({len(lig_bonds)}): {lig_bonds}")

        modified = content
        fixes    = 0

        # ── Corrige COM_OUT bonds para a parte RECEPTOR ──────────────────────
        # O bug original: COM_OUT usa índices REC sem o offset do complexo, ou
        # usa índice errado. Corrigimos os primeiros len(rec_bonds) COM bonds
        # para corresponder exatamente aos REC bonds.
        if rec_bonds and com_bonds:
            n_rec_bonds = len(rec_bonds)
            # Só corrigir COM bonds que claramente correspondem ao receptor
            for k, ((cw0, cw1), (rr0, rr1)) in enumerate(zip(com_bonds[:n_rec_bonds], rec_bonds)):
                old = f'bond COM_OUT.{cw0}.SG COM_OUT.{cw1}.SG'
                new = f'bond COM_OUT.{rr0}.SG COM_OUT.{rr1}.SG'
                if old != new:
                    modified = modified.replace(old, new, 1)
                    fixes += 1
                    wlog(f"[tleap-wrapper] FIX-REC: '{old}' → '{new}'")

        # ── Corrige LIG_OUT bonds para proteína-ligante com SS bonds ─────────
        # O bug para SKTI: LIG_OUT usa índices relativos ao COMPLEXO (incorreto)
        # em vez de índices relativos ao ligante. Detectamos contando resíduos.
        if lig_bonds:
            lig_pdb = None
            for fname in ('_GMXMMPBSA_LIG.pdb', '_GMXMMPBSA_ligand.pdb'):
                candidate = os.path.join(work_dir, fname)
                if os.path.exists(candidate):
                    lig_pdb = candidate
                    break
            n_lig_res = count_pdb_residues(lig_pdb) if lig_pdb else 0
            wlog(f"[tleap-wrapper] LIG PDB: {lig_pdb}  n_res={n_lig_res}")

            rec_pdb = None
            for fname in ('_GMXMMPBSA_REC.pdb', '_GMXMMPBSA_receptor.pdb'):
                candidate = os.path.join(work_dir, fname)
                if os.path.exists(candidate):
                    rec_pdb = candidate
                    break
            n_rec_res = count_pdb_residues(rec_pdb) if rec_pdb else 0
            wlog(f"[tleap-wrapper] REC PDB: {rec_pdb}  n_res={n_rec_res}")

            if n_lig_res > 0 and n_rec_res > 0:
                for (lb0, lb1) in lig_bonds:
                    idx0, idx1 = int(lb0), int(lb1)
                    # Se os índices são maiores que o tamanho do ligante,
                    # provavelmente estão com offset do receptor (bug)
                    if idx0 > n_lig_res or idx1 > n_lig_res:
                        corrected0 = idx0 - n_rec_res
                        corrected1 = idx1 - n_rec_res
                        if corrected0 > 0 and corrected1 > 0:
                            old = f'bond LIG_OUT.{lb0}.SG LIG_OUT.{lb1}.SG'
                            new = f'bond LIG_OUT.{corrected0}.SG LIG_OUT.{corrected1}.SG'
                            modified = modified.replace(old, new, 1)
                            fixes += 1
                            wlog(f"[tleap-wrapper] FIX-LIG: '{old}' → '{new}'")
                        else:
                            wlog(f"[tleap-wrapper] AVISO: índice LIG corrigido < 1 ({corrected0},{corrected1}) — ignorando")
                    else:
                        wlog(f"[tleap-wrapper] LIG bond {lb0}-{lb1} parece correto (dentro de n_lig={n_lig_res})")
            elif lig_bonds:
                wlog(f"[tleap-wrapper] AVISO: não foi possível contar resíduos (PDB não encontrado) — LIG bonds inalterados")

        # ── Corrige COM_OUT bonds para LIGANTE (bonds além dos do receptor) ───
        # Quando o ligante é proteína, gmx_MMPBSA às vezes escreve SS do ligante
        # como COM_OUT em vez de LIG_OUT, usando índices do ligante sem offset.
        extra_com_bonds = com_bonds[len(rec_bonds):]
        if extra_com_bonds and n_lig_res > 0 and n_rec_res > 0:
            wlog(f"[tleap-wrapper] Extra COM bonds (ligante): {extra_com_bonds}")
            for (cw0, cw1) in extra_com_bonds:
                idx0, idx1 = int(cw0), int(cw1)
                # Se os índices são menores que n_rec_res, são relativos ao ligante
                # sem o offset; adicionar o offset
                if idx0 <= n_lig_res and idx1 <= n_lig_res:
                    corrected0 = idx0 + n_rec_res
                    corrected1 = idx1 + n_rec_res
                    old = f'bond COM_OUT.{cw0}.SG COM_OUT.{cw1}.SG'
                    new = f'bond COM_OUT.{corrected0}.SG COM_OUT.{corrected1}.SG'
                    modified = modified.replace(old, new, 1)
                    fixes += 1
                    wlog(f"[tleap-wrapper] FIX-COM-LIG: '{old}' → '{new}'")

        if fixes > 0:
            open(fpath, 'w').write(modified)
            wlog(f"[tleap-wrapper] {fixes} correção(ões) aplicada(s) → {fpath}")
        elif not rec_bonds and not com_bonds and not lig_bonds:
            wlog("[tleap-wrapper] sem pontes SS — sem correção necessária")
        else:
            wlog("[tleap-wrapper] índices já corretos — sem modificação")

        # Log parcial do leap.in para diagnóstico
        bond_lines = [l.strip() for l in content.splitlines() if 'bond' in l.lower() or '.SG' in l]
        wlog(f"[tleap-wrapper] linhas SS no leap.in:\n" + '\n'.join(bond_lines[:30]))
        break

path_dirs = os.environ.get('PATH', '').split(':')
for d in path_dirs:
    if 'bin_patch' in d:
        continue
    for name in ('tleap', 'teLeap'):
        exe = os.path.join(d, name)
        if os.path.isfile(exe) and os.access(exe, os.X_OK):
            wlog(f"[tleap-wrapper] executando: {exe} {' '.join(args)}")
            ret = subprocess.run([exe] + args)
            wlog(f"[tleap-wrapper] retcode: {ret.returncode}")
            sys.exit(ret.returncode)

wlog("[tleap-wrapper] ERRO FATAL: tleap real não encontrado no PATH")
sys.exit(1)
WEOF
chmod +x bin_patch/tleap
echo "[OK] bin_patch/tleap criado"

# ── Executa gmx_MMPBSA no ambiente isolado ────────────────────────────────────
CURRENT_DIR="$PWD"
echo ""
echo "[MMGBSA] Iniciando gmx_MMPBSA (pode demorar 20-60 min)..."
echo "[MMGBSA] Saídas em: $CURRENT_DIR"
echo ""

mamba run -n mmgbsa-env bash -c "
export PATH=${CURRENT_DIR}/bin_patch:\$PATH
echo '[mmgbsa-env] PATH patch ativo'

gmx_MMPBSA -O \\
    -i mmgbsa.in \\
    -cs '${ABS_CS}' \\
    -ct '${ABS_XTC}' \\
    -ci '${ABS_NDX}' \\
    -cg Receptor Ligante \\
    -o  FINAL_RESULTS_MMGBSA.dat \\
    -eo mmgbsa_results.csv \\
    -deo decomp_results.csv \\
    -nogui \\
    2>&1
" 2>&1 | tee mmgbsa.log

# ── Verifica saídas ───────────────────────────────────────────────────────────
echo ""
if [[ -f "FINAL_RESULTS_MMGBSA.dat" ]]; then
    echo "=== RESULTADO MM-GBSA: ${SAMPLE_ID} ==="
    grep -A 5 "DELTA TOTAL" FINAL_RESULTS_MMGBSA.dat || true
    echo ""
    echo "[OK] FINAL_RESULTS_MMGBSA.dat gerado"
    [[ -f "decomp_results.csv" ]] || {
        echo "resid,resname,total" > decomp_results.csv
        echo "AVISO: decomp_results.csv não gerado — arquivo vazio criado"
    }
    echo "[OK] MM-GBSA concluído para ${SAMPLE_ID}"
    echo "[OK] Saídas em: ${OUT_DIR}"
else
    echo "ERRO: FINAL_RESULTS_MMGBSA.dat não gerado"
    echo "--- Últimas 50 linhas do log ---"
    tail -50 mmgbsa.log
    exit 1
fi
