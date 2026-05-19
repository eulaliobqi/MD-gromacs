#!/usr/bin/env bash
# ============================================================================
# Pipeline de DM: complexo proteina-peptideo via GROMACS
# Generico para qualquer tripsina + inibidor peptidico
#
# Uso:
#   ./dm_pipeline.sh -r receptor.pdb -l ligante.pdb [opcoes]
#
# Opcoes:
#   -r FILE    Receptor PDB (obrigatorio)
#   -l FILE    Ligante PDB (obrigatorio)
#   -p NAME    Nome do projeto (default: dm_complex)
#   -H FLOAT   pH (default: 7.4)
#   -t INT     Tempo de producao em ns (default: 100)
#   -F STR     Force field (default: amber99sb-ildn)
#   -W STR     Modelo de agua (default: tip3p)
#   -S FLOAT   Concentracao de NaCl em M (default: 0.15)
#   -B STR     Diretorio base (default: $HOME/gromacs/<projeto>)
#   -G INT     GPU id (default: 0)
#   -T INT     OMP threads (default: 8)
#   -h         Mostra ajuda
# ============================================================================
set -euo pipefail

# Defaults
PROJ_NAME="dm_complex"
PH=7.4
PROD_NS=100
FF="amber99sb-ildn"
WATER="tip3p"
ION_CONC=0.15
GPU_ID=0
NTOMP=8
BASE=""

# Parse args
while getopts "r:l:p:H:t:F:W:S:B:G:T:h" opt; do
    case $opt in
        r) RECEPTOR_PDB="$OPTARG" ;;
        l) LIGANTE_PDB="$OPTARG" ;;
        p) PROJ_NAME="$OPTARG" ;;
        H) PH="$OPTARG" ;;
        t) PROD_NS="$OPTARG" ;;
        F) FF="$OPTARG" ;;
        W) WATER="$OPTARG" ;;
        S) ION_CONC="$OPTARG" ;;
        B) BASE="$OPTARG" ;;
        G) GPU_ID="$OPTARG" ;;
        T) NTOMP="$OPTARG" ;;
        h) sed -n '2,20p' "$0"; exit 0 ;;
        *) echo "Opcao invalida"; exit 1 ;;
    esac
done

[[ -z "${RECEPTOR_PDB:-}" || -z "${LIGANTE_PDB:-}" ]] && {
    echo "ERRO: -r e -l sao obrigatorios. Use -h para ajuda."; exit 1; }
[[ -z "$BASE" ]] && BASE="${HOME}/gromacs/${PROJ_NAME}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RECEPTOR_PDB="$(realpath "$RECEPTOR_PDB")"
LIGANTE_PDB="$(realpath "$LIGANTE_PDB")"

# Calcula nsteps de producao (dt = 0.002 ps -> 500 steps/ps)
PROD_STEPS=$((PROD_NS * 500000))

# ----- Helpers -----
GREEN='\033[0;32m'; BLUE='\033[0;34m'; YELLOW='\033[1;33m'; RED='\033[0;31m'; NC='\033[0m'
banner() { echo -e "\n${BLUE}+--- $1 ---${NC}\n"; }
ok()     { echo -e "  ${GREEN}OK${NC}  $1"; }
warn()   { echo -e "  ${YELLOW}!!${NC}  $1"; }
err()    { echo -e "${RED}ERRO: $1${NC}" >&2; exit 1; }
check()  { [[ -s "$1" ]] || err "Arquivo nao encontrado: $1"; }

# Detecta GROMACS
GMX=$(command -v gmx_mpi 2>/dev/null || command -v gmx 2>/dev/null || echo "")
[[ -z "$GMX" ]] && err "gmx ou gmx_mpi nao encontrado."

# mdrun para producao/equilibracao (com PME GPU)
mdrun() {
    if [[ "$GMX" == *"gmx_mpi"* ]]; then
        mpirun --use-hwthread-cpus -np 1 "$GMX" mdrun \
            -ntomp "$NTOMP" -nb gpu -pme gpu -bonded gpu \
            -gpu_id "$GPU_ID" -pin on "$@"
    else
        "$GMX" mdrun -ntomp "$NTOMP" -nb gpu -pme gpu -bonded gpu \
            -gpu_id "$GPU_ID" -pin on "$@"
    fi
}

# mdrun para minimizacao de energia (steep nao suporta PME GPU)
em_mdrun() {
    if [[ "$GMX" == *"gmx_mpi"* ]]; then
        mpirun --use-hwthread-cpus -np 1 "$GMX" mdrun \
            -ntomp "$NTOMP" -nb gpu \
            -gpu_id "$GPU_ID" -pin on "$@"
    else
        "$GMX" mdrun -ntomp "$NTOMP" -nb gpu \
            -gpu_id "$GPU_ID" -pin on "$@"
    fi
}

# ============================================================
banner "ETAPA 0: Configuracao"
echo "  Projeto      : $PROJ_NAME"
echo "  Receptor     : $RECEPTOR_PDB"
echo "  Ligante      : $LIGANTE_PDB"
echo "  pH           : $PH"
echo "  Tempo        : $PROD_NS ns ($PROD_STEPS passos)"
echo "  Force field  : $FF / $WATER"
echo "  NaCl         : $ION_CONC M"
echo "  Diretorio    : $BASE"
echo "  GROMACS      : $GMX"

mkdir -p "$BASE"/{prep,topo,box,em,nvt,npt,prod,analise}
ok "Diretorios criados"

# ============================================================
banner "ETAPA 1: Preparacao do complexo (prepare_complex.py)"
cd "$BASE/prep"
python3 "$SCRIPT_DIR/prepare_complex.py" \
    --receptor "$RECEPTOR_PDB" \
    --ligante  "$LIGANTE_PDB" \
    --pH "$PH" \
    --out-dir .
check complexo.pdb
ok "Complexo preparado"

# ============================================================
banner "ETAPA 2: Topologia (pdb2gmx)"
cd "$BASE/topo"
if [[ ! -s complexo.gro ]]; then
    printf '0\n0\n0\n0\n' | "$GMX" pdb2gmx \
        -f "$BASE/prep/complexo.pdb" \
        -o complexo.gro -p topol.top -i posre.itp \
        -ff "$FF" -water "$WATER" \
        -ignh -ter -chainsep ter -merge no \
        2>&1 | tee pdb2gmx.log
fi
check complexo.gro
ok "Topologia gerada"

# ============================================================
banner "ETAPA 3: Caixa, solvatacao, ions"
cd "$BASE/box"
[[ ! -s box.gro ]] && "$GMX" editconf -f "$BASE/topo/complexo.gro" \
    -o box.gro -c -d 1.2 -bt dodecahedron && ok "Caixa dodecaedrica"
[[ ! -s solv.gro ]] && "$GMX" solvate -cp box.gro -cs spc216.gro \
    -p "$BASE/topo/topol.top" -o solv.gro && ok "Solvatado"

cat > ions.mdp <<EOF
integrator    = steep
nsteps        = 0
cutoff-scheme = Verlet
EOF

if [[ ! -s ions.gro ]]; then
    "$GMX" grompp -f ions.mdp -c solv.gro \
        -p "$BASE/topo/topol.top" -o ions.tpr -maxwarn 2
    echo "SOL" | "$GMX" genion -s ions.tpr -o ions.gro \
        -p "$BASE/topo/topol.top" -pname NA -nname CL \
        -neutral -conc "$ION_CONC"
    ok "Ions adicionados (NaCl $ION_CONC M)"
fi

# Copia para EM
cp ions.gro "$BASE/em/sistema.gro"
cp "$BASE/topo/topol.top" "$BASE/em/"
cp "$BASE/topo/"*.itp "$BASE/em/" 2>/dev/null || true

# ============================================================
banner "ETAPA 4: Minimizacao de energia"
cd "$BASE/em"
cat > em.mdp << 'EOF'
integrator      = steep
emtol           = 1000.0
emstep          = 0.01
nsteps          = 50000
cutoff-scheme   = Verlet
nstlist         = 10
coulombtype     = PME
rcoulomb        = 1.2
vdwtype         = Cut-off
rvdw            = 1.2
pbc             = xyz
EOF

if [[ ! -s em.gro ]]; then
    "$GMX" grompp -f em.mdp -c sistema.gro -p topol.top -o em.tpr -maxwarn 2
    em_mdrun -v -deffnm em
    ok "EM concluida"
fi
check em.gro

# ============================================================
banner "ETAPA 5: NVT (100 ps com restricoes)"
cd "$BASE/nvt"
cp "$BASE/em/topol.top" .
cp "$BASE/em/"*.itp . 2>/dev/null || true

cat > nvt.mdp << 'EOF'
define          = -DPOSRES
integrator      = md
dt              = 0.002
nsteps          = 50000
nstxout         = 500
nstvout         = 500
nstenergy       = 100
cutoff-scheme   = Verlet
nstlist         = 20
coulombtype     = PME
rcoulomb        = 1.2
vdwtype         = Cut-off
rvdw            = 1.2
constraints     = h-bonds
constraint-algorithm = LINCS
continuation    = no
gen-vel         = yes
gen-temp        = 300
gen-seed        = 42
tcoupl          = V-rescale
tc-grps         = Protein Non-Protein
tau-t           = 0.1 0.1
ref-t           = 300 300
pcoupl          = no
pbc             = xyz
EOF

if [[ ! -s nvt.gro ]]; then
    "$GMX" grompp -f nvt.mdp -c "$BASE/em/em.gro" -r "$BASE/em/em.gro" \
        -p topol.top -o nvt.tpr -maxwarn 2
    mdrun -v -deffnm nvt
    ok "NVT concluido"
fi

# ============================================================
banner "ETAPA 6: NPT (100 ps a 1 bar)"
cd "$BASE/npt"
cp "$BASE/nvt/topol.top" .
cp "$BASE/nvt/"*.itp . 2>/dev/null || true

cat > npt.mdp << 'EOF'
define          = -DPOSRES
integrator      = md
dt              = 0.002
nsteps          = 50000
nstxout         = 500
nstvout         = 500
nstenergy       = 100
cutoff-scheme   = Verlet
nstlist         = 20
coulombtype     = PME
rcoulomb        = 1.2
vdwtype         = Cut-off
rvdw            = 1.2
constraints     = h-bonds
constraint-algorithm = LINCS
continuation    = yes
tcoupl          = V-rescale
tc-grps         = Protein Non-Protein
tau-t           = 0.1 0.1
ref-t           = 300 300
pcoupl          = C-rescale
pcoupltype      = isotropic
tau-p           = 2.0
ref-p           = 1.0
compressibility = 4.5e-5
refcoord_scaling = com
pbc             = xyz
EOF

if [[ ! -s npt.gro ]]; then
    "$GMX" grompp -f npt.mdp -c "$BASE/nvt/nvt.gro" -r "$BASE/nvt/nvt.gro" \
        -t "$BASE/nvt/nvt.cpt" -p topol.top -o npt.tpr -maxwarn 2
    mdrun -v -deffnm npt
    ok "NPT concluido"
fi

# ============================================================
banner "ETAPA 7: Producao ${PROD_NS} ns"
cd "$BASE/prod"
cp "$BASE/npt/topol.top" .
cp "$BASE/npt/"*.itp . 2>/dev/null || true

cat > md.mdp << EOF
integrator      = md
dt              = 0.002
nsteps          = ${PROD_STEPS}
nstxout         = 0
nstvout         = 0
nstenergy       = 5000
nstlog          = 5000
nstxout-compressed = 5000
cutoff-scheme   = Verlet
nstlist         = 20
coulombtype     = PME
rcoulomb        = 1.2
vdwtype         = Cut-off
rvdw            = 1.2
constraints     = h-bonds
constraint-algorithm = LINCS
continuation    = yes
tcoupl          = V-rescale
tc-grps         = Protein Non-Protein
tau-t           = 0.1 0.1
ref-t           = 300 300
pcoupl          = Parrinello-Rahman
pcoupltype      = isotropic
tau-p           = 2.0
ref-p           = 1.0
compressibility = 4.5e-5
pbc             = xyz
comm-mode       = Linear
nstcomm         = 100
EOF

if [[ ! -s md.tpr ]]; then
    "$GMX" grompp -f md.mdp -c "$BASE/npt/npt.gro" -t "$BASE/npt/npt.cpt" \
        -p topol.top -o md.tpr -maxwarn 2
fi

if [[ ! -s md.gro ]]; then
    CPI=""; [[ -s md.cpt ]] && CPI="-cpi md.cpt"
    mdrun -v -deffnm md $CPI
    ok "Producao concluida"
fi

# ============================================================
banner "ETAPA 8: Pos-processamento (PBC + alinhamento)"
cd "$BASE/prod"

# nojump primeiro evita o spike de PBC
[[ ! -s md_nojump.xtc ]] && \
    echo 0 | "$GMX" trjconv -s md.tpr -f md.xtc \
        -o md_nojump.xtc -pbc nojump && ok "nojump OK"

[[ ! -s md_nopbc.xtc ]] && \
    echo -e '1\n0' | "$GMX" trjconv -s md.tpr -f md_nojump.xtc \
        -o md_nopbc.xtc -pbc mol -center -ur compact && ok "PBC corrigido"

[[ ! -s md_fit.xtc ]] && \
    echo -e '4\n0' | "$GMX" trjconv -s md.tpr -f md_nopbc.xtc \
        -o md_fit.xtc -fit rot+trans && ok "Alinhado pelo backbone"

# ============================================================
banner "ETAPA 9: Analises (delegado a run_analyses.sh)"
bash "$SCRIPT_DIR/run_analyses.sh" "$BASE"

# ============================================================
banner "ETAPA 10: Figuras (delegado a plot_results.py)"
python3 "$SCRIPT_DIR/plot_results.py" --analise-dir "$BASE/analise" \
    --titulo "$PROJ_NAME - DM ${PROD_NS} ns @ pH ${PH}"

banner "PIPELINE CONCLUIDO"
echo "  Trajetoria  : $BASE/prod/md_fit.xtc"
echo "  Analises    : $BASE/analise/"
echo "  Figura      : $BASE/analise/painel_completo.png"
