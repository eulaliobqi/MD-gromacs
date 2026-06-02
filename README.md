# MD-GROMACS — Pipeline Nextflow para Dinâmica Molecular

Pipeline automatizado para simulações de dinâmica molecular de complexos
**proteína–peptídeo** via GROMACS, usando Nextflow DSL2.

Desenvolvido e validado para tripsina digestiva (*Spodoptera* spp.) + inibidores
peptídicos (GORE4) a pH 8,2 por 100–200 ns.

---

## Índice

1. [Visão geral](#1-visão-geral)
2. [Estrutura do repositório](#2-estrutura-do-repositório)
3. [Fluxo do pipeline](#3-fluxo-do-pipeline)
4. [Pré-requisitos](#4-pré-requisitos)
5. [Passo a passo completo](#5-passo-a-passo-completo)
   - [5.1 Clonar o repositório](#51-clonar-o-repositório)
   - [5.2 Instalar dependências](#52-instalar-dependências)
   - [5.3 Preparar as entradas (PDBs)](#53-preparar-as-entradas-pdbs)
   - [5.4 Criar o samplesheet](#54-criar-o-samplesheet)
   - [5.5 Configurar execução](#55-configurar-execução)
   - [5.6 Executar o pipeline](#56-executar-o-pipeline)
   - [5.7 Monitorar e retomar](#57-monitorar-e-retomar)
   - [5.8 Ver os resultados](#58-ver-os-resultados)
6. [Exemplos práticos](#6-exemplos-práticos)
7. [Parâmetros](#7-parâmetros)
8. [Saídas detalhadas](#8-saídas-detalhadas)
9. [Solução de problemas](#9-solução-de-problemas)
10. [Notas técnicas](#10-notas-técnicas)

---

## 1. Visão geral

Este pipeline automatiza **todas as etapas** de uma DM de complexo proteína–peptídeo:
da preparação do PDB até os gráficos de análise prontos para publicação.

**O que o pipeline faz:**

- Prepara o complexo: protonação de His por pH, detecção automática de pontes
  dissulfeto (CYX), atribuição de cadeias A/B
- Gera topologia com force field AMBER99sb-ILDN e água TIP3P
- Constrói a caixa **cúbica** com margem de 2,0 nm, solvatação e adição de **KCl 0,10 M**
- Equilibração: minimização de energia → NVT 200 ps → NPT 500 ps (Berendsen)
- Produção: 100–200 ns NPT com Parrinello-Rahman e GPU
- Pós-processamento: correção de PBC e alinhamento da trajetória
- Análises: RMSD, RMSF, raio de giro, contatos, pontes de H, SASA,
  distâncias aos resíduos catalíticos/S1
- Figuras: painel completo e PNGs individuais prontos para artigo

**Parâmetros físicos otimizados para lepidópteros:**

| Parâmetro | Valor | Justificativa |
|-----------|-------|---------------|
| pH | 8,2 | Midgut alcalino de *Spodoptera*; His57 → HIE (forma ativa) |
| Cátion | K⁺ | Hemolínfa de insetos é K⁺-dominada (vs Na⁺ em mamíferos) |
| [KCl] | 0,10 M | Midgut hipotônico (literatura 0,10–0,15 M) |
| Temperatura | 300 K | Fisiológica larval (25–30 °C) |
| Force field | AMBER99sb-ILDN / TIP3P | Padrão para proteínas globulares |

**Receptores disponíveis neste repositório** (clusters HADDOCK vs GORE4):

| Receptor | Clusters |
|----------|----------|
| ACR157   | 10       |
| DN1441   | 10       |
| DN773    | 8        |
| DN1937   | 10       |
| QCL936   | 10       |
| XP273    | 10       |
| XP352    | 10       |

Cada cluster HADDOCK contém 4 poses (modelos 1–4).
Total: **68 pares** possíveis para simulação.

---

## 2. Estrutura do repositório

```
MD-gromacs/
│
├── main.nf                        # Workflow principal Nextflow
├── nextflow.config                # Parâmetros e perfis de execução
├── environment.yml                # Ambiente conda (gromacs, python, nextflow...)
│
├── modules/local/                 # Um módulo por etapa do pipeline
│   ├── prepare_complex/main.nf    # Preparação do complexo (CYX, HIS, chains)
│   ├── topology/main.nf           # pdb2gmx → topol.top + *.itp
│   ├── box_solvate_ions/main.nf   # editconf (cúbico 2 nm) + solvate + genion (KCl)
│   ├── minimization/main.nf       # Minimização de energia (steep)
│   ├── nvt/main.nf                # Equil. NVT 200 ps
│   ├── npt/main.nf                # Equil. NPT 500 ps (Berendsen)
│   ├── production/main.nf         # DM de produção (N ns, Parrinello-Rahman)
│   ├── postprocess/main.nf        # Correção PBC + alinhamento
│   ├── analyses/main.nf           # RMSD, RMSF, Rg, contatos, H-bonds, SASA
│   ├── analyses_triad/main.nf     # Distâncias aos resíduos catalíticos/S1
│   └── plot/main.nf               # Figuras PNG
│
├── conf/
│   ├── base.config                # Labels de recursos (CPU/GPU/memória/tempo)
│   ├── local.config               # Perfil para execução local
│   └── slurm.config               # Perfil para cluster SLURM (HPC)
│
├── bin/                           # Scripts auxiliares (adicionados ao PATH pelo Nextflow)
│   ├── prepare_complex.py         # Pré-processamento do PDB
│   ├── run_analyses.sh            # Re-executar análises manualmente
│   └── plot_results.py            # Gerar figuras manualmente
│
├── assets/
│   ├── samplesheet_example.csv    # Exemplo de samplesheet
│   ├── samplesheet_xp273_c1.csv   # XP273-GORE4 cluster 1
│   ├── samplesheet_xp352_c4.csv   # XP352-GORE4 cluster 4
│   └── samplesheet_acr157_c2.csv  # ACR157-GORE4 cluster 2
│
├── ACR157-GORE4-HADDOCK/          # Clusters HADDOCK: receptor ACR157 + peptídeo GORE4
├── DN1441-GORE4-HADDOCK/          #   cluster1_1.pdb ... clusterN_4.pdb
├── DN773-GORE4-HADDOCK/           #   (cadeias A=receptor, B=ligante)
├── DN1937-GORE4-HADDOCK/
├── QCL936-GORE4-HADDOCK/
├── XP273-GORE4-HADDOCK/
└── XP352-GORE4-HADDOCK/
```

---

## 3. Fluxo do pipeline

```
Entradas
   receptor.pdb  +  ligante.pdb
         │
         ▼
┌─────────────────────┐
│   PREPARE_COMPLEX   │  Protonação His por pH, CYX dissulfeto,
│                     │  chain A/B, renumeração do ligante
└──────────┬──────────┘
           │ complexo.pdb
           ▼
┌─────────────────────┐
│      TOPOLOGY       │  gmx pdb2gmx → complexo.gro + topol.top + *.itp
└──────────┬──────────┘
           │
           ▼
┌─────────────────────┐
│  BOX_SOLVATE_IONS   │  Cúbico 2,0 nm · TIP3P · KCl 0,10 M
└──────────┬──────────┘
           │ ions.gro
           ▼
┌─────────────────────┐
│    MINIMIZATION     │  Steep 50 000 passos  (sem PME GPU)
└──────────┬──────────┘
           │ em.gro
           ▼
┌─────────────────────┐
│        NVT          │  200 ps · 300 K · restrições de posição
└──────────┬──────────┘
           │ nvt.gro + nvt.cpt
           ▼
┌─────────────────────┐
│        NPT          │  500 ps · 1 bar · Berendsen
└──────────┬──────────┘
           │ npt.gro + npt.cpt
           ▼
┌─────────────────────┐
│     PRODUCTION      │  100–200 ns NPT · Parrinello-Rahman · GPU
└──────────┬──────────┘
           │ md.tpr + md.xtc
           ▼
┌─────────────────────┐
│     POSTPROCESS     │  nojump → pbc mol center → fit backbone
└──────────┬──────────┘
           │ md_fit.xtc
           ▼
┌─────────────────────┐
│      ANALYSES       │  RMSD backbone/ligante · RMSF · Rg
│                     │  contatos < 0,4 nm · pontes de H · SASA
└──────────┬──────────┘
           │ *.xvg
           ▼
┌─────────────────────┐
│   ANALYSES_TRIAD    │  Distâncias ligante → resíduos catalíticos/S1
│                     │  (triad_1..4 definidos no samplesheet)
└──────────┬──────────┘
           │ dist_triad*.xvg
           ▼
┌─────────────────────┐
│        PLOT         │  painel_completo.png + PNGs individuais
└─────────────────────┘
```

Cada etapa é um **processo Nextflow independente**: se o pipeline for
interrompido, `-resume` retoma do ponto exato onde parou (desde que nenhum
parâmetro de simulação tenha mudado).

---

## 4. Pré-requisitos

### Software

| Ferramenta  | Versão mínima | Observação                          |
|-------------|---------------|-------------------------------------|
| Nextflow    | 23.04         | `nextflow -version`                 |
| Conda/Mamba | qualquer      | recomendado: mamba (mais rápido)    |
| GROMACS     | 2022          | instalado via conda ou módulo HPC   |
| Python      | 3.10          | instalado via conda                 |
| Java        | 11            | necessário para o Nextflow          |

### Hardware

- **CPU**: mínimo 4 cores para equilibração; 8+ recomendado
- **GPU NVIDIA**: fortemente recomendado para produção (CUDA 11+)
- **Memória**: 16 GB RAM mínimo; 32 GB recomendado para sistemas grandes
- **Disco**: ~5–10 GB por simulação de 100 ns

---

## 5. Passo a passo completo

### 5.1 Clonar o repositório

```bash
git clone https://github.com/eulaliobqi/MD-gromacs.git
cd MD-gromacs
```

Verifique o conteúdo:

```bash
ls
# ACR157-GORE4-HADDOCK/  DN1441-GORE4-HADDOCK/  main.nf  nextflow.config ...
```

---

### 5.2 Instalar dependências

**Opção A — Ambiente conda (recomendada):**

```bash
# Criar o ambiente a partir do arquivo environment.yml
mamba env create -f environment.yml

# Ativar o ambiente
mamba activate md-gromacs

# Verificar instalações
gmx_mpi --version
nextflow -version
python --version
```

**Opção B — GROMACS já instalado no servidor (módulo HPC):**

```bash
# Carregue o módulo do cluster (exemplo)
module load gromacs/2022

# Instale apenas python e nextflow via conda
mamba create -n md-nextflow python=3.10 nextflow numpy matplotlib mdtraj
mamba activate md-nextflow
```

---

### 5.3 Preparar as entradas (PDBs)

O pipeline precisa de **dois PDBs separados** para cada simulação:
- `receptor.pdb` — proteína (cadeia A)
- `ligante.pdb`  — peptídeo (cadeia B)

#### Caso use os clusters HADDOCK do repositório

Os clusters HADDOCK (`clusterN_M.pdb`) contêm **complexo completo** com
cadeias A e B no mesmo arquivo. É preciso separá-las:

```bash
# Separar cadeia A (receptor) e cadeia B (ligante) de um cluster
CLUSTER="DN1441-GORE4-HADDOCK/cluster1_1.pdb"

awk '/^ATOM/ && substr($0,22,1)=="A"' $CLUSTER >  receptor_DN1441_c1.pdb
echo "TER" >> receptor_DN1441_c1.pdb
echo "END" >> receptor_DN1441_c1.pdb

awk '/^ATOM/ && substr($0,22,1)=="B"' $CLUSTER >  ligante_DN1441_c1.pdb
echo "TER" >> ligante_DN1441_c1.pdb
echo "END" >> ligante_DN1441_c1.pdb
```

#### Script para separar todos os clusters de um receptor de uma vez

```bash
RECEPTOR="XP273"
HADDOCK_DIR="${RECEPTOR}-GORE4-HADDOCK"
PREP_DIR="inputs/${RECEPTOR}"
mkdir -p "$PREP_DIR"

for pdb in ${HADDOCK_DIR}/cluster*_1.pdb; do
    base=$(basename $pdb .pdb)
    awk '/^ATOM/ && substr($0,22,1)=="A"' $pdb > ${PREP_DIR}/rec_${base}.pdb
    echo -e "TER\nEND" >> ${PREP_DIR}/rec_${base}.pdb

    awk '/^ATOM/ && substr($0,22,1)=="B"' $pdb > ${PREP_DIR}/lig_${base}.pdb
    echo -e "TER\nEND" >> ${PREP_DIR}/lig_${base}.pdb
done

ls ${PREP_DIR}/
```

---

### 5.4 Criar o samplesheet

O samplesheet é um arquivo CSV com **7 colunas**:

| Coluna      | Descrição                                              |
|-------------|--------------------------------------------------------|
| `sample_id` | Nome único para a simulação                            |
| `receptor`  | Caminho para o PDB do receptor                         |
| `ligand`    | Caminho para o PDB do ligante                          |
| `triad_1`   | Número do resíduo 1 para distâncias (ex: Tyr do sítio ativo) |
| `triad_2`   | Número do resíduo 2 (ex: Asp da tríade catalítica)    |
| `triad_3`   | Número do resíduo 3 (ex: Ser nucleofílica)             |
| `triad_4`   | Número do resíduo 4 (ex: Ile do bolsão S1)            |

> Os campos `triad_1..4` definem os resíduos-alvo para o módulo
> `ANALYSES_TRIAD` (distâncias ligante–resíduo ao longo da trajetória).
> Use `0` para desabilitar um resíduo.

**Exemplo — GORE4 tripsina (resíduos catalíticos reais):**

```bash
cat > samplesheet.csv << 'EOF'
sample_id,receptor,ligand,triad_1,triad_2,triad_3,triad_4
xp273-gore4-c1,inputs/XP273/rec_cluster1_1.pdb,inputs/XP273/lig_cluster1_1.pdb,83,132,234,229
acr157-gore4-c2,inputs/ACR157/rec_cluster2_1.pdb,inputs/ACR157/lig_cluster2_1.pdb,69,114,211,205
EOF
```

**Exemplo — samplesheets prontos no repositório:**

```bash
# XP273-GORE4 cluster 1
cat assets/samplesheet_xp273_c1.csv

# ACR157-GORE4 cluster 2
cat assets/samplesheet_acr157_c2.csv

# XP352-GORE4 cluster 4
cat assets/samplesheet_xp352_c4.csv
```

> Os caminhos podem ser relativos ao diretório de execução ou absolutos.

---

### 5.5 Configurar execução

#### Execução local (workstation com GPU)

Nenhuma configuração extra necessária. O perfil `local` usa os recursos
da máquina diretamente.

#### Execução em cluster SLURM (HPC)

Edite `conf/slurm.config` e ajuste os nomes das filas do seu cluster:

```bash
nano conf/slurm.config
```

```groovy
// Troque 'cpu' e 'gpu' pelos nomes reais das filas do seu cluster
withLabel: 'process_medium' {
    queue = 'cpu'      // ← altere aqui
    ...
}
withLabel: 'process_gpu' {
    queue = 'gpu'      // ← altere aqui
    clusterOptions = '--gres=gpu:1'   // ← ajuste conforme seu cluster
    ...
}
```

Para descobrir as filas disponíveis no seu cluster:

```bash
sinfo          # lista todas as filas e nós
squeue -u $USER   # jobs em execução/pendentes
```

---

### 5.6 Executar o pipeline

#### Execução local

```bash
nextflow run main.nf \
    --input assets/samplesheet_xp273_c1.csv \
    --outdir ~/gromacs/results/XP273-GORE4/MD \
    --pH 8.2 \
    --time_ns 100 \
    -profile conda,local
```

#### Execução em cluster SLURM

```bash
nextflow run main.nf \
    --input assets/samplesheet_xp273_c1.csv \
    --outdir ~/gromacs/results/XP273-GORE4/MD \
    --pH 8.2 \
    --time_ns 100 \
    --gmx_cmd gmx_mpi \
    --mpi_cmd "mpirun --use-hwthread-cpus -np 1" \
    -profile slurm,conda
```

#### Sem GPU (apenas CPU)

```bash
nextflow run main.nf \
    --input assets/samplesheet_xp273_c1.csv \
    --use_gpu false \
    -profile local,conda
```

> **Dica:** execute o Nextflow dentro de uma sessão `screen` para não perder
> a execução se a conexão SSH cair:
> ```bash
> screen -S md_pipeline
> nextflow run main.nf --input assets/samplesheet_xp273_c1.csv -profile conda,local
> # Ctrl+A, D para desanexar; screen -r md_pipeline para reconectar
> ```

> **Atenção `-resume`:** use `-resume` apenas para retomar um run interrompido
> com os **mesmos parâmetros**. Se alterar `box_dist`, `nacl_conc`, tempos de
> equilibração ou force field, rode sem `-resume` em um novo `--outdir`.

---

### 5.7 Monitorar e retomar

**Acompanhar em tempo real:**

```bash
# O Nextflow imprime o progresso no terminal:
# [d6/3a9bc1] PREPARE_COMPLEX (xp273-gore4-c1) [100%] 1 of 1 ✔
# [f2/1c88e2] TOPOLOGY (xp273-gore4-c1)        [ 50%] 0 of 1
```

**Ver log detalhado de um processo:**

```bash
# O ID do processo aparece no início de cada linha (ex: d6/3a9bc1)
cat work/d6/3a9bc1*/\.command.log
```

**Retomar execução interrompida (mesmos parâmetros):**

```bash
nextflow run main.nf --input assets/samplesheet_xp273_c1.csv -profile conda,local -resume
```

**Verificar status dos jobs no SLURM:**

```bash
squeue -u $USER
```

---

### 5.8 Ver os resultados

Os resultados ficam em `~/gromacs/results/<sample_id>/` (padrão)
ou no diretório definido por `--outdir`.

```bash
# Ver arquivos gerados para uma simulação
ls ~/gromacs/results/xp273-gore4-c1/

# Verificar se a produção terminou
ls ~/gromacs/results/xp273-gore4-c1/prod/
# md.tpr  md.xtc  md.gro  md.cpt  md_fit.xtc

# Abrir o painel de análises
ls ~/gromacs/results/xp273-gore4-c1/analise/
# rmsd_backbone.xvg  rmsd_ligante.xvg  rmsf_residuos.xvg
# gyrate.xvg  numcont.xvg  hbond.xvg  sasa_receptor.xvg  sasa_ligante.xvg
# dist_triad1.xvg  dist_triad2.xvg  dist_triad3.xvg  dist_triad4.xvg
# painel_completo.png  rmsd_bb.png  rmsd_lig.png  ...
```

Baixe o painel de figuras para visualizar:

```bash
# No seu computador local:
scp usuario@servidor:~/gromacs/results/xp273-gore4-c1/analise/painel_completo.png .
```

---

## 6. Exemplos práticos

### Exemplo 1 — Samplesheet de teste incluso no repositório

```bash
nextflow run main.nf \
    --input assets/samplesheet_example.csv \
    --outdir ~/gromacs/results/teste \
    --pH 8.2 \
    --time_ns 100 \
    -profile conda,local
```

### Exemplo 2 — Um único cluster HADDOCK

```bash
# 1. Separar as cadeias
awk '/^ATOM/ && substr($0,22,1)=="A"' XP273-GORE4-HADDOCK/cluster1_1.pdb > rec.pdb
echo -e "TER\nEND" >> rec.pdb
awk '/^ATOM/ && substr($0,22,1)=="B"' XP273-GORE4-HADDOCK/cluster1_1.pdb > lig.pdb
echo -e "TER\nEND" >> lig.pdb

# 2. Criar samplesheet (com resíduos catalíticos do GORE4)
cat > sheet.csv << 'EOF'
sample_id,receptor,ligand,triad_1,triad_2,triad_3,triad_4
xp273_c1,rec.pdb,lig.pdb,83,132,234,229
EOF

# 3. Executar
nextflow run main.nf \
    --input sheet.csv \
    --outdir ~/gromacs/results/XP273-GORE4/MD \
    --pH 8.2 --time_ns 100 \
    -profile conda,local
```

### Exemplo 3 — Todos os clusters top-1 de todos os receptores

```bash
# Gera inputs separados para todos os cluster*_1.pdb de todos os receptores
mkdir -p inputs
echo "sample_id,receptor,ligand,triad_1,triad_2,triad_3,triad_4" > samplesheet_all.csv

for dir in *-GORE4-HADDOCK; do
    receptor_name=$(echo $dir | sed 's/-GORE4-HADDOCK//')
    for pdb in ${dir}/cluster*_1.pdb; do
        base=$(basename $pdb .pdb)
        sample="${receptor_name}_${base}"

        rec="inputs/rec_${sample}.pdb"
        lig="inputs/lig_${sample}.pdb"

        awk '/^ATOM/ && substr($0,22,1)=="A"' $pdb > $rec; echo -e "TER\nEND" >> $rec
        awk '/^ATOM/ && substr($0,22,1)=="B"' $pdb > $lig; echo -e "TER\nEND" >> $lig

        # triad_1..4 = 0 para usar fallback global do nextflow.config
        echo "${sample},${rec},${lig},0,0,0,0" >> samplesheet_all.csv
    done
done

# Verifica: deve listar ~68 linhas de dados
wc -l samplesheet_all.csv

# Executa em paralelo no cluster (Nextflow gerencia a fila)
nextflow run main.nf \
    --input samplesheet_all.csv \
    --outdir ~/gromacs/results/ALL \
    --pH 8.2 \
    --time_ns 100 \
    --gmx_cmd gmx_mpi \
    --mpi_cmd "mpirun --use-hwthread-cpus -np 1" \
    -profile slurm,conda
```

### Exemplo 4 — Re-executar apenas análises (DM já concluída)

```bash
# Se a trajetória já existe e você quer refazer análises/figuras:
bash bin/run_analyses.sh ~/gromacs/results/xp273-gore4-c1

python3 bin/plot_results.py \
    --analise-dir ~/gromacs/results/xp273-gore4-c1/analise \
    --titulo "XP273-GORE4 cluster1 - 100 ns pH 8.2" \
    --window-ns 5
```

### Exemplo 5 — Retomar produção interrompida manualmente

```bash
cd ~/gromacs/results/xp273-gore4-c1/prod

# gmx_mpi com checkpoint
gmx_mpi mdrun -deffnm md -cpi md.cpt \
    -nb gpu -pme gpu -bonded gpu \
    -gpu_id 0 -ntomp 8 -pin on
```

---

## 7. Parâmetros

Todos os parâmetros têm valores padrão e podem ser sobrescritos na linha de comando
com `--nome_parametro valor`.

| Parâmetro       | Padrão               | Descrição                                                        |
|-----------------|----------------------|------------------------------------------------------------------|
| `--input`       | *obrigatório*        | Caminho para o samplesheet CSV                                   |
| `--outdir`      | `~/gromacs/results`  | Diretório de saída                                               |
| `--pH`          | `8.2`                | pH da simulação (define protonação da His)                       |
| `--time_ns`     | `200`                | Tempo de produção em nanossegundos                               |
| `--forcefield`  | `amber99sb-ildn`     | Force field GROMACS                                              |
| `--water`       | `tip3p`              | Modelo de água                                                   |
| `--cation`      | `K`                  | Cátion a adicionar (`K` para insetos, `NA` para mamíferos)       |
| `--nacl_conc`   | `0.10`               | Concentração de sal em mol/L                                     |
| `--temperature` | `300`                | Temperatura em Kelvin                                            |
| `--box_dist`    | `2.0`                | Distância mínima proteína–borda da caixa (nm)                    |
| `--box_type`    | `cubic`              | Tipo de caixa (`cubic`, `dodecahedron`, `triclinic`)             |
| `--gmx_cmd`     | `gmx_mpi`            | Binário do GROMACS (`gmx` ou `gmx_mpi`)                         |
| `--mpi_cmd`     | `mpirun -np 1`       | Prefixo MPI                                                      |
| `--use_gpu`     | `true`               | Habilita aceleração GPU no mdrun                                 |
| `--gpu_id`      | `0`                  | ID da GPU a usar                                                 |
| `--ntomp`       | `8`                  | Número de threads OpenMP                                         |
| `--maxwarn`     | `2`                  | Máximo de warnings aceitos pelo grompp                           |
| `--triad_1`     | `0`                  | Resíduo 1 para distâncias (fallback global; preferir samplesheet)|
| `--triad_2`     | `0`                  | Resíduo 2 para distâncias                                        |
| `--triad_3`     | `0`                  | Resíduo 3 para distâncias                                        |
| `--triad_4`     | `0`                  | Resíduo 4 para distâncias (bolsão S1)                            |

**Protonação da Histidina por pH:**

| pH              | Forma   | Descrição                             |
|-----------------|---------|---------------------------------------|
| < 6,5           | `HIP`   | Ambos os nitrogênios protonados       |
| 6,5 – 8,0       | `HID`   | H no Nδ (delta)                       |
| > 8,0 (padrão)  | `HIE`   | H no Nε (epsilon) — pH 8,2 usa HIE   |

**Resíduos catalíticos de GORE4 (tripsina *Spodoptera*):**

| triad_N | Resíduo | Papel                   |
|---------|---------|-------------------------|
| triad_1 | Tyr83   | Sítio ativo periférico  |
| triad_2 | Asp132  | Tríade catalítica       |
| triad_3 | Ser234  | Ser nucleofílica        |
| triad_4 | Ile229  | Bolsão S1 (especificidade) |

---

## 8. Saídas detalhadas

Para cada `sample_id` do samplesheet, o pipeline gera:

```
~/gromacs/results/<sample_id>/
│
├── prep/
│   ├── complexo.pdb          # Complexo após pré-processamento
│   ├── receptor_fixed.pdb    # Receptor com CYX e HIS corrigidos
│   └── ligante_fixed.pdb     # Ligante renumerado, chain B
│
├── topo/
│   ├── complexo.gro          # Estrutura GROMACS
│   ├── topol.top             # Topologia principal
│   ├── posre.itp             # Restrições de posição
│   └── topol_Protein_*.itp   # ITPs por cadeia
│
├── box/
│   └── ions.gro              # Sistema solvável + íons
│
├── em/
│   └── em.gro                # Estrutura minimizada
│
├── nvt/
│   ├── nvt.gro               # Estrutura após NVT
│   └── nvt.cpt               # Checkpoint NVT
│
├── npt/
│   ├── npt.gro               # Estrutura após NPT
│   └── npt.cpt               # Checkpoint NPT
│
├── prod/
│   ├── md.tpr                # Arquivo de run (binário)
│   ├── md.xtc                # Trajetória bruta
│   ├── md_fit.xtc            # Trajetória pós-processada (alinhada)
│   ├── md.gro                # Estrutura final
│   └── md.cpt                # Checkpoint (para retomar)
│
└── analise/
    ├── lig.ndx               # Arquivo de índice (grupos Ligante/Receptor)
    ├── rmsd_backbone.xvg     # RMSD do backbone do receptor (nm vs tempo)
    ├── rmsd_ligante.xvg      # RMSD do ligante (nm vs tempo)
    ├── rmsf_residuos.xvg     # RMSF por resíduo (flexibilidade)
    ├── gyrate.xvg            # Raio de giro da proteína (nm vs tempo)
    ├── numcont.xvg           # N.° de contatos receptor–ligante < 0,4 nm
    ├── hbond.xvg             # N.° de pontes de H receptor–ligante
    ├── sasa_receptor.xvg     # SASA do receptor ao longo da trajetória
    ├── sasa_ligante.xvg      # SASA do ligante ao longo da trajetória
    ├── dist_triad1.xvg       # Distância ligante → resíduo triad_1
    ├── dist_triad2.xvg       # Distância ligante → resíduo triad_2
    ├── dist_triad3.xvg       # Distância ligante → resíduo triad_3
    ├── dist_triad4.xvg       # Distância ligante → resíduo triad_4
    ├── painel_completo.png   # Figura completa com todas as métricas ← principal
    ├── rmsd_bb.png
    ├── rmsd_lig.png
    ├── rmsf.png
    ├── rg.png
    ├── ncont.png
    └── hbond.png
```

---

## 9. Solução de problemas

### GROMACS não encontrado

```
ERRO: gmx ou gmx_mpi não encontrado
```

```bash
# Verifique se o ambiente está ativo
mamba activate md-gromacs
which gmx_mpi

# Ou carregue o módulo do cluster
module load gromacs
```

### pdb2gmx falha com "residue not found"

Normalmente indica resíduos não padrão no PDB. Verifique:

```bash
grep "^ATOM\|^HETATM" receptor.pdb | awk '{print substr($0,18,3)}' | sort -u
```

Se houver `MSE`, `SEP`, `TPO` ou outros resíduos modificados, eles precisam
ser convertidos para os equivalentes padrão antes de entrar no pipeline.

### Spike no RMSD (pico abrupto na trajetória)

O passo `POSTPROCESS` já aplica `nojump` antes de `mol -center` para
evitar esse problema. Se persistir, refaça manualmente com `-pbc cluster`:

```bash
cd ~/gromacs/results/<sample_id>/prod

printf '1\n0\n' | gmx trjconv -s md.tpr -f md.xtc \
    -o md_cluster.xtc -pbc cluster -ur compact -center

printf '4\n0\n' | gmx trjconv -s md.tpr -f md_cluster.xtc \
    -o md_fit2.xtc -fit rot+trans
```

### Produção interrompida antes de completar N ns

```bash
# Verifique quantos ns foram simulados
cd ~/gromacs/results/<sample_id>/prod
gmx check -f md.xtc 2>&1 | grep "Last frame"

# Retome com checkpoint (via Nextflow — recomendado)
nextflow run main.nf --input assets/samplesheet.csv -profile conda,local -resume

# Ou retome manualmente
gmx_mpi mdrun -deffnm md -cpi md.cpt \
    -nb gpu -pme gpu -bonded gpu -gpu_id 0 -ntomp 8 -pin on
```

### cadeia B não detectada nas análises

```
ERRO: cadeia B não detectada em complexo.pdb
```

O `prepare_complex.py` não encontrou átomos com `chain == B` no PDB do ligante.
Verifique se o arquivo de entrada tem o campo de cadeia preenchido:

```bash
awk '{print substr($0,22,1)}' ligante.pdb | sort -u
```

Se estiver vazio, o `prepare_complex.py` atribui automaticamente chain B.
Certifique-se de que o arquivo começa com `ATOM` (não `HETATM`).

### Artefatos de PBC (complexo atravessa a caixa)

A caixa cúbica com 2,0 nm de margem é o padrão atual. Se ainda houver
artefatos (complexo muito elongado), aumente via linha de comando:

```bash
nextflow run main.nf --input samplesheet.csv --box_dist 2.5 -profile conda,local
```

---

## 10. Notas técnicas

Correções específicas implementadas e validadas:

| Problema                          | Solução                                                                |
|-----------------------------------|------------------------------------------------------------------------|
| HIS com protonação errada         | Forma automática por pH: HIP/HID/HIE conforme pKa da His (~6)         |
| CYS em pontes dissulfeto          | Auto-detecção por distância SG–SG < 2,5 Å → renomeada para CYX        |
| Ligante transladado para centro   | Pose dockada **preservada** (não transloca para o centróide)           |
| Conflito de numeração de resíduos | Ligante renumerado para vir após o receptor                            |
| Chain IDs ausentes                | Cadeia A=receptor, B=ligante, TER entre as cadeias                    |
| Spike no RMSD por PBC             | `nojump` aplicado antes de `mol -center`                               |
| PME GPU em minimização            | Desabilitado para o integrador `steep` (não suportado)                 |
| tc-grps incorreto                 | Usa `Protein Non-Protein` (peptídeo já pertence ao grupo Protein)      |
| Propagação de meta entre módulos  | Dados de resíduos (triad_1..4) injetados via canal separado (`ch_triad_input`), não no mapa `meta` global |

---

## Citação / Referência

Se usar este pipeline em publicações, cite o GROMACS e os force fields:

- **GROMACS**: Abraham et al., *SoftwareX* 1–2 (2015) 19–25
- **AMBER99sb-ILDN**: Lindorff-Larsen et al., *Proteins* 78 (2010) 1950–1958
- **TIP3P**: Jorgensen et al., *J. Chem. Phys.* 79 (1983) 926–935
- **Nextflow**: Di Tommaso et al., *Nat. Biotechnol.* 35 (2017) 316–319
