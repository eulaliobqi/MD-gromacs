# Pipeline de Dinâmica Molecular: Complexo Proteína-Peptídeo via GROMACS

Documentação técnica do pipeline `dm_pipeline.sh` e seus scripts auxiliares. O objetivo é automatizar todo o processo de simulação de dinâmica molecular (DM) de um complexo proteína-peptídeo, desde a preparação dos arquivos PDB até a geração de figuras de análise e o cálculo de energia de ligação (MM-GBSA).

---

## Visão Geral

O pipeline é composto por quatro scripts principais que se chamam em sequência:

```
dm_pipeline.sh
├── prepare_complex.py       ← Etapa 1: prepara o complexo
├── [GROMACS embutido]       ← Etapas 2–8: simulação
├── run_analyses.sh          ← Etapa 9: análises
│   └── [GROMACS embutido]
├── plot_results.py          ← Etapa 10: figuras
└── scripts/run_mmgbsa.sh   ← Opcional: energia MM-GBSA
```

Existe também uma versão Nextflow (`main.nf`) que executa os mesmos passos de forma paralela e reprodutível em ambientes de HPC/cloud.

---

## Como Usar

```bash
./dm_pipeline.sh \
  -r receptor.pdb \
  -l ligante.pdb \
  -p nome_projeto \
  -H 7.4 \
  -t 100
```

**Parâmetros principais:**

| Flag | Significado | Padrão |
|------|-------------|--------|
| `-r` | PDB do receptor (obrigatório) | — |
| `-l` | PDB do ligante peptídico (obrigatório) | — |
| `-p` | Nome do projeto | `dm_complex` |
| `-H` | pH da simulação | `7.4` |
| `-t` | Tempo de produção em ns | `100` |
| `-F` | Campo de força | `amber99sb-ildn` |
| `-W` | Modelo de água | `tip3p` |
| `-S` | Concentração de NaCl em M | `0.15` |
| `-B` | Diretório base de saída | `$HOME/gromacs/<projeto>` |
| `-G` | ID da GPU | `0` |
| `-T` | Threads OMP | `8` |

---

## Estrutura de Diretórios Criada

Após execução, a seguinte hierarquia é criada em `$BASE`:

```
$BASE/
├── prep/         ← arquivos PDB preparados
├── topo/         ← topologia GROMACS (.top, .itp, .gro)
├── box/          ← caixa, solvente e íons
├── em/           ← minimização de energia
├── nvt/          ← equilíbrio canônico (NVT)
├── npt/          ← equilíbrio isobárico-isotérmico (NPT)
├── prod/         ← simulação de produção + trajetórias
└── analise/      ← arquivos .xvg + figuras .png
```

---

## Etapas Detalhadas

### Etapa 0 — Configuração e Validação

**Script:** `dm_pipeline.sh` (linhas 37–115)

O script inicia parseando os argumentos via `getopts`. Em seguida:

1. Valida que `-r` e `-l` foram fornecidos.
2. Converte os caminhos para absolutos com `realpath`.
3. Calcula o número de passos de produção: `PROD_STEPS = PROD_NS × 500.000` (passo de integração `dt = 0.002 ps` → 500 passos por picossegundo).
4. Detecta o binário do GROMACS (`gmx_mpi` tem prioridade sobre `gmx`).
5. Define duas funções internas:
   - `mdrun()` — para equilíbrio e produção, com `PME` descarregado na GPU.
   - `em_mdrun()` — para minimização, sem PME na GPU (o integrador `steep` não suporta PME GPU).
6. Cria os diretórios `prep`, `topo`, `box`, `em`, `nvt`, `npt`, `prod`, `analise`.

---

### Etapa 1 — Preparação do Complexo (`prepare_complex.py`)

**Script:** `prepare_complex.py`

Converte dois PDBs separados (receptor + ligante) em um único `complexo.pdb` pronto para o GROMACS. As operações são:

#### 1a. Leitura dos átomos
Lê apenas linhas `ATOM` e `HETATM` de cada arquivo. Ignora cabeçalhos e registros `CONECT`.

#### 1b. Detecção de pontes dissulfeto
Percorre todos os átomos `SG` de resíduos `CYS`. Para cada par de cisteínas, calcula a distância euclidiana. Se a distância for menor que 2,5 Å (padrão), as cisteínas são marcadas como formando ponte dissulfeto.

#### 1c. Protonação da histidina por pH
Determina a forma da histidina com base no pH informado:
- pH < 6,5 → `HIP` (ambos os nitrogênios protonados, carga +1)
- 6,5 ≤ pH < 8,0 → `HID` (próton no N-delta, neutro)
- pH ≥ 8,0 → `HIE` (próton no N-epsilon, neutro)

#### 1d. Pré-processamento do receptor
- Cisteínas em pontes S-S: renomeadas de `CYS` para `CYX` e o hidrogênio `HG` é removido.
- Histidinas: renomeadas para a forma determinada pelo pH.

#### 1e. Renumeração e atribuição de cadeias
- Receptor → Cadeia A, numeração começa em 1.
- Ligante → Cadeia B, numeração continua após o último resíduo do receptor (evita conflito no `pdb2gmx`).
- A pose do ligante (coordenadas XYZ) é preservada integralmente, sem translação.
- Numeração serial atômica é contínua entre as duas cadeias.

#### 1f. Saídas geradas em `prep/`
- `receptor_fixed.pdb` — receptor corrigido isolado.
- `ligante_fixed.pdb` — ligante renumerado isolado.
- `complexo.pdb` — arquivo final com as duas cadeias, separadas por `TER`.

---

### Etapa 2 — Geração de Topologia (`pdb2gmx`)

**Script:** `dm_pipeline.sh`, diretório `topo/`

Executa `gmx pdb2gmx` sobre o `complexo.pdb`:

```
printf '0\n0\n0\n0\n' | gmx pdb2gmx \
    -f complexo.pdb \
    -o complexo.gro \
    -p topol.top \
    -i posre.itp \
    -ff amber99sb-ildn \
    -water tip3p \
    -ignh -ter -chainsep ter -merge no
```

- `-ignh`: ignora hidrogênios do PDB e os regenera com base no campo de força.
- `-ter`: processa terminais (N-terminal/C-terminal) interativamente (quatro `0` entram via `printf` para aceitar os terminais padrão das quatro cadeias presentes).
- `-merge no`: mantém as cadeias separadas na topologia.
- Gera `topol.top` (topologia principal) e `posre.itp` (restrições de posição).

A etapa é pulada se `complexo.gro` já existir (idempotência).

---

### Etapa 3 — Caixa, Solvatação e Íons

**Script:** `dm_pipeline.sh`, diretório `box/`

#### 3a. Caixa dodecaédrica
```bash
gmx editconf -f complexo.gro -o box.gro -c -d 1.2 -bt dodecahedron
```
Cria uma caixa dodecaédrica com pelo menos 1,2 nm de distância entre o complexo e as bordas. O dodecaedro minimiza o volume de solvente em comparação com caixas cúbicas.

#### 3b. Solvatação
```bash
gmx solvate -cp box.gro -cs spc216.gro -p topol.top -o solv.gro
```
Preenche a caixa com moléculas de água TIP3P (a geometria inicial vem do arquivo `spc216.gro` incluso no GROMACS). A topologia é atualizada automaticamente.

#### 3c. Adição de íons
Um arquivo `ions.mdp` mínimo é gerado (apenas para criar o `.tpr` necessário pelo `genion`). Em seguida:
```bash
gmx grompp -f ions.mdp -c solv.gro -p topol.top -o ions.tpr
echo "SOL" | gmx genion -s ions.tpr -o ions.gro \
    -pname NA -nname CL -neutral -conc 0.15
```
Substitui moléculas de água aleatórias por íons Na⁺ e Cl⁻ para neutralizar a carga do sistema e atingir a concentração fisiológica de NaCl (0,15 M por padrão).

Os arquivos de topologia (`.top` e `.itp`) são copiados para `em/`.

---

### Etapa 4 — Minimização de Energia (EM)

**Script:** `dm_pipeline.sh`, diretório `em/`

Parâmetros do `em.mdp`:
- `integrator = steep` — método de descida mais íngreme (steepest descent).
- `emtol = 1000.0` kJ·mol⁻¹·nm⁻¹ — critério de convergência da força máxima.
- `nsteps = 50000` — máximo de passos.
- PME para eletrostática longa distância com `rcoulomb = 1.2 nm`.

Executa via `em_mdrun` (sem PME GPU, pois `steep` não suporta). Verifica que `em.gro` foi gerado antes de prosseguir.

---

### Etapa 5 — Equilíbrio NVT (100 ps)

**Script:** `dm_pipeline.sh`, diretório `nvt/`

Ensemble canônico (N, V, T constantes) com restrições de posição (`-DPOSRES`) nos átomos pesados do complexo.

Parâmetros relevantes do `nvt.mdp`:
- `integrator = md`, `dt = 0.002 ps`, `nsteps = 50000` → 100 ps.
- `gen-vel = yes`, `gen-temp = 300 K` — velocidades geradas da distribuição de Maxwell-Boltzmann na semente 42.
- `tcoupl = V-rescale` — termostato para dois grupos: `Protein` e `Non-Protein` (solvente + íons), ambos a 300 K.
- `pcoupl = no` — sem barostato nesta etapa.
- `constraints = h-bonds` com `LINCS` — restringe ligações com hidrogênio, permitindo `dt = 2 fs`.

---

### Etapa 6 — Equilíbrio NPT (100 ps)

**Script:** `dm_pipeline.sh`, diretório `npt/`

Ensemble isobárico-isotérmico (N, P, T constantes), ainda com restrições de posição.

Diferenças em relação ao NVT:
- `continuation = yes` — velocidades lidas do checkpoint NVT (`nvt.cpt`).
- `gen-vel = no` — não regera velocidades.
- `pcoupl = C-rescale` — barostato isotrópico com `tau-p = 2.0 ps` e `ref-p = 1.0 bar`.
- `compressibility = 4.5e-5 bar⁻¹` — compressibilidade da água à temperatura ambiente.
- `refcoord_scaling = com` — escala as coordenadas de referência das restrições.

---

### Etapa 7 — Simulação de Produção

**Script:** `dm_pipeline.sh`, diretório `prod/`

A simulação livre (sem restrições de posição) pelo tempo definido em `-t`.

Parâmetros do `md.mdp`:
- `nsteps = PROD_STEPS` — calculado como `PROD_NS × 500.000`.
- `nstxout-compressed = 5000` → frame a cada 10 ps gravado na trajetória `.xtc`.
- `nstenergy = 5000`, `nstlog = 5000` — energia e log a cada 10 ps.
- `pcoupl = Parrinello-Rahman` — barostato mais rigoroso para produção.
- `comm-mode = Linear` com `nstcomm = 100` — remove translação do centro de massa do sistema.

Se já existir um `md.cpt` (checkpoint), a corrida é retomada automaticamente (`-cpi md.cpt`), permitindo recuperação após interrupção.

---

### Etapa 8 — Pós-processamento da Trajetória (PBC + Alinhamento)

**Scripts:** `dm_pipeline.sh` e `run_analyses.sh`, diretório `prod/`

Três passos sequenciais para corrigir artefatos de condições periódicas de contorno (PBC):

#### 8a. nojump
```bash
echo 0 | gmx trjconv -s md.tpr -f md.xtc -o md_nojump.xtc -pbc nojump
```
Evita "saltos" de moléculas que cruzam a borda da caixa. Processa o sistema inteiro (grupo `0 = System`).

#### 8b. Correção de PBC
```bash
echo -e '1\n0' | gmx trjconv -s md.tpr -f md_nojump.xtc \
    -o md_nopbc.xtc -pbc mol -center -ur compact
```
- Centraliza a proteína (grupo `1 = Protein`) na caixa.
- `-ur compact` — empacota moléculas de solvente quebradas na borda.

#### 8c. Alinhamento
```bash
echo -e '4\n0' | gmx trjconv -s md.tpr -f md_nopbc.xtc \
    -o md_fit.xtc -fit rot+trans
```
Alinha todos os frames pelo backbone (grupo `4 = Backbone`) removendo rotação e translação do corpo rígido. Esta trajetória (`md_fit.xtc`) é usada em todas as análises subsequentes.

---

### Etapa 9 — Análises (`run_analyses.sh`)

**Script:** `run_analyses.sh`

Todos os arquivos são gerados em `$BASE/analise/`. Links simbólicos apontam para `md.tpr` e `md_fit.xtc` em `prod/`.

#### 9a. Detecção automática do ligante
O script lê `complexo.pdb` e identifica o primeiro e último número de resíduo na Cadeia B usando `awk`. Esta faixa é usada para construir o índice.

#### 9b. Construção do índice de grupos (`make_ndx`)
Detecta quantos grupos padrão já existem no `.tpr` e cria dois grupos adicionais:
- **Ligante** — resíduos da Cadeia B.
- **Receptor** — todos os átomos da Proteína que não são Ligante.

O índice é salvo em `lig.ndx`.

#### 9c. RMSD do backbone
```bash
printf 'Backbone\nBackbone\n' | gmx rms -s sistema.tpr -f trajetoria.xtc \
    -n lig.ndx -o rmsd_backbone.xvg -tu ns
```
Calcula o RMSD de todos os frames em relação à estrutura de referência (frame 0). Mede a estabilidade global do complexo ao longo do tempo em nanossegundos.

#### 9d. RMSD do ligante
```bash
printf 'Ligante\nLigante\n' | gmx rms ...
```
Mesmo cálculo, mas restrito aos átomos do ligante peptídico. Indica se o peptídeo permanece no sítio de ligação.

#### 9e. RMSF por resíduo
```bash
printf 'Backbone\n' | gmx rmsf -s sistema.tpr -f trajetoria.xtc \
    -n lig.ndx -o rmsf_residuos.xvg -res -fit
```
Calcula a flutuação quadrática média (Root Mean Square Fluctuation) por resíduo, medindo regiões flexíveis vs. rígidas da proteína.

#### 9f. Raio de giro
```bash
printf 'Protein\n' | gmx gyrate -o gyrate.xvg -tu ns
```
Mede a compacidade global do complexo proteico ao longo da simulação.

#### 9g. Contatos receptor-ligante
```bash
printf 'Receptor\nLigante\n' | gmx mindist \
    -od mindist.xvg -on numcont.xvg -d 0.4 -tu ns
```
Conta o número de pares de átomos receptor-ligante com distância < 0,4 nm a cada frame. Indica a qualidade e estabilidade dos contatos de ligação.

#### 9h. Pontes de hidrogênio
```bash
printf 'Receptor\nLigante\n' | gmx hbond -num hbond.xvg -tu ns
```
Conta o número de pontes de hidrogênio entre receptor e ligante ao longo do tempo.

---

### Etapa 10 — Geração de Figuras (`plot_results.py`)

**Script:** `plot_results.py`

Lê os seis arquivos `.xvg` gerados na etapa anterior e produz:

#### Painel completo (`painel_completo.png`)
Grid 3×2 com os seguintes subgráficos:
1. **RMSD do Backbone** — linha azul (`navy`), com linha tracejada vermelha na média. Título inclui média ± desvio padrão.
2. **RMSD do Ligante** — linha laranja (`darkorange`).
3. **RMSF por Resíduo** — gráfico de barras verde (`seagreen`), com os 5 últimos resíduos destacados em vermelho (aproximação dos resíduos do ligante).
4. **Raio de Giro** — linha roxa (`purple`).
5. **Contatos Receptor-Ligante** — linha verde-azulado (`teal`).
6. **Pontes de Hidrogênio** — linha vermelha-índia (`indianred`).

Se um arquivo `.xvg` não for encontrado, o painel exibe `"(arquivo não encontrado)"` naquele subgráfico, sem travar a execução.

#### PNGs individuais
Um arquivo `.png` separado é gerado para cada métrica (`rmsd_bb.png`, `rmsd_lig.png`, `rmsf.png`, `rg.png`, `ncont.png`, `hbond.png`).

Todos os arquivos são salvos em `$BASE/analise/` com resolução de 150 DPI.

---

## Script Auxiliar: MM-GBSA (`scripts/run_mmgbsa.sh`)

**Uso:**
```bash
bash scripts/run_mmgbsa.sh <work_dir> <sys_name> <time_ns> <saltcon>
```

Calcula a energia livre de ligação pelo método MM-GBSA usando `gmx_MMPBSA`. Destinado a execução manual fora do Nextflow.

### Passo a Passo

#### 1. Preparação do diretório de saída
Cria `mmgbsa_<SYS_NAME>/` como irmão do diretório de trabalho Nextflow. Calcula o intervalo entre frames: `interval = max(1, TOTAL_FRAMES ÷ 200)`, mantendo no máximo 200 frames para análise.

#### 2. Vinculação de arquivos de entrada
Cria links simbólicos para `md.tpr`, `md_fit.xtc` e `lig.ndx` do diretório de trabalho Nextflow.

#### 3. Arquivo de configuração `mmgbsa.in`
Define os parâmetros para `gmx_MMPBSA`:
- `igb=2` — modelo GB de Onufriev-Bashford-Case (GBOBCII).
- `saltcon` — concentração iônica.
- `startframe/endframe` — intervalo completo da simulação.

#### 4. Wrapper `tleap` (correção de bug)
Cria um script Python em `bin_patch/tleap` que intercepta a chamada real ao `tleap`. Esse wrapper corrige um bug do `gmx_MMPBSA 1.6.5`: os índices de pontes dissulfeto (`SS bonds`) no complexo (`COM_OUT`) são gerados com offset errado (deslocados pelo tamanho do ligante). O fix consiste em copiar os índices corretos do receptor (`REC_OUT`) para o complexo, linha a linha via regex.

#### 5. Execução do `gmx_MMPBSA`
Executa dentro do ambiente conda `mmgbsa-env` (isolado do ambiente GROMACS principal para evitar conflitos com `cpptraj`). O `bin_patch/` é colocado no início do `PATH` para que o wrapper intercepte o `tleap` antes do binário real.

#### 6. Geração do painel de resultados
Se `FINAL_RESULTS_MMGBSA.dat` for criado com sucesso, chama `plot_results.py` com a flag `--mmgbsa-csv` para incluir os dados de energia livre no painel.

---

## Versão Nextflow (`main.nf`)

O arquivo `main.nf` executa o mesmo fluxo de forma paralela e reprodutível. Lê um `samplesheet.csv` com colunas `sample_id`, `receptor`, `ligand` e processa múltiplos sistemas simultaneamente.

**Ordem de execução dos processos:**

```
PREPARE_PH → PREPARE_COMPLEX → TOPOLOGY → BOX_SOLVATE_IONS
→ MINIMIZATION → NVT → NPT → PRODUCTION → POSTPROCESS
→ ANALYSES → PLOT
→ MMGBSA  (consome POSTPROCESS + ANALYSES em paralelo)
```

O canal `ch_mmgbsa` é construído juntando a trajetória pós-processada com os arquivos de análise:
```groovy
ch_mmgbsa = POSTPROCESS.out.fit
    .join(ANALYSES.out.xvg, by: [0])
    .map { meta, tpr, xtc, xvgs, ndx -> tuple(meta, tpr, xtc, ndx, xvgs) }
```

---

## Saídas Finais

| Arquivo | Localização | Descrição |
|---------|-------------|-----------|
| `complexo.pdb` | `prep/` | Complexo preparado para GROMACS |
| `topol.top` | `topo/` | Topologia AMBER99SB-ILDN |
| `md.tpr` | `prod/` | Arquivo de execução da simulação |
| `md_fit.xtc` | `prod/` | Trajetória pós-processada e alinhada |
| `rmsd_backbone.xvg` | `analise/` | RMSD do backbone vs. tempo |
| `rmsd_ligante.xvg` | `analise/` | RMSD do ligante vs. tempo |
| `rmsf_residuos.xvg` | `analise/` | RMSF por resíduo |
| `gyrate.xvg` | `analise/` | Raio de giro vs. tempo |
| `numcont.xvg` | `analise/` | Número de contatos vs. tempo |
| `hbond.xvg` | `analise/` | Pontes de hidrogênio vs. tempo |
| `painel_completo.png` | `analise/` | Painel 3×2 com todas as métricas |
| `FINAL_RESULTS_MMGBSA.dat` | `mmgbsa_*/` | Energia livre de ligação (MM-GBSA) |
