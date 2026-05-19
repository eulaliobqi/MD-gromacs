# Pipeline MD-GROMACS (Nextflow)

Pipeline automatizado para Dinâmica Molecular de complexos proteína-peptídeo via GROMACS.
Desenvolvido e validado com tripsina + inibidores peptídicos (DN1973+GORE4/LALAK).

## Estrutura do Pipeline

```
PREPARE_COMPLEX → TOPOLOGY → BOX_SOLVATE_IONS → MINIMIZATION → NVT → NPT → PRODUCTION → POSTPROCESS → ANALYSES → PLOT
```

| Etapa            | Descrição                                              | Recurso     |
|------------------|--------------------------------------------------------|-------------|
| PREPARE_COMPLEX  | CYX, HIS por pH, chain A/B, renumeração do ligante    | CPU low     |
| TOPOLOGY         | `pdb2gmx` amber99sb-ildn + TIP3P                      | CPU medium  |
| BOX_SOLVATE_IONS | Dodecaedro 1.2 nm, TIP3P, NaCl 0.15 M                 | CPU medium  |
| MINIMIZATION     | Steep 50k steps (sem PME GPU)                          | GPU         |
| NVT              | 100 ps, 300 K, com PR                                  | GPU         |
| NPT              | 100 ps, 1 bar, com PR                                  | GPU         |
| PRODUCTION       | **200 ns** NPT, Parrinello-Rahman                      | GPU (long)  |
| POSTPROCESS      | nojump → pbc mol center → fit backbone                 | CPU medium  |
| ANALYSES         | RMSD bb, RMSD lig, RMSF, Rg, contatos, H-bonds        | CPU medium  |
| PLOT             | Painel 3×2 + PNGs individuais                          | CPU low     |

## Pré-requisitos

```bash
mamba env create -f environment.yml
mamba activate md-gromacs
```

Ou instale manualmente:
```bash
mamba install -c conda-forge gromacs python>=3.10 numpy matplotlib mdtraj nextflow
```

## Uso Rápido (caso teste receptor_ph82 + ligante_ph82)

```bash
# 1. Crie o samplesheet
cat > samplesheet.csv << 'EOF'
sample_id,receptor,ligand
receptor_ph82_test,receptor_ph82.pdb,ligante_ph82.pdb
EOF

# 2. Execute localmente
nextflow run main.nf \
    --input samplesheet.csv \
    --pH 8.2 \
    --time_ns 200 \
    -profile local,conda

# 3. Retomar execução interrompida
nextflow run main.nf --input samplesheet.csv -resume
```

## Múltiplos Complexos em Paralelo

```bash
# samplesheet com vários pares receptor/ligante
cat > samplesheet.csv << 'EOF'
sample_id,receptor,ligand
XP273_c1,XP273-GORE4-HADDOCK/receptor.pdb,XP273-GORE4-HADDOCK/cluster1_1.pdb
DN773_c1,DN773-GORE4-HADDOCK/receptor.pdb,DN773-GORE4-HADDOCK/cluster1_1.pdb
EOF

nextflow run main.nf --input samplesheet.csv --pH 8.2 --time_ns 200 -profile slurm,conda
```

> **Nota HADDOCK**: os PDBs dos clusters HADDOCK contêm o complexo completo (cadeia A=receptor, B=ligante). Se precisar separar, use:
> ```bash
> awk '/^ATOM/ && substr($0,22,1)=="A"' cluster1_1.pdb > receptor.pdb; echo "TER" >> receptor.pdb
> awk '/^ATOM/ && substr($0,22,1)=="B"' cluster1_1.pdb > ligante.pdb;  echo "TER" >> ligante.pdb
> ```

## Parâmetros

| Parâmetro      | Default            | Descrição                                     |
|----------------|--------------------|-----------------------------------------------|
| `--input`      | obrigatório        | CSV com colunas: sample_id, receptor, ligand  |
| `--outdir`     | `~/gromacs/results`| Diretório de saída                            |
| `--pH`         | `8.2`              | pH (afeta protonação His: <6.5=HIP, 6.5-8=HID, >8=HIE) |
| `--time_ns`    | `200`              | Tempo de produção em ns                       |
| `--forcefield` | `amber99sb-ildn`   | Force field GROMACS                           |
| `--water`      | `tip3p`            | Modelo de água                                |
| `--nacl_conc`  | `0.15`             | Concentração NaCl (M)                         |
| `--temperature`| `300`              | Temperatura (K)                               |
| `--gmx_cmd`    | `gmx`              | Binário GROMACS (`gmx` ou `gmx_mpi`)          |
| `--mpi_cmd`    | `''`               | Prefixo MPI (ex: `mpirun --use-hwthread-cpus -np 1`) |
| `--use_gpu`    | `true`             | Habilita flags de GPU no mdrun                |
| `--gpu_id`     | `0`                | ID da GPU                                     |
| `--ntomp`      | `8`                | Threads OpenMP                                |

## Execução no Servidor (SLURM)

```bash
# Configura filas em conf/slurm.config (ajuste 'cpu' e 'gpu' para nomes do seu cluster)
nextflow run main.nf \
    --input samplesheet.csv \
    --pH 8.2 \
    --time_ns 200 \
    --gmx_cmd gmx_mpi \
    --mpi_cmd "mpirun --use-hwthread-cpus -np 1" \
    -profile slurm,conda \
    -resume
```

## Re-executar Apenas as Análises

```bash
# Se a DM já terminou e você só quer refazer análises/figuras:
bash bin/run_analyses.sh ~/gromacs/results/receptor_ph82_test
python3 bin/plot_results.py --analise-dir ~/gromacs/results/receptor_ph82_test/analise \
    --titulo "receptor_ph82_test - DM 200ns pH 8.2"
```

## Saída por Amostra

```
results/<sample_id>/
  prep/    complexo.pdb, receptor_fixed.pdb, ligante_fixed.pdb
  topo/    complexo.gro, topol.top, *.itp
  box/     ions.gro, topol.top
  em/      em.gro
  nvt/     nvt.gro, nvt.cpt
  npt/     npt.gro, npt.cpt
  prod/    md.tpr, md.xtc, md.gro, md.cpt, md_fit.xtc
  analise/ *.xvg, lig.ndx, painel_completo.png, *.png
```

## Correções Aplicadas (validadas em DN1973+GORE4)

- **HIS por pH**: HIP <6.5 / HID 6.5–8 / HIE >8 (pH 8.2 → HIE)
- **CYX dissulfeto**: auto-detecção SG-SG <2.5 Å, sem hardcode
- **Pose docada preservada**: ligante NÃO é transladado para o centróide
- **Numeração do ligante**: sempre renumerado após o receptor (evita conflito pdb2gmx)
- **Chain IDs**: A=receptor, B=ligante, com TER entre cadeias
- **PBC spike**: `nojump` antes de `mol -center` elimina saltos atômicos no RMSD
- **PME GPU**: desabilitado para minimização (`steep` não suporta)
- **tc-grps**: usa `Protein Non-Protein` (peptídeo ligante já é parte do grupo Protein)

## Retomar Produção Interrompida (manual)

```bash
cd ~/gromacs/results/<sample_id>/prod
gmx_mpi mdrun -deffnm md -cpi md.cpt -nb gpu -pme gpu -bonded gpu \
    -gpu_id 0 -ntomp 8 -pin on
```

Ou simplesmente use `-resume` no Nextflow — o pipeline detecta `md.cpt` e continua.
