#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { PREPARE_PH       } from './modules/local/prepare_ph/main.nf'
include { PREPARE_COMPLEX  } from './modules/local/prepare_complex/main.nf'
include { TOPOLOGY         } from './modules/local/topology/main.nf'
include { BOX_SOLVATE_IONS } from './modules/local/box_solvate_ions/main.nf'
include { MINIMIZATION     } from './modules/local/minimization/main.nf'
include { NVT              } from './modules/local/nvt/main.nf'
include { NPT              } from './modules/local/npt/main.nf'
include { PRODUCTION       } from './modules/local/production/main.nf'
include { POSTPROCESS      } from './modules/local/postprocess/main.nf'
include { ANALYSES         } from './modules/local/analyses/main.nf'
include { PLOT             } from './modules/local/plot/main.nf'
include { STABILITY_FILTER } from './modules/local/stability_filter/main.nf'
include { CLUSTERING       } from './modules/local/clustering/main.nf'
include { MMGBSA_ROBUST    } from './modules/local/mmgbsa_robust/main.nf'
include { MMGBSA_INTERPRET } from './modules/local/mmgbsa_interpret/main.nf'

workflow {
    if (!params.input) {
        error "Informe o samplesheet: --input samplesheet.csv"
    }

    Channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true, strip: true)
        .map { row ->
            def meta     = [id: row.sample_id]
            def receptor = file(row.receptor, checkIfExists: true)
            def ligand   = file(row.ligand,   checkIfExists: true)
            tuple(meta, receptor, ligand)
        }
        .set { ch_input }

    PREPARE_PH(ch_input)
    PREPARE_COMPLEX(PREPARE_PH.out.protonated)
    TOPOLOGY(PREPARE_COMPLEX.out.complexo)
    BOX_SOLVATE_IONS(TOPOLOGY.out.topology)
    MINIMIZATION(BOX_SOLVATE_IONS.out.system)
    NVT(MINIMIZATION.out.system)
    NPT(NVT.out.system)
    PRODUCTION(NPT.out.system)
    POSTPROCESS(PRODUCTION.out.traj)

    // Junta complexo.pdb (preparação) com trajetória pós-processada
    ch_analyses = PREPARE_COMPLEX.out.complexo
        .join(POSTPROCESS.out.fit, by: [0])

    ANALYSES(ch_analyses)
    PLOT(ANALYSES.out.xvg)

    // ── Pipeline MM-GBSA robusto ─────────────────────────────────────────────
    // Etapa 1: junta POSTPROCESS + ANALYSES e filtra a fase estável
    ch_mmgbsa_prep = POSTPROCESS.out.fit
        .join(ANALYSES.out.xvg, by: [0])
        .map { meta, tpr, xtc, xvgs, ndx -> tuple(meta, tpr, xtc, xvgs, ndx) }

    STABILITY_FILTER(ch_mmgbsa_prep)

    // Etapa 2: clustering conformacional + subsampling
    CLUSTERING(STABILITY_FILTER.out.stable)

    // Etapa 3: validação + gmx_MMPBSA + decomposição por resíduo
    MMGBSA_ROBUST(CLUSTERING.out.for_mmgbsa)

    // Etapa 4: interpretação automática + painel de figuras
    ch_interpret = MMGBSA_ROBUST.out.results
        .join(STABILITY_FILTER.out.report, by: [0])
        .map { meta, csv, dat, decomp, stab -> tuple(meta, csv, dat, decomp, stab) }
    MMGBSA_INTERPRET(ch_interpret)
}
