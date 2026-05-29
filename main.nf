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
include { ANALYSES_SASA   } from './modules/local/analyses_sasa/main.nf'
include { ANALYSES_TRIAD  } from './modules/local/analyses_triad/main.nf'

workflow {
    if (!params.input) {
        error "Informe o samplesheet: --input samplesheet.csv"
    }

    // Canal principal — meta sem campos de tríade para preservar cache upstream
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

    // Canal de resíduos de interesse — injetado apenas em ANALYSES_TRIAD
    Channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true, strip: true)
        .map { row ->
            tuple(
                [id: row.sample_id],
                row.triad_1 ?: params.triad_1,
                row.triad_2 ?: params.triad_2,
                row.triad_3 ?: params.triad_3
            )
        }
        .set { ch_triad_params }

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

    // ── Análises estendidas ────────────────────────────────────────────────────
    ch_extended = POSTPROCESS.out.fit
        .join(ANALYSES.out.xvg, by: [0])
        .map { meta, tpr, xtc, xvgs, ndx -> tuple(meta, tpr, xtc, ndx) }

    ANALYSES_SASA(ch_extended)

    // Injeta resíduos de interesse no canal só para ANALYSES_TRIAD
    ch_triad_input = ch_extended
        .join(ch_triad_params, by: [0])
        .map { meta, tpr, xtc, ndx, t1, t2, t3 ->
            tuple(meta + [triad_1: t1, triad_2: t2, triad_3: t3], tpr, xtc, ndx)
        }

    ANALYSES_TRIAD(ch_triad_input)
}
