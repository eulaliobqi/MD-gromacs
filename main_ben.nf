#!/usr/bin/env nextflow
// Pipeline de DM para benzamidina (BEN) contra tripsinas Spodoptera.
// Usa TOPOLOGY_SMALL_MOLECULE (pdb2gmx + ACPYPE/GAFF2) em vez de TOPOLOGY.
// Receptores já protonados em pH 8.2 → PREPARE_PH omitido.
// Downstream (BOX_SOLVATE_IONS em diante) idêntico ao main.nf.
nextflow.enable.dsl = 2

include { PREPARE_COMPLEX        } from './modules/local/prepare_complex/main.nf'
include { TOPOLOGY_SMALL_MOLECULE} from './modules/local/topology_small_molecule/main.nf'
include { BOX_SOLVATE_IONS       } from './modules/local/box_solvate_ions/main.nf'
include { MINIMIZATION           } from './modules/local/minimization/main.nf'
include { NVT                    } from './modules/local/nvt/main.nf'
include { NPT                    } from './modules/local/npt/main.nf'
include { PRODUCTION             } from './modules/local/production/main.nf'
include { POSTPROCESS            } from './modules/local/postprocess/main.nf'
include { ANALYSES_BEN           } from './modules/local/analyses_ben/main.nf'
include { PLOT                   } from './modules/local/plot/main.nf'
include { ANALYSES_TRIAD         } from './modules/local/analyses_triad/main.nf'

workflow {
    if (!params.input) {
        error "Informe o samplesheet: --input samplesheet.csv"
    }

    // Canal principal: receptor (já protonado pH 8.2) + pose BEN (do Vina)
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

    // Canal de resíduos catalíticos — injetado apenas em ANALYSES_TRIAD
    Channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true, strip: true)
        .map { row ->
            tuple(
                [id: row.sample_id],
                row.triad_1 ?: params.triad_1,
                row.triad_2 ?: params.triad_2,
                row.triad_3 ?: params.triad_3,
                row.triad_4 ?: params.triad_4
            )
        }
        .set { ch_triad_params }

    // Receptor já protonado → vai direto para PREPARE_COMPLEX (sem PREPARE_PH)
    PREPARE_COMPLEX(ch_input)

    // TOPOLOGY_SMALL_MOLECULE: pdb2gmx (receptor) + ACPYPE GAFF2 (BEN) + merge
    TOPOLOGY_SMALL_MOLECULE(PREPARE_COMPLEX.out.complexo)

    BOX_SOLVATE_IONS(TOPOLOGY_SMALL_MOLECULE.out.topology)
    MINIMIZATION(BOX_SOLVATE_IONS.out.system)
    NVT(MINIMIZATION.out.system)
    NPT(NVT.out.system)
    PRODUCTION(NPT.out.system)
    POSTPROCESS(PRODUCTION.out.traj)

    // ANALYSES_BEN detecta BEN via HETATM cadeia B (não ATOM como ANALYSES padrão)
    ch_analyses = PREPARE_COMPLEX.out.complexo
        .join(POSTPROCESS.out.fit, by: [0])

    ANALYSES_BEN(ch_analyses)

    // ANALYSES_TRIAD recebe lig.ndx de ANALYSES_BEN — funciona sem alteração
    ch_extended = POSTPROCESS.out.fit
        .join(ANALYSES_BEN.out.xvg, by: [0])
        .map { meta, tpr, xtc, xvgs, ndx -> tuple(meta, tpr, xtc, ndx) }

    ch_triad_input = ch_extended
        .join(ch_triad_params, by: [0])
        .map { meta, tpr, xtc, ndx, t1, t2, t3, t4 ->
            tuple(meta + [triad_1: t1, triad_2: t2, triad_3: t3, triad_4: t4], tpr, xtc, ndx)
        }

    ANALYSES_TRIAD(ch_triad_input)

    ch_plot_input = ANALYSES_BEN.out.xvg
        .join(ANALYSES_TRIAD.out.triad, by: [0])
        .join(ANALYSES_TRIAD.out.info,  by: [0])
        .map { meta, xvgs, ndx, d1, d2, d3, d4, s1, s2, s3, s4, info ->
            def all_files = (xvgs instanceof List ? xvgs : [xvgs]) +
                            [d1, d2, d3, d4, s1, s2, s3, s4, info]
            tuple(meta, all_files, ndx)
        }

    PLOT(ch_plot_input)
}
