process ANALYSES_SASA {
    tag "${meta.id}"
    label 'process_medium'

    publishDir { "${params.outdir}/${meta.id}/analise" }, mode: 'copy'

    input:
    tuple val(meta), path(md_tpr), path(md_fit_xtc), path(lig_ndx)

    output:
    tuple val(meta), path("sasa_protein.xvg"), path("sasa_ligante.xvg"), emit: sasa

    script:
    """
    echo "=== ANALYSES_SASA: ${meta.id} ===" >&2

    # SASA total do receptor ao longo da trajetória
    ${params.gmx_cmd} sasa \\
        -s ${md_tpr} \\
        -f ${md_fit_xtc} \\
        -n ${lig_ndx} \\
        -surface 'Protein' \\
        -output  'Protein' \\
        -o sasa_protein.xvg \\
        -tu ns

    # SASA do ligante no complexo (valores baixos = peptídeo enterrado)
    ${params.gmx_cmd} sasa \\
        -s ${md_tpr} \\
        -f ${md_fit_xtc} \\
        -n ${lig_ndx} \\
        -surface 'Ligante' \\
        -output  'Ligante' \\
        -o sasa_ligante.xvg \\
        -tu ns

    echo "[OK] SASA concluída para ${meta.id}" >&2
    """
}
