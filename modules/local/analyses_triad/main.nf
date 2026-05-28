process ANALYSES_TRIAD {
    tag "${meta.id}"
    label 'process_low'

    publishDir { "${params.outdir}/${meta.id}/analise" }, mode: 'copy'

    input:
    tuple val(meta), path(md_tpr), path(md_fit_xtc), path(lig_ndx)

    output:
    tuple val(meta), path("dist_ser.xvg"), path("dist_his.xvg"), path("dist_asp.xvg"), emit: triad

    script:
    def triad_ser = params.triad_ser
    def triad_his = params.triad_his
    def triad_asp = params.triad_asp
    """
    echo "=== ANALYSES_TRIAD: ${meta.id} ===" >&2
    echo "Tríade catalítica: Ser${triad_ser} His${triad_his} Asp${triad_asp}" >&2

    # Conta grupos existentes no lig.ndx para determinar próximo índice disponível
    N_CURR=\$(echo q | ${params.gmx_cmd} make_ndx \\
        -f ${md_tpr} -n ${lig_ndx} -o _tmp_count.ndx 2>&1 \\
        | grep -cE "^ *[0-9]+ +[A-Za-z]")
    rm -f _tmp_count.ndx
    SER_IDX=\${N_CURR}
    HIS_IDX=\$((N_CURR + 1))
    ASP_IDX=\$((N_CURR + 2))

    # Adiciona grupos da tríade ao ndx
    ${params.gmx_cmd} make_ndx -f ${md_tpr} -n ${lig_ndx} -o triad.ndx << MNDX
r ${triad_ser}
name \${SER_IDX} Ser_Cat
r ${triad_his}
name \${HIS_IDX} His_Cat
r ${triad_asp}
name \${ASP_IDX} Asp_Cat
q
MNDX

    # Distância mínima ligante–Serina catalítica
    printf 'Ligante\\nSer_Cat\\n' | ${params.gmx_cmd} mindist \\
        -s ${md_tpr} -f ${md_fit_xtc} \\
        -n triad.ndx \\
        -od dist_ser.xvg \\
        -tu ns

    # Distância mínima ligante–Histidina catalítica
    printf 'Ligante\\nHis_Cat\\n' | ${params.gmx_cmd} mindist \\
        -s ${md_tpr} -f ${md_fit_xtc} \\
        -n triad.ndx \\
        -od dist_his.xvg \\
        -tu ns

    # Distância mínima ligante–Aspartato catalítico
    printf 'Ligante\\nAsp_Cat\\n' | ${params.gmx_cmd} mindist \\
        -s ${md_tpr} -f ${md_fit_xtc} \\
        -n triad.ndx \\
        -od dist_asp.xvg \\
        -tu ns

    echo "[OK] Distâncias tríade concluídas para ${meta.id}" >&2
    """
}
