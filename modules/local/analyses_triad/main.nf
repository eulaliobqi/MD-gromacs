process ANALYSES_TRIAD {
    tag "${meta.id}"
    label 'process_low'

    publishDir { "${params.outdir}/${meta.id}/analise" }, mode: 'copy'

    input:
    tuple val(meta), path(md_tpr), path(md_fit_xtc), path(lig_ndx)

    output:
    tuple val(meta), path("dist_ser.xvg"), path("dist_his.xvg"), path("dist_asp.xvg"), emit: triad
    tuple val(meta), path("triad_info.txt"), emit: info

    script:
    def r1 = meta.triad_1
    def r2 = meta.triad_2
    def r3 = meta.triad_3
    """
    echo "=== ANALYSES_TRIAD: ${meta.id} ===" >&2
    echo "Resíduos de interesse: ${r1} ${r2} ${r3}" >&2

    # Grava rótulos para o plot_results.py
    printf '${r1}\\n${r2}\\n${r3}\\n' > triad_info.txt

    # Conta grupos existentes no lig.ndx para determinar próximo índice disponível
    N_CURR=\$(echo q | ${params.gmx_cmd} make_ndx \\
        -f ${md_tpr} -n ${lig_ndx} -o _tmp_count.ndx 2>&1 \\
        | grep -cE "^ *[0-9]+ +[A-Za-z]")
    rm -f _tmp_count.ndx
    RES1_IDX=\${N_CURR}
    RES2_IDX=\$((N_CURR + 1))
    RES3_IDX=\$((N_CURR + 2))

    # Adiciona grupos dos resíduos de interesse ao ndx
    ${params.gmx_cmd} make_ndx -f ${md_tpr} -n ${lig_ndx} -o triad.ndx << MNDX
r ${r1}
name \${RES1_IDX} Res1_Cat
r ${r2}
name \${RES2_IDX} Res2_Cat
r ${r3}
name \${RES3_IDX} Res3_Cat
q
MNDX

    # Distância mínima ligante–resíduo 1
    printf 'Ligante\\nRes1_Cat\\n' | ${params.gmx_cmd} mindist \\
        -s ${md_tpr} -f ${md_fit_xtc} \\
        -n triad.ndx \\
        -od dist_ser.xvg \\
        -tu ns

    # Distância mínima ligante–resíduo 2
    printf 'Ligante\\nRes2_Cat\\n' | ${params.gmx_cmd} mindist \\
        -s ${md_tpr} -f ${md_fit_xtc} \\
        -n triad.ndx \\
        -od dist_his.xvg \\
        -tu ns

    # Distância mínima ligante–resíduo 3
    printf 'Ligante\\nRes3_Cat\\n' | ${params.gmx_cmd} mindist \\
        -s ${md_tpr} -f ${md_fit_xtc} \\
        -n triad.ndx \\
        -od dist_asp.xvg \\
        -tu ns

    echo "[OK] Distâncias resíduos-chave concluídas para ${meta.id}" >&2
    """
}
