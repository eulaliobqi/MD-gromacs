process MMGBSA {
    tag "${meta.id}"
    label 'process_medium'
    errorStrategy 'ignore'

    publishDir { "${params.outdir}/${meta.id}/mmgbsa" }, mode: 'copy'

    input:
    tuple val(meta), path(md_tpr), path(md_xtc), path(lig_ndx), path(xvg_files)

    output:
    tuple val(meta), path("mmgbsa_results.csv"), path("FINAL_RESULTS_MMGBSA.dat"), emit: results
    tuple val(meta), path("painel_final.png"), emit: panel

    script:
    def sys_name     = meta.id
    def total_frames = (params.time_ns as int) * 100
    def interval     = Math.max(1, total_frames.intdiv(200))
    def titulo       = "${meta.id} — DM ${params.time_ns} ns @ pH ${params.pH} [MM-GBSA]"
    """
    cat > mmgbsa.in << 'MDP_EOF'
&general
sys_name="${sys_name}",
startframe=1,
endframe=${total_frames},
interval=${interval},
verbose=0,
/
&gb
igb=2,
saltcon=${params.nacl_conc},
/
MDP_EOF

    mamba run -n mmgbsa-env gmx_MMPBSA -O \\
        -i mmgbsa.in \\
        -cs ${md_tpr} \\
        -ct ${md_xtc} \\
        -ci ${lig_ndx} \\
        -cg Receptor Ligante \\
        -o FINAL_RESULTS_MMGBSA.dat \\
        -eo mmgbsa_results.csv \\
        -nogui 2>&1 | tee mmgbsa.log

    plot_results.py \\
        --analise-dir . \\
        --titulo "${titulo}" \\
        --mmgbsa-csv mmgbsa_results.csv \\
        --output painel_final.png
    """
}
