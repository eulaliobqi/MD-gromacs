process BOX_SOLVATE_IONS {
    tag "${meta.id}"
    label 'process_medium'

    publishDir { "${params.outdir}/${meta.id}/box" }, mode: 'copy'

    input:
    tuple val(meta), path(gro), path(top), path(itps)

    output:
    tuple val(meta), path("ions.gro"), path("topol.top"), path("*.itp"), emit: system

    script:
    """
    # Copia topology (solvate e genion modificam o arquivo)
    cp ${top} topol.top

    # Caixa dodecaédrica com margem de ${params.box_dist} nm
    ${params.gmx_cmd} editconf \\
        -f ${gro} -o box.gro \\
        -c -d ${params.box_dist} -bt ${params.box_type}

    # Solvatação TIP3P
    ${params.gmx_cmd} solvate \\
        -cp box.gro -cs spc216.gro \\
        -p topol.top -o solv.gro

    # MDP mínimo para grompp (apenas para gerar ions.tpr)
    cat > ions.mdp << 'MDP_EOF'
integrator    = steep
nsteps        = 0
cutoff-scheme = Verlet
MDP_EOF

    ${params.gmx_cmd} grompp \\
        -f ions.mdp -c solv.gro \\
        -p topol.top -o ions.tpr \\
        -maxwarn ${params.maxwarn}

    # Adiciona NaCl ${params.nacl_conc} M + neutralização
    echo "SOL" | ${params.gmx_cmd} genion \\
        -s ions.tpr -o ions.gro \\
        -p topol.top -pname NA -nname CL \\
        -neutral -conc ${params.nacl_conc}
    """
}
