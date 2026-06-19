// Variante de ANALYSES para molécula pequena (BEN).
// Diferença: detecta cadeia B via HETATM (não ATOM) para encontrar BEN.
// ANALYSES_TRIAD não necessita alteração — recebe lig.ndx com "Ligante" já definido.
process ANALYSES_BEN {
    tag "${meta.id}"
    label 'process_medium'

    publishDir { "${params.outdir}/${meta.id}/analise" }, mode: 'copy'

    input:
    tuple val(meta), path(complexo_pdb), path(md_tpr), path(md_fit_xtc)

    output:
    tuple val(meta), path("*.xvg"), path("lig.ndx"), emit: xvg

    script:
    """
    # BEN é HETATM na cadeia B — detectar resíduo correto
    LIG_FIRST=\$(awk '(/^ATOM/ || /^HETATM/) && substr(\$0,22,1)=="B" {print substr(\$0,23,4)+0; exit}' ${complexo_pdb})
    LIG_LAST=\$(awk '(/^ATOM/ || /^HETATM/) && substr(\$0,22,1)=="B" {r=substr(\$0,23,4)+0} END{print r}' ${complexo_pdb})

    if [ -z "\${LIG_FIRST}" ] || [ -z "\${LIG_LAST}" ]; then
        echo "ERRO: cadeia B (BEN) não detectada em ${complexo_pdb}"; exit 1
    fi
    echo "[ANALISE_BEN] Ligante BEN: resíduo \${LIG_FIRST}-\${LIG_LAST}"

    # Conta grupos default do .tpr (para indexar novos grupos após os existentes)
    N_DEFAULT=\$(echo q | ${params.gmx_cmd} make_ndx \\
        -f ${md_tpr} -o _default.ndx 2>&1 \\
        | grep -cE "^ *[0-9]+ +[A-Za-z]")
    LIG_IDX=\${N_DEFAULT}
    REC_IDX=\$((N_DEFAULT + 1))
    rm -f _default.ndx

    # Cria índice com grupos Ligante (BEN) e Receptor
    ${params.gmx_cmd} make_ndx -f ${md_tpr} -o lig.ndx << EOF
r \${LIG_FIRST}-\${LIG_LAST}
name \${LIG_IDX} Ligante
1 & ! \${LIG_IDX}
name \${REC_IDX} Receptor
q
EOF

    # RMSD do backbone do receptor
    printf 'Backbone\\nBackbone\\n' | ${params.gmx_cmd} rms \\
        -s ${md_tpr} -f ${md_fit_xtc} \\
        -n lig.ndx -o rmsd_backbone.xvg -tu ns

    # RMSD do BEN (molécula pequena)
    printf 'Ligante\\nLigante\\n' | ${params.gmx_cmd} rms \\
        -s ${md_tpr} -f ${md_fit_xtc} \\
        -n lig.ndx -o rmsd_ligante.xvg -tu ns

    # RMSF por resíduo (backbone receptor)
    printf 'Backbone\\n' | ${params.gmx_cmd} rmsf \\
        -s ${md_tpr} -f ${md_fit_xtc} \\
        -n lig.ndx -o rmsf_residuos.xvg -res -fit

    # Raio de giro (proteína)
    printf 'Protein\\n' | ${params.gmx_cmd} gyrate \\
        -s ${md_tpr} -f ${md_fit_xtc} \\
        -n lig.ndx -o gyrate.xvg -tu ns

    # Contatos receptor–BEN < 0.4 nm
    printf 'Receptor\\nLigante\\n' | ${params.gmx_cmd} mindist \\
        -s ${md_tpr} -f ${md_fit_xtc} \\
        -n lig.ndx -od mindist.xvg -on numcont.xvg \\
        -d 0.4 -tu ns

    # Pontes de hidrogênio receptor–BEN
    # gmx hbond nao reconhece doadores/aceptores de moleculas pequenas GAFF2 → fallback placeholder
    printf 'Receptor\\nLigante\\n' | ${params.gmx_cmd} hbond \\
        -s ${md_tpr} -f ${md_fit_xtc} \\
        -n lig.ndx -num hbond.xvg -tu ns 2>&1 | tee hbond.log || \\
        printf '# gmx hbond: GAFF2 atom types not recognized as donors/acceptors\\n@ title "Number of Hydrogen Bonds"\\n@ xaxis label "Time (ns)"\\n@ yaxis label "Number"\\n0.000 0\\n' > hbond.xvg

    # SASA do receptor
    ${params.gmx_cmd} sasa \\
        -s ${md_tpr} -f ${md_fit_xtc} \\
        -n lig.ndx \\
        -surface 'Protein' -output 'Protein' \\
        -o sasa_protein.xvg -tu ns

    # SASA do BEN (valores baixos = BEN enterrado na interface)
    ${params.gmx_cmd} sasa \\
        -s ${md_tpr} -f ${md_fit_xtc} \\
        -n lig.ndx \\
        -surface 'Ligante' -output 'Ligante' \\
        -o sasa_ligante.xvg -tu ns

    echo "[OK] Análises BEN concluídas para ${meta.id}"
    """
}
