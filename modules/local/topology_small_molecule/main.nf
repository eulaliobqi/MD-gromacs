process TOPOLOGY_SMALL_MOLECULE {
    tag "${meta.id}"
    label 'process_medium'

    publishDir { "${params.outdir}/${meta.id}/topo" }, mode: 'copy'

    input:
    tuple val(meta), path(complexo_pdb)

    output:
    tuple val(meta), path("complexo.gro"), path("topol.top"), path("*.itp"), emit: topology

    script:
    """
    echo "=== TOPOLOGY_SMALL_MOLECULE: ${meta.id} ==="

    # 1. Separar cadeia A (receptor, ATOM) e cadeia B (BEN, ATOM ou HETATM)
    awk '/^ATOM/ && substr(\$0,22,1)=="A" {print}' ${complexo_pdb} > receptor.pdb
    echo "TER" >> receptor.pdb
    echo "END" >> receptor.pdb

    awk '(/^ATOM/ || /^HETATM/) && substr(\$0,22,1)=="B" {print}' ${complexo_pdb} > ben.pdb
    echo "END" >> ben.pdb

    NATOM_REC=\$(grep -c "^ATOM" receptor.pdb || echo 0)
    NATOM_BEN=\$(grep -cE "^ATOM|^HETATM" ben.pdb || echo 0)
    echo "  Receptor: \${NATOM_REC} átomos ATOM"
    echo "  BEN:      \${NATOM_BEN} átomos"

    # 2. Topologia do receptor com AMBER99SB-ILDN (opcao 0 N/C terminal)
    printf '0\\n0\\n' | ${params.gmx_cmd} pdb2gmx \\
        -f receptor.pdb \\
        -o receptor.gro \\
        -p receptor.top \\
        -i posre_Protein.itp \\
        -ff ${params.forcefield} \\
        -water ${params.water} \\
        -ignh -ter \\
        2>&1 | tee pdb2gmx.log

    if [ ! -f receptor.gro ]; then
        echo "ERRO: pdb2gmx falhou — ver pdb2gmx.log"; exit 1
    fi

    # 3. Topologia do BEN com ACPYPE (GAFF2, carga +1 = amidínio pH 8.2)
    #    -c bcc: cargas AM1-BCC (antechamber)
    #    -n 1: carga total +1 (protonado)
    #    -a gaff2: campo de força GAFF2
    #    -b BEN: nome base da molécula
    acpype -i ben.pdb -c bcc -n 1 -a gaff2 -b BEN -o gmx \\
        2>&1 | tee acpype.log

    # Verificar se ACPYPE gerou os arquivos necessários
    if [ ! -f "BEN.acpype/BEN_GMX.gro" ] || [ ! -f "BEN.acpype/BEN_GMX.itp" ]; then
        echo "AVISO: AM1-BCC falhou — tentando cargas Gasteiger (-c gas)"
        acpype -i ben.pdb -c gas -n 1 -a gaff2 -b BEN -o gmx \\
            2>&1 | tee acpype_gas.log
    fi

    if [ ! -f "BEN.acpype/BEN_GMX.gro" ]; then
        echo "ERRO: ACPYPE falhou com BCC e GAS — ver acpype.log"; exit 1
    fi

    # 4. Mesclar GROs e patchar topologia
    python3 ${projectDir}/bin/merge_small_molecule_topology.py \\
        --protein-gro receptor.gro \\
        --ligand-gro  BEN.acpype/BEN_GMX.gro \\
        --protein-top receptor.top \\
        --ligand-itp  BEN.acpype/BEN_GMX.itp \\
        --ligand-mol  BEN \\
        --out-gro     complexo.gro \\
        --out-top     topol.top

    if [ ! -f complexo.gro ]; then
        echo "ERRO: merge_small_molecule_topology.py falhou"; exit 1
    fi

    # 5. Copiar ITPs para publishDir
    cp posre_Protein.itp .
    cp BEN.acpype/BEN_GMX.itp .
    cp BEN.acpype/posre_BEN.itp .

    # Verificar resultado
    NTOTAL=\$(awk 'NR==2 {print \$1}' complexo.gro)
    echo "[OK] complexo.gro: \${NTOTAL} átomos totais"
    echo "[OK] Topologia mesclada: topol.top"
    """
}
