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

    # Wrapper Python que intercepta tleap e corrige o bug SS bonds do gmx_MMPBSA.
    # O bug: indices COM_OUT gerados com offset errado (N_ligante a menos).
    # Fix: copiar indices dos bonds do REC_OUT para o COM_OUT.
    mkdir -p bin_patch
    cat > bin_patch/tleap << 'WEOF'
#!/usr/bin/env python3
import sys, os, re, subprocess
args = sys.argv[1:]
for i, a in enumerate(args):
    if a == '-f' and i+1 < len(args):
        f = args[i+1]
        if os.path.exists(f):
            t = open(f).read()
            rec = re.findall(r'bond REC_OUT[.](\d+)[.]SG REC_OUT[.](\d+)[.]SG', t)
            com = re.findall(r'bond COM_OUT[.](\d+)[.]SG COM_OUT[.](\d+)[.]SG', t)
            for w, r in zip(com, rec):
                t = t.replace('bond COM_OUT.{}.SG COM_OUT.{}.SG'.format(w[0],w[1]),
                               'bond COM_OUT.{}.SG COM_OUT.{}.SG'.format(r[0],r[1]),1)
            open(f,'w').write(t)
        break
for d in os.environ.get('PATH','').split(':'):
    if 'bin_patch' in d: continue
    for name in ('tleap','teLeap'):
        c = os.path.join(d,name)
        if os.path.isfile(c) and os.access(c,os.X_OK):
            sys.exit(subprocess.run([c]+args).returncode)
sys.exit(1)
WEOF
    chmod +x bin_patch/tleap

    # Executa gmx_MMPBSA com bin_patch no PATH para que o wrapper intercepte tleap
    mamba run -n mmgbsa-env bash -c "
export PATH=\$PWD/bin_patch:\\\$PATH
gmx_MMPBSA -O \\
    -i mmgbsa.in \\
    -cs ${md_tpr} \\
    -ct ${md_xtc} \\
    -ci ${lig_ndx} \\
    -cg Receptor Ligante \\
    -o FINAL_RESULTS_MMGBSA.dat \\
    -eo mmgbsa_results.csv \\
    -nogui
" 2>&1 | tee mmgbsa.log
    test -f FINAL_RESULTS_MMGBSA.dat || { echo "gmx_MMPBSA falhou — veja mmgbsa.log"; exit 1; }

    plot_results.py \\
        --analise-dir . \\
        --titulo "${titulo}" \\
        --mmgbsa-csv mmgbsa_results.csv \\
        --output painel_final.png
    """
}
