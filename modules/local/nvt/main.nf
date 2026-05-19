process NVT {
    tag "${meta.id}"
    label 'process_gpu'

    publishDir { "${params.outdir}/${meta.id}/nvt" }, mode: 'copy'

    input:
    tuple val(meta), path(em_gro), path(top, stageAs: 'input.top'), path(itps)

    output:
    tuple val(meta), path("nvt.gro"), path("nvt.cpt"), path("topol.top"), path("*.itp"), emit: system

    script:
    def gpu_flags = params.use_gpu ? "-nb gpu -pme gpu -bonded gpu -gpu_id ${params.gpu_id}" : ""
    def mpi       = params.mpi_cmd ?: ""
    def temp      = params.temperature
    """
    cp ${top} topol.top

    cat > nvt.mdp << MDP_EOF
define          = -DPOSRES
integrator      = md
dt              = 0.002
nsteps          = 50000
nstxout         = 500
nstvout         = 500
nstenergy       = 100
cutoff-scheme   = Verlet
nstlist         = 20
coulombtype     = PME
rcoulomb        = 1.2
vdwtype         = Cut-off
rvdw            = 1.2
constraints     = h-bonds
constraint-algorithm = LINCS
continuation    = no
gen-vel         = yes
gen-temp        = ${temp}
gen-seed        = 42
tcoupl          = V-rescale
tc-grps         = Protein Non-Protein
tau-t           = 0.1 0.1
ref-t           = ${temp} ${temp}
pcoupl          = no
pbc             = xyz
MDP_EOF

    ${params.gmx_cmd} grompp \\
        -f nvt.mdp \\
        -c ${em_gro} -r ${em_gro} \\
        -p topol.top -o nvt.tpr \\
        -maxwarn ${params.maxwarn}

    ${mpi} ${params.gmx_cmd} mdrun \\
        -v -deffnm nvt \\
        -ntomp ${params.ntomp} \\
        -pin on ${gpu_flags}
    """
}
