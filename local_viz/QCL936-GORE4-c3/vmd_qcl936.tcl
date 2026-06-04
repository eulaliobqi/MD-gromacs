# VMD startup script — QCL936-GORE4 cluster3 (100 ns, pH 8.2)
# Uso: vmd -e vmd_qcl936.tcl
# Ou abrir VMD e: source vmd_qcl936.tcl

# ── Carregar estrutura + trajetória ────────────────────────────────────────
mol new estrutura_qcl936_c3.pdb type pdb waitfor all
mol addfile traj_qcl936_c3.xtc type xtc waitfor all

set mol [molinfo top get id]

# ── Limpar representação padrão ────────────────────────────────────────────
mol delrep 0 $mol

# ── Receptor (chain A) — NewCartoon azul ──────────────────────────────────
mol addrep $mol
mol modstyle  0 $mol NewCartoon
mol modselect 0 $mol "chain A"
mol modcolor  0 $mol ColorID 0
mol modmaterial 0 $mol AOChalky

# ── Ligante (chain B) — Licorice laranja ──────────────────────────────────
mol addrep $mol
mol modstyle  1 $mol Licorice 0.3 12 12
mol modselect 1 $mol "chain B"
mol modcolor  1 $mol ColorID 3
mol modmaterial 1 $mol AOShiny

# ── Superfície de contato do ligante (transparente) ───────────────────────
mol addrep $mol
mol modstyle  2 $mol QuickSurf 1.5 0.5 1.0 1
mol modselect 2 $mol "chain B"
mol modcolor  2 $mol ColorID 3
mol modmaterial 2 $mol Transparent

# ── Resíduos catalíticos (His92, Asp114, Ser211, Asp241) — VDW vermelho ───
mol addrep $mol
mol modstyle  3 $mol Licorice 0.4 12 12
mol modselect 3 $mol "resid 92 114 211 241 and chain A"
mol modcolor  3 $mol ColorID 1
mol modmaterial 3 $mol AOShiny

# ── Labels nos Cα dos resíduos catalíticos ────────────────────────────────
foreach {resid label} {92 "His92" 114 "Asp114" 211 "Ser211" 241 "Asp241(S1)"} {
    set sel [atomselect $mol "resid $resid and name CA and chain A"]
    if {[$sel num] > 0} {
        set idx [lindex [$sel get index] 0]
        label add Atoms $mol/$idx
    }
    $sel delete
}

# ── Vista inicial centrada no sítio ativo ─────────────────────────────────
set site [atomselect $mol "resid 92 114 211 241 and chain A"]
set center [measure center $site weight mass]
$site delete

display resetview
molinfo $mol set center_matrix [list [transoffset [vecscale -1 $center]]]

# ── Estilo de fundo e iluminação ──────────────────────────────────────────
color Display Background white
display shadows on
display ambientocclusion on
display aoambient 0.8
display aodirect 0.3

puts "QCL936-GORE4 c3 carregado. [molinfo $mol get numframes] frames."
puts "Tríade catalítica: His92 · Asp114 · Ser211 | Bolsão S1: Asp241"
puts "Animação: play/pause com tecla L ou botão play no rodapé."
