# ============================================================
# Visualizacao interativa no PyMOL (local)
# Uso:  pymol visualizar.pml
# ============================================================

# Carrega estrutura + trajetoria
load estrutura.pdb, complex
load_traj trajetoria.xtc, complex

# Alinha pelo backbone para reduzir movimento global
intra_fit complex and chain A and name CA

# === Estilo visual ===
hide everything
bg_color white

# Receptor (chain A): cartoon azul
show cartoon, chain A
color marine, chain A
set cartoon_transparency, 0.0

# Ligante (chain B): sticks laranja + superficie semitransparente
show sticks, chain B
show surface, chain B
color orange, chain B
set transparency, 0.4, chain B
set surface_color, lightorange, chain B

# Highlight residuos do sitio ativo (His-Asp-Ser typical da serina-protease)
# AJUSTE conforme sua tripsina! Para DN773/DN1973 sao geralmente:
select sitio_ativo, chain A and resi 57+102+195
show sticks, sitio_ativo
color red, sitio_ativo

# Camera centralizada no ligante
orient chain B
zoom complex, 5

# Labels uteis
set label_size, 14
set label_color, black

# Velocidade de animacao
set movie_loop, 1
set movie_fps, 15

# Mostra controles
print "================================================"
print " Trajetoria carregada: %d frames" % cmd.count_states()
print " Use ▶ no canto inferior direito para animar"
print " ScrollMouse: zoom | Botao esquerdo: rotacao"
print "================================================"
