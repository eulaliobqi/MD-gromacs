# ============================================================
# Gera frames PNG para virar GIF/MP4 - PyMOL local com GUI
# Uso:  pymol gerar_filme.pml
# Apos: ffmpeg para montar o filme (veja README_LOCAL.md)
# ============================================================

load estrutura.pdb, complex
load_traj trajetoria.xtc, complex
intra_fit complex and chain A and name CA

# Estilo
hide everything
bg_color white
set ray_shadows, 0
set ray_opaque_background, 1
set antialias, 2

show cartoon, chain A
color marine, chain A

show sticks, chain B
show surface, chain B
color orange, chain B
set transparency, 0.3, chain B
set surface_color, lightorange, chain B

orient chain B
zoom complex, 5

# === Renderiza cada frame com ray-tracing ===
# Como roda local com GUI, ray-trace funciona perfeito
import os
os.makedirs('frames', exist_ok=True)

n = cmd.count_states('complex')
print(f"Renderizando {n} frames com ray-tracing (~30s a 2min cada)...")

for i in range(1, n + 1):
    cmd.frame(i)
    cmd.set('state', i)
    cmd.ray(1280, 720)
    cmd.png(f'frames/frame_{i:04d}.png', dpi=150)
    if i % 10 == 0:
        print(f"  ... frame {i}/{n}")

print(f"Pronto! {n} PNGs em frames/")
print("Para montar o GIF/MP4, veja README_LOCAL.md")
