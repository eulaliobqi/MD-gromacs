# PyMOL publication figure — xp273-skti-c2  (t = 85.8 ns)
load /home/eulalio/gromacs/MD-gromacs/poses/xp273-skti-c2/receptor.pdb, receptor
load /home/eulalio/gromacs/MD-gromacs/poses/xp273-skti-c2/ligand.pdb, ligand
bg_color white

hide everything, receptor
show cartoon, receptor
color gray75, receptor
set cartoon_fancy_helices, 1

hide everything, ligand
show cartoon, ligand
color cyan, ligand
set cartoon_transparency, 0.0, ligand

# (nenhum resíduo-chave disponível)
distance hbonds, receptor, ligand, 3.5, mode=2
set dash_color, red
set dash_width, 2.5
set dash_gap, 0.35
set dash_radius, 0.06
hide labels, hbonds

orient ligand
zoom ligand, 15
rotate y, 15
rotate x, -10

set ray_shadows, 1
set ray_opaque_background, 1
set antialias, 2
set ambient, 0.35
set direct, 0.65
set specular, 0.20

ray 2400, 1800
png /home/eulalio/gromacs/MD-gromacs/poses/xp273-skti-c2/pose.png, dpi=300, ray=1
quit
