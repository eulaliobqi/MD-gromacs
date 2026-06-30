# PyMOL publication figure — qcl936-gore4-c3  (t = 80.2 ns)
load /home/eulalio/gromacs/MD-gromacs/poses/qcl936-gore4-c3/receptor.pdb, receptor
load /home/eulalio/gromacs/MD-gromacs/poses/qcl936-gore4-c3/ligand.pdb, ligand
bg_color white

hide everything, receptor
show cartoon, receptor
color gray75, receptor
set cartoon_fancy_helices, 1
set cartoon_rect_length, 1.2

hide everything, ligand
show sticks, ligand
color cyan, ligand
set stick_radius, 0.18
show spheres, ligand and name CA
set sphere_scale, 0.25, ligand and name CA

# (nenhum resíduo-chave disponível)
distance hbonds, receptor, ligand, 3.5, mode=2
set dash_color, red
set dash_width, 2.5
set dash_gap, 0.35
set dash_radius, 0.06
hide labels, hbonds

create pocket, byres receptor within 5 of ligand
show surface, pocket
color tv_blue, pocket
set transparency, 0.60, pocket

orient ligand
zoom ligand, 10
rotate y, 15
rotate x, -10

set ray_shadows, 1
set ray_opaque_background, 1
set antialias, 2
set ambient, 0.35
set direct, 0.65
set specular, 0.20

ray 2400, 1800
png /home/eulalio/gromacs/MD-gromacs/poses/qcl936-gore4-c3/pose.png, dpi=300, ray=1
quit
