# Visualizacao local da DM (Tripsina + LALAK)

## Pre-requisitos no seu PC

### Windows
1. PyMOL: baixe de https://pymol.org/2/ ou via conda:
   ```
   conda install -c conda-forge pymol-open-source
   ```
2. ffmpeg: https://www.gyan.dev/ffmpeg/builds/ (adicione ao PATH)

### Mac
```bash
brew install pymol ffmpeg
```

### Linux (Ubuntu/Debian)
```bash
sudo apt install pymol ffmpeg imagemagick
# ou via conda:
mamba install -c conda-forge pymol-open-source ffmpeg
```

## Arquivos baixados do servidor

```
visualizacao_local/
  estrutura.pdb          # Frame 0 (tempo zero)
  estrutura_final.pdb    # Frame final (100 ns)
  trajetoria.xtc         # ~200 frames reduzidos
  visualizar.pml         # Script para visualizacao interativa
  gerar_filme.pml        # Script para gerar PNGs
  README_LOCAL.md        # Este arquivo
```

## 1) Visualizacao interativa (mais simples)

```bash
cd visualizacao_local/
pymol visualizar.pml
```

Vai abrir a GUI do PyMOL com:
- Receptor em cartoon azul
- Ligante (LALAK) em sticks + superficie laranja
- Triade catalitica destacada em vermelho

**Controles:**
- `▶` no canto inferior direito = anima
- Scroll do mouse = zoom
- Botao esquerdo + arrastar = rotacao
- `Ctrl+Click` em atomo = mostra info

## 2) Snapshot estatico (alta resolucao para artigo)

Na linha de comando do PyMOL (apos abrir):
```
ray 1920, 1080
png snapshot.png, dpi=300
```

## 3) GIF / MP4 da trajetoria

### Passo 1: gerar frames (PyMOL)
```bash
pymol -cq gerar_filme.pml
# Demora 5-30 min dependendo da CPU (200 frames ray-traced)
```

Saida: `frames/frame_0001.png` a `frame_0200.png`

### Passo 2: montar GIF (ffmpeg)
```bash
# GIF leve (~5-15 MB)
ffmpeg -framerate 25 -i frames/frame_%04d.png \
    -vf "fps=15,scale=600:-1:flags=lanczos,palettegen" palette.png

ffmpeg -framerate 25 -i frames/frame_%04d.png -i palette.png \
    -lavfi "fps=15,scale=600:-1[x];[x][1:v]paletteuse" \
    trajetoria.gif
```

### Passo 3: montar MP4 (alta qualidade para apresentacao)
```bash
ffmpeg -framerate 25 -i frames/frame_%04d.png \
    -c:v libx264 -pix_fmt yuv420p -crf 18 \
    trajetoria.mp4
```

## 4) Visualizacao 3D online (sem instalar nada)

Arraste `estrutura_final.pdb` para:
- https://molstar.org/viewer/ (Mol*, melhor opcao)
- https://nglviewer.org/ngl/

Ambos sao 100% web, funcionam em qualquer navegador.

## 5) Sessao PyMOL salvavel (.pse)

Apos abrir com visualizar.pml, na linha de comando:
```
save sessao.pse
```

Depois e so abrir `sessao.pse` direto no PyMOL para retomar.

## Dicas

**Reduzir mais a trajetoria** (se estiver pesada):
```bash
gmx trjconv -s sistema.tpr -f trajetoria.xtc \
    -o traj_leve.xtc -skip 5     # 1 a cada 5 frames
```

**Rotacao 360 da estrutura** (no PyMOL):
```
mset 1 x180
mview store, 1
turn y, 360
mview store, 180
mplay
```

**Mudar cores** na linha de comando do PyMOL:
```
color cyan, chain A          # receptor cian
color magenta, chain B       # ligante magenta
util.cnc                     # carbon green/N blue/O red
```
