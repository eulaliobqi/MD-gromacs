from docx import Document
from docx.shared import Pt, Inches, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml.ns import qn
from docx.oxml import OxmlElement

doc = Document()

# ── Estilos gerais ──────────────────────────────────────────────────────────
style_normal = doc.styles['Normal']
style_normal.font.name = 'Times New Roman'
style_normal.font.size = Pt(12)

for s in ['Heading 1', 'Heading 2', 'Heading 3']:
    h = doc.styles[s]
    h.font.name = 'Times New Roman'
    h.font.color.rgb = RGBColor(0, 0, 0)

doc.styles['Heading 1'].font.size = Pt(14)
doc.styles['Heading 2'].font.size = Pt(13)
doc.styles['Heading 3'].font.size = Pt(12)

# margens 2,5 cm
for section in doc.sections:
    section.top_margin    = Inches(0.984)
    section.bottom_margin = Inches(0.984)
    section.left_margin   = Inches(1.181)
    section.right_margin  = Inches(0.984)

def add_heading(text, level=1):
    h = doc.add_heading(text, level=level)
    h.alignment = WD_ALIGN_PARAGRAPH.LEFT
    return h

def add_para(text, justify=True, bold_parts=None):
    p = doc.add_paragraph()
    if justify:
        p.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
    p.paragraph_format.first_line_indent = Pt(36)
    p.paragraph_format.space_after = Pt(0)
    run = p.add_run(text)
    run.font.name = 'Times New Roman'
    run.font.size = Pt(12)
    return p

def add_para_noindent(text):
    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
    p.paragraph_format.space_after = Pt(0)
    run = p.add_run(text)
    run.font.name = 'Times New Roman'
    run.font.size = Pt(12)
    return p

# ── Título ──────────────────────────────────────────────────────────────────
titulo = doc.add_paragraph()
titulo.alignment = WD_ALIGN_PARAGRAPH.CENTER
r = titulo.add_run('METODOLOGIA')
r.bold = True
r.font.name = 'Times New Roman'
r.font.size = Pt(14)
doc.add_paragraph()

# ════════════════════════════════════════════════════════════════════════════
# 1. DOCKING MOLECULAR
# ════════════════════════════════════════════════════════════════════════════
add_heading('1. Docking Molecular', level=1)

add_para(
    'As estruturas tridimensionais dos peptídeos candidatos a inibidores foram geradas '
    'por modelagem de novo utilizando o servidor PEPstrMOD (Singh et al., 2015), '
    'considerando o pH fisiológico do intestino médio de Spodoptera spp. (pH 8,2). '
    'A estrutura cristalográfica da tripsina digestiva de Spodoptera frugiperda (GORE4) '
    'foi obtida a partir do Protein Data Bank (PDB) ou modelada por homologia utilizando '
    'o servidor AlphaFold2, seguida de refinamento estrutural com o servidor ModRefiner.'
)

add_para(
    'O preparo dos receptores para o docking foi realizado com o software PDB2PQR v3.6 '
    '(Dolinsky et al., 2007), utilizando o campo de força AMBER e atribuição de estados '
    'de protonação pelo método PROPKA em pH 8,2, compatível com o ambiente intestinal '
    'alcalino de lepidópteros. Moléculas de água e ligantes heterólogos foram removidos '
    'da estrutura do receptor antes da análise.'
)

add_para(
    'O docking proteína–peptídeo foi conduzido no servidor HADDOCK 2.4 '
    '(High Ambiguity Driven DOCKing) (van Zundert et al., 2016; Dominguez et al., 2003). '
    'As restrições ambíguas de interação (AIRs) foram definidas com base nos resíduos '
    'do sítio ativo da tripsina — a tríade catalítica formada por Serina (Ser195), '
    'Histidina (His57) e Aspartato (Asp102) — como resíduos "ativos" do receptor, '
    'e nos resíduos básicos dos peptídeos (lisina e arginina) como resíduos ativos '
    'dos ligantes, em conformidade com a especificidade de clivagem de tripsinas por '
    'resíduos básicos.'
)

add_para(
    'O protocolo padrão do HADDOCK foi aplicado em três etapas: (i) minimização de '
    'energia de corpo rígido (it0), gerando 1.000 estruturas; (ii) refinamento '
    'semirígido em espaço de torção (it1), retendo as 200 melhores estruturas; e '
    '(iii) refinamento final em solvente explícito (water refinement), com seleção '
    'das 200 melhores poses. As estruturas foram agrupadas utilizando RMSD de 7,5 Å '
    'sobre os átomos da interface proteína–peptídeo. Os clusters foram classificados '
    'pelo HADDOCK score, calculado como combinação linear das energias de van der Waals '
    '(peso = 1,0), eletrostática (peso = 0,2), dessolvatação (peso = 1,0) e restrições '
    'ambíguas (peso = 0,1). A pose representativa de cada cluster (estrutura de menor '
    'energia dentro do cluster) foi selecionada para as simulações de dinâmica molecular.'
)

# ════════════════════════════════════════════════════════════════════════════
# 2. DINÂMICA MOLECULAR
# ════════════════════════════════════════════════════════════════════════════
add_heading('2. Dinâmica Molecular', level=1)

# 2.1
add_heading('2.1 Preparo do sistema', level=2)

add_para(
    'Os complexos tripsina–peptídeo obtidos pelo docking foram preparados para '
    'as simulações de dinâmica molecular utilizando um pipeline automatizado '
    'desenvolvido em Nextflow DSL2 (Di Tommaso et al., 2017), integrando ferramentas '
    'do pacote GROMACS 2026 (Abraham et al., 2015) e PDB2PQR v3.6. '
    'Todos os sistemas foram tratados em pH 8,2 para reproduzir o ambiente fisiológico '
    'do intestino médio de Spodoptera spp.'
)

add_para(
    'O preparo dos complexos compreendeu as seguintes etapas: (i) atribuição de estados '
    'de protonação de todos os resíduos ionizáveis pelo método PROPKA em pH 8,2, '
    'com posterior conversão para o formato GROMACS; (ii) fusão das coordenadas do '
    'receptor e do peptídeo em uma única estrutura do complexo; e (iii) geração da '
    'topologia com o campo de força AMBER99SB-ILDN (Lindorff-Larsen et al., 2010), '
    'adequado para simulações de proteínas e peptídeos.'
)

add_para(
    'O complexo foi centralizado em uma caixa dodecaédrica com margem mínima de '
    '1,2 nm entre o soluto e as faces da caixa, garantindo ausência de interações '
    'entre imagens periódicas. A caixa foi preenchida com moléculas de água no '
    'modelo TIP3P (Jorgensen et al., 1983). A neutralização do sistema e o '
    'ajuste da concentração de íons a 0,15 M de NaCl foram realizados pela '
    'substituição aleatória de moléculas de água por íons Na⁺ e Cl⁻, '
    'mimetizando a força iônica fisiológica intestinal.'
)

# 2.2
add_heading('2.2 Minimização de energia e equilibração', level=2)

add_para(
    'A minimização de energia foi realizada pelo método do gradiente conjugado '
    'com convergência definida quando a força máxima sobre qualquer átomo foi '
    'inferior a 1.000 kJ mol⁻¹ nm⁻¹. Após a minimização, o sistema foi '
    'submetido a duas etapas de equilibração termodinâmica com restrições de '
    'posição aplicadas aos átomos pesados do soluto (constante de força: '
    '1.000 kJ mol⁻¹ nm⁻²).'
)

add_para(
    'Na primeira etapa de equilibração (ensemble NVT), o sistema foi aquecido '
    'gradualmente até 300 K durante 100 ps, utilizando o termostato de velocidades '
    'reescalonadas (v-rescale) com constante de acoplamento de 0,1 ps. '
    'Na segunda etapa (ensemble NPT), a pressão foi ajustada a 1 bar durante '
    '100 ps, com o barostato de Parrinello-Rahman (constante de acoplamento de '
    '2,0 ps), mantendo a temperatura de 300 K.'
)

# 2.3
add_heading('2.3 Simulação de produção', level=2)

add_para(
    'Após a equilibração, foram realizadas simulações de dinâmica molecular de '
    'produção com duração de 100 ns para cada complexo tripsina–peptídeo, '
    'no ensemble NPT (300 K, 1 bar). O passo de integração foi de 2 fs, '
    'com as ligações envolvendo hidrogênio restringidas pelo algoritmo LINCS '
    '(Hess et al., 1997). As interações eletrostáticas de longo alcance foram '
    'tratadas pelo método Particle Mesh Ewald (PME) (Darden et al., 1993), '
    'com frequência de atualização da lista de vizinhos a cada 10 passos. '
    'As coordenadas foram salvas a cada 10 ps, gerando trajetórias de 10.000 frames '
    'para análise. As simulações foram executadas em servidor Linux Debian com GPU '
    'NVIDIA, utilizando o binário gmx_mpi compilado com suporte CUDA para '
    'aceleração de cálculos de interações de curto alcance e PME na GPU.'
)

# 2.4
add_heading('2.4 Análise das trajetórias', level=2)

add_para(
    'As trajetórias de produção foram pós-processadas com GROMACS para correção '
    'de condições periódicas de contorno (PBC), centralização do complexo e '
    'alinhamento pelo backbone do receptor. As seguintes análises foram realizadas '
    'para cada sistema:'
)

# Lista de análises
analises = [
    ('Estabilidade estrutural (RMSD):',
     'o desvio quadrático médio (RMSD) foi calculado para o backbone do receptor '
     '(Cα, C, N, O) e separadamente para os átomos pesados do peptídeo, tomando '
     'a estrutura do início da produção como referência, utilizando gmx rms.'),
    ('Flexibilidade por resíduo (RMSF):',
     'a flutuação quadrática média das posições atômicas (RMSF) foi calculada '
     'por resíduo sobre o backbone do receptor ao longo de toda a trajetória '
     'de produção, com gmx rmsf.'),
    ('Raio de giro (Rg):',
     'a compactação estrutural do complexo foi monitorada pelo raio de giro '
     'calculado sobre todos os átomos da proteína com gmx gyrate.'),
    ('Contatos e distância mínima:',
     'o número de contatos intermoleculares receptor–peptídeo (distância < 0,4 nm) '
     'e a distância mínima entre os grupos foram calculados com gmx mindist.'),
    ('Pontes de hidrogênio:',
     'as pontes de hidrogênio intermoleculares entre receptor e peptídeo foram '
     'identificadas e contadas ao longo da trajetória com gmx hbond, '
     'utilizando os critérios geométricos padrão '
     '(distância doador–aceptor ≤ 0,35 nm; ângulo ≤ 30°).'),
    ('Área de superfície acessível ao solvente (SASA):',
     'a SASA total do receptor e do peptídeo no complexo foi calculada com '
     'gmx sasa, utilizando o algoritmo de Connolly com raio de sonda de 0,14 nm. '
     'Valores de SASA reduzidos do peptídeo ao longo do tempo indicam '
     'enterramento progressivo na interface de ligação.'),
    ('Distâncias à tríade catalítica:',
     'a distância mínima entre os átomos pesados do peptídeo e cada resíduo '
     'da tríade catalítica (Ser195, His57, Asp102) foi monitorada ao longo '
     'da trajetória com gmx mindist, como indicador de ocupação e bloqueio '
     'do sítio ativo da tripsina.'),
]

for titulo_item, descricao in analises:
    p = doc.add_paragraph(style='List Bullet')
    p.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
    p.paragraph_format.space_after = Pt(2)
    run_bold = p.add_run(titulo_item + ' ')
    run_bold.bold = True
    run_bold.font.name = 'Times New Roman'
    run_bold.font.size = Pt(12)
    run_normal = p.add_run(descricao)
    run_normal.font.name = 'Times New Roman'
    run_normal.font.size = Pt(12)

doc.add_paragraph()

add_para(
    'A fase estável de cada trajetória foi determinada pela identificação do '
    'plateau no perfil de RMSD do backbone, descartando a fase de equilibração '
    'inicial. Os critérios de inibição peptídica foram avaliados de forma '
    'integrada, considerando: (i) RMSD do peptídeo < 0,3 nm (estabilidade '
    'conformacional na interface); (ii) presença de pelo menos uma ponte de '
    'hidrogênio persistente com o receptor; (iii) distância à tríade catalítica '
    '< 0,5 nm; (iv) redução da SASA do peptídeo (enterramento); e (v) número '
    'de contatos intermoleculares estável ao longo da fase estável da simulação.'
)

# ════════════════════════════════════════════════════════════════════════════
# REFERÊNCIAS
# ════════════════════════════════════════════════════════════════════════════
add_heading('Referências', level=1)

refs = [
    'Abraham, M. J. et al. GROMACS: High performance molecular simulations through '
    'multi-level parallelism from laptops to supercomputers. SoftwareX, v. 1–2, '
    'p. 19–25, 2015.',

    'Darden, T.; York, D.; Pedersen, L. Particle mesh Ewald: An N·log(N) method '
    'for Ewald sums in large systems. Journal of Chemical Physics, v. 98, '
    'p. 10089–10092, 1993.',

    'Di Tommaso, P. et al. Nextflow enables reproducible computational workflows. '
    'Nature Biotechnology, v. 35, p. 316–319, 2017.',

    'Dolinsky, T. J. et al. PDB2PQR: expanding and upgrading automated preparation '
    'of biomolecular structures for molecular simulations. Nucleic Acids Research, '
    'v. 35, p. W522–W525, 2007.',

    'Dominguez, C.; Boelens, R.; Bonvin, A. M. J. J. HADDOCK: A protein–protein '
    'docking approach based on biochemical or biophysical information. Journal of '
    'the American Chemical Society, v. 125, p. 1731–1737, 2003.',

    'Hess, B. et al. LINCS: A linear constraint solver for molecular simulations. '
    'Journal of Computational Chemistry, v. 18, p. 1463–1472, 1997.',

    'Jorgensen, W. L. et al. Comparison of simple potential functions for simulating '
    'liquid water. Journal of Chemical Physics, v. 79, p. 926–935, 1983.',

    'Lindorff-Larsen, K. et al. Improved side-chain torsion potentials for the Amber '
    'ff99SB protein force field. Proteins, v. 78, p. 1950–1958, 2010.',

    'Singh, S. et al. PEPstrMOD: structure prediction of peptides containing natural, '
    'non-natural and modified residues. Biology Direct, v. 10, p. 73, 2015.',

    'van Zundert, G. C. P. et al. The HADDOCK2.2 Web Server: User-Friendly '
    'Integrative Modeling of Biomolecular Complexes. Journal of Molecular Biology, '
    'v. 428, p. 720–725, 2016.',
]

for ref in refs:
    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
    p.paragraph_format.space_after = Pt(6)
    p.paragraph_format.left_indent = Pt(36)
    p.paragraph_format.first_line_indent = Pt(-36)
    run = p.add_run(ref)
    run.font.name = 'Times New Roman'
    run.font.size = Pt(12)

doc.save('metodologia_artigo.docx')
print('Arquivo gerado: metodologia_artigo.docx')
