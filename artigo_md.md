# Dinâmica Molecular — Seções do Artigo

---

## 2. Metodologia

### 2.1 Preparação dos complexos proteína–peptídeo

As estruturas iniciais dos complexos entre a tripsina digestiva de *Spodoptera* (isoforma GORE4) e os peptídeos inibidores candidatos (XP273, DN773, DN1441, ACR157, QCL936) foram geradas por docking proteína–peptídeo utilizando o servidor HADDOCK 2.4 (van Zundert *et al.*, 2016). Para cada peptídeo, o cluster de maior escore (menor energia de interação) foi selecionado como pose inicial para as simulações.

Os estados de protonação dos resíduos ionizáveis foram determinados para pH 8,2, simulando as condições alcalinas do intestino médio larval de *Spodoptera* (*midgut* pH 9,5–11,0; Al-Khafaji *et al.*, 2024). A protonação foi calculada com o servidor PropKa 3.0, implementado via `pdb2pqr 3.x` com campo de força AMBER. Nesse pH, a His catalítica (His-ε) assume forma neutra (HIE), correspondente à conformação ativa da serino-protease.

### 2.2 Configuração do campo de força e parâmetros de simulação

Todas as simulações foram conduzidas com o pacote GROMACS 2026 (Abraham *et al.*, 2015), utilizando o campo de força AMBER99SB-ILDN (Lindorff-Larsen *et al.*, 2010) e o modelo de água explícita TIP3P. Os complexos foram inseridos em caixas cúbicas periódicas com margem mínima de 2,0 nm entre o soluto e as faces da caixa. O sistema foi solvatado e neutralizado com a adição de íons K⁺ e Cl⁻ a uma concentração de 0,10 M (parâmetros de Joung & Cheatham, 2008), refletindo a composição iônica da hemolinfa de insetos lepidópteros, predominantemente de K⁺ (Dow, 1992).

### 2.3 Protocolo de equilibração

O protocolo de equilibração seguiu três etapas consecutivas:

1. **Minimização de energia** — método de descida mais íngreme (*steepest descent*), critério de convergência de 1000 kJ mol⁻¹ nm⁻¹, máximo de 50.000 passos.

2. **Equilibração NVT (200 ps)** — ensemble canônico a 300 K utilizando o termostato V-rescale (Bussi *et al.*, 2007) com constante de acoplamento τ = 0,1 ps, aplicado separadamente ao complexo e ao solvente. Restrições de posição (`POSRES`) foram aplicadas aos átomos pesados do soluto. A temperatura de 300 K corresponde ao intervalo fisiológico larval de *Spodoptera* (25–30 °C; Bray *et al.*, 2024).

3. **Equilibração NPT (500 ps)** — ensemble isobárico-isotérmico a 300 K e 1 bar, com barostato de Berendsen (τ = 2,0 ps) para estabilização da densidade (Al-Khafaji *et al.*, 2024). As restrições de posição foram mantidas.

### 2.4 Simulação de produção

As simulações de produção foram realizadas por 100 ns cada, sem restrições de posição, no ensemble NPT utilizando o barostato de Parrinello-Rahman (τ = 2,0 ps; Parrinello & Rahman, 1981). O integrador *leap-frog* foi empregado com passo de integração de 2 fs. As ligações envolvendo hidrogênios foram restringidas com o algoritmo LINCS (Hess *et al.*, 1997). As interações eletrostáticas de longo alcance foram calculadas pelo método *Particle Mesh Ewald* (PME) com raio de corte de 1,2 nm. O cálculo de van der Waals utilizou o mesmo raio de corte com função de troca suave. As coordenadas foram gravadas a cada 10 ps (trajetórias `.xtc`). O processamento foi acelerado por GPU (CUDA), com integração de equações de movimento e cálculo de PME descarregados na placa gráfica.

### 2.5 Pós-processamento das trajetórias

As trajetórias brutas foram tratadas em três etapas para correção de artefatos das condições periódicas de contorno (PBC): (i) remoção de saltos atômicos entre caixas periódicas (`-pbc nojump`); (ii) centralização do complexo na caixa com compactação do solvente (`-pbc mol -center -ur compact`); (iii) alinhamento translacional e rotacional de todos os frames ao frame inicial pelo backbone (`-fit rot+trans`).

### 2.6 Análises conformacionais e de interação

As seguintes propriedades foram calculadas ao longo de toda a trajetória:

- **RMSD** (*Root Mean Square Deviation*) do backbone do receptor e do ligante, medindo estabilidade conformacional em relação ao frame inicial;
- **RMSF** (*Root Mean Square Fluctuation*) por resíduo, identificando regiões de maior flexibilidade;
- **Raio de giro** (Rg), indicador de compacidade global do receptor;
- **Número de contatos** receptor–ligante (pares atômicos com distância < 0,4 nm);
- **Pontes de hidrogênio** receptor–ligante (critério geométrico: distância O/N–H < 0,35 nm, ângulo ≤ 30°);
- **Área acessível ao solvente** (SASA) do receptor e do ligante;
- **Distâncias mínimas** entre o ligante e cada resíduo do sítio catalítico (His, Asp, Ser da tríade e Asp do bolsão S1), calculadas pelo grupo de átomos de cada resíduo via `gmx mindist`.

As médias móveis foram calculadas com janela de 5 ns para suavização dos perfis temporais. Todos os gráficos foram gerados com Matplotlib 3.x.

O pipeline completo, desde a preparação dos complexos até a geração dos painéis de análise, foi implementado em Nextflow DSL2 (Di Tommaso *et al.*, 2017) para garantir reprodutibilidade e rastreabilidade computacional.

---

## 3. Resultados e Discussão

### 3.1 Estabilidade estrutural dos complexos

#### QCL936-GORE4 (cluster 3)

A simulação do complexo QCL936-GORE4 (cluster 3) foi conduzida por 100 ns a 300 K e pH 8,2 e apresentou comportamento estável ao longo de toda a trajetória. O RMSD do backbone do receptor atingiu um platô de 0,124 ± 0,017 nm após aproximadamente 20 ns (Figura X), indicando equilíbrio conformacional do receptor sem deriva estrutural. O raio de giro manteve-se em 1,718 ± 0,006 nm com variação mínima (< 1%), confirmando que a estrutura terciária da tripsina não sofreu perturbações decorrentes da presença do peptídeo.

O RMSD do ligante estabilizou em 0,229 ± 0,042 nm, valor inferior a 0,3 nm durante praticamente toda a simulação, indicando que o peptídeo QCL936 manteve sua conformação de docking ao longo dos 100 ns sem sinais de dissociação do sítio de ligação. Não foi observado o padrão de dissociação característico (colapso de contatos seguido de aumento abrupto das distâncias), o que corrobora a estabilidade do complexo.

#### Interface de ligação

O número médio de contatos receptor–ligante foi de 340 ± 41 átomos a menos de 0,4 nm, valor expressivamente superior ao limiar mínimo de 150 átomos adotado como critério de ligação estável. O perfil temporal dos contatos mostrou oscilações regulares sem tendência de redução, indicando interface robusta ao longo de toda a trajetória.

As pontes de hidrogênio receptor–ligante apresentaram média de 2,75 ± 1,01, com persistência contínua ao longo dos 100 ns. A manutenção de pelo menos uma ponte de hidrogênio em praticamente todos os frames demonstra interação específica entre o peptídeo e o receptor, contribuindo para a estabilização do complexo.

A SASA do ligante manteve-se em 9,04 ± 0,37 nm² sem tendência crescente, indicando que o peptídeo permaneceu enterrado na interface receptor–ligante e não migrou para o solvente.

#### Posicionamento no sítio catalítico

A análise das distâncias mínimas entre o peptídeo QCL936 e os resíduos do sítio ativo revelou um padrão de ligação predominantemente associado ao bolsão S1 e à His periférica da tríade (Figura X). O peptídeo manteve distância média de 0,242 ± 0,05 nm da His92 e de 0,227 ± 0,04 nm do Asp241 (bolsão S1) durante toda a simulação, ambas claramente inferiores ao limiar de 0,5 nm adotado para caracterização de contato catalítico direto.

Por outro lado, as distâncias ao Asp114 (1,719 ± 0,18 nm) e à Ser211 (1,195 ± 0,14 nm) permaneceram superiores a 0,5 nm ao longo da maior parte da trajetória, indicando que o QCL936 não interage diretamente com o núcleo Asp-Ser da tríade catalítica. Uma leve aproximação a esses resíduos foi observada após os 80 ns, sugerindo reorientação conformacional tardia do peptídeo, porém sem cruzar o limiar de contato catalítico.

Esse padrão de ligação — ancoragem no bolsão S1 e na His da tríade, com ausência de contato com Asp e Ser — caracteriza o QCL936 como um inibidor com mecanismo predominantemente de bloqueio do reconhecimento de substrato, ocupando o bolsão de especificidade (S1) e impedindo o acesso da região de clivagem do substrato ao sítio ativo, sem necessariamente inibir diretamente o mecanismo catalítico nucleofílico. Este perfil é biologicamente relevante, dado que o bolsão S1 de serino-proteases determina a especificidade pelo substrato (Perona & Craik, 1995), e inibidores que o ocupam de forma estável constituem candidatos potentes e seletivos.

---

<!-- PLACEHOLDER: adicionar resultados de XP273, DN773, DN1441, ACR157 seguindo o mesmo formato -->

---

## Referências (a inserir)

- Abraham, M.J. *et al.* (2015) GROMACS: A high performance molecular dynamics package. *SoftwareX*, 1–2, 19–25.
- Al-Khafaji, K. *et al.* (2024) Molecular dynamics of serine protease–peptide complexes. *Biomolecules*, 14, XXX.
- Bray, S.A. *et al.* (2024) Best practices for MD equilibration of protein–peptide complexes. *J. Phys. Chem. B*, 128, XXX.
- Bussi, G., Donadio, D. & Parrinello, M. (2007) Canonical sampling through velocity rescaling. *J. Chem. Phys.*, 126, 014101.
- Di Tommaso, P. *et al.* (2017) Nextflow enables reproducible computational workflows. *Nat. Biotechnol.*, 35, 316–319.
- Dow, J.A.T. (1992) pH gradients in lepidopteran midgut. *J. Exp. Biol.*, 172, 355–375.
- Hess, B. *et al.* (1997) LINCS: A linear constraint solver for molecular simulations. *J. Comput. Chem.*, 18, 1463–1472.
- Joung, I.S. & Cheatham, T.E. (2008) Determination of alkali and halide monovalent ion parameters for use in explicitly solvated biomolecular simulations. *J. Phys. Chem. B*, 112, 9020–9041.
- Lindorff-Larsen, K. *et al.* (2010) Improved side-chain torsion potentials for the AMBER ff99SB protein force field. *Proteins*, 78, 1950–1958.
- Parrinello, M. & Rahman, A. (1981) Polymorphic transitions in single crystals. *J. Appl. Phys.*, 52, 7182–7190.
- Perona, J.J. & Craik, C.S. (1995) Structural basis of substrate specificity in the serine proteases. *Protein Sci.*, 4, 337–360.
- van Zundert, G.C.P. *et al.* (2016) The HADDOCK2.2 web server. *J. Mol. Biol.*, 428, 720–725.
