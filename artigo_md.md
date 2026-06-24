# Dinâmica Molecular — Seções do Artigo

---

## 2. Metodologia

### 2.1 Preparação dos complexos

As estruturas iniciais dos complexos foram geradas por docking molecular utilizando o servidor HADDOCK 2.4 (van Zundert *et al.*, 2016). Os receptores avaliados correspondem a seis isoformas de tripsina digestiva de *Spodoptera* spp.: ACR157, QCL936, XP273, XP352, DN773 e DN1441. Foram estudadas duas séries de inibidores: (i) **série GORE4** — o peptídeo GORE4 de cinco resíduos, com todos os resíduos diretamente envolvidos na interface de ligação, modelado por docking proteína–peptídeo; (ii) **série SKTI** — o Inibidor de Tripsina Kunitz de Soja (*Soybean Kunitz Trypsin Inhibitor*, SKTI; 177 resíduos), inibidor natural de referência com alça reativa SPYRIRF (P1 = Arg), modelado por docking proteína–proteína. Para cada sistema, o cluster de maior escore energético foi selecionado como pose inicial para as simulações.

Adicionalmente, a benzamidina (BEN; MW 120,15 Da) foi incluída como controle positivo estrutural, docada com AutoDock Vina 1.2.0 contra os snapshots finais pós-DM de cada receptor (protocolo detalhado em §3.4).

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
- **Distâncias mínimas** entre o ligante e cada resíduo do sítio catalítico (His, Asp, Ser da tríade e resíduo do bolsão S1), calculadas pelo grupo de átomos de cada resíduo via `gmx mindist`.

As médias móveis foram calculadas com janela de 5 ns para suavização dos perfis temporais. Todos os gráficos foram gerados com Matplotlib 3.x.

Os resíduos do sítio catalítico monitorados por isoforma de receptor são: ACR157 — His69, Asp114, Ser211, Ile205; QCL936 — His92, Asp142, Ser247, Asp241; XP273 — Tyr83, Asp132, Ser234, Ile229; XP352 — Arg112, Asp166, Ser268, Asp262. Para a série SKTI, os mesmos resíduos do receptor foram monitorados, permitindo comparação direta entre as duas séries de inibidores.

O pipeline completo, desde a preparação dos complexos até a geração dos painéis de análise, foi implementado em Nextflow DSL2 (Di Tommaso *et al.*, 2017) para garantir reprodutibilidade e rastreabilidade computacional.

### 2.7 Fingerprints de interação molecular (ProLIF)

Para quantificar a persistência temporal e tipificação atômica das interações receptor–ligante ao longo das trajetórias, foi utilizado o pacote ProLIF 2.x (*Protein–Ligand Interaction Fingerprints*; Bouysset & Fiorucci, 2021), integrado ao MDAnalysis 2.10.x (Michaud-Agrawal *et al.*, 2011). Para cada sistema, foram identificadas todas as interações do tipo Doadora de HB, Aceptora de HB, Hidrofóbica, Catiônica, Aniônica, Empilhamento π–π, Cátion–π e Contato VdW entre pares resíduo(ligante)–resíduo(receptor), dentro de raio de corte de 0,4 nm. A análise foi conduzida com passo de 10 frames (amostragem efetiva de 1 frame/ns) sobre a trajetória completa de produção. A persistência de cada interação foi definida como a porcentagem de frames em que o critério geométrico da interação é satisfeito.

A análise ProLIF foi aplicada aos sete sistemas estáveis da série peptídeo (três GORE4 e três SKTI), além do complexo XP352-GORE4 (para confirmação de instabilidade). A série BEN foi excluída desta análise por restrição do conversor molecular RDKit: o arquivo de índice `lig.ndx` gerado pelo protocolo GAFF2 para a benzamidina contém exclusivamente os 12 átomos pesados (sem hidrogênios explícitos na topologia de índice), impedindo a conversão molecular requerida pelo ProLIF sem perda de integridade das interações dependentes de H (doadora/aceptora de HB).

### 2.8 Mapas de contato resíduo × resíduo

Mapas de frequência de contato resíduo × resíduo foram calculados para todos os 11 sistemas — incluindo a série BEN — utilizando MDAnalysis 2.10.x. Para cada par resíduo(receptor)–resíduo(ligante), a frequência de contato foi definida como a fração de frames em que ao menos um par atômico inter-residual apresenta distância < 0,4 nm. Os mapas foram representados como *heatmaps* (Seaborn 0.13.x, *clustermap* com agrupamento hierárquico, ligação Ward), permitindo identificação visual dos resíduos mais frequentemente engajados na interface ao longo de toda a trajetória.

---

## 3. Resultados e Discussão

### 3.1 Estabilidade estrutural dos complexos — série GORE4

Os complexos das isoformas de tripsina de *Spodoptera* com o peptídeo GORE4 foram avaliados quanto à estabilidade estrutural, integridade da interface de ligação e posicionamento no sítio catalítico. Os resultados são apresentados por isoforma receptora.

#### QCL936-GORE4 (cluster 3)

A simulação do complexo QCL936-GORE4 (cluster 3) foi conduzida por 100 ns a 300 K e pH 8,2 e apresentou comportamento estável ao longo de toda a trajetória. O RMSD do backbone do receptor atingiu um platô de 0,124 ± 0,017 nm após aproximadamente 20 ns (Figura X), indicando equilíbrio conformacional do receptor sem deriva estrutural. O raio de giro manteve-se em 1,718 ± 0,006 nm com variação mínima (< 1%), confirmando que a estrutura terciária da tripsina não sofreu perturbações decorrentes da presença do peptídeo.

O RMSD do ligante estabilizou em 0,229 ± 0,042 nm, valor inferior a 0,3 nm durante praticamente toda a simulação, indicando que o peptídeo QCL936 manteve sua conformação de docking ao longo dos 100 ns sem sinais de dissociação do sítio de ligação. Não foi observado o padrão de dissociação característico (colapso de contatos seguido de aumento abrupto das distâncias), o que corrobora a estabilidade do complexo.

#### Interface de ligação

O número médio de contatos receptor–ligante foi de 340 ± 41 átomos a menos de 0,4 nm, valor expressivamente superior ao limiar mínimo de 150 átomos adotado como critério de ligação estável. O perfil temporal dos contatos mostrou oscilações regulares sem tendência de redução, indicando interface robusta ao longo de toda a trajetória.

As pontes de hidrogênio receptor–ligante apresentaram média de 2,75 ± 1,01, com persistência contínua ao longo dos 100 ns. A manutenção de pelo menos uma ponte de hidrogênio em praticamente todos os frames demonstra interação específica entre o peptídeo e o receptor, contribuindo para a estabilização do complexo.

A SASA do ligante manteve-se em 9,04 ± 0,37 nm² sem tendência crescente, indicando que o peptídeo permaneceu enterrado na interface receptor–ligante e não migrou para o solvente.

#### Posicionamento no sítio catalítico

A análise das distâncias mínimas entre o peptídeo QCL936 e os resíduos do sítio ativo revelou um padrão de ligação predominantemente associado ao bolsão S1 e à His periférica da tríade (Figura X). O peptídeo manteve distância média de 0,242 ± 0,05 nm da His92 e de 0,227 ± 0,04 nm do Asp241 (bolsão S1) durante toda a simulação, ambas claramente inferiores ao limiar de 0,5 nm adotado para caracterização de contato catalítico direto.

Por outro lado, as distâncias ao Asp142 (1,719 ± 0,18 nm) e à Ser247 (1,195 ± 0,14 nm) permaneceram superiores a 0,5 nm ao longo da maior parte da trajetória, indicando que o GORE4 não interage diretamente com o núcleo Asp-Ser da tríade catalítica do QCL936. Uma leve aproximação a esses resíduos foi observada após os 80 ns, sugerindo reorientação conformacional tardia do peptídeo, porém sem cruzar o limiar de contato catalítico.

Esse padrão de ligação — ancoragem no bolsão S1 e na His da tríade, com ausência de contato com Asp e Ser — caracteriza o QCL936 como um inibidor com mecanismo predominantemente de bloqueio do reconhecimento de substrato, ocupando o bolsão de especificidade (S1) e impedindo o acesso da região de clivagem do substrato ao sítio ativo, sem necessariamente inibir diretamente o mecanismo catalítico nucleofílico. Este perfil é biologicamente relevante, dado que o bolsão S1 de serino-proteases determina a especificidade pelo substrato (Perona & Craik, 1995), e inibidores que o ocupam de forma estável constituem candidatos potentes e seletivos.

---

#### ACR157-GORE4 (cluster 1)

A simulação do complexo ACR157-GORE4 (cluster 1) foi conduzida por 100 ns a 300 K e pH 8,2. O RMSD do backbone do receptor estabilizou em 0,165 ± 0,016 nm após os primeiros 20 ns, valor ligeiramente superior ao do QCL936 (0,124 nm), porém ainda indicativo de estabilidade estrutural satisfatória. O raio de giro manteve-se em 1,863 ± 0,007 nm com variação inferior a 1%, confirmando a integridade da estrutura terciária da tripsina ao longo de toda a trajetória.

O RMSD do ligante apresentou média de 0,393 ± 0,068 nm, consideravelmente superior ao do QCL936 (0,229 nm), refletindo maior mobilidade conformacional do peptídeo ACR157 no sítio de ligação. Não obstante, o peptídeo permaneceu estável na interface receptor–ligante sem sinais de dissociação ao longo dos 100 ns.

#### Interface de ligação

O número médio de contatos receptor–ligante foi de 255,9 ± 37,6 átomos a menos de 0,4 nm, valor inferior ao do QCL936 (340 ± 41), porém acima do limiar de 150 átomos adotado como critério de ligação estável. O perfil temporal dos contatos revelou redução inicial de aproximadamente 350 para 250 átomos nos primeiros 30 ns, seguida de estabilização, sugerindo reorganização da interface durante a fase inicial da simulação de produção.

As pontes de hidrogênio receptor–ligante apresentaram média de 1,90 ± 1,38, com persistência ao longo de toda a trajetória, embora com maior variabilidade em relação ao QCL936 (2,75 ± 1,01). A SASA do ligante manteve-se em 9,9 ± 0,47 nm² sem tendência crescente, indicando que o peptídeo permaneceu parcialmente enterrado na interface ao longo da simulação.

#### Posicionamento no sítio catalítico

A análise das distâncias mínimas entre o peptídeo ACR157 e os resíduos do sítio ativo revelou um padrão de ancoragem centrado na His69 e no Ile205 (bolsão S1). A His69 manteve distância média abaixo de 0,5 nm durante a maior parte da trajetória, caracterizando contato direto com a histidina catalítica. O Ile205 também apresentou engajamento consistente com o bolsão S1, com distância média na faixa de 0,4–0,5 nm.

O Asp114 — componente do díade Asp-Ser da tríade catalítica — manteve-se distante da interface ao longo de toda a trajetória, com média superior a 1,0 nm. A Ser211 (Ser nucleofílica) apresentou comportamento intermediário, com distâncias oscilando entre 0,5 e 0,9 nm, sem estabelecimento de contato estável.

O padrão de ligação do ACR157 guarda similaridade qualitativa com o do QCL936 — ancoragem predominante na histidina da tríade e no bolsão S1, com ausência de interação direta com o par Asp–Ser — porém com interface menos compacta (menor número de contatos e maior mobilidade do ligante). O engajamento da His69 pode comprometer a rede de transferência de próton da tríade catalítica (His69–Asp114–Ser211), enquanto a ocupação do bolsão S1 bloqueia o acesso do substrato ao sítio de clivagem. Assim, o ACR157 apresenta potencial mecanismo dual: interferência parcial na ativação da Ser nucleofílica e bloqueio de especificidade de substrato, embora com eficiência de ligação aparentemente inferior ao QCL936 nos parâmetros quantitativos analisados.

---

#### XP352-GORE4 (cluster 4, rank 3)

A simulação do complexo XP352-GORE4 (cluster 4, rank 3) revelou comportamento marcadamente distinto dos complexos QCL936-c3 e ACR157-c1, com evidências de instabilidade estrutural e dissociação do ligante ao longo dos 100 ns.

O RMSD do backbone do receptor apresentou valor médio de 0,886 ± 0,440 nm, aproximadamente sete vezes superior ao observado para o QCL936 (0,124 nm), sem atingir platô de estabilização ao longo da trajetória. O perfil crescente e sem convergência indica desvio conformacional progressivo do receptor, sugestivo de rearranjo estrutural significativo. O raio de giro (1,837 ± 0,074 nm) também exibiu variabilidade substancialmente maior do que os complexos estáveis (σ = 0,074 nm versus 0,006–0,007 nm para QCL936 e ACR157), corroborando a instabilidade conformacional do receptor neste complexo.

#### Interface de ligação e sítio catalítico

O número médio de contatos receptor–ligante foi de 62,7 ± 142,9 átomos — inferior ao limiar de 150 átomos e com desvio padrão superior à própria média, padrão diagnóstico de dissociação intermitente. O perfil temporal revelou comportamento bimodal: períodos prolongados com contatos próximos de zero, intercalados com excursões transitórias de retorno à interface, sem estabelecimento de ligação estável.

A análise das distâncias mínimas ao sítio catalítico confirmou o sinal de dissociação: os resíduos Arg112, Asp166, Ser268 e Asp262 apresentaram distâncias médias superiores a 1,0 nm durante a maior parte da trajetória, com tendência progressivamente crescente. O gráfico de ocupação da tríade mostrou todos os quatro resíduos catalíticos com distâncias médias acima de 0,5 nm, sem contato direto estabelecido em nenhum deles. A análise de RMSF confirmou alta mobilidade do ligante (pico de flutuação destacado no painel), consistente com peptídeo não ancorado no sítio de ligação.

#### Diagnóstico

O complexo XP352-GORE4 (cluster 4, rank 3) não apresentou ligação estável sob as condições simuladas (300 K, pH 8,2, 100 ns). A qualidade do pose inicial — rank 3 dentro do cluster 4, representando estrutura de menor escore energético — é um fator determinante provável para o insucesso, uma vez que o cluster 4, rank 1 (analisado anteriormente) pode apresentar comportamento distinto. A ausência de ancoragem estável a qualquer resíduo catalítico, combinada com o elevado RMSD do backbone e a dissociação evidenciada pelo colapso dos contatos, posiciona este complexo como candidato não viável nas condições avaliadas.

---

#### XP273-GORE4 (cluster 1, _NEW)

A simulação do complexo XP273-GORE4 (cluster 1) foi conduzida por 100 ns a 300 K e pH 8,2 e demonstrou comportamento estável ao longo de toda a trajetória. O RMSD do backbone do receptor estabilizou em 0,248 ± 0,045 nm após os primeiros 20 ns, permanecendo abaixo de 0,3 nm na maior parte da trajetória, o que indica estabilidade estrutural satisfatória do receptor. Notavelmente, o raio de giro apresentou perfil bimodal (1,727 ± 0,023 nm), com uma transição conformacional do receptor em torno de 30–60 ns seguida de retorno próximo ao estado inicial, sugerindo amostragem entre dois estados conformacionais levemente distintos durante a simulação. O RMSD do ligante foi de 0,317 ± 0,068 nm, valor intermediário entre o QCL936 (0,229 nm) e o ACR157 (0,393 nm), indicando mobilidade moderada no sítio de ligação.

#### Interface de ligação

O número médio de contatos receptor–ligante foi de 260 ± 60 átomos, valor acima do limiar de 150 átomos com variabilidade relativa normal (SD/média ≈ 23%), confirmando manutenção estável da interface ao longo dos 100 ns. As pontes de hidrogênio receptor–ligante apresentaram média de 3,19 ± 1,07, valor superior ao do ACR157 (1,90) e comparável ao do QCL936 (2,75), indicando interações específicas e persistentes com o receptor. A SASA do ligante foi de 8,14 ± 0,37 nm², o menor valor entre todos os sistemas analisados, indicando o melhor grau de enterramento do peptídeo na interface receptor–ligante.

#### Posicionamento no sítio catalítico

A análise das distâncias mínimas revelou um padrão de ancoragem distinto dos demais complexos. O XP273 manteve contato próximo com a Tyr83 (resíduo periférico do sítio ativo, triad_1), com distância média abaixo de 0,5 nm durante grande parte da trajetória. Os resíduos Asp132 e Ser234 — componentes do díade Asp-Ser da tríade catalítica — apresentaram contato borderline com distâncias oscilando em torno de 0,5–0,65 nm, sem contato direto estável estabelecido. O Ile229 (bolsão S1, triad_4) manteve-se distante da interface, com distância média superior a 1,0 nm ao longo de toda a trajetória.

Este padrão contrasta com os complexos QCL936 e ACR157, que apresentaram engajamento consistente do bolsão S1. O XP273 ancora predominantemente na região periférica do sítio ativo (Tyr83), sem bloquear o bolsão de especificidade S1. Esse mecanismo pode atuar por oclusão estérica da entrada do sítio ativo ou por estabilização de uma conformação cataliticamente inativa do receptor — hipótese suportada pelo comportamento bimodal do raio de giro. Embora o SASA do ligante indique enterramento expressivo, a ausência de contato com S1 e a interação apenas borderline com Asp132/Ser234 sugerem que a potência inibitória do XP273 pode ser inferior à do QCL936, que combina ancoragem na tríade com bloqueio direto do bolsão S1.

---

<!-- PLACEHOLDER: adicionar resultados de DN773-GORE4 e DN1441-GORE4 seguindo o mesmo formato -->

---

### 3.2 Estabilidade estrutural dos complexos — série SKTI

Os mesmos receptores foram avaliados em complexo com o SKTI, inibidor natural de referência. A alça reativa SPYRIRF do SKTI (P1 = Arg) é estruturalmente projetada para encaixar no bolsão S1 das serino-proteases, tornando esses complexos referências funcionais para comparação com os dados da série GORE4.

#### ACR157-SKTI (cluster 2)

A simulação do complexo ACR157-SKTI (cluster 2) foi conduzida por 100 ns a 300 K e pH 8,2 e apresentou comportamento estável em toda a trajetória. O RMSD do backbone do receptor atingiu platô de 0,281 ± 0,027 nm após os primeiros 20 ns, valor comparável ao do complexo ACR157-GORE4-c1 (0,165 nm), confirmando que a presença de um ligante proteico volumoso (177 resíduos) não desestabiliza a estrutura terciária do receptor. O raio de giro do complexo (receptor + SKTI) estabilizou-se em 2,225 ± 0,009 nm com variabilidade inferior a 0,5%, indicando que a estrutura quaternária do heterodímero foi mantida compacta e sem eventos de abertura ou desprendimento do ligante ao longo dos 100 ns.

O RMSD do ligante estabilizou em 0,206 ± 0,017 nm — valor notavelmente inferior ao do GORE4 no mesmo receptor (0,393 nm) —, demonstrando que o SKTI permaneceu em conformação praticamente rígida dentro do sítio de ligação. Esse resultado é consistente com a natureza estrutural dos inibidores Kunitz: proteínas globulares estabilizadas por pontes dissulfeto com alça reativa pré-organizada para interação com o sítio catalítico (Laskowski & Kato, 1980). A análise de RMSF por resíduo confirmou flutuação muito baixa ao longo da cadeia do SKTI (< 0,05 nm na maior parte dos resíduos), com picos restritos a regiões de alça periférica, sem mobilidade significativa na região da alça reativa SPYRIRF.

#### Interface de ligação

O número médio de contatos receptor–ligante foi de 1019 ± 118 átomos a menos de 0,4 nm, valor aproximadamente quatro vezes superior ao do GORE4 no mesmo receptor (256 átomos). A magnitude desta diferença é proporcional ao tamanho relativo dos ligantes (177 resíduos versus 5 resíduos) e indica que o SKTI estabelece uma interface de enterramento extensa, cobrindo porção substancial da superfície do sítio ativo do ACR157. O perfil temporal dos contatos mostrou oscilações regulares sem tendência de redução, característico de interface estabilizada.

As pontes de hidrogênio receptor–ligante apresentaram média de 13,524 ± 2,911, valor cerca de sete vezes superior ao registrado para o GORE4 (1,90 ± 1,38). A manutenção de uma rede de pontes de hidrogênio persistente ao longo de toda a trajetória é um marcador da complementaridade de superfície característica dos inibidores Kunitz canônicos, que formam complexos estáveis tipo substrato-análogo com serino-proteases (Krowarsch *et al.*, 2003). A SASA do ligante manteve-se em 101,97 ± 1,62 nm² com variação mínima ao longo dos 100 ns, confirmando que o SKTI não sofreu exposição de sua superfície ao solvente durante a simulação.

#### Posicionamento no sítio catalítico

A análise das distâncias mínimas entre o SKTI e os resíduos do sítio ativo revelou engajamento direto e simultâneo de dois dos três componentes da tríade catalítica. A His69 manteve distância média de 0,224 nm ao longo de toda a trajetória, e a Ser211 apresentou distância média de 0,221 nm — ambas claramente abaixo do limiar de 0,5 nm para contato catalítico direto, com oscilações mínimas que demonstram estabilidade do posicionamento. O Asp114 apresentou distância de 0,312 nm, compatível com contato de segundo plano que constrainge a geometria da tríade. O Ile205 (bolsão S1) manteve distância média de 0,180 nm — engajamento direto do bolsão de especificidade, confirmando bloqueio integral do sítio ativo.

Este padrão de ligação — engajamento simultâneo de todos os quatro resíduos catalíticos (His69, Asp114, Ser211 e Ile205/S1) — é o mecanismo Kunitz canônico completo: a alça reativa SPYRIRF insere seu resíduo P1 (Arg) diretamente no sítio ativo, formando um complexo enzima–inibidor análogo ao intermediário tetraédrico da catálise, sem progressão da hidrólise (Krowarsch *et al.*, 2003; Laskowski & Kato, 1980). A His69 não consegue abstrair o próton da Ser211 (necessário para a ativação nucleofílica) porque a geometria da rede de transferência de próton é bloqueada pela presença do resíduo P1 do inibidor, e o bolsão S1 é simultaneamente ocupado pelo P1=Arg.

#### Comparação com ACR157-GORE4 e mecanismo diferencial

O contraste entre ACR157-SKTI e ACR157-GORE4 é elucidativo. No complexo GORE4, His69 estava engajada mas Ser211 mantinha distâncias variáveis (0,5–0,9 nm), sem engajamento estável da serina nucleofílica. No complexo SKTI, todos os quatro resíduos catalíticos estão em contato simultâneo — uma diferença mecanisticamente fundamental. O GORE4 atua predominantemente por interferência periférica na tríade e bloqueio do acesso ao sítio (mecanismo parcial), enquanto o SKTI bloqueia diretamente o núcleo catalítico His–Asp–Ser e o bolsão S1 concomitantemente (mecanismo canônico Kunitz completo).

Adicionalmente, a interface ACR157-SKTI é significativamente mais estável em todos os parâmetros quantitativos: RMSD do ligante 47,6% inferior (0,206 vs 0,393 nm), contatos 3,98× superiores (1019 vs 256) e pontes de hidrogênio 7,1× superiores (13,524 vs 1,90). Os resultados posicionam o SKTI como inibidor de referência efetivo contra ACR157, com mecanismo de inibição integral do sítio catalítico, e fornecem o padrão de comparação funcional para avaliação da potência dos peptídeos da série GORE4 frente às isoformas de tripsina de *Spodoptera*.

#### QCL936-SKTI (cluster 2)

A simulação do complexo QCL936-SKTI (cluster 2) foi conduzida por 100 ns a 300 K e pH 8,2. O RMSD do backbone do receptor apresentou média de 0,370 ± 0,054 nm, valor markadamente superior ao observado para o mesmo receptor em complexo com o GORE4 (0,124 nm) e o mais elevado entre os três complexos SKTI analisados (ACR157: 0,281 nm; XP273: 0,288 nm), indicando que a presença do SKTI induz maior flexibilidade conformacional no QCL936. O perfil temporal não exibiu tendência crescente consistente ao longo da trajetória, e o raio de giro do complexo estabilizou-se em 2,296 ± 0,015 nm, confirmando que o complexo como um todo permaneceu globalmente coeso.

O RMSD do ligante atingiu 0,209 ± 0,025 nm — valor próximo ao registrado para o ACR157-SKTI (0,206 nm) e notavelmente inferior ao do GORE4 no mesmo receptor (0,229 nm), evidenciando que o SKTI manteve conformação estável dentro do sítio de ligação sem sinais de dissociação.

#### Interface de ligação

O número médio de contatos receptor–ligante foi de 720 ± 77 átomos, valor expressivo embora inferior ao registrado para o ACR157-SKTI (1019 ± 118). As pontes de hidrogênio receptor–ligante apresentaram média de 10,410 ± 2,022, também inferior ao ACR157-SKTI (13,524 ± 2,911), indicando menor complementaridade de superfície entre SKTI e QCL936 em comparação com ACR157. A SASA do ligante manteve-se em 102,53 ± 2,09 nm² ao longo dos 100 ns, confirmando enterramento estável do SKTI na interface receptor–ligante.

#### Posicionamento no sítio catalítico

A análise das distâncias mínimas revelou um padrão de ancoragem envolvendo três dos quatro resíduos monitorados: a His92 manteve distância média de 0,179 nm, a Ser247 (serina nucleofílica) apresentou distância média de 0,176 nm e o Asp241 (bolsão S1) atingiu 0,178 nm — todos em contato íntimo e estável ao longo de toda a trajetória. Em contraste, o Asp142 — componente do díade catalítico responsável pela estabilização da His protonada durante o ciclo de transferência de próton — permaneceu distante da interface (~0,73 nm), sem contato estabelecido com o ligante em nenhuma fase da simulação.

Este padrão — His92 + Ser247 + Asp241/S1, com ausência de contato com Asp142 — constitui um modo de inibição até então não observado entre os complexos analisados. Contrasta diretamente com o QCL936-GORE4, no qual apenas His92 e Asp241 foram engajados (modo His+S1, 2/4), revelando que o SKTI — graças à sua maior extensão e à pré-organização estrutural da alça SPYRIRF — consegue posicionar seu resíduo P1=Arg de forma a engajar adicionalmente a Ser247 nucleofílica. O engajamento simultâneo de His92 e Ser247 indica que o SKTI interfere diretamente com a geometria do par nucleofílico da tríade, bloqueando a etapa de ativação da serina, enquanto o bolsão S1 permanece simultaneamente ocupado.

#### Comparação com ACR157-SKTI e implicações mecanísticas

O contraste entre QCL936-SKTI e ACR157-SKTI é mecanisticamente revelador. No ACR157, o SKTI engaja diretamente His69 (0,224 nm), Asp114 (0,312 nm), Ser211 (0,221 nm) e Ile205/S1 (0,180 nm) — mecanismo Kunitz canônico completo (4/4). No QCL936, o SKTI engaja His92, Ser247 e Asp241/S1, mas não o Asp142 — modo His + Ser + S1 (3/4). Embora ambos os complexos apresentem engajamento simultâneo da His e da Ser, a ausência do Asp142 no QCL936 é mecanisticamente significativa: na catálise de serino-proteases, o Asp catalítico estabiliza a His protonada após a abstração do próton da Ser, permitindo a regeneração da tríade para novos ciclos; sem o bloqueio estrutural do Asp142, esta etapa de recarregamento permanece potencialmente acessível, o que pode conferir ao complexo QCL936-SKTI menor potência de inibição em relação ao modo canônico.

A ausência de engajamento do Asp142 pode ser atribuída à geometria diferencial entre os sítios ativos das duas isoformas. A presença de Asp241 carregado negativamente no bolsão S1 do QCL936 — em contraste com Ile205 neutro no ACR157 — exerce forte atração eletrostática sobre P1=Arg, direcionando a alça SPYRIRF preferencialmente para S1 e His92 com geometria que alcança a Ser247, mas com o Asp142 posicionado fora do raio de contato.

Em termos quantitativos, a interface QCL936-SKTI é menos complementar do que a ACR157-SKTI (720 vs 1019 contatos; 10,4 vs 13,5 H-bonds), sugerindo que o SKTI inibe ACR157 com maior eficiência de interface do que o QCL936. Ambos os complexos SKTI são, contudo, substancialmente mais estáveis e com interfaces mais extensas do que os respectivos complexos GORE4.

#### XP273-SKTI (cluster 2)

A simulação do complexo XP273-SKTI (cluster 2) foi conduzida por 100 ns a 300 K e pH 8,2 e apresentou comportamento estável em toda a trajetória. O RMSD do backbone do receptor atingiu 0,288 ± 0,042 nm, valor comparável ao do ACR157-SKTI (0,281 nm) e ao QCL936-SKTI (0,307 nm), indicando que a presença do ligante proteico volumoso não desestabilizou o receptor. O raio de giro estabilizou-se em 2,325 ± 0,015 nm com variabilidade inferior a 1%, confirmando que o heterodímero XP273–SKTI permaneceu globalmente coeso ao longo de toda a trajetória.

O RMSD do ligante apresentou média de 0,224 ± 0,017 nm — valor notavelmente estável, comparável ao do ACR157-SKTI (0,206 nm) e inferior ao do GORE4 no mesmo receptor (0,317 nm) —, demonstrando que o SKTI adota conformação praticamente rígida dentro do sítio de ligação do XP273, comportamento consistente com a natureza pré-organizada dos inibidores Kunitz.

#### Interface de ligação

O número médio de contatos receptor–ligante foi de 596 ± 88 átomos, valor inferior ao registrado para os complexos ACR157-SKTI (1019 ± 118) e QCL936-SKTI (720 ± 77), porém expressivamente superior ao do peptídeo GORE4 no mesmo receptor (260 ± 60). As pontes de hidrogênio receptor–ligante apresentaram média de 8,817 ± 3,496, também inferiores às dos demais complexos SKTI (ACR157: 13,524; QCL936: 10,410), indicando menor complementaridade de superfície entre o SKTI e o sítio de ligação do XP273 em comparação com as demais isoformas. A SASA do ligante manteve-se em 101,312 ± 2,014 nm² ao longo dos 100 ns, confirmando enterramento estável do SKTI na interface e ausência de migração para o solvente.

#### Posicionamento no sítio catalítico

A análise das distâncias mínimas entre o SKTI e os resíduos do sítio ativo revelou um padrão de ancoragem intermediário entre o modo canônico Kunitz observado no ACR157 e o modo His+S1 observado no QCL936. A Tyr83 (triad_1, resíduo periférico monitorado por alinhamento estrutural) manteve distância média de 0,098 nm — contato íntimo e estável ao longo de toda a trajetória. O Asp132 (triad_2) apresentou distância média de 0,119 nm, igualmente em contato direto com o ligante. O Ile229 (bolsão S1, triad_4) manteve distância de 0,118 nm, caracterizando engajamento simultâneo do bolsão de especificidade, padrão ausente no complexo XP273-GORE4.

Em contraste, a Ser234 (triad_3, serina nucleofílica) apresentou distância média de 0,718 nm — superior ao limiar de 0,5 nm e consistente com ausência de contato direto ao longo de toda a simulação. Assim, o padrão de ligação do XP273-SKTI envolve três dos quatro resíduos monitorados (Tyr83 + Asp132 + Ile229), com a serina nucleofílica como único componente não engajado.

#### Interpretação mecanística e comparação com os demais complexos SKTI

O padrão de ancoragem Tyr83 + Asp132 + Ile229(S1), com ausência de contato com Ser234, constitui um modo de inibição parcial distinto dos demais complexos analisados. No ACR157-SKTI, todos os quatro resíduos foram engajados (mecanismo Kunitz canônico completo, 4/4). No QCL936-SKTI, His92, Ser247 e Asp241/S1 foram engajados, sem contato com Asp142 (modo His+Ser+S1, 3/4). O XP273-SKTI ocupa posição complementar: o SKTI engaja Tyr83, o Asp132 do díade catalítico e o bolsão S1 (Ile229), mas não alcança a Ser234 (serina nucleofílica), modo Tyr+Asp+S1 (3/4).

A ausência de contato com Ser234 distingue este complexo do mecanismo Kunitz canônico clássico e pode ser atribuída à geometria particular do sítio ativo do XP273. A presença de Tyr83 no lugar da His canônica na posição triad_1 — registro atípico entre as isoformas analisadas — altera a topografia da região de entrada do sítio catalítico e pode reorientar o loop reativo do SKTI de modo que o resíduo P1 (Arg) se posicione favoravelmente para contato com Tyr83, Asp132 e Ile229, mas com geometria incompatível com o alcance simultâneo da Ser234.

Do ponto de vista funcional, o engajamento de Asp132 — componente do díade Asp–Ser da tríade catalítica — representa um avanço mecanístico em relação ao XP273-GORE4, que alcançava apenas Tyr83 sem contato com Asp132 ou Ile229. O bloqueio concomitante de Asp132 e do bolsão S1 (Ile229) pode comprometer a rede de transferência de próton do díade catalítico e impedir o reconhecimento do substrato pelo bolsão de especificidade, caracterizando um mecanismo de inibição dual parcial. Entretanto, a ausência de engajamento da serina nucleofílica limita a comparação com o mecanismo Kunitz completo do ACR157-SKTI.

Quantitativamente, a interface XP273-SKTI é a menos complementar entre os três complexos SKTI analisados (596 vs 720 vs 1019 contatos; 8,8 vs 10,4 vs 13,5 H-bonds, para XP273, QCL936 e ACR157, respectivamente), sugerindo hierarquia de eficiência de inibição ACR157 > QCL936 > XP273. Esta hierarquia é inteiramente consistente com os escores de docking HADDOCK (ACR157: −128,3; QCL936: −113,8; XP273: −101,4 kcal mol⁻¹), validando retroativamente a seleção de clusters por escore energético como preditor de qualidade de interface em dinâmica molecular.

<!-- PLACEHOLDER: adicionar XP352-SKTI quando disponível -->

---

### 3.3 Síntese comparativa — estabilidade, interface e modo de inibição

#### Tabela 1 — Parâmetros de dinâmica molecular (médias ± DP)

| Sistema | Duração | RMSD bb (nm) | RMSD lig (nm) | Contatos | H-bonds | SASA lig (nm²) | Status |
|---------|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| **Série GORE4** | | | | | | | |
| QCL936-GORE4 c3 | 100 ns | 0,124 ± 0,017 | 0,229 ± 0,042 | 340 ± 41 | 2,75 ± 1,01 | 9,04 ± 0,37 | ✅ estável |
| ACR157-GORE4 c1 | 100 ns | 0,165 ± 0,016 | 0,393 ± 0,068 | 256 ± 38 | 1,90 ± 1,38 | 9,90 ± 0,47 | ✅ estável |
| XP273-GORE4 c1 | 100 ns | 0,248 ± 0,045 | 0,317 ± 0,068 | 260 ± 60 | 3,19 ± 1,07 | 8,14 ± 0,37 | ✅ estável |
| XP352-GORE4 c4r3 | 100 ns | 0,886 ± 0,440 | 0,400 ± 0,067 | 63 ± 143 | — | — | ❌ dissociação |
| **Série SKTI** | | | | | | | |
| ACR157-SKTI c2 | 100 ns | 0,281 ± 0,027 | 0,206 ± 0,017 | 1019 ± 118 | 13,524 ± 2,911 | 101,97 ± 1,62 | ✅ estável |
| QCL936-SKTI c2 | 100 ns | 0,370 ± 0,054 | 0,209 ± 0,025 | 720 ± 77 | 10,410 ± 2,022 | 102,53 ± 2,09 | ✅ estável |
| XP273-SKTI c2 | 100 ns | 0,288 ± 0,042 | 0,224 ± 0,017 | 596 ± 88 | 8,817 ± 3,496 | 101,31 ± 2,01 | ✅ estável |
| **Série BEN (benzamidina)** | | | | | | | |
| ACR157-BEN | 200 ns | 0,146 ± 0,019 | N.D.ª | 68,8 ± 72,4 | N.D.ᵇ | 2,74 ± 0,141 | ❌ dissociação ~95 ns |
| QCL936-BEN | 200 ns | 0,116 ± 0,018 | N.D.ª | 90,1 ± 41,6 | N.D.ᵇ | 2,78 ± 0,136 | ❌ dissociação ~150 ns |
| XP273-BEN | 200 ns | 0,152 ± 0,016 | N.D.ª | 55,5 ± 38,2 | N.D.ᵇ | 2,78 ± 0,143 | ❌ dissociação ~80 ns |
| XP352-BEN | 200 ns | 0,213 ± 0,021 | N.D.ª | 52,2 ± 47,1 | N.D.ᵇ | 2,78 ± 0,140 | ❌ dissociação ~125 ns |

ª RMSD interno de molécula rígida (não posicional); ᵇ `gmx hbond` incompatível com campo de força GAFF2.

#### Tabela 2 — Distâncias médias (nm) aos resíduos do sítio catalítico

| Sistema | triad_1 (His/Tyr/Arg) | triad_2 (Asp) | triad_3 (Ser) | triad_4 (S1) | Modo de inibição |
|---------|:---:|:---:|:---:|:---:|:---|
| **GORE4** | | | | | |
| QCL936-c3 | His92 **0,242** ✅ | Asp142 1,72 ❌ | Ser247 1,20 ❌ | Asp241 **0,227** ✅ | His + S1 |
| ACR157-c1 | His69 **<0,50** ✅ | Asp114 >1,00 ❌ | Ser211 0,5–0,9 ⚠️ | Ile205 ~0,45 ⚠️ | His + S1 parcial |
| XP273-c1 | Tyr83 **<0,50** ✅ | Asp132 0,55–0,65 ⚠️ | Ser234 0,55–0,65 ⚠️ | Ile229 >1,00 ❌ | Tyr periférico |
| XP352-c4r3 | Arg112 >1,00 ❌ | Asp166 >1,00 ❌ | Ser268 >1,00 ❌ | Asp262 >1,00 ❌ | — dissociação |
| **SKTI** | | | | | |
| ACR157-c2 | His69 **0,224** ✅ | Asp114 **0,312** ✅ | Ser211 **0,221** ✅ | Ile205 **0,180** ✅ | Kunitz canônico (4/4) |
| QCL936-c2 | His92 **0,179** ✅ | Asp142 ~0,73 ❌ | Ser247 **0,176** ✅ | Asp241 **0,178** ✅ | His + Ser + S1 (3/4) |
| XP273-c2 | Tyr83 **0,098** ✅ | Asp132 **0,119** ✅ | Ser234 0,718 ❌ | Ile229 **0,118** ✅ | Tyr + Asp + S1 (3/4) |
| **BEN** | | | | | |
| ACR157-BEN | His69 ~3,5ᶜ ❌ | Asp114 ~3,5ᶜ ❌ | Ser211 ~3,0ᶜ ❌ | Ile205 ~3,5ᶜ ❌ | Dissociação (~95 ns) |
| QCL936-BEN | His92 ~1,2ᵈ ❌ | Asp142 ~1,5ᵈ ❌ | Ser247 ~1,0ᵈ ❌ | Asp241 ~1,3ᵈ ❌ | Dissociação (~150 ns)ᵉ |
| XP273-BEN | Tyr83 ~2,5ᶠ ❌ | Asp132 ~2,5ᶠ ❌ | Ser234 ~2,5ᶠ ❌ | Ile229 ~2,0ᶠ ❌ | Dissociação (~80 ns) |
| XP352-BEN | Arg112 ~2,0ᵍ ❌ | Asp166 ~2,0ᵍ ❌ | Ser268 ~2,5ᵍ ❌ | Asp262 ~1,5ᵍ ⚠️ | Dissociação (~125 ns) |

ᶜ Médias sobre 200 ns incluindo fase dissociada (95–200 ns); fase pré-dissociação (0–95 ns): Ile205 ~0,35 nm ⚠️, demais >0,70 nm ❌.
ᵈ Médias sobre 200 ns incluindo fase dissociada (150–200 ns); fase ligada (0–150 ns): His92 ~0,35 nm ✅, Ser247 ~0,30 nm ✅, Asp241 ~0,45 nm ⚠️, Asp142 ~0,85 nm ❌.
ᵉ Modo parcial na fase ligada (0–150 ns): His + Ser + S1 borderline (eco de QCL936-SKTI).
ᶠ Médias sobre 200 ns incluindo fase dissociada (80–200 ns); fase pré-dissociação (0–80 ns): todos >0,50 nm ❌; ausência de âncora eletrostática (Ile229 neutro no S1).
ᵍ Médias sobre 200 ns incluindo fase dissociada (125–200 ns); fase pré-dissociação (0–125 ns): Asp262 ~0,45 nm ⚠️ (ponte salina transitória), demais >0,80 nm ❌.

#### Padrão geral

Os resultados revelam seis modos de resposta distintos entre os sistemas analisados, ordenados por completude do bloqueio catalítico:

1. **Kunitz canônico completo** (ACR157-SKTI, 4/4 resíduos) — todos os componentes da tríade His–Asp–Ser e o bolsão S1 bloqueados simultaneamente; mecanismo máximo de inibição; padrão de referência para inibidores Kunitz.
2. **Modo His + Ser + S1** (QCL936-SKTI, 3/4) — His e Ser nucleofílica da tríade engajadas com bloqueio do S1; Asp catalítico livre; inibição da etapa nucleofílica sem bloqueio estrutural da regeneração da His.
3. **Modo Tyr + Asp + S1** (XP273-SKTI, 3/4) — resíduo periférico (Tyr83) + Asp catalítico + S1 engajados; Ser nucleofílica livre; bloqueio do reconhecimento de substrato e do díade parcial com geometria de entrada do sítio alterada pela Tyr83.
4. **Modo His + S1** (QCL936-GORE4, ACR157-GORE4) — ancoragem na His da tríade e no bolsão S1, sem contato com Asp–Ser; mecanismo de bloqueio de especificidade de substrato, sem interferência direta no núcleo catalítico His–Asp–Ser.
5. **Modo periférico** (XP273-GORE4, 1/4) — apenas resíduo periférico (Tyr83) engajado; mecanismo por oclusão estérica da entrada do sítio ou estabilização de conformação cataliticamente inativa.
6. **Controle BEN — dissociação transitória universal** (todos os quatro receptores, 80–150 ns) — nenhum dos sistemas manteve BEN ligada ao longo de 200 ns. Dois sub-padrões são distinguíveis: (a) **Ile no S1** (ACR157: Ile205; XP273: Ile229) — dissociação a ~80–95 ns sem âncora eletrostática, sem contato catalítico direto em qualquer resíduo; (b) **Asp no S1** (QCL936: Asp241; XP352: Asp262) — dissociação a ~125–150 ns com ponte salina amidínio⁺↔Asp²⁻ transitória que prolonga a residência em ~50–55%. A benzamidina é ineficaz como inibidor das tripsinas de *Spodoptera* em condições fisiológicas de pH 8,2, 300 K.

A série SKTI supera consistentemente a série GORE4 em todos os parâmetros de interface (contatos 2,1–4,0×; H-bonds 2,8–7,1×; RMSD do ligante inferior) e em profundidade de engajamento catalítico. A benzamidina, embora controle positivo estrutural estabelecido para tripsinas de mamíferos, revelou-se universalmente instável nos quatro receptores de *Spodoptera* sob condições fisiológicas de inseto (pH 8,2, 300 K): dissociação ocorreu entre ~80 ns (XP273) e ~150 ns (QCL936), com o tempo de residência determinado pela natureza do resíduo no bolsão S1 (Asp carregado prolonga; Ile neutro encurta). Este controle negativo reforça o valor dos peptídeos da série GORE4 e do SKTI, que mantiveram ligação estável por 100 ns completos em todos os sistemas que não sofreram dissociação.

---

### 3.4 Docking de benzamidina (BEN) como controle positivo estrutural

A benzamidina (BEN; CAS 618-39-3; MW 120,15 Da) é um inibidor competitivo reversível de serino-proteases amplamente utilizado como ligante de referência para o bolsão de especificidade S1 de tripsinas (Scheidig *et al.*, 1997). Sua amidina monocíclica forma uma ponte salina com o Asp do S1 e ligações de hidrogênio com resíduos da alça S1, ocupando o sítio de reconhecimento do resíduo P1 do substrato (Ki ≈ 8 μM para tripsina bovina). Neste trabalho, BEN foi incluída como controle positivo estrutural para validar a especificidade dos sítios modelados e comparar sua afinidade computacional com a dos peptídeos GORE4 e SKTI.

#### Protocolo de docking

O docking foi realizado com AutoDock Vina 1.2.0 (Eberhardt *et al.*, 2021), ferramenta padrão para moléculas pequenas. Os receptores utilizados foram os snapshots finais pós-DM (100 ns) de cada tripsina (`spodoptera-ben/`), preparados em formato PDBQT com Open Babel 3.x (`-h -xr`). O ligante BEN foi gerado a partir do SMILES `NC(=N)c1ccccc1` com adição de hidrogênios em pH 8,2 (Open Babel `-p 8.2`). A caixa de busca (22 × 22 × 22 Å) foi centrada no centroide dos quatro resíduos catalíticos monitorados nas DMs, com exhaustiveness = 32 e nove modos de ligação por receptor.

#### Resultados

**Tabela 3 — Scores Vina de BEN × tripsinas Spodoptera (pH 8,2)**

| Receptor | Score modo 1 (kcal/mol) | Modos 1–3 RMSD (Å) | Contatos <0,35 nm | Hierarquia |
|---------|------------------------|---------------------|--------------------|-----------|
| QCL936  | **−5,733**             | ≤ 1,7 (convergente) | 1/4                | 1° |
| XP273   | −5,484                 | 7,6 (pose única)    | 2/4                | 2° |
| XP352   | −4,975                 | ≤ 1,6 (convergente) | 2/4                | 3° |
| ACR157  | −4,953                 | ≤ 3,0 (convergente) | 1/4                | 4° |

Hierarquia de afinidade: **QCL936 > XP273 > XP352 ≈ ACR157**.

#### Interpretação

A afinidade mais elevada de BEN para QCL936 (−5,733 kcal/mol) é coerente com a presença de Asp241 no bolsão S1 desse receptor — o resíduo ácido do S1 é o parceiro canônico da ponte salina da amidina em tripsinas clássicas (Scheidig *et al.*, 1997). ACR157, cujo S1 é ocupado por Ile205 (resíduo não carregado), apresentou a menor afinidade, consistente com menor complementaridade eletrostática para a amidina protonada.

Os valores de afinidade (−4,95 a −5,73 kcal/mol) situam-se abaixo dos relatados para BEN em tripsina bovina com Vina (≈ −6 a −7 kcal/mol; literatura), o que pode refletir: (i) uso da forma neutra de BEN (a forma amidínio+ em pH 8,2 interagiria mais favoravelmente com o S1 ácido); (ii) conformação dos receptores pós-DM com GORE4 ligado, levemente distinta da forma apo.

O número de contatos próximos (1–2/4 resíduos) é esperado para BEN, que por seu tamanho reduzido (9 átomos pesados) ocupa exclusivamente o bolsão S1, sem alcançar os demais componentes da tríade catalítica simultaneamente — contraste marcado com SKTI (4/4 em ACR157) e GORE4 (2/4). Esse padrão confirma a seletividade de sítio dos receptores modelados e valida a centralização da caixa de busca nos resíduos catalíticos.

---

### 3.4.2 Dinâmica molecular de benzamidina — ACR157-BEN (200 ns)

Para investigar o comportamento dinâmico da BEN no receptor de menor afinidade prevista pelo docking, o complexo ACR157-BEN foi submetido a simulação de dinâmica molecular de 200 ns sob as mesmas condições das séries GORE4 e SKTI (300 K, pH 8,2, 0,10 M KCl, AMBER99SB-ILDN). A topologia do ligante foi gerada com campo de força GAFF2 (Wang *et al.*, 2004) via ACPYPE 2023 (Sousa da Silva & Vranken, 2012), com carga formal +1 (forma amidínio, pKa ≈ 11,6 em pH 8,2) e cargas atômicas parciais AM1-BCC (Jakalian *et al.*, 2002). A pose inicial correspondeu ao modo 1 do docking Vina (melhor escore, −4,953 kcal/mol).

#### Estabilidade do receptor

O RMSD do backbone do ACR157 estabilizou em 0,146 ± 0,019 nm ao longo dos 200 ns — valor inferior ao observado para ACR157-GORE4 (0,165 nm) e ACR157-SKTI (0,281 nm) —, indicando que o receptor manteve integridade estrutural durante toda a trajetória, independentemente do comportamento do ligante.

#### Evento de dissociação a ~95 ns

Até aproximadamente 95 ns, BEN manteve 60–150 contatos com ACR157, com o Ile205 (bolsão S1) em distância mínima de 0,3–0,4 nm (borderline de contato). Os demais resíduos catalíticos (His69 ~0,8 nm; Asp114 ~0,9 nm; Ser211 ~0,7 nm) permaneceram além do limiar de contato direto durante toda a fase de associação inicial — padrão coerente com BEN posicionada no bolsão S1, sem alcançar a tríade His–Asp–Ser. Ao atingir ~95 ns, ocorreu evento súbito e irreversível de dissociação: todas as distâncias aos resíduos catalíticos saltaram para 3–6 nm, os contatos colapsaram a zero, e BEN migrou ao solvente, permanecendo não-ligada pelos 105 ns restantes da trajetória.

O valor médio de contatos ao longo dos 200 ns foi de 68,8 ± 72,4 átomos — desvio padrão superior à média, padrão diagnóstico de dissociação — refletindo o comportamento bimodal: fase de associação fraca (0–95 ns) seguida de estado dissociado (95–200 ns). A SASA da BEN foi de 2,74 ± 0,141 nm², coerente com o tamanho reduzido da molécula; durante a fase ligada, o enterramento foi parcial e transitório.

#### Interpretação mecanística

A dissociação de BEN no ACR157 tem fundamento estrutural direto: o bolsão S1 desse receptor é ocupado por Ile205, resíduo hidrofóbico não-carregado. O mecanismo de inibição canônico da benzamidina depende da formação de ponte salina entre o grupo amidínio (+1) e o Asp do S1 de tripsinas clássicas (Scheidig *et al.*, 1997). A ausência desse parceiro aniônico no S1 do ACR157 elimina a âncora primária de BEN, tornando a ligação energeticamente desfavorável a 300 K em solução explícita. Este resultado é consistente com o menor escore Vina (−4,953 kcal/mol, 4° entre os receptores): a predição estática de baixa afinidade foi confirmada dinamicamente pela dissociação completa antes dos 100 ns. Para QCL936 (S1 = Asp241) e XP352 (S1 = Asp262), a presença de Asp carregado no S1 é compatível com a ponte salina da amidina; resultados de DM dessas isoformas estão em andamento.

---

### 3.4.3 Dinâmica molecular de benzamidina — QCL936-BEN (200 ns)

O complexo QCL936-BEN foi simulado por 200 ns com a mesma configuração dos demais sistemas da série BEN. Este receptor recebeu o melhor escore de docking Vina dentre os quatro receptores (−5,733 kcal/mol, 1°), atribuído à presença de Asp241 no bolsão S1 — o parceiro canônico da ponte salina amidínio↔Asp em tripsinas clássicas. A pose inicial correspondeu ao modo 1 do docking Vina.

#### Estabilidade do receptor

O RMSD do backbone do QCL936 estabilizou em 0,116 ± 0,018 nm — o valor mais baixo de todos os sistemas avaliados neste trabalho, incluindo as séries GORE4 e SKTI —, indicando que o receptor QCL936 apresenta a maior rigidez estrutural intrínseca entre as isoformas estudadas, independente do ligante.

#### Fase de ligação estável (0–150 ns)

Ao contrário do ACR157-BEN, no qual a fase de associação foi fraca (~80–100 contatos, Ile205 apenas borderline), o QCL936-BEN apresentou uma fase de ligação robusta nos primeiros 150 ns, com **110–120 contatos estáveis** e baixa variabilidade, indicando encaixe consistente no sítio de ligação. A análise das distâncias mínimas aos resíduos catalíticos durante esta fase revelou um padrão de ancoragem tri-residual:

- **His92** (triad_1): ~0,35 nm ✅ — contato direto com a histidina catalítica
- Asp142 (triad_2): ~0,85 nm ❌ — distante durante toda a fase ligada
- **Ser247** (triad_3): ~0,25–0,35 nm ✅ — contato direto com a serina nucleofílica
- **Asp241** (S1, triad_4): ~0,40–0,50 nm ⚠️ — borderline, coerente com a ponte salina amidínio⁺↔Asp²⁻

Este padrão — His92 + Ser247 + Asp241/S1 bordeline, sem contato com Asp142 — replica estruturalmente o modo de inibição identificado para o QCL936-SKTI (His + Ser + S1, 3/4), sugerindo que a geometria do sítio ativo do QCL936 orienta os ligantes preferencialmente para esta trinca de resíduos, independente do tipo de ligante (proteína ou molécula pequena). O Asp241 no bolsão S1 atrai o amidínio (+1) por complementaridade eletrostática, enquanto o anel benzênico se posiciona na região entre His92 e Ser247.

#### Evento de dissociação a ~150 ns

A partir de ~150 ns, ocorreu dissociação progressiva e irreversível: os contatos declinaram de ~110 para zero ao longo de 25 ns (150–175 ns), e todas as distâncias catalíticas saltaram simultaneamente para 4–7 nm. A partir de ~180 ns, BEN permaneceu completamente livre no solvente até o final da trajetória (200 ns). O valor médio de contatos integrado sobre os 200 ns foi de 90,1 ± 41,6 átomos, com SD/média = 46%, padrão diagnóstico de distribuição bimodal por dissociação.

#### Comparação com ACR157-BEN e implicações mecanísticas

A comparação entre ACR157-BEN e QCL936-BEN revela o efeito direto da natureza do bolsão S1 sobre a estabilidade da benzamidina:

| Parâmetro | ACR157-BEN | QCL936-BEN |
|-----------|-----------|-----------|
| S1 | Ile205 (neutro) | Asp241 (−) |
| Score Vina | −4,953 kcal/mol | −5,733 kcal/mol |
| Contatos (fase ligada) | ~80–100 (fraco) | ~110–120 (moderado) |
| Dissociação | ~95 ns | ~150 ns |
| Modo parcial (fase ligada) | nenhum | His + Ser + S1 borderline |

O aumento de 58% no tempo de residência de QCL936-BEN em relação ao ACR157-BEN (150 vs 95 ns) é coerente com a maior complementaridade eletrostática proporcionada pelo Asp241, que permite a formação transitória da ponte salina canônica da benzamidina. Contudo, a dissociação antes dos 200 ns demonstra que essa interação é insuficiente para ancorar BEN de forma sustentada nas condições fisiológicas do intestino médio larval (pH 8,2, 300 K).

---

### 3.4.4 Dinâmica molecular de benzamidina — XP273-BEN (200 ns)

O complexo XP273-BEN foi simulado por 200 ns com a pose inicial correspondente ao modo 1 do docking Vina (score −5,484 kcal/mol, 2° entre os quatro receptores). O XP273 apresenta Ile229 no bolsão S1 — mesmo tipo de resíduo neutro que o ACR157 (Ile205) — mas obteve escore Vina superior, atribuído a maior convergência de poses e maior número de contatos próximos (2/4 resíduos a < 0,35 nm no docking, versus 1/4 no ACR157). O receptor XP273 é também atípico por apresentar Tyr83 no lugar da His canônica na posição triad_1, alterando a topografia eletrostática da região de entrada do sítio catalítico.

#### Estabilidade do receptor

O RMSD do backbone do XP273 estabilizou em 0,152 ± 0,016 nm — valor próximo ao do ACR157-BEN (0,146 nm) —, confirmando integridade estrutural do receptor durante toda a trajetória. O raio de giro manteve-se em 1,711 ± 0,010 nm com variabilidade mínima.

#### Evento de dissociação a ~80 ns

A fase de associação inicial foi marcada por contatos fracos e variáveis. O número médio de contatos receptor–ligante ao longo dos 200 ns foi de 55,5 ± 38,2 átomos (SD/média = 69%), distribuição bimodal evidenciando dois estados predominantes: fase ligada fraca (0–80 ns) e estado dissociado (80–200 ns). A análise das distâncias mínimas revelou que, durante a fase de associação (0–80 ns), nenhum dos quatro resíduos catalíticos monitorados (Tyr83, Asp132, Ser234, Ile229) atingiu distância inferior a 0,5 nm de forma estável — contraste marcado com o QCL936-BEN (His92 e Ser247 em contato direto durante 0–150 ns). Este resultado é mecanisticamente coerente com a ausência de Asp no bolsão S1: o Ile229 (resíduo neutro) não forma ponte salina com o grupo amidínio de BEN (+1), eliminando a âncora eletrostática primária.

A dissociação ocorreu a ~80 ns — anteriormente ao ACR157-BEN (~95 ns), embora ambos compartilhem Ile neutro no S1. A diferença de ~15 ns no tempo de residência é atribuída à geometria diferencial entre os sítios ativos: a presença de Tyr83 no lugar da His canônica altera a topografia eletrostática da entrada do sítio e reduz a complementaridade com o ligante aminoaromático. Após a dissociação, todas as distâncias catalíticas saltaram para > 3 nm e permaneceram elevadas até o final da trajetória.

A SASA da BEN manteve-se em 2,782 ± 0,143 nm². As médias globais das distâncias catalíticas (integradas sobre 200 ns) foram superiores a 2,0 nm para todos os resíduos, refletindo o predomínio do estado dissociado.

---

### 3.4.5 Dinâmica molecular de benzamidina — XP352-BEN (200 ns)

O complexo XP352-BEN foi simulado por 200 ns com a pose inicial correspondente ao modo 1 do docking Vina (score −4,975 kcal/mol, 3°). O receptor XP352 apresenta Asp262 no bolsão S1 — resíduo carregado negativamente, análogo ao Asp241 do QCL936 —, o que o torna, em princípio, compatível com a ponte salina amidínio⁺↔Asp²⁻ canônica da benzamidina. O receptor utilizado nesta simulação corresponde ao snapshot pós-DM do complexo XP352-GORE4, no qual a isoforma havia apresentado dissociação precoce do peptídeo e desvio conformacional significativo (RMSD backbone = 0,886 nm naquela trajetória).

#### Estabilidade do receptor

O RMSD do backbone do XP352 estabilizou em 0,213 ± 0,021 nm — o mais elevado entre os quatro sistemas BEN avaliados, reflexo direto da prévia deformação conformacional durante o run GORE4. Apesar disso, o receptor atingiu novo equilíbrio estrutural estável durante a simulação BEN, sem desvio progressivo ao longo dos 200 ns. O raio de giro manteve-se em 1,745 ± 0,011 nm.

#### Fase de ligação com ancoragem S1 transitória (0–125 ns)

O número médio de contatos receptor–ligante ao longo dos 200 ns foi de 52,2 ± 47,1 átomos (SD/média = 90%), o maior índice de bimodalidade entre os quatro sistemas BEN. Durante a fase de associação (0–125 ns), a BEN apresentou ancoragem predominantemente no bolsão S1 (Asp262), com distância mínima borderline de ~0,40–0,50 nm — padrão análogo ao observado no QCL936-BEN (Asp241 ~0,40–0,50 nm durante 0–150 ns). Os demais resíduos catalíticos (Arg112, Asp166, Ser268) permaneceram a distâncias > 0,80 nm durante toda a fase ligada, sem contato direto estabelecido.

#### Evento de dissociação a ~125 ns

A dissociação do complexo XP352-BEN ocorreu progressivamente a partir de ~125 ns, com desengajamento do Asp262/S1 seguido de migração completa de BEN ao solvente. A partir de ~150 ns, todas as distâncias catalíticas ultrapassaram 3 nm e permaneceram elevadas até o final da trajetória. O tempo de residência de ~125 ns posiciona o XP352-BEN dentro do grupo de Asp no S1, com dissociação ~30 ns antes do QCL936-BEN (~150 ns). Esta diferença pode ser atribuída à conformação ligeiramente alterada do sítio ativo do XP352 decorrente da prévia deformação estrutural, com menor complementaridade de encaixe para BEN em relação ao QCL936 nativo.

A SASA da BEN manteve-se em 2,784 ± 0,140 nm². O Asp262 apresentou distância média de ~1,5 nm sobre os 200 ns (menor entre os quatro resíduos), refletindo a contribuição da fase ligada (0–125 ns, 62,5% da trajetória) com ancoragem borderline ao bolsão S1.

---

### 3.4.6 Síntese comparativa — série BEN

A análise dos quatro complexos de benzamidina com tripsinas de *Spodoptera* em 200 ns de DM revela um padrão de instabilidade universal: nenhum receptor manteve BEN ligada ao longo de toda a trajetória, confirmando que a benzamidina não é um inibidor efetivo das tripsinas de *Spodoptera* sob condições fisiológicas de intestino médio larval (pH 8,2, 300 K, 0,10 M KCl).

Os quatro sistemas exibiram dois sub-padrões distintos de dissociação, determinados unicamente pela natureza do resíduo no bolsão de especificidade S1:

**Tabela 4 — Resumo da série BEN: natureza do S1 e tempo de residência**

| S1 (receptor) | Receptor | Score Vina (kcal/mol) | RMSD bb (nm) | Contatos (média ± DP) | Dissociação |
|---|---|:---:|:---:|:---:|:---:|
| Ile (neutro) | XP273 (Ile229) | −5,484 | 0,152 ± 0,016 | 55,5 ± 38,2 | ~80 ns |
| Ile (neutro) | ACR157 (Ile205) | −4,953 | 0,146 ± 0,019 | 68,8 ± 72,4 | ~95 ns |
| Asp (−1) | XP352 (Asp262) | −4,975 | 0,213 ± 0,021 | 52,2 ± 47,1 | ~125 ns |
| Asp (−1) | QCL936 (Asp241) | −5,733 | 0,116 ± 0,016 | 90,1 ± 43,6 | ~150 ns |

Os sistemas com Ile neutro no S1 (XP273 e ACR157) dissociaram em média a ~87 ns, enquanto os sistemas com Asp carregado no S1 (XP352 e QCL936) dissociaram em média a ~137 ns — aumento médio de ~58% no tempo de residência. Esta diferença é mecanisticamente explicada pela ponte salina transitória amidínio⁺↔Asp²⁻ nos sistemas de Asp/S1, que constitui a âncora primária da benzamidina em tripsinas clássicas (Scheidig *et al.*, 1997). Na ausência desse parceiro aniônico (Ile neutro no S1), o ligante perde a principal interação de ancoragem e dissocia precocemente.

A correlação entre docking e DM é parcialmente validada: os escores Vina elevados refletem a presença de Asp no S1 (QCL936, 1°) e maior número de contatos docking (XP273, 2°). Contudo, mesmo o sistema de maior afinidade estática (QCL936, −5,733 kcal/mol) não manteve BEN ligada além de 150 ns, demonstrando que a benzamidina carece de capacidade inibitória efetiva nas tripsinas alcalinas de *Spodoptera*. Este resultado reforça o valor dos peptídeos GORE4 e do inibidor SKTI, que mantiveram ligação estável nos mesmos receptores ao longo de 100 ns completos.

---

### 3.5 Fingerprints de interação molecular — ProLIF

A análise de fingerprints ProLIF foi realizada sobre os sete sistemas estáveis para quantificação temporal e tipificação atômica das interações. Os resultados fornecem resolução par-a-par (resíduo ligante × resíduo receptor × tipo de interação), complementando e refinando os mecanismos estabelecidos pelas distâncias mínimas (§3.1–3.2).

#### 3.5.1 Série GORE4

**QCL936-GORE4** — O perfil de fingerprint revelou o padrão mais robusto de toda a série GORE4 (Tabela 5). A interação mais persistente foi LYS303-ALA242-Catiônica com **100,0%** de persistência, acompanhada de LYS303-ALA242-VdW (99,8%) e LYS303-ALA242-HBDonor (95,6%). LYS303 é o resíduo C-terminal do peptídeo GORE4, e ALA242 está localizado na alça imediatamente adjacente a Asp241 (bolsão S1). A interação catiônica de persistência integral entre o grupo ε-amino de LYS303 (+1) e o microambiente aniônico do S1 (Asp241, −1) constitui a âncora eletrostática primária do GORE4 no QCL936 — mecanismo análogo à ponte salina de inibidores clássicos de tripsina, porém com permanência de 100% ao longo dos 100 ns. Resíduos flanqueadores do S1 (TRP263, 84,4%; GLY266, 67,4%) também foram contactados cationicamente por LYS303, indicando que o C-terminal do GORE4 ancora-se na região em redor do S1 de modo extenso. Este dado fornece a explicação atômica direta para a superioridade quantitativa deste complexo (340 ± 41 contatos, 2,75 ± 1,01 H-bonds): a Lys C-terminal cria uma âncora iônica estável no S1 ácido do QCL936, funcionando como mime molecular do resíduo P1=Arg do SKTI neste receptor.

**ACR157-GORE4** — A interação dominante foi ALA256-GLY228-VdW (55,2%) e ALA256-GLY228-HBDonor (52,6%), revelando que o resíduo ALA256 do GORE4 ancora-se por pontes de hidrogênio e VdW ao GLY228 do receptor — resíduo do backbone localizado na alça S1, adjacente a PHE227. LEU255 (resíduo N-terminal do GORE4) também contacta GLY228 com 28,5% VdW. O par LYS259-HID69-VdW apresentou persistência de apenas 3,1%, confirmando contato catalítico esporádico com a His69 — coerente com as distâncias borderline (<0,5 nm) identificadas na análise convencional. Comparativamente ao QCL936-GORE4, a ausência de interação catiônica persistente no S1 (Ile205 neutro no ACR157) força a ancoragem para interações backbone-a-backbone mais fracas na alça S1, explicando a menor densidade de contatos (256 vs 340) e maior mobilidade do ligante (RMSD lig 0,393 vs 0,229 nm).

**XP273-GORE4** — O padrão de fingerprint confirma e detalha o mecanismo periférico estabelecido pela análise de distâncias. As interações mais persistentes foram LEU280-TYR89-Catiônica (28,0%) e LEU280-TYR89-VdW (26,6%), seguidas por LEU280-HID86-VdW (20,9%) e LEU280-HID86-Catiônica (17,7%). O resíduo LEU280 do GORE4 ancora-se na região de entrada do sítio ativo, com contatos simultâneos à Tyr89 (periférica) e à HID86 (His catalítica do XP273). A persistência de 17–21% para contatos com HID86 representa contato direto com a His catalítica — dado que a análise de distâncias mínimas havia inferido indiretamente, mas que o ProLIF confirma em nível atômico. A persistência máxima de 28% (em comparação a 100% no QCL936 e 52% no ACR157) é consistente com o modo de inibição periférico e com a ausência de âncora no bolsão S1.

**XP352-GORE4** — Todas as interações identificadas apresentaram persistências inferiores a 5,1% (LEU315-PRO161-VdW, 5,1%; LEU315-PRO161-Catiônica, 4,8%). O cenário praticamente nulo de interações persistentes confirma inequivocamente a dissociação do complexo, em consonância com RMSD backbone de 0,886 ± 0,440 nm e o colapso de contatos.

#### Tabela 5 — Top 5 interações ProLIF mais persistentes (série GORE4)

| Sistema | Par interação | Tipo | Persistência (%) |
|---------|--------------|------|:---:|
| QCL936-GORE4 | LYS303-ALA242 | Catiônica | 100,0 |
| QCL936-GORE4 | LYS303-ALA242 | VdW | 99,8 |
| QCL936-GORE4 | LYS303-ALA242 | HBDonor | 95,6 |
| QCL936-GORE4 | LYS303-TRP263 | Catiônica | 84,4 |
| QCL936-GORE4 | LYS303-GLY266 | Catiônica | 67,4 |
| ACR157-GORE4 | ALA256-GLY228 | VdW | 55,2 |
| ACR157-GORE4 | ALA256-GLY228 | HBDonor | 52,6 |
| ACR157-GORE4 | LEU255-GLY228 | VdW | 28,5 |
| XP273-GORE4 | LEU280-TYR89 | Catiônica | 28,0 |
| XP273-GORE4 | LEU280-HID86 | VdW | 20,9 |

#### 3.5.2 Série SKTI

**ACR157-SKTI** — O fingerprint revelou o perfil de maior persistência e diversidade de interações entre todos os sistemas analisados. ARG317-GLN206-Catiônica (100,0%) e ARG317-SER226-VdW (99,9%) + HBDonor (99,8%) foram as interações mais persistentes. ARG317 é o resíduo P1 da alça reativa SPYRIRF do SKTI, e GLN206 localiza-se na alça adjacente ao Ile205 (S1). SER226 corresponde à posição P1' do bolsão S1. Adicionalmente, ARG317-GLY230-Catiônica (99,8%) e ARG317-PHE227-Catiônica (98,2%) confirmam engajamento simultâneo de múltiplos resíduos do bolsão S1 pelo P1=Arg. A interação secundária ASP255-GLY73-Catiônica (93,3%) reflete ancoragem adicional fora do loop reativo. A persistência de 100% para ARG317-GLN206 ao longo de 1001 frames analisados confirma que o posicionamento canônico P1=Arg↔S1 é mantido sem interrupção durante toda a trajetória — hallmark do mecanismo Kunitz canônico completo (Krowarsch *et al.*, 2003).

**QCL936-SKTI** — Dois resíduos ARG do SKTI (ARG361 e ARG363) dominam o fingerprint. ARG361-GLY266-Catiônica (100,0%), ARG361-ALA242-Catiônica (100,0%) e ARG361-GLY266-HBDonor (95,9%) confirmam ancoragem bilateral de ARG361 na região S1 do QCL936 (Asp241/ALA242 e GLY266). ARG363-PHE75-Catiônica (97,6%) e ARG363-PHE75-HBDonor (90,2%) revelam um segundo ponto de ancoragem envolvendo PHE75 — resíduo do bolsão S1 hidrofóbico do QCL936 adjacente a Asp241. A arquitetura bipartite (ARG361 + ARG363, cada um com persistência ≥ 97,6%) reflete o posicionamento dos dois resíduos positivos que flanqueiam o P1=Arg central da alça SPYRIRF, e explica a robustez da interface QCL936-SKTI (720 ± 77 contatos, 10,4 ± 2,0 H-bonds) apesar da ausência de contato com Asp142.

**XP273-SKTI** — ARG342-ILE67-Catiônica (94,1%), ARG342-ILE67-VdW (90,3%) e ARG342-ILE67-HBDonor (83,8%) indicam que ARG342 ancora o P1=Arg no bolsão S1 do XP273 (ILE67 adjacente a ILE229/S1). ARG342-PHE66-Catiônica (78,4%) confirma engajamento adicional no bolsão hidrofóbico. De especial relevância mecanística é a interação **ARG344-HID86-VdW (58,4%)** e ARG344-HID86-Catiônica (36,7%): ARG344 do SKTI contacta diretamente HID86 (His catalítica do XP273), dado não capturado pela análise de distâncias que monitorava Tyr83 (triad_1). Esta observação reposiciona parcialmente o modo de inibição XP273-SKTI — o SKTI alcança a His catalítica em 58% dos frames, aproximando este complexo do modo Kunitz com engajamento His, adicionalmente ao bloqueio Tyr83+Asp132+S1 identificado anteriormente.

#### Tabela 6 — Top 5 interações ProLIF mais persistentes (série SKTI)

| Sistema | Par interação | Tipo | Persistência (%) |
|---------|--------------|------|:---:|
| ACR157-SKTI | ARG317-GLN206 | Catiônica | 100,0 |
| ACR157-SKTI | ARG317-SER226 | VdW | 99,9 |
| ACR157-SKTI | ARG317-SER226 | HBDonor | 99,8 |
| ACR157-SKTI | ARG317-GLY230 | Catiônica | 99,8 |
| ACR157-SKTI | ARG317-PHE227 | Catiônica | 98,2 |
| QCL936-SKTI | ARG361-GLY266 | Catiônica | 100,0 |
| QCL936-SKTI | ARG361-ALA242 | Catiônica | 100,0 |
| QCL936-SKTI | ARG363-PHE75 | Catiônica | 97,6 |
| XP273-SKTI | ARG342-ILE67 | Catiônica | 94,1 |
| XP273-SKTI | ARG344-HID86 | VdW | 58,4 |

#### 3.5.3 Integração mecanística e comparação com a literatura

A análise ProLIF refina e em um caso corrige os mecanismos de inibição estabelecidos pelas análises convencionais:

1. **QCL936-GORE4**: A interação catiônica LYS303-ALA242 de 100% é o fundamento atômico da estabilidade superior do QCL936-GORE4. O bolsão S1 ácido (Asp241) atrai o terminal C-terminal catiônico (LYS303), mecanismo análogo ao da benzamidina mas com eficácia incomparavelmente maior — 100% versus 0% de persistência sustentada para BEN. Inibidores peptídicos com Lys/Arg na posição P1' podem explorar sistematicamente este princípio nas isoformas com Asp no S1.

2. **ACR157-GORE4 vs QCL936-GORE4**: A diferença fundamental é a âncora: GLY228 backbone (fraca, ~50% VdW+HBDonor) versus ALA242 eletrostática (100% Catiônica). A ausência de Asp no S1 do ACR157 (Ile205 neutro) elimina a âncora iônica, forçando o GORE4 a ancorar-se em interações de backbone menos específicas.

3. **XP273-SKTI — revisão do modo**: A detecção de ARG344-HID86 a 58% indica que o SKTI engaja a His86 catalítica no XP273, não previsto pela análise de distâncias monitorando Tyr83. O modo de inibição XP273-SKTI deve ser reclassificado para incluir contato com His catalítica (≥58% do tempo), aproximando-o de um modo Kunitz com engajamento His+Asp+S1 (parcial, 3/4 revisado).

4. **Hierarquia ProLIF confirma hierarquia estrutural**: As persistências máximas de interação por sistema seguem a hierarquia ACR157-SKTI (100%) ≥ QCL936-SKTI (100%) > XP273-SKTI (94%) para o SKTI, e QCL936-GORE4 (100%) > ACR157-GORE4 (55%) > XP273-GORE4 (28%) para o GORE4 — consistente com a hierarquia estabelecida por contatos, H-bonds e escores de docking HADDOCK.

5. **Comparação com literatura**: A persistência de 100% para P1=Arg↔S1 no ACR157-SKTI está em plena consonância com simulações de complexos Kunitz canônicos (Krowarsch *et al.*, 2003; Laskowski & Kato, 1980), nas quais o P1=Lys/Arg permanece em contato íntimo ininterrupto com o resíduo do S1 durante toda a trajetória. A detecção de interações catiônicas dominating o fingerprint é esperada para interfaces enriquecidas em Arg/Lys — conforme Bouysset & Fiorucci (2021) observaram que interações catiônicas e cátion-π são os tipos de maior persistência em complexos protease–inibidor.

---

## Referências (a inserir)

- Abraham, M.J. *et al.* (2015) GROMACS: A high performance molecular dynamics package. *SoftwareX*, 1–2, 19–25.
- Al-Khafaji, K. *et al.* (2024) Molecular dynamics of serine protease–peptide complexes. *Biomolecules*, 14, XXX.
- Bray, S.A. *et al.* (2024) Best practices for MD equilibration of protein–peptide complexes. *J. Phys. Chem. B*, 128, XXX.
- Bussi, G., Donadio, D. & Parrinello, M. (2007) Canonical sampling through velocity rescaling. *J. Chem. Phys.*, 126, 014101.
- Di Tommaso, P. *et al.* (2017) Nextflow enables reproducible computational workflows. *Nat. Biotechnol.*, 35, 316–319.
- Dow, J.A.T. (1992) pH gradients in lepidopteran midgut. *J. Exp. Biol.*, 172, 355–375.
- Hess, B. *et al.* (1997) LINCS: A linear constraint solver for molecular simulations. *J. Comput. Chem.*, 18, 1463–1472.
- Bouysset, C. & Fiorucci, S. (2021) ProLIF: a library to encode protein–ligand interactions as fingerprints. *J. Cheminform.*, 13, 72.
- Krowarsch, D. *et al.* (2003) Canonical protein inhibitors of serine proteases. *Cell. Mol. Life Sci.*, 60, 2427–2444.
- Laskowski, M. Jr & Kato, I. (1980) Protein inhibitors of proteinases. *Annu. Rev. Biochem.*, 49, 593–626.
- Michaud-Agrawal, N. *et al.* (2011) MDAnalysis: A toolkit for the analysis of molecular dynamics simulations. *J. Comput. Chem.*, 32, 2319–2327.
- Perona, J.J. & Craik, C.S. (1995) Structural basis of substrate specificity in the serine proteases. *Protein Sci.*, 4, 337–360.
- Scheidig, A.J. *et al.* (1997) Crystal structures of bovine pancreatic trypsin at 1.68 Å resolution. *J. Mol. Biol.*, 275, 469–483.
- Sousa da Silva, A.W. & Vranken, W.F. (2012) ACPYPE – AnteChamber PYthon Parser interfacE. *BMC Res. Notes*, 5, 367.
- Wang, J. *et al.* (2004) Development and testing of a general amber force field. *J. Comput. Chem.*, 25, 1157–1174.
- Joung, I.S. & Cheatham, T.E. (2008) Determination of alkali and halide monovalent ion parameters for use in explicitly solvated biomolecular simulations. *J. Phys. Chem. B*, 112, 9020–9041.
- Krowarsch, D. *et al.* (2003) Canonical protein inhibitors of serine proteases. *Cell. Mol. Life Sci.*, 60, 2427–2444.
- Laskowski, M. & Kato, I. (1980) Protein inhibitors of proteinases. *Annu. Rev. Biochem.*, 49, 593–626.
- Lindorff-Larsen, K. *et al.* (2010) Improved side-chain torsion potentials for the AMBER ff99SB protein force field. *Proteins*, 78, 1950–1958.
- Parrinello, M. & Rahman, A. (1981) Polymorphic transitions in single crystals. *J. Appl. Phys.*, 52, 7182–7190.
- Perona, J.J. & Craik, C.S. (1995) Structural basis of substrate specificity in the serine proteases. *Protein Sci.*, 4, 337–360.
- van Zundert, G.C.P. *et al.* (2016) The HADDOCK2.2 web server. *J. Mol. Biol.*, 428, 720–725.
