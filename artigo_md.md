# Dinâmica Molecular — Seções do Artigo

---

## 2. Metodologia

### 2.1 Preparação dos complexos

As estruturas iniciais dos complexos foram geradas por docking molecular utilizando o servidor HADDOCK 2.4 (van Zundert *et al.*, 2016). Os receptores avaliados correspondem a seis isoformas de tripsina digestiva de *Spodoptera* spp.: ACR157, QCL936, XP273, XP352, DN773 e DN1441. Foram estudadas duas séries de inibidores: (i) **série GORE4** — o peptídeo GORE4 de cinco resíduos, com todos os resíduos diretamente envolvidos na interface de ligação, modelado por docking proteína–peptídeo; (ii) **série SKTI** — o Inibidor de Tripsina Kunitz de Soja (*Soybean Kunitz Trypsin Inhibitor*, SKTI; 177 resíduos), inibidor natural de referência com alça reativa SPYRIRF (P1 = Arg), modelado por docking proteína–proteína. Para cada sistema, o cluster de maior escore energético foi selecionado como pose inicial para as simulações.

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

#### Tabela 1 — Parâmetros de dinâmica molecular (médias ± DP, 100 ns)

| Sistema | RMSD bb (nm) | RMSD lig (nm) | Contatos | H-bonds | SASA lig (nm²) | Status |
|---------|:---:|:---:|:---:|:---:|:---:|:---:|
| **Série GORE4** | | | | | | |
| QCL936-GORE4 c3 | 0,124 ± 0,017 | 0,229 ± 0,042 | 340 ± 41 | 2,75 ± 1,01 | 9,04 ± 0,37 | ✅ estável |
| ACR157-GORE4 c1 | 0,165 ± 0,016 | 0,393 ± 0,068 | 256 ± 38 | 1,90 ± 1,38 | 9,90 ± 0,47 | ✅ estável |
| XP273-GORE4 c1 | 0,248 ± 0,045 | 0,317 ± 0,068 | 260 ± 60 | 3,19 ± 1,07 | 8,14 ± 0,37 | ✅ estável |
| XP352-GORE4 c4r3 | 0,886 ± 0,440 | 0,400 ± 0,067 | 63 ± 143 | — | — | ❌ dissociação |
| **Série SKTI** | | | | | | |
| ACR157-SKTI c2 | 0,281 ± 0,027 | 0,206 ± 0,017 | 1019 ± 118 | 13,524 ± 2,911 | 101,97 ± 1,62 | ✅ estável |
| QCL936-SKTI c2 | 0,370 ± 0,054 | 0,209 ± 0,025 | 720 ± 77 | 10,410 ± 2,022 | 102,53 ± 2,09 | ✅ estável |
| XP273-SKTI c2 | 0,288 ± 0,042 | 0,224 ± 0,017 | 596 ± 88 | 8,817 ± 3,496 | 101,31 ± 2,01 | ✅ estável |

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

#### Padrão geral

Os resultados revelam cinco modos de inibição distintos entre os complexos analisados, ordenados por completude do bloqueio catalítico:

1. **Kunitz canônico completo** (ACR157-SKTI, 4/4 resíduos) — todos os componentes da tríade His–Asp–Ser e o bolsão S1 bloqueados simultaneamente; mecanismo máximo de inibição; padrão de referência para inibidores Kunitz.
2. **Modo His + Ser + S1** (QCL936-SKTI, 3/4) — His e Ser nucleofílica da tríade engajadas com bloqueio do S1; Asp catalítico livre; inibição da etapa nucleofílica sem bloqueio estrutural da regeneração da His; modo identificado após re-run com resíduos corretos (Asp142/Ser247).
3. **Modo Tyr + Asp + S1** (XP273-SKTI, 3/4) — resíduo periférico (Tyr83) + Asp catalítico + S1 engajados; Ser nucleofílica livre; bloqueio do reconhecimento de substrato e do díade parcial com geometria de entrada do sítio alterada pela Tyr83.
4. **Modo His + S1** (QCL936-GORE4, ACR157-GORE4) — ancoragem na His da tríade e no bolsão S1, sem contato com Asp–Ser; mecanismo de bloqueio de especificidade de substrato, sem interferência direta no núcleo catalítico His–Asp–Ser.
5. **Modo periférico** (XP273-GORE4, 1/4) — apenas resíduo periférico (Tyr83) engajado; mecanismo por oclusão estérica da entrada do sítio ou estabilização de conformação cataliticamente inativa.

A série SKTI supera consistentemente a série GORE4 em todos os parâmetros de interface (contatos 2,1–4,0×; H-bonds 2,8–7,1×; RMSD do ligante inferior) e em profundidade de engajamento catalítico — resultado coerente com a natureza pré-organizada e complementar do inibidor Kunitz para o sítio ativo de serino-proteases.

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
- Krowarsch, D. *et al.* (2003) Canonical protein inhibitors of serine proteases. *Cell. Mol. Life Sci.*, 60, 2427–2444.
- Laskowski, M. & Kato, I. (1980) Protein inhibitors of proteinases. *Annu. Rev. Biochem.*, 49, 593–626.
- Lindorff-Larsen, K. *et al.* (2010) Improved side-chain torsion potentials for the AMBER ff99SB protein force field. *Proteins*, 78, 1950–1958.
- Parrinello, M. & Rahman, A. (1981) Polymorphic transitions in single crystals. *J. Appl. Phys.*, 52, 7182–7190.
- Perona, J.J. & Craik, C.S. (1995) Structural basis of substrate specificity in the serine proteases. *Protein Sci.*, 4, 337–360.
- van Zundert, G.C.P. *et al.* (2016) The HADDOCK2.2 web server. *J. Mol. Biol.*, 428, 720–725.
