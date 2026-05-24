# TVGD — Limitações Conhecidas

Este documento resume as limitações atuais da **Teoria do Vácuo Gravitacional Dispersivo (TVGD)** e do pacote público de validação.

O objetivo não é enfraquecer a teoria, mas deixar claro para avaliadores externos quais partes já estão bem estruturadas e quais ainda precisam de derivação, testes independentes ou likelihoods formais.

---

## 1. Status geral

A TVGD já possui uma arquitetura efetiva avançada, com testes numéricos em múltiplos regimes:

- galáxias e relação de aceleração radial;
- validação externa com LITTLE THINGS;
- buracos negros com saturação estrutural;
- cosmologia de fundo;
- tensão de Hubble;
- crescimento de estruturas via RSD / fsigma8;
- setor covariante efetivo com escalar, vetor, PPN e condição c_T = 1;
- critérios explícitos de falsificabilidade.

Mesmo assim, a teoria ainda deve ser apresentada como um **programa efetivo em desenvolvimento**, não como teoria final fechada.

---

## 2. Ação covariante final

A TVGD já avançou além de uma simples fenomenologia. Existe uma arquitetura covariante efetiva com:

- campo escalar estrutural;
- setor vetorial tipo æther;
- ramo com velocidade tensorial c_T = 1;
- auditoria PPN forte no ramo testado;
- checagens de estabilidade escalar, vetorial e tensorial;
- recuperação efetiva do regime galáctico;
- compatibilidade preliminar com crescimento de estruturas;
- setor de buracos negros por saturação estrutural.

A limitação atual não é ausência de covariantização. A limitação é mais específica:

> transformar a arquitetura covariante efetiva já construída em uma ação final única, com derivação variacional completa e demonstração explícita de que os regimes galáctico, cosmológico, PPN e de buracos negros emergem como limites consistentes do mesmo núcleo teórico.

---

## 3. Derivação de a0

A escala `a0` aparece como parâmetro central no regime galáctico:

```text
u = sqrt(g_bar / a0)
```

e

```text
P(u) = 1 - exp(-u)
```

A TVGD já mostra como essa escala funciona no regime efetivo, mas ainda falta uma derivação final mostrando se `a0` emerge de:

- constantes fundamentais;
- dinâmica do campo estrutural;
- estado cosmológico de fundo;
- escala de relaxação do vácuo;
- ou combinação desses fatores.

Enquanto isso não for fechado, `a0` deve ser tratado como escala fenomenológica efetiva.

---

## 4. Likelihoods completas

Vários testes públicos usam estatísticas efetivas, diagonais ou simplificadas.

Isso é aceitável para validação inicial, mas não substitui análises formais com:

- matrizes de covariância completas;
- likelihoods oficiais;
- priors transparentes;
- comparação contra ΛCDM, MOND e modelos de halos;
- tratamento de sistemáticos observacionais;
- validação independente por terceiros.

Setores onde isso é especialmente importante:

- Pantheon+SH0ES;
- CMB;
- BAO;
- RSD / fsigma8;
- lentes gravitacionais;
- SPARC com tratamento completo de M*/L e incertezas.

---

## 5. Comparação Bayesiana

A TVGD possui vários testes de χ², resíduos, AIC/BIC aproximado e robustez.

Ainda falta uma comparação Bayesiana formal contra:

- ΛCDM;
- MOND;
- AQUAL / TeVeS-like;
- modelos com halos NFW;
- modelos empíricos da relação de aceleração radial;
- extensões de energia escura dinâmica.

Essa comparação é importante porque uma teoria alternativa deve demonstrar não apenas bom ajuste, mas também boa relação entre ajuste, número de parâmetros e poder preditivo.

---

## 6. Dados externos e reprodutibilidade

O repositório prioriza código e documentação. Alguns dados externos podem não estar incluídos por tamanho, licença ou praticidade.

Isso significa que alguns scripts podem exigir:

- download manual de dados;
- colocação do arquivo na pasta esperada;
- adaptação de caminho local;
- instalação de pacotes adicionais.

Para publicação formal, cada script deve indicar claramente:

- fonte dos dados;
- referência científica;
- arquivo esperado;
- colunas usadas;
- cortes aplicados;
- saída esperada.

---

## 7. Setor galáctico

O setor galáctico é um dos mais fortes da TVGD, mas ainda possui limitações:

- dependência de tratamento correto de massa-luminosidade estelar;
- sensibilidade a outliers observacionais;
- necessidade de validação externa além do SPARC;
- necessidade de comparação formal contra halos de matéria escura;
- necessidade de separar testes exploratórios de testes pré-registrados.

A função central

```text
P(u) = 1 - exp(-u)
```

é falsificável. Amostras limpas que violem sistematicamente essa forma enfraqueceriam a teoria.

---

## 8. LITTLE THINGS

O teste com LITTLE THINGS é uma validação externa importante, mas ainda deve ser descrito como teste efetivo.

Limitações:

- reconstrução de componentes bariônicas pode depender da forma dos dados disponíveis;
- amostra menor que SPARC;
- necessidade de pipeline independente;
- necessidade de comparação formal com modelos de halo.

Mesmo assim, é um bloco favorável porque não vem da mesma amostra principal do SPARC.

---

## 9. Buracos negros

O setor de buracos negros preserva benchmarks externos importantes:

```text
r_h / r_s ≈ 1
r_ph / r_s ≈ 1.5
b_shadow / r_s ≈ 2.598076
```

A limitação atual é que a saturação estrutural ainda precisa ser ligada de forma variacional completa à ação covariante final.

Também faltam testes mais formais com:

- Kerr completo;
- modos quase-normais;
- estabilidade dinâmica;
- acreção;
- comparação sistemática com EHT e LIGO/Virgo/KAGRA.

---

## 10. Cosmologia e tensão de Hubble

O setor cosmológico da TVGD é promissor porque busca preservar o regime primordial e permitir uma taxa local efetiva maior.

Limitações atuais:

- necessidade de likelihood conjunta CMB+BAO+SN+RSD;
- controle de degenerescências com parâmetros cosmológicos;
- tratamento de sistemáticos em supernovas;
- validação com Planck/ACT/SPT;
- comparação com modelos de energia escura dinâmica;
- teste robusto de crescimento de estruturas.

Esse setor deve ser apresentado como promissor, não como fechado.

---

## 11. RSD / crescimento de estruturas

O teste RSD / fsigma8 mostra compatibilidade efetiva com ΛCDM em compilação diagonal.

Limitações:

- necessidade de covariâncias completas;
- correções fiduciais;
- likelihoods oficiais;
- comparação com diferentes compilações;
- validação contra CMB lensing e clustering.

---

## 12. Lentes gravitacionais

A TVGD prevê amplificação estrutural via:

```text
alpha(u) = 1 / P(u)
```

Limitações atuais:

- necessidade de mapas de lenteamento com máscaras bem controladas;
- tratamento de massa bariônica projetada;
- comparação com halos NFW;
- controle de seleção de amostra;
- validação em lentes fortes e fracas separadamente.

---

## 13. Setor ambiental

Alguns testes ambientais ou de vazios apresentaram resultados mistos ou falhos dependendo da máscara, seleção e overlap angular.

Esse setor não deve ser usado como evidência principal neste estágio.

Ele pode permanecer como linha exploratória, mas exige:

- catálogo de ambiente dedicado;
- overlap real com supernovas;
- máscara angular controlada;
- pré-registro dos cortes;
- validação independente.

---

## 14. Testes exploratórios

A pasta de testes exploratórios não deve ser interpretada como evidência consolidada.

Testes exploratórios servem para:

- gerar hipóteses;
- diagnosticar problemas;
- testar novas assinaturas;
- preparar previsões futuras.

A evidência principal deve vir apenas dos scripts marcados como validação principal.

---

## 15. Linguagem científica recomendada

Ao apresentar a TVGD publicamente, recomenda-se usar frases como:

- “teoria efetiva em desenvolvimento”;
- “arquitetura covariante efetiva avançada”;
- “validação pública inicial”;
- “compatível em testes efetivos”;
- “promissora, mas ainda dependente de likelihoods completas”;
- “falsificável por dados independentes”.

Evitar frases como:

- “prova final”;
- “teoria já comprovada”;
- “derrubou definitivamente ΛCDM”;
- “resolveu todos os problemas”;
- “não pode ser falsificada”.

---

## 16. Próximos passos prioritários

Os próximos passos mais importantes são:

1. consolidar a ação covariante final;
2. derivar variacionalmente os limites galáctico, cosmológico, PPN e de buracos negros;
3. rodar likelihoods completas com covariâncias;
4. fazer comparação Bayesiana formal;
5. validar SPARC e LITTLE THINGS com pipelines independentes;
6. ampliar testes de lenteamento;
7. preparar release pública versionada do repositório;
8. publicar notebooks reproduzíveis com instruções claras.

---

## 17. Conclusão

A TVGD já possui uma cadeia de testes e uma arquitetura covariante efetiva mais avançada do que uma simples parametrização fenomenológica.

A limitação principal atual não é “falta total de covariante”, mas sim a necessidade de transformar a estrutura efetiva já testada em uma formulação final, variacional, única e publicável, acompanhada por likelihoods completas e reprodução independente.
