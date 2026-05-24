# TVGD — Limitações Conhecidas e Status da Covariante

Este documento registra, de forma curta e objetiva, o estado atual da formulação covariante da **Teoria do Vácuo Gravitacional Dispersivo (TVGD)**.

A TVGD não deve mais ser apresentada como uma simples fenomenologia. O estágio atual é o de uma **teoria efetiva covariante avançada**, com ação estruturada, recuperação do limite galáctico, mecanismo de screening local, setor vetorial compatível com ondas gravitacionais e extensão efetiva para buracos negros por saturação estrutural.

O ponto que ainda falta não é “ter uma covariante”. A covariante efetiva já existe. O que falta é o fechamento final da teoria fundamental.

---

## 1. Status atual

A formulação atual pode ser resumida assim:

| Nível | Situação |
|---|---|
| Fenomenologia galáctica isolada | Superada |
| Teoria efetiva covariante | Forte e em desenvolvimento avançado |
| Teoria fundamental final | Ainda não fechada |

A TVGD já possui uma arquitetura covariante efetiva capaz de conectar os setores principais da teoria. Ainda assim, por rigor científico, ela deve ser descrita como **teoria efetiva em desenvolvimento**, e não como teoria final comprovada.

---

## 2. Estrutura covariante efetiva

A Camada II da TVGD organiza a teoria por meio de uma ação efetiva com setores separados:

```math
S_{TVGD}
=
S_g
+
S_{proj}
+
S_{theta}
+
S_U
+
S_{lambda}
+
S_m
```

onde:

| Termo | Papel físico |
|---|---|
| \(S_g\) | gravidade métrica |
| \(S_{proj}\) | resposta estrutural projetada do vácuo |
| \(S_{theta}\) | setor cosmológico homogêneo |
| \(S_U\) | setor vetorial tipo Einstein-Aether |
| \(S_{lambda}\) | normalização do vetor temporal |
| \(S_m\) | matéria bariônica |

Essa estrutura já é suficiente para tratar a TVGD como uma teoria efetiva covariante testável.

---

## 3. Recuperação do regime galáctico

No limite galáctico quase-estático, a TVGD recupera a relação central entre aceleração bariônica e aceleração efetiva:

```math
g_{TVGD}
=
\frac{g_{bar}}{P(u)}
```

com

```math
P(u)
=
1-e^{-u}
```

e

```math
u
=
\sqrt{\frac{g_{bar}}{a_0}} .
```

No regime de baixa aceleração, \(u \ll 1\), temos \(P(u) \simeq u\). Assim:

```math
g_{TVGD}
\simeq
\sqrt{g_{bar}a_0}.
```

Esse resultado é um dos pontos fortes da TVGD: o comportamento galáctico do tipo MOND-like surge como limite efetivo da resposta estrutural do vácuo, e não apenas como uma regra empírica isolada.

---

## 4. Screening local

Em regiões de alta aceleração, \(u \gg 1\), a função de resposta tende a:

```math
P(u) \to 1.
```

Portanto:

```math
g_{TVGD}
\to
g_{bar}.
```

Esse é o mecanismo de screening local da TVGD. Ele permite que a teoria recupere aproximadamente a gravidade usual no Sistema Solar e em regimes de alta aceleração, preservando testes como GPS, atraso de Shapiro, órbitas planetárias e limites PPN.

---

## 5. Setor vetorial e compatibilidade com ondas gravitacionais

A formulação covariante inclui um vetor temporal unitário \(U^\mu\), com restrição:

```math
U^\mu U_\mu = -1.
```

O setor vetorial pode ser escolhido em um ramo compatível com a velocidade das ondas gravitacionais. A condição principal é:

```math
c_1+c_3=0.
```

Nesse ramo:

```math
c_T=1.
```

Isso é essencial para compatibilidade com as restrições observacionais de ondas gravitacionais.

Os testes efetivos já realizados também indicam uma região de parâmetros fracos compatível com limites PPN fortes. Portanto, o setor vetorial não é uma fragilidade imediata da teoria no ramo testado.

---

## 6. Buracos negros e saturação estrutural

No setor de buracos negros, a TVGD usa uma saturação estrutural do vácuo. Uma forma efetiva é:

```math
Q_\psi(x)
=
1-
\exp\left[
-\left(\frac{x}{x_c}\right)^n
\right],
\qquad
x=\frac{r}{r_s}.
```

A métrica efetiva pode ser escrita como:

```math
A(x)
=
1-
\frac{Q_\psi(x)}{x}.
```

Esse setor preserva os principais benchmarks externos:

```math
r_h/r_s \simeq 1,
\qquad
r_{ph}/r_s \simeq 1.5,
\qquad
b_{shadow}/r_s \simeq 2.598076.
```

A interpretação é que a resposta estrutural do vácuo pode regularizar o núcleo efetivo sem destruir o exterior observado. O ponto ainda em aberto é derivar essa saturação diretamente da ação covariante final.

---

## 7. O que ainda falta

As limitações essenciais da TVGD são poucas e bem definidas:

1. **Fechamento da ação final:** transformar a Camada II efetiva em uma formulação final única.
2. **Equações acopladas completas:** derivar explicitamente o sistema completo envolvendo métrica, campo estrutural, vetor, setor cosmológico e matéria.
3. **Conservação covariante:** demonstrar de forma completa a conservação do tensor energia-momento total.
4. **Perturbações cosmológicas:** implementar o sistema completo em códigos como CLASS ou CAMB.
5. **Origem de \(a_0\):** derivar rigorosamente a escala \(a_0\), possivelmente ligada a escalas cosmológicas como \(cH_0\).
6. **Likelihoods completas:** substituir testes efetivos por análises formais com matrizes de covariância, priors e comparação Bayesiana.
7. **Validação independente:** permitir reprodução externa dos principais resultados.

Esses pontos não indicam ausência de base teórica. Eles indicam o caminho necessário para transformar a TVGD de uma teoria efetiva covariante avançada em uma teoria fundamental final.

---

## 8. Formulação curta recomendada

A forma mais correta de apresentar o estágio atual é:

> A TVGD já possui uma formulação covariante efetiva, com recuperação do regime galáctico, screening local, setor vetorial compatível com \(c_T=1\) e saturação estrutural em buracos negros. A limitação atual não é a ausência de covariante, mas o fechamento final da ação fundamental, das equações acopladas, das perturbações cosmológicas e da validação estatística completa.

---

## 9. Conclusão

A TVGD avançou além de uma fenomenologia isolada. A teoria já possui uma base covariante efetiva capaz de organizar seus principais setores físicos.

O status correto é:

> **A TVGD possui uma ação covariante efetiva testável; o que ainda falta é o fechamento fundamental e a validação estatística completa.**
