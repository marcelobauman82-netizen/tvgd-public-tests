# TVGD — Limitações Conhecidas e Status da Covariante

Este documento resume o estado atual da **Teoria do Vácuo Gravitacional Dispersivo (TVGD)** após o avanço da formulação covariante efetiva.

A conclusão principal é simples: a TVGD **não deve mais ser descrita como uma fenomenologia sem covariante**. O estágio atual é mais forte. A teoria já possui uma **Camada II covariante efetiva**, com ação estruturada, campo escalar, vetor temporal unitário, projetor espacial, setor vetorial tipo Einstein–Æther, restrição de normalização, mecanismo de screening e recuperação do limite galáctico.

O que ainda falta não é “encontrar uma covariante”, mas sim fechar a formulação fundamental final: derivar todas as equações acopladas, demonstrar conservação total, controlar perturbações cosmológicas completas, provar estabilidade global e deduzir rigorosamente a origem de \(a_0\).

---

## 1. Status em uma frase

A formulação atual da TVGD pode ser descrita como:

> **Uma teoria gravitacional efetiva covariante, com núcleo galáctico recuperado, setor vetorial controlado, limite local preservado e extensão forte por saturação estrutural, ainda em processo de fechamento fundamental.**

Portanto, a classificação honesta é:

| Nível | Status |
|---|---|
| Fenomenologia simples | Superada |
| Teoria efetiva covariante | Forte e avançada |
| Teoria fundamental final | Ainda em desenvolvimento |

---

## 2. Arquitetura covariante efetiva

A Camada II da TVGD organiza a teoria como uma ação efetiva composta por setores físicos distintos:

$$
S_{\mathrm{TVGD}}
=
S_g
+
S_{\mathrm{proj}}
+
S_{\theta}
+
S_U
+
S_{\lambda}
+
S_m .
$$

Aqui:

| Termo | Interpretação |
|---|---|
| \(S_g\) | setor gravitacional tipo Einstein–Hilbert |
| \(S_{\mathrm{proj}}\) | setor estrutural projetado do vácuo |
| \(S_{\theta}\) | setor homogêneo/cosmológico |
| \(S_U\) | setor vetorial tipo Einstein–Æther |
| \(S_{\lambda}\) | restrição de normalização do vetor |
| \(S_m\) | ação da matéria |

Essa estrutura já permite tratar a TVGD como uma teoria efetiva covariante testável, não apenas como uma relação empírica ajustada a galáxias.

---

## 3. Campo estrutural e projetor espacial

A covariante usa um campo escalar estrutural \(\psi\) e um vetor temporal unitário \(U^\mu\). O vetor define o referencial físico local do vácuo.

A métrica espacial projetada é:

$$
h^{\mu\nu}
=
g^{\mu\nu}
+
U^\mu U^\nu .
$$

A variável estrutural projetada é:

$$
Y_s
=
\frac{
h^{\mu\nu}
\nabla_{\mu}\psi
\nabla_{\nu}\psi
}{
a_0^2
}.
$$

Esse ponto é importante porque o setor galáctico da TVGD não depende simplesmente de um gradiente escalar arbitrário. Ele depende da parte espacial do gradiente no referencial definido por \(U^\mu\).

---

## 4. Recuperação do regime galáctico

No limite quase-estático galáctico, a dinâmica temporal do campo estrutural fica desprezível:

$$
U^\mu \nabla_\mu \psi \simeq 0 .
$$

Nesse regime, o setor projetado recupera a equação efetiva:

$$
\nabla \cdot
\left[
P(u)\nabla\psi
\right]
=
4\pi G\rho_{\mathrm{bar}} .
$$

Em simetria esférica, isso leva à relação central:

$$
P(u)\,g_{\mathrm{TVGD}}
=
g_{\mathrm{bar}} .
$$

Logo:

$$
g_{\mathrm{TVGD}}
=
\frac{g_{\mathrm{bar}}}{P(u)} .
$$

A resposta estrutural usada no regime galáctico é:

$$
P(u)
=
1-e^{-u},
\qquad
u
=
\sqrt{\frac{g_{\mathrm{bar}}}{a_0}} .
$$

No limite profundo, \(u\ll 1\), temos:

$$
P(u)\simeq u .
$$

Portanto:

$$
g_{\mathrm{TVGD}}
\simeq
\sqrt{g_{\mathrm{bar}}a_0}.
$$

Esse é um dos pontos mais fortes da formulação atual: o limite galáctico não aparece isolado. Ele pode ser recuperado como limite efetivo de uma estrutura covariante.

---

## 5. A função \(P(u)\) não é apenas um ajuste

A forma

$$
P(u)=1-e^{-u}
$$

pode ser interpretada como resposta constitutiva do vácuo gravitacional.

Na formulação covariante, essa resposta pode ser associada a uma função efetiva de ação. Em termos qualitativos, a teoria não precisa tratar \(P(u)\) apenas como uma curva empírica. Ela pode ser vista como derivada de uma função cinética ou de uma ação projetada dependente de \(Y_s\).

Isso fortalece a posição teórica da TVGD: a função central do regime galáctico tem um caminho covariante para ser incorporada à ação.

---

## 6. Estabilidade efetiva

O setor escalar possui uma estrutura de estabilidade controlada por uma combinação efetiva do tipo:

$$
D(u)
=
P(u)
+
u^2 e^{-u}.
$$

A velocidade efetiva de propagação pode ser escrita como:

$$
c_s^2
=
\frac{P(u)}{D(u)} .
$$

Nos ramos auditados, a teoria satisfaz:

$$
P(u)>0,
\qquad
D(u)>0,
\qquad
0<c_s^2\leq 1 .
$$

Isso indica ausência de ghost escalar, ausência de instabilidade de gradiente e propagação causal no domínio efetivo testado.

Esse resultado não substitui uma análise Hamiltoniana global do sistema completo, mas remove um risco imediato importante.

---

## 7. Screening local

Em regimes de alta aceleração:

$$
u\gg 1 .
$$

Nesse limite:

$$
P(u)\to 1 .
$$

Assim:

$$
g_{\mathrm{TVGD}}
\to
g_{\mathrm{bar}} .
$$

Esse é o mecanismo de screening local da TVGD. Ele permite que a teoria modifique regimes galácticos de baixa aceleração sem destruir automaticamente os testes do Sistema Solar, como órbitas planetárias, GPS, atraso de Shapiro e deflexão solar da luz.

---

## 8. Setor vetorial e ondas gravitacionais

A formulação inclui um setor vetorial tipo Einstein–Æther, com acoplamentos \(c_1,c_2,c_3,c_4\).

O ramo mais importante impõe:

$$
c_1+c_3=0 .
$$

Com essa escolha, a velocidade tensorial fica preservada:

$$
c_T=1 .
$$

Esse ponto é essencial para compatibilidade com restrições de ondas gravitacionais.

Nos testes efetivos já realizados, há também uma região de acoplamento fraco compatível com limites PPN fortes. Isso torna o setor vetorial viável no ramo testado.

---

## 9. Restrição de normalização

O vetor temporal satisfaz:

$$
U^\mu U_\mu=-1 .
$$

Essa condição é imposta por um multiplicador de Lagrange.

No limite galáctico quase-estático, a restrição não introduz uma fonte espúria dominante. Isso é importante porque o setor vetorial deve organizar a estrutura causal da teoria sem destruir a recuperação do regime galáctico.

---

## 10. Cosmologia

A formulação covariante separa o setor galáctico projetado do setor cosmológico homogêneo.

Em um fundo FLRW homogêneo, a parte espacial projetada pode se anular ou ficar subdominante. Assim, o mesmo mecanismo que explica o regime galáctico não precisa dominar automaticamente a expansão homogênea do Universo.

O setor cosmológico pode ser carregado por uma variável homogênea, frequentemente representada por \(\theta\), ou por uma densidade efetiva estrutural \(\Omega_\psi(a)\).

Essa separação é uma vantagem teórica: ela permite que a TVGD preserve o regime primordial e, ao mesmo tempo, modele transições tardias associadas à tensão de Hubble.

A limitação atual é que ainda faltam perturbações cosmológicas completas, implementação formal em CLASS ou CAMB e likelihoods conjuntas com CMB, BAO, supernovas e crescimento de estruturas.

---

## 11. Buracos negros e saturação estrutural

O setor de buracos negros usa uma saturação estrutural do vácuo. Uma forma efetiva estudada é:

$$
Q_\psi(x)
=
1-
\exp\left[
-\left(
\frac{x}{x_c}
\right)^n
\right],
\qquad
x=\frac{r}{r_s}.
$$

A métrica efetiva pode ser escrita como:

$$
A(x)
=
1-
\frac{Q_\psi(x)}{x}.
$$

Para \(x\gg x_c\), a solução recupera o exterior semelhante a Schwarzschild.

Os benchmarks externos preservados são:

$$
r_h/r_s\simeq 1,
\qquad
r_{\mathrm{ph}}/r_s\simeq 1.5,
\qquad
b_{\mathrm{shadow}}/r_s\simeq 2.598076 .
$$

Esse setor é forte porque preserva observáveis externos enquanto regulariza o núcleo efetivo.

A limitação ainda existente é que a saturação estrutural de buracos negros precisa ser derivada de forma completa a partir da ação covariante final, incluindo estabilidade dinâmica, modos quase-normais, acreção e comparação sistemática com EHT e LIGO/Virgo/KAGRA.

---

## 12. Origem de \(a_0\)

A escala \(a_0\) é central para o regime galáctico.

A TVGD sugere uma conexão natural entre \(a_0\) e escalas cosmológicas, da ordem de:

$$
a_0
\sim
\frac{cH_0}{2\pi}.
$$

Isso reduz a impressão de que \(a_0\) seja apenas um parâmetro arbitrário. No entanto, a dedução fundamental ainda não está fechada.

Ainda falta demonstrar rigorosamente se \(a_0\) emerge de um modo estrutural do vácuo, de uma condição de horizonte, de uma escala de relaxação cosmológica ou de uma combinação desses fatores.

Portanto, o status correto é: \(a_0\) é uma escala efetiva bem motivada, mas ainda não uma constante fundamental deduzida de primeiros princípios.

---

## 13. Limitações reais que permanecem

As limitações principais da TVGD, no estágio atual, são:

1. derivar todas as equações acopladas da ação covariante completa;
2. demonstrar explicitamente a conservação covariante total;
3. calcular os potenciais fracos \(\Phi\) e \(\Psi\);
4. derivar o slip gravitacional e a combinação de lenteamento \(\Phi+\Psi\);
5. implementar perturbações cosmológicas completas;
6. provar estabilidade global do sistema acoplado;
7. deduzir rigorosamente a origem de \(a_0\);
8. formalizar a variável universal de saturação;
9. substituir testes efetivos por likelihoods completas;
10. realizar comparação Bayesiana formal contra \(\Lambda\)CDM, MOND e modelos com halos;
11. obter reprodução independente por terceiros.

Esses pontos não anulam os avanços já obtidos. Eles definem o caminho para transformar a TVGD de teoria efetiva covariante avançada em uma teoria fundamental final.

---

## 14. Setores fortes no estágio atual

Os seguintes blocos podem ser apresentados como fortes dentro da formulação efetiva atual:

| Setor | Estado |
|---|---|
| Limite galáctico | Recuperado |
| Limite profundo \(g\simeq\sqrt{g_{\mathrm{bar}}a_0}\) | Recuperado |
| Screening local | Preservado |
| Ação covariante efetiva | Construída em Camada II |
| Setor vetorial com \(c_T=1\) | Compatível no ramo testado |
| PPN | Compatível no ramo fraco auditado |
| Estabilidade escalar/vetorial | Favorável nos domínios testados |
| Fundo cosmológico | Separado do setor projetado galáctico |
| Buracos negros | Benchmarks externos preservados |
| Falsificabilidade | Critérios explícitos definidos |

---

## 15. Setores promissores, mas ainda não fechados

Os seguintes pontos devem ser tratados com cautela:

| Setor | O que falta |
|---|---|
| Origem de \(a_0\) | Derivação fundamental |
| Setor cosmológico homogêneo | Potencial e dinâmica final |
| Perturbações cosmológicas | Implementação completa |
| Lenteamento relativístico | Derivar \(\Phi\), \(\Psi\) e slip |
| Estabilidade global | Hamiltoniano completo |
| Saturação universal | Variável covariante única |
| Buracos negros | Derivação variacional final |
| Estatística | Likelihoods completas e evidência Bayesiana |
| Reprodutibilidade | Reprodução independente |

---

## 16. Linguagem científica recomendada

A formulação mais justa é:

> A TVGD possui uma formulação covariante efetiva avançada, capaz de recuperar o regime galáctico, preservar o limite local, incorporar um setor vetorial compatível com \(c_T=1\) e modelar saturação estrutural em buracos negros. O que ainda falta é o fechamento fundamental completo, com derivação acoplada, conservação total, perturbações cosmológicas, estabilidade global e validação estatística formal.

Essa formulação é mais correta do que dizer que a TVGD “ainda não tem covariante”. Ela tem uma covariante efetiva. O que falta é transformá-la em uma teoria final.

---

## 17. Conclusão

O desenvolvimento covariante muda a forma correta de apresentar as limitações da TVGD.

Antes, a limitação parecia ser encontrar uma formulação covariante. Agora, o status é mais forte: a TVGD já possui uma Camada II covariante efetiva, com recuperação dos principais limites físicos e auditorias preliminares de estabilidade.

O que permanece em aberto está em um nível mais alto: ação fundamental única, equações acopladas completas, conservação covariante total, estabilidade global, perturbações cosmológicas completas, origem rigorosa de \(a_0\), variável universal de saturação e validação estatística independente.

Em resumo:

> **A TVGD já possui uma ação covariante efetiva testável; o que falta é o fechamento fundamental e a validação estatística completa.**
