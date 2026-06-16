# TVGD–MEM — Testes de Bancada

## Explicação, objetivo físico e relação com a hipótese de estrutura efetiva do vácuo

## Objetivo geral

Os testes de bancada não provam diretamente a TVGD, mas transformam a hipótese em grandezas mensuráveis.

A ideia central é verificar se uma possível estrutura efetiva do vácuo pode se manifestar experimentalmente como:

```text
estrutura efetiva do vácuo
        ↓
modulação de fase
        ↓
memória ou perda de coerência
        ↓
assinatura angular
        ↓
ponte entre micro-ondas e óptica
```

Na formulação operacional da TVGD–MEM, a ponte matemática central é:

$$
V_{\rm corr}^{\rm optico}
\approx
M_{\rm struct}
\approx
C_{\rm corr}^{\rm micro}.
$$

---

# Teste 1 — Calibração de referência da bancada

## O que faz

Mede a resposta do sistema sem perturbação, definindo a linha de base da bancada.

A transmissão medida pelo VNA é modelada por:

$$
S_{21}^{\rm ref}(f)
=
|S_{21}(f)|e^{i\phi_{\rm ref}(f)}.
$$

Também podem ser extraídos a fase de referência:

$$
\phi_{\rm ref}(f)
$$

e o atraso de grupo:

$$
\tau_g(f)
=
-\frac{1}{2\pi}
\frac{d\phi}{df}.
$$

## O que ajuda a sustentar

Este teste sustenta a confiabilidade instrumental.

Antes de qualquer interpretação física, é necessário mostrar que a bancada é suficientemente estável para servir como referência.

---

# Teste 2 — Memória de fase estática em micro-ondas

## O que faz

Compara diferentes condições experimentais:

- referência;
- fraca;
- média;
- forte;
- controle.

Em cada condição, mede-se a fase relativa em cada frequência:

$$
\Delta\phi_j(f_k)
=
\phi_j(f_k)
-
\phi_{\rm ref}(f_k).
$$

A coerência fasorial é então calculada por:

$$
C_{\rm corr}(j)
=
\left|
\frac{1}{N}
\sum_{k=1}^{N}
e^{i\Delta\phi_j(f_k)}
\right|.
$$

## Assinatura esperada

$$
C_{\rm ref}
>
C_{\rm fraca}
>
C_{\rm media}
>
C_{\rm forte},
\qquad
C_{\rm controle}
\approx
0.
$$

## O que ajuda a sustentar

Este teste sustenta a existência de uma memória de fase efetiva mensurável.

Em linguagem TVGD–MEM, esse teste avalia a hipótese operacional:

$$
C_{\rm corr}^{\rm micro}
\approx
M_{\rm struct}.
$$

---

# Teste 3 — Robustez temporal

## O que faz

Repete o Teste 2 em instantes diferentes:

$$
T_0,\ T_1,\ T_2,\ldots
$$

e verifica se a hierarquia de coerência se mantém no tempo.

## Assinatura esperada

$$
C_{\rm ref}
>
C_{\rm fraca}
>
C_{\rm media}
>
C_{\rm forte}
$$

em várias repetições temporais.

## O que ajuda a sustentar

Este teste sustenta a reprodutibilidade do efeito.

Se o comportamento aparece de modo consistente ao longo do tempo, ele deixa de parecer uma flutuação isolada ou um artefato numérico.

---

# Teste 4 — Resposta angular diferencial

## O que faz

Gira a bancada ou altera sua orientação espacial, medindo a resposta em função do ângulo:

$$
\chi.
$$

Na TVGD–MEM, a previsão efetiva de bancada é uma assinatura angular dipolar do tipo:

$$
\epsilon_{\rm mem}(u,\chi)
=
\epsilon_0 e^{-u}\cos\chi.
$$

## Assinatura esperada

A assinatura esperada é uma variação aproximadamente cossenoidal, com:

- máximos e mínimos em orientações opostas;
- valor intermediário próximo de 90 graus;
- repetibilidade quando a orientação é retomada.

## O que ajuda a sustentar

Este teste sustenta a hipótese de uma estrutura efetiva projetada ou direcional.

Ele é um dos testes mais importantes porque um ruído puramente aleatório não tende a reproduzir uma assinatura angular limpa e repetível.

---

# Teste 5 — Blindagem estrutural ou variação do parâmetro u

## O que faz

Altera uma condição associada ao parâmetro efetivo:

$$
u,
$$

que controla o fator de blindagem.

Na formulação efetiva, a amplitude residual esperada segue:

$$
\epsilon_{\rm mem}(u)
=
\epsilon_0 e^{-u}.
$$

## Assinatura esperada

$$
u \uparrow
\quad\Rightarrow\quad
e^{-u}\downarrow
\quad\Rightarrow\quad
\epsilon_{\rm mem}\downarrow.
$$

## O que ajuda a sustentar

Este teste sustenta a ideia de screening ou blindagem estrutural.

Se a resposta experimental segue aproximadamente uma lei exponencial em `u`, isso fortalece a consistência interna da hipótese TVGD–MEM.

---

# Teste 6 — Dupla fenda de micro-ondas

## O que faz

Usa uma fonte de micro-ondas e uma geometria de duas fendas, gerando um padrão interferométrico de intensidade.

Um modelo efetivo para o padrão na tela é:

$$
I(x)
=
I_{\rm env}(x)
\left[
1
+
V_0 M_{\rm struct}
\cos(\Delta\phi(x))
\right].
$$

## Assinaturas esperadas

Se:

$$
M_{\rm struct}
\approx
1,
$$

as franjas são fortes.

Se:

$$
0
<
M_{\rm struct}
<
1,
$$

as franjas enfraquecem.

Se:

$$
M_{\rm struct}
\rightarrow
0,
$$

o padrão interferométrico praticamente desaparece.

## O que ajuda a sustentar

Este teste sustenta a ponte direta entre geometria de dupla fenda, fase relativa e memória estrutural.

Ele conecta micro-ondas à lógica clássica da interferência de duas trajetórias.

---

# Teste 7 — Comparação micro-ondas × luz

## O que faz

Compara uma medida interferométrica em micro-ondas com uma medida de visibilidade de franjas em óptica.

No óptico, usa-se a visibilidade corrigida:

$$
V_{\rm corr}^{\rm optico}.
$$

No micro-ondas, usa-se a coerência fasorial:

$$
C_{\rm corr}^{\rm micro}.
$$

A hipótese operacional é:

$$
V_{\rm corr}^{\rm optico}
\approx
M_{\rm struct}
\approx
C_{\rm corr}^{\rm micro}.
$$

## O que ajuda a sustentar

Este teste sustenta a universalidade operacional da métrica de memória de fase.

Se a mesma métrica descreve tanto a perda de visibilidade óptica quanto a perda de coerência fasorial em micro-ondas, a hipótese deixa de parecer um mero ajuste específico de um único sistema.

---

# Teste 8 — Controle nulo

## O que faz

Cria uma condição onde não deveria haver assinatura TVGD–MEM.

Exemplos:

- fase embaralhada artificialmente;
- orientação aleatória;
- canal sem interferência;
- carga casada;
- dados permutados.

## Assinaturas esperadas

$$
C_{\rm controle}
\approx
0,
$$

e:

$$
A_{\cos\chi}
\approx
0.
$$

## O que ajuda a sustentar

Este teste sustenta a ausência de falso positivo.

Ele é essencial porque mostra que o método não encontra efeito onde não deveria existir sinal estruturado.

---

# Teste 9 — Teste adversarial instrumental

## O que faz

Introduz deliberadamente problemas de bancada, como:

- maus contatos;
- reflexões;
- altas perdas;
- ripple;
- ruído;
- falhas de de-embedding;
- desalinhamento;
- deriva térmica.

## O que ajuda a sustentar

Este teste sustenta a separação entre efeito físico e artefato instrumental.

Se o algoritmo começar a acusar “sinal TVGD” em uma bancada claramente defeituosa, então o protocolo ainda não é confiável.

---

# Conjunto mínimo forte de bancada

Um conjunto enxuto e forte de testes seria:

1. calibração/referência;
2. memória estática;
3. robustez temporal;
4. resposta angular;
5. controle nulo;
6. comparação com óptico.

Esse conjunto já oferece uma linha experimental sólida para avaliar a hipótese TVGD–MEM.

---

# O que cada teste ajuda a sustentar em teoria

| Teste | O que mede | O que sustenta |
|---|---|---|
| Calibração | estabilidade da bancada | confiabilidade instrumental |
| Memória estática | `C_corr` | memória de fase efetiva |
| Robustez temporal | repetição do efeito | reprodutibilidade |
| Angular | dependência `cos χ` | estrutura projetada/direcional |
| Blindagem `u` | queda `e^{-u}` | screening estrutural |
| Dupla fenda micro-ondas | franjas eletromagnéticas | fase relativa interferométrica |
| Micro-ondas × óptico | `C_corr ≈ V_corr` | universalidade operacional |
| Controle nulo | ausência de sinal | não é falso positivo |
| Adversarial | falhas instrumentais | separa sinal de artefato |

---

# Conclusão

Se os testes de bancada passarem, a conclusão correta e cientificamente prudente é:

> Os testes de bancada não provam definitivamente a TVGD nem a estrutura fundamental do vácuo. Porém, se passarem com controles adequados, eles fornecem evidência experimental compatível com a hipótese de uma estrutura efetiva do vácuo, manifestada como memória de fase, perda controlada de coerência fasorial, resposta angular e equivalência entre visibilidade óptica e coerência de micro-ondas.

Em resumo, a TVGD deixa de ser apenas uma ideia matemática abstrata e passa a ter assinaturas experimentais mensuráveis.
