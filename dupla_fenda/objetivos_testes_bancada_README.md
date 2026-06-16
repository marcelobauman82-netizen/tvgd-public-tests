[testes_bancada_TVGD_MEM_README_equacoes_corrigidas.md](https://github.com/user-attachments/files/29005504/testes_bancada_TVGD_MEM_README_equacoes_corrigidas.md)
# Testes de Bancada para TVGD-MEM

**Explicação, objetivo físico e relação com a hipótese de estrutura efetiva do vácuo**

Este documento organiza os testes de bancada propostos para a TVGD-MEM.

A ideia central não é afirmar que os testes provam diretamente a TVGD, mas transformar a hipótese em grandezas mensuráveis.

## Objetivo geral

A cadeia física investigada é:

```text
estrutura efetiva do vácuo
-> modulação de fase
-> memória ou perda de coerência
-> assinatura angular
-> ponte entre micro-ondas e óptica
```

Na formulação operacional da TVGD-MEM, a ponte matemática central é:

$$
V_{\mathrm{corr}}^{\mathrm{optico}} \approx M_{\mathrm{struct}} \approx C_{\mathrm{corr}}^{\mathrm{micro}}
$$

Onde:

- $V_{\mathrm{corr}}^{\mathrm{optico}}$ é a visibilidade corrigida de franjas ópticas;
- $M_{\mathrm{struct}}$ é a métrica efetiva de memória estrutural;
- $C_{\mathrm{corr}}^{\mathrm{micro}}$ é a coerência fasorial medida em micro-ondas.

---

## Teste 1 — Calibração de referência da bancada

### O que faz

Mede a resposta do sistema sem perturbação, definindo a linha de base instrumental da bancada.

A transmissão medida pelo VNA pode ser representada por:

$$
S_{21}^{\mathrm{ref}}(f)=|S_{21}(f)|e^{i\phi_{\mathrm{ref}}(f)}
$$

Também podem ser extraídos a fase de referência $\phi_{\mathrm{ref}}(f)$ e o atraso de grupo:

$$
\tau_g(f)=-\frac{1}{2\pi}\frac{d\phi}{df}
$$

### O que ajuda a sustentar

Confiabilidade instrumental.

Antes de qualquer interpretação física, é necessário mostrar que a bancada é suficientemente estável para servir como referência.

---

## Teste 2 — Memória de fase estática em micro-ondas

### O que faz

Compara diferentes condições experimentais:

```text
referência
fraca
média
forte
controle
```

Para cada condição, mede-se a fase relativa em cada frequência:

$$
\Delta\phi_j(f_k)=\phi_j(f_k)-\phi_{\mathrm{ref}}(f_k)
$$

Depois calcula-se a coerência fasorial:

$$
C_{\mathrm{corr}}(j)=\left|\frac{1}{N}\sum_{k=1}^{N}e^{i\Delta\phi_j(f_k)}\right|
$$

### Assinatura esperada

$$
C_{\mathrm{ref}}>C_{\mathrm{fraca}}>C_{\mathrm{media}}>C_{\mathrm{forte}}
$$

e

$$
C_{\mathrm{controle}}\approx 0
$$

### O que ajuda a sustentar

A existência de uma memória de fase efetiva mensurável.

Em linguagem TVGD-MEM, esse teste avalia a hipótese operacional:

$$
C_{\mathrm{corr}}^{\mathrm{micro}}\approx M_{\mathrm{struct}}
$$

---

## Teste 3 — Robustez temporal

### O que faz

Repete o Teste 2 em instantes diferentes:

$$
T_0,\ T_1,\ T_2,\ldots
$$

e verifica se a hierarquia se mantém ao longo do tempo.

### Assinatura esperada

$$
C_{\mathrm{ref}}>C_{\mathrm{fraca}}>C_{\mathrm{media}}>C_{\mathrm{forte}}
$$

em várias repetições temporais.

### O que ajuda a sustentar

Reprodutibilidade do efeito.

Se o comportamento aparece de modo consistente ao longo do tempo, ele deixa de parecer uma flutuação isolada ou um artefato numérico.

---

## Teste 4 — Resposta angular diferencial

### O que faz

Gira a bancada ou altera sua orientação espacial, medindo a resposta em função do ângulo efetivo $\chi$.

Na TVGD-MEM, a previsão efetiva de bancada é uma assinatura angular dipolar:

$$
\epsilon_{\mathrm{mem}}(u,\chi)=\epsilon_0 e^{-u}\cos\chi
$$

### Assinatura esperada

Uma variação aproximadamente cossenoidal:

```text
máximo em uma orientação
mínimo na orientação oposta
valor intermediário próximo de 90 graus
```

### O que ajuda a sustentar

Uma estrutura efetiva projetada ou direcional.

Este é um dos testes mais importantes, porque um ruído puramente aleatório não tende a reproduzir uma assinatura angular limpa e repetível.

---

## Teste 5 — Blindagem estrutural ou variação do parâmetro u

### O que faz

Altera uma condição associada ao parâmetro efetivo $u$, que controla o fator de blindagem.

Na formulação efetiva, a amplitude residual esperada segue:

$$
\epsilon_{\mathrm{mem}}(u)=\epsilon_0 e^{-u}
$$

### Assinatura esperada

$$
u\uparrow \Rightarrow e^{-u}\downarrow \Rightarrow \epsilon_{\mathrm{mem}}\downarrow
$$

### O que ajuda a sustentar

A ideia de screening ou blindagem estrutural.

Se a resposta experimental segue aproximadamente uma lei exponencial em $u$, isso fortalece a consistência interna da hipótese TVGD-MEM.

---

## Teste 6 — Dupla fenda de micro-ondas

### O que faz

Usa uma fonte de micro-ondas e uma geometria de duas fendas, gerando um padrão interferométrico de intensidade.

Um modelo efetivo para o padrão na tela é:

$$
I(x)=I_{\mathrm{env}}(x)\left[1+V_0M_{\mathrm{struct}}\cos(\Delta\phi(x))\right]
$$

### Assinaturas esperadas

- se $M_{\mathrm{struct}}\approx 1$, as franjas são fortes;
- se $0<M_{\mathrm{struct}}<1$, as franjas enfraquecem;
- se $M_{\mathrm{struct}}\to 0$, o padrão interferométrico praticamente desaparece.

### O que ajuda a sustentar

A ponte direta entre geometria de dupla fenda, fase relativa e memória estrutural.

Esse teste conecta micro-ondas à lógica clássica da interferência de duas trajetórias.

---

## Teste 7 — Comparação micro-ondas x luz

### O que faz

Compara uma medida interferométrica em micro-ondas com uma medida de visibilidade de franjas em óptica.

No óptico, usa-se:

$$
V_{\mathrm{corr}}^{\mathrm{optico}}
$$

No micro-ondas, usa-se:

$$
C_{\mathrm{corr}}^{\mathrm{micro}}
$$

A hipótese operacional é:

$$
V_{\mathrm{corr}}^{\mathrm{optico}}\approx M_{\mathrm{struct}}\approx C_{\mathrm{corr}}^{\mathrm{micro}}
$$

### O que ajuda a sustentar

Universalidade operacional da métrica de memória de fase.

Se a mesma métrica descreve tanto a perda de visibilidade óptica quanto a perda de coerência fasorial em micro-ondas, a hipótese deixa de parecer um mero ajuste específico de um único sistema.

---

## Teste 8 — Controle nulo

### O que faz

Cria uma condição onde não deveria haver assinatura TVGD-MEM.

Exemplos:

- fase embaralhada artificialmente;
- orientação aleatória;
- canal sem interferência;
- carga casada;
- dados permutados.

### Assinaturas esperadas

$$
C_{\mathrm{controle}}\approx 0
$$

e

$$
A_{\cos\chi}\approx 0
$$

### O que ajuda a sustentar

Ausência de falso positivo.

Esse teste é essencial porque mostra que o método não encontra sinal onde não deveria existir sinal estruturado.

---

## Teste 9 — Teste adversarial instrumental

### O que faz

Introduz deliberadamente problemas de bancada, como:

- maus contatos;
- reflexões;
- altas perdas;
- ripple;
- ruído;
- falhas de de-embedding;
- desalinhamento;
- deriva térmica.

### O que ajuda a sustentar

Separação entre efeito físico e artefato instrumental.

Se o algoritmo começar a acusar "sinal TVGD" em uma bancada claramente defeituosa, então o protocolo ainda não é confiável.

---

## Conjunto mínimo forte de bancada

Um conjunto enxuto e forte seria:

1. Calibração/referência;
2. Memória estática;
3. Robustez temporal;
4. Resposta angular;
5. Controle nulo;
6. Comparação com óptico.

Esse conjunto já oferece uma linha experimental sólida para avaliar a hipótese TVGD-MEM.

---

## Resumo dos testes

| Teste | O que mede | O que sustenta |
|---|---|---|
| Calibração | Estabilidade da bancada | Confiabilidade instrumental |
| Memória estática | $C_{\mathrm{corr}}$ | Memória de fase efetiva |
| Robustez temporal | Repetição do efeito | Reprodutibilidade |
| Angular | Dependência $\cos\chi$ | Estrutura projetada/direcional |
| Blindagem $u$ | Queda $e^{-u}$ | Screening estrutural |
| Dupla fenda micro-ondas | Franjas eletromagnéticas | Fase relativa interferométrica |
| Micro x óptico | $C_{\mathrm{corr}}\approx V_{\mathrm{corr}}$ | Universalidade operacional |
| Controle nulo | Ausência de sinal | Não é falso positivo |
| Adversarial | Falhas instrumentais | Separação entre sinal e artefato |

---

## Conclusão

Se os testes de bancada passarem, a conclusão correta e cientificamente prudente é:

> Os testes de bancada não provam definitivamente a TVGD nem a estrutura fundamental do vácuo. Porém, se passarem com controles adequados, eles fornecem evidência experimental compatível com a hipótese de uma estrutura efetiva do vácuo, manifestada como memória de fase, perda controlada de coerência fasorial, resposta angular e equivalência entre visibilidade óptica e coerência de micro-ondas.

Em resumo:

> A TVGD deixa de ser apenas uma ideia matemática abstrata e passa a ter assinaturas experimentais mensuráveis.

---

## Observação metodológica

Este material deve ser lido como documentação de apoio para testes de bancada.

Ele não apresenta uma prova definitiva de nova física.

Ele organiza um protocolo de exclusão:

```text
1. Primeiro, tenta-se explicar tudo pela física convencional.
2. Depois, removem-se artefatos instrumentais e falsos positivos.
3. Só então se avalia se resta um excedente compatível com TVGD-MEM.
```

---

## Nota sobre renderização no GitHub

Este arquivo usa equações em formato GitHub Markdown com `$$ ... $$`.

Se a visualização de equações não aparecer corretamente em algum ambiente, as mesmas equações continuam legíveis no texto e podem ser copiadas diretamente para LaTeX ou Overleaf.
