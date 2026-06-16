[testes_bancada_TVGD_MEM_README.md](https://github.com/user-attachments/files/29005332/testes_bancada_TVGD_MEM_README.md)
# Testes de Bancada para TVGD-MEM

**Explicação, objetivo físico e relação com a hipótese de estrutura efetiva do vácuo**

Este documento organiza os testes de bancada propostos para a TVGD-MEM.  
A ideia central não é afirmar que os testes provam diretamente a TVGD, mas transformar a hipótese em grandezas mensuráveis.

A cadeia física investigada é:

```
estrutura efetiva do vácuo
-> modulação de fase
-> memória ou perda de coerência
-> assinatura angular
-> ponte entre micro-ondas e óptica
```

Na formulação operacional da TVGD-MEM, a ponte matemática central é:

```
V_corr_optico ≈ M_struct ≈ C_corr_micro
```

Onde:

- `V_corr_optico` é a visibilidade corrigida de franjas ópticas;
- `M_struct` é a métrica efetiva de memória estrutural;
- `C_corr_micro` é a coerência fasorial medida em micro-ondas.

---

## Teste 1 — Calibração de referência da bancada

### O que faz

Mede a resposta do sistema sem perturbação, definindo a linha de base instrumental da bancada.

A transmissão medida pelo VNA pode ser representada por:

```
S21_ref(f) = |S21(f)| * exp(i * phi_ref(f))
```

Também podem ser extraídos:

```
phi_ref(f)
```

e o atraso de grupo:

```
tau_g(f) = -(1 / (2*pi)) * dphi/df
```

### O que ajuda a sustentar

Confiabilidade instrumental.

Antes de qualquer interpretação física, é necessário mostrar que a bancada é suficientemente estável para servir como referência.

---

## Teste 2 — Memória de fase estática em micro-ondas

### O que faz

Compara diferentes condições experimentais:

```
referência
fraca
média
forte
controle
```

Para cada condição, mede-se a fase relativa em cada frequência:

```
Delta_phi_j(f_k) = phi_j(f_k) - phi_ref(f_k)
```

Depois calcula-se a coerência fasorial:

```
C_corr(j) = | (1/N) * sum_{k=1}^{N} exp(i * Delta_phi_j(f_k)) |
```

### Assinatura esperada

```
C_ref > C_fraca > C_media > C_forte
C_controle ≈ 0
```

### O que ajuda a sustentar

A existência de uma memória de fase efetiva mensurável.

Em linguagem TVGD-MEM, esse teste avalia a hipótese operacional:

```
C_corr_micro ≈ M_struct
```

---

## Teste 3 — Robustez temporal

### O que faz

Repete o Teste 2 em instantes diferentes:

```
T0, T1, T2, ...
```

e verifica se a hierarquia se mantém ao longo do tempo.

### Assinatura esperada

```
C_ref > C_fraca > C_media > C_forte
```

em várias repetições temporais.

### O que ajuda a sustentar

Reprodutibilidade do efeito.

Se o comportamento aparece de modo consistente ao longo do tempo, ele deixa de parecer uma flutuação isolada ou um artefato numérico.

---

## Teste 4 — Resposta angular diferencial

### O que faz

Gira a bancada ou altera sua orientação espacial, medindo a resposta em função do ângulo efetivo `chi`.

Na TVGD-MEM, a previsão efetiva de bancada é uma assinatura angular dipolar:

```
epsilon_mem(u, chi) = epsilon_0 * exp(-u) * cos(chi)
```

### Assinatura esperada

Uma variação aproximadamente cossenoidal:

```
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

Altera uma condição associada ao parâmetro efetivo `u`, que controla o fator de blindagem.

Na formulação efetiva, a amplitude residual esperada segue:

```
epsilon_mem(u) = epsilon_0 * exp(-u)
```

### Assinatura esperada

```
u aumenta -> exp(-u) diminui -> epsilon_mem diminui
```

### O que ajuda a sustentar

A ideia de screening ou blindagem estrutural.

Se a resposta experimental segue aproximadamente uma lei exponencial em `u`, isso fortalece a consistência interna da hipótese TVGD-MEM.

---

## Teste 6 — Dupla fenda de micro-ondas

### O que faz

Usa uma fonte de micro-ondas e uma geometria de duas fendas, gerando um padrão interferométrico de intensidade.

Um modelo efetivo para o padrão na tela é:

```
I(x) = I_env(x) * [1 + V0 * M_struct * cos(Delta_phi(x))]
```

### Assinaturas esperadas

- se `M_struct ≈ 1`, as franjas são fortes;
- se `0 < M_struct < 1`, as franjas enfraquecem;
- se `M_struct -> 0`, o padrão interferométrico praticamente desaparece.

### O que ajuda a sustentar

A ponte direta entre geometria de dupla fenda, fase relativa e memória estrutural.

Esse teste conecta micro-ondas à lógica clássica da interferência de duas trajetórias.

---

## Teste 7 — Comparação micro-ondas x luz

### O que faz

Compara uma medida interferométrica em micro-ondas com uma medida de visibilidade de franjas em óptica.

No óptico, usa-se:

```
V_corr_optico
```

No micro-ondas, usa-se:

```
C_corr_micro
```

A hipótese operacional é:

```
V_corr_optico ≈ M_struct ≈ C_corr_micro
```

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

```
C_controle ≈ 0
A_cos_chi ≈ 0
```

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
| Memória estática | `C_corr` | Memória de fase efetiva |
| Robustez temporal | Repetição do efeito | Reprodutibilidade |
| Angular | Dependência `cos(chi)` | Estrutura projetada/direcional |
| Blindagem `u` | Queda `exp(-u)` | Screening estrutural |
| Dupla fenda micro-ondas | Franjas eletromagnéticas | Fase relativa interferométrica |
| Micro x óptico | `C_corr ≈ V_corr` | Universalidade operacional |
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

```
1. Primeiro, tenta-se explicar tudo pela física convencional.
2. Depois, removem-se artefatos instrumentais e falsos positivos.
3. Só então se avalia se resta um excedente compatível com TVGD-MEM.
```
