# TVGD — Limitações Conhecidas e Status da Covariante

Este documento atualiza as limitações conhecidas da **Teoria do Vácuo Gravitacional Dispersivo (TVGD)** à luz do documento `covariante.pdf`.

A conclusão principal é que a TVGD **não está mais no estágio de mera fenomenologia sem base covariante**. O estado atual é mais forte: existe uma **Camada II covariante efetiva**, com ação explícita, campo escalar estrutural, vetor temporal unitário, projetor espacial, setor tipo Einstein-Aether, restrição de normalização, recuperação do limite galáctico, estabilidade efetiva, screening local, setor cosmológico separado e setor forte por saturação estrutural.

Ainda assim, a teoria deve ser apresentada como **formulação covariante efetiva avançada**, não como teoria fundamental final.

---

## 1. Resumo honesto do status atual

A TVGD deve ser descrita neste estágio como:

```text
Teoria efetiva covariante avançada, com ação de Camada II explícita,
limites físicos recuperados e auditorias de estabilidade preliminares,
mas ainda sem fechamento fundamental completo.
```

O avanço importante é que a covariante atual já contém uma ação efetiva consolidada do tipo:

```text
S_TVGD = S_g + S_proj + S_theta + S_U + S_lambda + S_m
```

ou, na versão projetada mínima:

```text
S_TVGD = S_g + S_proj + S_lambda + S_U + S_m
```

onde:

```text
S_g       = setor Einstein-Hilbert
S_proj    = setor estrutural projetado do vácuo
S_theta   = setor homogêneo/cosmológico
S_U       = setor vetorial tipo Einstein-Aether
S_lambda  = restrição de normalização do vetor U
S_m       = ação da matéria
```

Portanto, a limitação atual não é mais “falta de covariante”. A limitação correta é:

```text
falta transformar a Camada II covariante efetiva em uma teoria fundamental final,
com derivação completa das equações acopladas, conservação total, perturbações
cosmológicas, estabilidade Hamiltoniana completa e origem rigorosa de a0.
```

---

## 2. O que já está forte na covariante

### 2.1 Ação projetada do vácuo

A TVGD usa um campo escalar estrutural `psi` e um vetor temporal unitário `U^mu`.

O projetor espacial é:

```text
h^{mu nu} = g^{mu nu} + U^mu U^nu
```

A variável projetada estrutural é:

```text
Y_s = h^{mu nu} grad_mu psi grad_nu psi / a0^2
```

Essa construção permite separar a resposta espacial do vácuo no referencial físico definido por `U^mu`.

### 2.2 Recuperação do regime galáctico

No limite quase-estático galáctico:

```text
U^mu grad_mu psi ≈ 0
h^{ij} -> delta^{ij}
```

a equação estrutural reduz para:

```text
div[ P(u) grad psi ] = 4 pi G rho_bar
```

Em simetria esférica:

```text
P(u) g_TVGD = g_bar
```

logo:

```text
g_TVGD = g_bar / P(u)
```

com:

```text
P(u) = 1 - exp(-u)
u = sqrt(g_bar / a0)
```

No limite profundo:

```text
g_TVGD ≈ sqrt(g_bar * a0)
```

Assim, a Camada II covariante recupera o núcleo galáctico da TVGD.

### 2.3 A função P(u) não fica apenas como ajuste

O documento covariante registra duas rotas complementares:

1. Uma função cinética covariante `Q_cov(u)` cuja derivada gera:

```text
P(u) = 1 - exp(-u)
```

2. Uma ação projetada `H(Y_s)` definida para gerar o fator constitutivo no limite fraco.

Isso fortalece a interpretação de que `P(u)` não é apenas uma curva ajustada, mas uma resposta estrutural derivável de uma função de ação efetiva.

### 2.4 Estabilidade e causalidade efetiva

O núcleo escalar foi auditado com:

```text
D(u) = P(u) + u^2 exp(-u)
```

e velocidade efetiva:

```text
c_s^2 = P(u) / D(u)
```

O domínio testado satisfaz:

```text
P(u) > 0
D(u) > 0
0 < c_s^2 <= 1
```

No limite profundo, a velocidade efetiva tende para valor positivo da ordem de `2/3`, e no limite de alta aceleração tende para `1`.

Isso remove, no setor testado, risco imediato de ghost escalar, instabilidade de gradiente e propagação superluminal.

### 2.5 Screening local

Para alta aceleração:

```text
u >> 1
P(u) -> 1
alpha(u) = 1/P(u) -> 1
```

Logo:

```text
g_TVGD -> g_bar
```

Esse mecanismo fornece blindagem local automática para Terra, Sistema Solar, GPS, Shapiro, deflexão solar e órbitas planetárias, dentro do ramo efetivo testado.

### 2.6 Setor vetorial

O setor vetorial é do tipo Einstein-Aether, com coeficientes `c1`, `c2`, `c3`, `c4`.

A condição central usada é:

```text
c1 + c3 = 0
```

Isso preserva:

```text
c_T = 1
```

ou seja, velocidade tensorial igual à velocidade da luz.

O documento registra ramo de acoplamento fraco com modos positivos e parâmetros PPN preferenciais pequenos. No limite quase-estático galáctico:

```text
grad U ≈ 0
L_U ≈ 0
T_U ≈ 0
```

Portanto, o setor vetorial não destrói a fonte galáctica recuperada pelo setor projetado.

### 2.7 Restrição de normalização do vetor

A condição:

```text
U^mu U_mu = -1
```

é imposta por multiplicador de Lagrange `lambda_U`.

No limite galáctico quase-estático:

```text
U dot q ≈ 0
lambda_U ≈ 0
```

Logo, a restrição não introduz fonte espúria na dinâmica galáctica.

### 2.8 Tensor energia-momento projetado

O setor projetado já possui tensor energia-momento explícito da forma efetiva:

```text
T_proj = termo isotrópico proporcional a H(Y_s)
       + termo anisotrópico proporcional a H_Y q_mu q_nu
```

Nos testes quasi-estáticos, a fonte projetada reproduz a fonte gravitacional efetiva esperada e os termos locais de pressão/isotropia/anisotropia aparecem subdominantes no regime auditado.

Isso é um avanço relevante: o setor já não está apenas como equação de força, mas possui objeto tensorial efetivo.

### 2.9 Cosmologia homogênea

A parte projetada espacial não domina automaticamente o fundo FLRW, pois para campo homogêneo:

```text
Y_s = 0
```

Isso é positivo: separa o setor galáctico projetado do setor cosmológico homogêneo.

A evolução de fundo fica carregada pelo setor `theta` ou por uma densidade efetiva do tipo:

```text
Omega_psi(a)
```

com parametrização efetiva testável.

Esse ponto evita que o mesmo termo galáctico gere automaticamente uma cosmologia errada.

### 2.10 Buracos negros e saturação estrutural

O setor forte usa uma lei de saturação estrutural do tipo:

```text
S(x) = 1 - exp[-(x/x_c)^n]
```

com métrica efetiva:

```text
A(x) = 1 - S(x)/x
```

Para `x >> x_c`, recupera-se o exterior de Schwarzschild.

Os benchmarks preservados são:

```text
x_h = 1
x_ph = 1.5
b_shadow / r_s = 2.598076
x_ISCO ≈ 3
```

A condição conservadora de regularidade direta do núcleo é:

```text
n >= 6
```

Assim, a TVGD possui uma ponte forte entre saturação estrutural local de buracos negros e saturação universal do vácuo, ainda que o acoplamento final desse setor à ação completa permaneça em desenvolvimento.

---

## 3. O que mudou em relação às limitações antigas

A limitação antiga dizia, aproximadamente:

```text
derivação final que conecte regime galáctico, cosmológico e buracos negros
dentro da mesma ação efetiva.
```

Essa frase ficou fraca demais para o estágio atual.

A versão corrigida deve ser:

```text
A TVGD já possui uma Camada II covariante efetiva que conecta, em nível
variacional e de limites, o regime galáctico, o setor cosmológico homogêneo,
o setor vetorial/PPN e o regime forte por saturação. O que ainda falta é fechar
a teoria fundamental completa: equações acopladas finais, conservação total,
perturbações cosmológicas, estabilidade Hamiltoniana global e origem rigorosa
de a0 e dos setores de saturação.
```

---

## 4. Limitações reais que permanecem

### 4.1 Equações acopladas completas

Ainda falta derivar e auditar integralmente o sistema acoplado:

```text
g_{mu nu} + psi + U^mu + theta + lambda_U + matéria
```

A ação existe como Camada II efetiva, mas a teoria final precisa mostrar explicitamente todas as equações de campo e seus limites.

### 4.2 Conservação covariante total

É necessário verificar explicitamente:

```text
nabla_mu T_total^{mu nu} = 0
```

incluindo:

```text
T_m
T_proj
T_theta
T_U
T_lambda
```

Essa checagem é essencial para garantir consistência relativística completa.

### 4.3 Potenciais fracos Phi e Psi

Ainda falta derivar completamente os potenciais métricos fracos:

```text
Phi
Psi
```

e a combinação de lenteamento:

```text
Phi + Psi
```

bem como o slip gravitacional:

```text
Phi - Psi
```

A TVGD já tem estimativas e fonte anisotrópica controlada, mas falta a derivação relativística final e o confronto sistemático com lentes fracas e fortes.

### 4.4 Perturbações cosmológicas completas

O setor FLRW de fundo é promissor, mas ainda falta implementar o sistema completo em código cosmológico do tipo CLASS/CAMB com:

```text
delta_psi
delta_theta
delta_U
perturbações métricas
crescimento de estruturas
lenteamento CMB
```

A compatibilidade com CMB de fundo e RSD diagonal é forte preliminarmente, mas a análise final exige likelihoods oficiais.

### 4.5 Estabilidade além do regime quase-estático

A estabilidade escalar, vetorial e PPN foi auditada em ramos efetivos.

Ainda falta uma análise global de:

```text
Hamiltoniano quadrático completo
ghosts no sistema acoplado
instabilidades de gradiente fora do limite quasi-estático
modos mistos psi-U-g
estabilidade cosmológica em alto redshift
```

Esse é um dos pontos mais importantes antes de tratar a covariante como final.

### 4.6 Origem rigorosa de a0

O documento covariante fortalece a interpretação:

```text
a0 ≈ c H0 / (2 pi)
```

e mostra que o valor sai na ordem correta de grandeza.

Mas ainda falta derivar rigorosamente o fator `2 pi` a partir de:

```text
modo fundamental do campo estrutural
termodinâmica de horizonte
condição global de relaxamento
ou princípio de saturação do vácuo
```

Portanto, `a0` não deve mais ser tratado como parâmetro arbitrário puro, mas ainda não está fundamentalmente deduzido de forma final.

### 4.7 Setor theta

O setor `theta` carrega a evolução homogênea cosmológica e pode representar saturação/relaxamento temporal.

Ainda falta derivar de forma final:

```text
V_theta(theta)
Z_theta(theta)
Omega_psi(a)
w(a)
```

a partir de princípio estrutural mais profundo.

### 4.8 Saturação universal

A lei exponencial de saturação está conceitualmente muito forte:

```text
S(Theta) = 1 - exp(-Theta)
```

e gera reduções compatíveis com:

```text
P(u) no regime galáctico
Q_psi(x) no regime de buracos negros
Xi(z) no setor temporal/cosmológico
```

Mas ainda falta formalizar a variável covariante universal `Theta` em todos os regimes.

### 4.9 Comparação Bayesiana e likelihoods

Ainda faltam análises formais com:

```text
likelihoods oficiais
matrizes de covariância completas
priors declarados
Bayesian evidence
comparação com Lambda-CDM
comparação com MOND/RMOND/TeVeS/AQUAL
comparação com halos NFW
```

O repositório público atual é uma base de validação e auditoria, não substitui uma análise estatística cosmológica final.

### 4.10 Reprodutibilidade independente

Os scripts públicos são um passo essencial, mas ainda falta:

```text
reprodução por terceiros
ambientes congelados
testes automatizados
dados versionados
notebooks com checksums
relatórios de saída esperada
```

Esse ponto é importante para responder críticas de overfitting.

---

## 5. Setores que devem ser tratados como fortes

Os seguintes blocos podem ser apresentados como fortes no estágio efetivo:

```text
1. Recuperação do limite galáctico g_TVGD = g_bar / P(u)
2. Limite profundo g_TVGD ≈ sqrt(g_bar * a0)
3. Screening local para u >> 1
4. Estabilidade escalar efetiva: P > 0, D > 0, 0 < c_s^2 <= 1
5. Setor vetorial com c_T = 1 no ramo c1 + c3 = 0
6. Restrição U^mu U_mu = -1 sem fonte espúria galáctica
7. Tensor energia-momento projetado explícito
8. Separação entre setor projetado galáctico e fundo FLRW
9. Recuperação dos benchmarks externos de buracos negros
10. Estrutura de ação covariante de Camada II
```

---

## 6. Setores que devem ser tratados como promissores, mas não fechados

Os seguintes blocos devem ser descritos com cautela:

```text
1. Origem fundamental de a0
2. Fator 2 pi em a0 ≈ cH0/(2 pi)
3. Potenciais Phi e Psi completos
4. Slip gravitacional em lentes fracas externas
5. Perturbações cosmológicas completas
6. Hamiltoniano quadrático do sistema acoplado
7. Setor theta fundamental
8. Variável universal de saturação Theta
9. Conexão final entre saturação primordial e buracos negros
10. Bayesian evidence contra Lambda-CDM e MOND
```

---

## 7. Linguagem científica recomendada

Use:

```text
A TVGD possui uma formulação covariante efetiva de Camada II.
```

Use:

```text
A covariante já recupera os limites galáctico, local, cosmológico de fundo
e forte efetivo, mas ainda não é uma teoria fundamental final.
```

Use:

```text
O setor vetorial possui ramo compatível com c_T = 1 e PPN pequeno.
```

Use:

```text
A função P(u) é derivável de uma estrutura de ação efetiva projetada/cinética,
não apenas de um ajuste fenomenológico.
```

Evite:

```text
A covariante está totalmente fechada.
```

Evite:

```text
A TVGD já substitui definitivamente Lambda-CDM.
```

Evite:

```text
A origem de a0 está provada de forma fundamental.
```

Evite:

```text
Todos os modos cosmológicos já foram demonstrados estáveis.
```

---

## 8. Texto curto para substituir no README

Sugestão para a seção de limitações do README:

```text
A TVGD já possui uma arquitetura covariante efetiva avançada, com ação de
Camada II contendo setor estrutural projetado, vetor temporal unitário, setor
tipo Einstein-Aether, restrição de normalização, setor cosmológico homogêneo
e recuperação explícita do limite galáctico. A limitação atual não é ausência
de covariante, mas o fechamento fundamental: derivar todas as equações
acopladas, verificar conservação covariante total, calcular Phi e Psi,
implementar perturbações cosmológicas completas, provar estabilidade global
do sistema acoplado e deduzir rigorosamente a origem de a0 e da variável
universal de saturação.
```

---

## 9. Classificação atual

A classificação honesta após o documento covariante é:

```text
TVGD fenomenológica simples: superada.
TVGD covariante efetiva: forte e avançada.
TVGD fundamental final: ainda não fechada.
```

Ou, em uma frase:

```text
A TVGD já possui uma ação covariante efetiva testável; o que falta é o
fechamento fundamental e a validação estatística completa.
```

---

## 10. Conclusão

O documento `covariante.pdf` muda a forma correta de apresentar as limitações.

Antes, parecia que a TVGD ainda precisava encontrar uma covariante. Agora, o status correto é diferente:

```text
A TVGD já tem uma Camada II covariante efetiva matematicamente explícita,
com recuperação dos limites físicos principais e auditorias de estabilidade.
```

O que permanece em aberto é mais alto nível:

```text
fechar a ação fundamental única,
derivar o sistema acoplado completo,
provar conservação e estabilidade global,
implementar perturbações cosmológicas completas,
deduzir a origem rigorosa de a0,
e validar tudo com likelihoods independentes.
```

Esse é um status muito mais forte e mais justo para a teoria.
