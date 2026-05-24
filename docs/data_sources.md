[data_sources.md](https://github.com/user-attachments/files/28193869/data_sources.md)
# TVGD — Fontes de Dados

Este documento lista as principais bases públicas e classes de dados usadas ou previstas nos testes da TVGD.

---

## 1. Objetivo

A TVGD depende de validação em múltiplas escalas. Para isso, os testes usam ou preveem o uso de dados de:

- galáxias;
- curvas de rotação;
- lentes gravitacionais;
- supernovas;
- crescimento de estruturas;
- cosmologia de fundo;
- buracos negros;
- ondas gravitacionais.

Nem todos os dados estão incluídos diretamente no repositório. Em muitos casos, o ideal é baixar os catálogos a partir das fontes oficiais.

---

## 2. SPARC

**Uso:** curvas de rotação galácticas e relação de aceleração radial.

**Setor:** galáxias / regime de baixa aceleração.

**Objetivo na TVGD:** testar a função

```math
P(u)=1-e^{-u}
```

e a previsão

```math
g_{\rm TVGD}=\frac{g_{\rm bar}}{P(u)}.
```

**Observação:** SPARC é central para o setor galáctico e deve ser tratado com controle cuidadoso de outliers, incertezas e massa-luminosidade estelar.

---

## 3. LITTLE THINGS

**Uso:** validação externa em galáxias anãs.

**Setor:** galáxias / teste externo.

**Objetivo na TVGD:** verificar se a resposta estrutural permanece competitiva fora da amostra principal do SPARC.

**Status:** usado como teste externo efetivo no pacote público.

---

## 4. Pantheon+SH0ES

**Uso:** supernovas, expansão local, tensão de Hubble e testes de transição estrutural.

**Setor:** cosmologia / \(H_0\).

**Objetivo na TVGD:** investigar se uma transição estrutural pode preservar o regime primordial próximo de \(H_0 \approx 67.4\) enquanto permite valores locais maiores, próximos de \(H_0 \approx 72-73\).

**Cuidado:** testes formais devem usar matriz de covariância completa e controle de calibradores, fluxo local, peculiar velocities e seleção de amostra.

---

## 5. RSD / \(f\sigma_8\)

**Uso:** crescimento de estruturas.

**Setor:** cosmologia perturbativa.

**Objetivo na TVGD:** verificar se a teoria preserva o crescimento observado de estruturas.

**Status:** em teste diagonal efetivo, a TVGD V1.1 ficou praticamente indistinguível de ΛCDM.

**Cuidado:** para publicação, o ideal é usar likelihoods oficiais e correções fiduciais.

---

## 6. BAO

**Uso:** distâncias cosmológicas e escala acústica.

**Setor:** cosmologia de fundo.

**Objetivo na TVGD:** verificar se a expansão efetiva preserva distâncias compatíveis com ΛCDM.

**Cuidado:** BAO formal exige tratamento conjunto com CMB, SN e matriz de covariância.

---

## 7. CMB

**Uso:** geometria primordial, espectro de anisotropias, lenteamento CMB e escala acústica.

**Setor:** cosmologia primordial.

**Objetivo na TVGD:** garantir que modificações tardias não destruam o sucesso do regime primordial.

**Ferramentas previstas:** CLASS, CAMB, Cobaya e likelihoods Planck/ACT/SPT quando disponíveis.

---

## 8. HSC e lenteamento fraco

**Uso:** comparação entre massa bariônica projetada, resposta TVGD e sinal de lenteamento fraco.

**Setor:** lentes gravitacionais.

**Objetivo na TVGD:** testar se a amplificação estrutural

```math
\alpha(u)=\frac{1}{P(u)}
```

é compatível com sinais de lenteamento.

---

## 9. COSMOS-Web / JWST

**Uso:** mapas de massa projetada, catálogos fotométricos e possíveis comparações com mapas de convergência.

**Setor:** lentes / estrutura em grande escala.

**Status:** setor exploratório. Deve ser tratado como previsão e preparação para testes futuros, não como evidência principal enquanto os mapas e máscaras não estiverem formalmente controlados.

---

## 10. Buracos negros e EHT/LIGO

**Uso:** sombra, fótonsfera, horizonte, ringdown e compatibilidade com observações fortes.

**Setor:** campo forte.

**Objetivo na TVGD:** preservar os benchmarks externos de Schwarzschild/Kerr enquanto regulariza a estrutura interna por saturação efetiva.

Benchmarks principais:

```math
r_h/r_s \approx 1,
```

```math
r_{\rm ph}/r_s \approx 1.5,
```

```math
b_{\rm shadow}/r_s \approx 2.598076.
```

---

## 11. Política de dados do repositório

Este repositório prioriza código, documentação e reprodutibilidade.

Dados muito grandes ou com licenças específicas não devem ser redistribuídos sem permissão. O ideal é que cada script indique:

- nome da base;
- fonte oficial;
- arquivo esperado;
- pasta esperada;
- instrução de download;
- referência científica.

---

## 12. Estrutura futura sugerida

```text
data/
├── raw/
├── processed/
├── external/
└── README.md
```

A pasta `data/` pode ser mantida fora do GitHub se os arquivos forem grandes.

---

## 13. Resumo

As fontes de dados da TVGD cobrem regimes independentes. Isso é essencial para evitar que a teoria seja ajustada apenas a um tipo de observação.

A força do programa está em sobreviver simultaneamente a dados de galáxias, cosmologia, crescimento de estruturas, lentes, buracos negros e limites locais.
