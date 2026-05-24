# TVGD — Teoria do Vácuo Gravitacional Dispersivo

**Autor:** Marcelo Bauman  
**Status:** teoria efetiva covariante em desenvolvimento, com pacote público inicial de validação numérica, reprodutibilidade e falsificabilidade.

---

## 1. Visão geral

A **Teoria do Vácuo Gravitacional Dispersivo (TVGD)** é uma proposta efetiva para descrever fenômenos gravitacionais em múltiplas escalas sem assumir, como ponto de partida, que matéria escura e energia escura sejam componentes fundamentais.

A ideia central é que o vácuo gravitacional possui uma resposta estrutural dependente do regime físico. Em regiões de alta aceleração ou alta curvatura, a teoria deve recuperar aproximadamente a gravidade usual. Em regimes de baixa aceleração, grandes escalas cosmológicas ou núcleos compactos, essa resposta pode modificar a relação efetiva entre matéria bariônica, geometria e observáveis.

O estado atual da TVGD é mais forte do que uma fenomenologia isolada: a teoria já possui uma **Camada II covariante efetiva**, com ação estruturada, campo escalar estrutural, vetor temporal unitário, setor vetorial compatível com ondas gravitacionais, mecanismo de screening local e extensão efetiva para buracos negros por saturação estrutural.

Este repositório disponibiliza códigos públicos para auditoria, reprodução, comparação e crítica independente.

---

## 2. Equação fenomenológica central

No regime galáctico, a resposta principal usada nos testes públicos é:

**P(u) = 1 − exp(−u)**

com:

**u = sqrt(g_bar / a0)**

A aceleração efetiva é escrita como:

**g_TVGD = g_bar / P(u)**

No limite de baixa aceleração, quando **u << 1**, temos:

**P(u) ≈ u**

Portanto:

**g_TVGD ≈ sqrt(g_bar · a0)**

Esse limite reproduz uma relação de aceleração radial profunda do tipo MOND-like, mas dentro da interpretação da TVGD como resposta estrutural do vácuo.

---

## 3. Interpretação física

Na TVGD, a função **P(u)** representa a fração efetiva de resposta gravitacional estrutural do vácuo.

O fator de amplificação efetiva é:

**alpha(u) = 1 / P(u)**

De maneira qualitativa:

- em regimes solares e laboratoriais, **u >> 1**, então **P(u) → 1**;
- em galáxias e regiões de baixa aceleração, **P(u) < 1**, produzindo amplificação efetiva;
- em buracos negros, a teoria usa saturação estrutural do vácuo para preservar o exterior e regularizar o núcleo efetivo;
- em cosmologia, a resposta estrutural pode aparecer como uma transição temporal efetiva;
- na formulação covariante efetiva, a TVGD conecta campo escalar, vetor temporal, resposta projetada do vácuo e setores de saturação.

---

## 4. Status covariante

A TVGD não deve mais ser descrita como uma teoria sem covariante. O estágio atual é o de uma **teoria efetiva covariante avançada**.

A Camada II pode ser resumida por uma ação efetiva com setores:

**S_TVGD = S_g + S_proj + S_theta + S_U + S_lambda + S_m**

onde:

| Termo | Papel físico |
|---|---|
| **S_g** | setor gravitacional métrico |
| **S_proj** | resposta estrutural projetada do vácuo |
| **S_theta** | setor cosmológico homogêneo |
| **S_U** | setor vetorial tipo Einstein-Aether |
| **S_lambda** | restrição de normalização do vetor temporal |
| **S_m** | matéria bariônica |

O setor vetorial admite um ramo compatível com velocidade tensorial igual à da luz:

**c1 + c3 = 0 ⇒ c_T = 1**

Isso é importante para compatibilidade com restrições de ondas gravitacionais.

O ponto ainda em aberto não é “ter uma covariante”, mas sim transformar essa Camada II efetiva em uma formulação fundamental final, com equações acopladas completas, conservação covariante total, perturbações cosmológicas e validação estatística formal.

---

## 5. Setores testados

Este repositório contém códigos públicos de validação para os seguintes setores:

| Setor | Objetivo | Status atual |
|---|---|---|
| Galáxias / SPARC | Testar a relação g_bar → g_TVGD | Evidência robusta/promissora |
| LITTLE THINGS | Validação externa em galáxias anãs | Passou forte como teste externo efetivo |
| Buracos negros | Preservar horizonte, fótonsfera e sombra com núcleo regular | Passou forte nos benchmarks efetivos |
| Cosmologia de fundo | Checar compatibilidade com expansão tipo Lambda-CDM | Compatível em testes efetivos |
| Tensão de Hubble | Modelar transição estrutural entre H0 primordial e local | Promissor, ainda em desenvolvimento |
| RSD / fsigma8 | Verificar crescimento de estruturas | Compatível com Lambda-CDM em teste diagonal |
| PPN / covariante efetiva | Checar riscos solares, vetoriais e velocidade das ondas gravitacionais | Passou forte no ramo efetivo testado |
| Falsificabilidade | Definir previsões que podem derrubar a teoria | Checklist público incluído |

---

## 6. Estrutura do repositório

Estrutura atual/recomendada:

```text
TVGD-Exames-Publicos/
├── README.md
├── requirements.txt
├── LICENSE
├── CITATION.cff
├── docs/
│   ├── README.md
│   ├── theory_summary.md
│   ├── validation_map.md
│   ├── known_limitations.md
│   ├── reproducibility_guide.md
│   ├── falsifiability.md
│   └── data_sources.md
├── notebooks/
│   ├── main_validation/
│   ├── black_holes/
│   ├── hubble/
│   └── exploratory/
├── paper/
├── results/
├── scripts/
└── src/
```

Observação: se o navegador traduzir automaticamente a interface do GitHub, use os nomes reais que aparecem na URL do arquivo.

---

## 7. Instalação

Clone o repositório:

```bash
git clone https://github.com/marcelobauman82-netizen/TVGD-Exames-Publicos.git
cd TVGD-Exames-Publicos
```

Crie um ambiente Python, se desejar:

```bash
python -m venv .venv
```

Ative o ambiente.

No Windows:

```bash
.venv\Scripts\activate
```

No Linux/macOS:

```bash
source .venv/bin/activate
```

Instale as dependências:

```bash
pip install -r requirements.txt
```

---

## 8. Teste rápido de reprodutibilidade

Para verificar a estrutura do repositório, dependências mínimas e sintaxe dos scripts principais, rode:

```bash
python notebooks/main_validation/23_reproducibility_quickstart.py
```

Para também executar um subconjunto leve de scripts:

```bash
python notebooks/main_validation/23_reproducibility_quickstart.py --run-basic
```

As saídas serão salvas em:

```text
results/23_reproducibility_quickstart/
```

Esse script é o ponto de entrada recomendado para avaliadores externos.

---

## 9. Painel estatístico dos testes

Para gerar um resumo estatístico dos testes principais, rode:

```bash
python notebooks/main_validation/22_model_comparison_aic_bic_summary.py
```

As saídas serão salvas em:

```text
results/22_model_comparison_aic_bic_summary/
```

Esse painel reúne os principais resultados, setores, status, valores de chi2, chi2_red, AIC/BIC aproximado e notas de interpretação.

---

## 10. Documentação complementar

A pasta [`docs/`](docs/) contém documentação de apoio:

- [`docs/theory_summary.md`](docs/theory_summary.md): resumo teórico da TVGD;
- [`docs/validation_map.md`](docs/validation_map.md): mapa dos testes públicos;
- [`docs/known_limitations.md`](docs/known_limitations.md): limitações conhecidas;
- [`docs/reproducibility_guide.md`](docs/reproducibility_guide.md): guia de reprodutibilidade;
- [`docs/data_sources.md`](docs/data_sources.md): fontes de dados usadas ou previstas;
- [`docs/falsifiability.md`](docs/falsifiability.md): critérios de falsificação.

---

## 11. Dados externos

Alguns testes dependem de bases públicas externas. Entre elas:

- SPARC rotation curves;
- LITTLE THINGS;
- Pantheon+SH0ES;
- compilações de fsigma8;
- dados de lenteamento fraco;
- catálogos cosmológicos e astrofísicos auxiliares.

Nem todos os dados pesados são incluídos diretamente no repositório. Quando necessário, os scripts indicam a fonte ou esperam que o usuário coloque os arquivos na pasta correta.

A separação entre código e dados pesados foi mantida para facilitar a publicação pública do repositório e evitar problemas de tamanho, licença ou redistribuição de catálogos.

---

## 12. Resultados principais resumidos

### 12.1 Regime galáctico

A função central é:

**P(u) = 1 − exp(−u)**

Ela reproduz o limite profundo:

**g_TVGD ≈ sqrt(g_bar · a0)**

Em amostras limpas do SPARC, a TVGD apresenta resíduos centrados próximos de zero e erro médio reduzido após remoção de outliers observacionais.

### 12.2 Validação externa com LITTLE THINGS

Em teste externo efetivo com LITTLE THINGS, a família exponencial próxima de **beta ≈ 1.1** supera o modelo bariônico puro.

Resultado registrado:

```text
Status: PASSOU_FORTE_EXTRISK1B_LITTLE_THINGS_EXTERNO
```

### 12.3 Buracos negros

O setor de buracos negros preserva os benchmarks externos de Schwarzschild:

```text
r_h / r_s ≈ 1
r_ph / r_s ≈ 1.5
b_shadow / r_s ≈ 2.598076
```

A estrutura interna é regularizada por saturação efetiva do vácuo, evitando a singularidade no modelo efetivo testado.

### 12.4 Covariante efetiva e PPN

O setor covariante efetivo possui ramos que preservam:

- velocidade tensorial c_T = 1;
- estabilidade dos modos no ramo auditado;
- limites PPN fortes;
- compatibilidade solar efetiva.

Resultado registrado:

```text
Status: PASSOU_FORTE_COVRISK_PPN
```

### 12.5 Crescimento de estruturas

Em RSD / fsigma8, a TVGD V1.1 fica praticamente indistinguível de Lambda-CDM em teste diagonal efetivo.

Resultado registrado:

```text
Status: PASSOU_FORTE_COVRISK8_RSD_FSIGMA8
```

### 12.6 Tensão de Hubble

A TVGD possui um setor efetivo de transição estrutural capaz de conectar um valor primordial próximo de:

```text
H0 ≈ 67.4
```

a um valor local efetivo maior, próximo de:

```text
H0 ≈ 72 - 73
```

Esse setor é promissor, mas ainda exige validação formal com likelihoods completas, CMB, BAO, supernovas, covariâncias e controle de efeitos sistemáticos.

---

## 13. Falsificabilidade

A TVGD pode ser enfraquecida ou derrubada por observações como:

1. curvas de rotação galácticas limpas que violem sistematicamente a forma **P(u) = 1 − exp(−u)**;
2. lentes gravitacionais que exijam amplificações incompatíveis com **alpha(u) = 1/P(u)**;
3. dados PPN ou ondas gravitacionais que excluam o ramo covariante efetivo saudável;
4. medidas de fsigma8, CMB, BAO ou supernovas que gerem tensão estatística forte contra a transição cosmológica proposta;
5. sombras de buracos negros ou ringdowns que rejeitem a saturação estrutural regularizada;
6. catálogos externos independentes que mostrem que a função **P(u)** falha de forma sistemática em diferentes classes de galáxias.

A teoria é construída para ser testável. O objetivo do repositório é permitir que outros pesquisadores rodem, critiquem e tentem falsificar os resultados.

---

## 14. Limitações atuais

A TVGD está em desenvolvimento avançado. As limitações essenciais são:

- fechar a ação fundamental final a partir da Camada II covariante efetiva;
- derivar explicitamente todas as equações acopladas;
- demonstrar conservação covariante total do tensor energia-momento;
- implementar perturbações cosmológicas completas;
- deduzir de forma rigorosa a origem de a0;
- substituir testes efetivos por likelihoods completas com covariâncias;
- realizar comparação Bayesiana formal contra Lambda-CDM, MOND e modelos com halos;
- obter reprodução independente por terceiros.

Essas limitações não significam ausência de base covariante. Significam que a TVGD já tem uma formulação efetiva testável, mas ainda não deve ser apresentada como teoria fundamental final.

---

## 15. Como rodar os principais scripts

### Quickstart

```bash
python notebooks/main_validation/23_reproducibility_quickstart.py
```

### Quickstart com execução leve

```bash
python notebooks/main_validation/23_reproducibility_quickstart.py --run-basic
```

### Painel estatístico

```bash
python notebooks/main_validation/22_model_comparison_aic_bic_summary.py
```

### Teste de buracos negros

```bash
python notebooks/main_validation/19_black_hole_structural_saturation.py
```

### Checklist de falsificabilidade

```bash
python notebooks/main_validation/20_falsifiability_prediction_checklist.py
```

### RSD / crescimento de estruturas

```bash
python notebooks/main_validation/18_covrisk_rsd_fsigma8.py
```

### PPN / risco covariante

```bash
python notebooks/main_validation/17_covrisk_scalar_vector_ppn.py
```

---

## 16. Como interpretar os resultados

Os resultados devem ser interpretados como validação efetiva inicial, não como prova final.

Uma saída do tipo **PASSOU_FORTE** significa que o teste específico passou nos critérios definidos pelo próprio script. Isso não significa que a teoria esteja definitivamente comprovada.

Uma saída do tipo **PASSOU_MODERADO** significa compatibilidade parcial ou dependente de hipóteses, cortes, amostra ou aproximações.

Uma saída do tipo **FALHOU** indica que o setor, modelo ou hipótese testada não passou no critério definido.

A TVGD deve ser avaliada pela cadeia completa de testes, não por um único resultado isolado.

---

## 17. Sobre matéria escura e energia escura

A TVGD não parte da hipótese de que matéria escura e energia escura sejam componentes fundamentais. Em vez disso, tenta reinterpretar parte dos fenômenos atribuídos a esses setores como efeitos de resposta estrutural do vácuo gravitacional.

Isso não significa ignorar as evidências observacionais. Pelo contrário: curvas de rotação, lentes, CMB, BAO, supernovas, crescimento de estruturas e buracos negros são tratados como restrições obrigatórias.

A pergunta central da TVGD é:

> Os fenômenos atualmente atribuídos aos setores escuros podem emergir, total ou parcialmente, de uma resposta efetiva do vácuo gravitacional?

---

## 18. Filosofia do projeto

Este repositório segue três princípios:

1. **Reprodutibilidade:** os códigos devem ser públicos e auditáveis.
2. **Falsificabilidade:** a teoria deve declarar claramente onde pode falhar.
3. **Cautela científica:** resultados promissores não devem ser apresentados como prova final antes de validação independente.

---

## 19. Como citar

Enquanto o artigo formal estiver em preparação, cite este repositório como:

```text
Bauman, M. TVGD — Teoria do Vácuo Gravitacional Dispersivo:
public validation repository. GitHub repository, 2026.
```

O arquivo [`CITATION.cff`](CITATION.cff) também fornece metadados de citação para o GitHub.

---

## 20. Aviso científico

Este projeto apresenta uma teoria efetiva alternativa/em desenvolvimento.

Os códigos e resultados devem ser tratados como material científico aberto para auditoria, reprodução, crítica e tentativa de falsificação.

A TVGD não é apresentada como teoria final consolidada, mas como um programa de pesquisa quantitativo, falsificável e progressivamente covariante.

---

## 21. Contato

**Autor:** Marcelo Bauman  
**GitHub:** https://github.com/marcelobauman82-netizen  
**Repositório:** https://github.com/marcelobauman82-netizen/TVGD-Exames-Publicos
