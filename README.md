[README_corrigido.md](https://github.com/user-attachments/files/28194645/README_corrigido.md)
# TVGD — Teoria do Vácuo Gravitacional Dispersivo

**Autor:** Marcelo Bauman  
**Status:** teoria efetiva em desenvolvimento, com pacote público inicial de validação numérica, reprodutibilidade e falsificabilidade.

---

## 1. Visão geral

A **Teoria do Vácuo Gravitacional Dispersivo (TVGD)** é uma proposta efetiva para descrever fenômenos gravitacionais em múltiplas escalas sem assumir, como ponto de partida, que matéria escura e energia escura sejam componentes fundamentais.

A ideia central é que o vácuo gravitacional possui uma resposta estrutural dependente do regime físico. Em regiões de alta aceleração ou alta curvatura, a teoria deve recuperar aproximadamente a gravidade usual. Em regimes de baixa aceleração, grandes escalas cosmológicas ou núcleos compactos, a resposta estrutural pode modificar a relação efetiva entre matéria bariônica, geometria e observáveis.

A TVGD é apresentada como uma teoria efetiva, falsificável e progressivamente covariantizada. O objetivo deste repositório é disponibilizar códigos públicos para auditoria, reprodução e crítica independente.

---

## 2. Equação fenomenológica central

No regime galáctico, a resposta principal usada nos testes públicos é

```math
P(u) = 1 - e^{-u}
```

com

```math
u = \sqrt{\frac{g_{\rm bar}}{a_0}}.
```

A aceleração efetiva é escrita como

```math
g_{\rm TVGD} = \frac{g_{\rm bar}}{P(u)}.
```

No limite de baixa aceleração, quando \(u \ll 1\),

```math
P(u) \approx u,
```

então

```math
g_{\rm TVGD}
\approx
\frac{g_{\rm bar}}{\sqrt{g_{\rm bar}/a_0}}
=
\sqrt{g_{\rm bar}a_0}.
```

Esse limite reproduz uma relação de aceleração radial profunda do tipo MOND-like, mas dentro da interpretação da TVGD como resposta estrutural do vácuo.

---

## 3. Interpretação física

Na TVGD, a função \(P(u)\) representa a fração efetiva de resposta gravitacional estrutural do vácuo.

O fator

```math
\alpha(u) = \frac{1}{P(u)}
```

funciona como uma amplificação efetiva da gravidade bariônica em regimes onde \(P(u)<1\).

De maneira qualitativa:

- em regimes solares e laboratoriais, \(u\gg1\), então \(P(u)\to1\);
- em galáxias e regiões de baixa aceleração, \(P(u)<1\), produzindo amplificação efetiva;
- em buracos negros, a teoria usa uma saturação estrutural regularizada;
- em cosmologia, a resposta estrutural pode aparecer como uma transição temporal efetiva;
- em uma formulação covariante efetiva, a TVGD busca conectar campos escalares, vetoriais e respostas estruturais sem violar limites solares, PPN e ondas gravitacionais.

---

## 4. Setores testados

Este repositório contém códigos públicos de validação para os seguintes setores:

| Setor | Objetivo | Status atual |
|---|---|---|
| Galáxias / SPARC | Testar a relação \(g_{\rm bar}\rightarrow g_{\rm TVGD}\) | Evidência robusta/promissora |
| LITTLE THINGS | Validação externa em galáxias anãs | Passou forte como teste externo efetivo |
| Buracos negros | Preservar horizonte, fótonsfera e sombra de Schwarzschild com núcleo regular | Passou forte nos benchmarks efetivos |
| Cosmologia de fundo | Checar compatibilidade com expansão tipo ΛCDM | Compatível em testes efetivos |
| Tensão de Hubble | Modelar transição estrutural entre \(H_0\) primordial e local | Promissor, ainda em desenvolvimento |
| RSD / \(f\sigma_8\) | Verificar crescimento de estruturas | Compatível com ΛCDM em teste diagonal |
| PPN / covariante efetiva | Checar riscos solares, vetoriais e velocidade das ondas gravitacionais | Passou forte no ramo efetivo testado |
| Falsificabilidade | Definir previsões que podem derrubar a teoria | Checklist público incluído |

---

## 5. Estrutura do repositório

Estrutura atual/recomendada:

```text
TVGD-Exames-Publicos/
├── README.md
├── requirements.txt
├── LICENSE
├── CITATION.cff
├── Docs/
│   ├── README.md
│   ├── theory_summary.md
│   ├── validation_map.md
│   ├── known_limitations.md
│   ├── reproducibility_guide.md
│   ├── falsifiability.md
│   └── data_sources.md
├── notebooks/
│   └── main_validation/
│       ├── 00_baseline_consistency_macro_micro.py
│       ├── 01_sparc_rar_tvgd.py
│       ├── 01B_sparc_ml_robustness.py
│       ├── 01C_sparc_rar_comparison.py
│       ├── 02_black_hole_regular_core.py
│       ├── 03_cosmology_background_bao.py
│       ├── 04B_controlled_strong_lensing_selected.py
│       ├── 05_stability_causality_screening.py
│       ├── 06_nonphantom_cosmology_background.py
│       ├── 07B_hubble_structural_phase_transition.py
│       ├── 07E_pantheon_full_covariance_hubble.py
│       ├── 07G_pantheon_noncalibrator_shape.py
│       ├── 07H_hubble_structural_phase_robustness.py
│       ├── 08_black_hole_structural_saturation.py
│       ├── 08B_rotating_black_hole_structural_saturation.py
│       ├── 09_falsifiability_prediction_checklist.py
│       ├── 10_euclid_forecast_structural_response.py
│       ├── 11_rmax_sparc_universality.py
│       ├── 12_rmax_transition_beta_robustness.py
│       ├── 15_extrisk_little_things_external.py
│       ├── 17_covrisk_scalar_vector_ppn.py
│       ├── 18_covrisk_rsd_fsigma8.py
│       ├── 19_black_hole_structural_saturation.py
│       ├── 20_falsifiability_prediction_checklist.py
│       ├── 22_model_comparison_aic_bic_summary.py
│       └── 23_reproducibility_quickstart.py
├── paper/
├── results/
├── scripts/
└── src/
```

> Observação: se o navegador traduzir automaticamente os nomes das pastas na interface do GitHub, use os nomes reais que aparecem na URL do arquivo. Os comandos abaixo assumem a estrutura padrão em inglês: `notebooks/`, `results/`, `paper/`, `scripts/` e `Docs/`.

---

## 6. Instalação

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

## 7. Teste rápido de reprodutibilidade

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

## 8. Painel estatístico dos testes

Para gerar um resumo estatístico dos testes principais, rode:

```bash
python notebooks/main_validation/22_model_comparison_aic_bic_summary.py
```

As saídas serão salvas em:

```text
results/22_model_comparison_aic_bic_summary/
```

Esse script gera:

```text
tvgd_model_comparison_summary.csv
tvgd_model_comparison_summary.json
tvgd_model_comparison_report.md
```

O objetivo desse painel é reunir, em um só lugar, os principais resultados, setores, status, valores de \(\chi^2\), \(\chi^2_{\rm red}\), AIC/BIC aproximado e notas de interpretação.

---

## 9. Documentação complementar

A pasta [`Docs/`](Docs/) contém documentação de apoio:

- [`Docs/theory_summary.md`](Docs/theory_summary.md): resumo teórico da TVGD;
- [`Docs/validation_map.md`](Docs/validation_map.md): mapa dos testes públicos;
- [`Docs/known_limitations.md`](Docs/known_limitations.md): limitações conhecidas;
- [`Docs/reproducibility_guide.md`](Docs/reproducibility_guide.md): guia de reprodutibilidade;
- [`Docs/data_sources.md`](Docs/data_sources.md): fontes de dados usadas ou previstas;
- [`Docs/falsifiability.md`](Docs/falsifiability.md): critérios de falsificação.

---

## 10. Dados externos

Alguns testes dependem de bases públicas externas. Entre elas:

- SPARC rotation curves;
- LITTLE THINGS;
- Pantheon+SH0ES;
- compilações de \(f\sigma_8\);
- dados de lenteamento fraco;
- catálogos cosmológicos e astrofísicos auxiliares.

Nem todos os dados pesados são incluídos diretamente no repositório. Quando necessário, os scripts indicam a fonte ou esperam que o usuário coloque os arquivos na pasta correta.

A separação entre código e dados pesados foi mantida para facilitar a publicação pública do repositório e evitar problemas de tamanho, licença ou redistribuição de catálogos.

---

## 11. Resultados principais resumidos

Alguns resultados consolidados nos testes internos e nos scripts públicos:

### 11.1 Regime galáctico

A função

```math
P(u)=1-e^{-u}
```

reproduz o limite profundo

```math
g_{\rm TVGD}\approx\sqrt{g_{\rm bar}a_0}.
```

Em amostras limpas do SPARC, a TVGD apresenta resíduos centrados próximos de zero e erro médio reduzido após remoção de outliers observacionais.

### 11.2 Validação externa com LITTLE THINGS

Em teste externo efetivo com LITTLE THINGS, a família exponencial próxima de \(\beta\approx1.1\) supera o modelo bariônico puro.

Resultado registrado:

```text
Status: PASSOU_FORTE_EXTRISK1B_LITTLE_THINGS_EXTERNO
```

### 11.3 Buracos negros

O setor de buracos negros preserva os benchmarks externos de Schwarzschild:

```math
r_h/r_s \approx 1,
```

```math
r_{\rm ph}/r_s \approx 1.5,
```

```math
b_{\rm shadow}/r_s \approx 2.598076.
```

A estrutura interna é regularizada por uma saturação efetiva do vácuo, evitando a singularidade no modelo efetivo testado.

### 11.4 Covariante efetiva e PPN

O setor covariante efetivo possui ramos que preservam:

- velocidade tensorial \(c_T=1\);
- estabilidade dos modos;
- limites PPN fortes;
- compatibilidade solar efetiva.

Resultado registrado:

```text
Status: PASSOU_FORTE_COVRISK_PPN
```

### 11.5 Crescimento de estruturas

Em RSD/\(f\sigma_8\), a TVGD V1.1 fica praticamente indistinguível de ΛCDM em teste diagonal efetivo.

Resultado registrado:

```text
Status: PASSOU_FORTE_COVRISK8_RSD_FSIGMA8
```

### 11.6 Tensão de Hubble

A TVGD possui um setor efetivo de transição estrutural capaz de conectar um valor primordial próximo de

```math
H_0 \approx 67.4
```

a um valor local efetivo maior, próximo de

```math
H_0 \approx 72 - 73.
```

Esse setor é promissor, mas ainda exige validação formal com likelihoods completas, CMB, BAO, supernovas, covariâncias e controle de efeitos sistemáticos.

---

## 12. Falsificabilidade

A TVGD pode ser enfraquecida ou derrubada por observações como:

1. curvas de rotação galácticas limpas que violem sistematicamente a forma \(P(u)=1-e^{-u}\);
2. lentes gravitacionais que exijam amplificações incompatíveis com \(\alpha=1/P(u)\);
3. dados PPN ou ondas gravitacionais que excluam o ramo covariante efetivo saudável;
4. medidas de \(f\sigma_8\), CMB, BAO ou supernovas que gerem tensão estatística forte contra a transição cosmológica proposta;
5. sombras de buracos negros ou ringdowns que rejeitem a saturação estrutural regularizada;
6. catálogos externos independentes que mostrem que a função \(P(u)\) falha de forma sistemática em diferentes classes de galáxias.

A teoria é construída para ser testável. O objetivo do repositório é permitir que outros pesquisadores rodem, critiquem e tentem falsificar os resultados.

---

## 13. Limitações atuais

A TVGD ainda está em desenvolvimento. Entre os pontos em aberto:

- derivação covariante completa a partir de uma ação final única;
- conexão formal definitiva entre \(a_0\), campos estruturais e constantes fundamentais;
- likelihoods completas com matrizes de covariância para todos os setores;
- comparação Bayesiana formal contra ΛCDM, MOND e modelos com halos de matéria escura;
- validação independente com pipelines externos;
- expansão dos testes para Euclid, Rubin/LSST, JWST, SKA e catálogos futuros;
- derivação final que conecte de forma rigorosa o regime galáctico, cosmológico e de buracos negros dentro da mesma ação efetiva.

---

## 14. Como rodar os principais scripts

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

## 15. Como interpretar os resultados

Os resultados devem ser interpretados como validação efetiva inicial, não como prova final.

Uma saída do tipo

```text
PASSOU_FORTE
```

significa que o teste específico passou nos critérios definidos pelo próprio script. Isso não significa que a teoria esteja definitivamente comprovada.

Uma saída do tipo

```text
PASSOU_MODERADO
```

significa compatibilidade parcial ou dependente de hipóteses, cortes, amostra ou aproximações.

Uma saída do tipo

```text
FALHOU
```

indica que o setor, modelo ou hipótese testada não passou no critério definido.

A TVGD deve ser avaliada pela cadeia completa de testes, não por um único resultado isolado.

---

## 16. Sobre matéria escura e energia escura

A TVGD não parte da hipótese de que matéria escura e energia escura sejam componentes fundamentais. Em vez disso, tenta reinterpretar parte dos fenômenos atribuídos a esses setores como efeitos de resposta estrutural do vácuo gravitacional.

Isso não significa ignorar as evidências observacionais. Pelo contrário: curvas de rotação, lentes, CMB, BAO, supernovas, crescimento de estruturas e buracos negros são tratados como restrições obrigatórias.

A pergunta central da TVGD é:

> os fenômenos atualmente atribuídos aos setores escuros podem emergir, total ou parcialmente, de uma resposta efetiva do vácuo gravitacional?

---

## 17. Filosofia do projeto

Este repositório segue três princípios:

1. **Reprodutibilidade:** os códigos devem ser públicos e auditáveis.
2. **Falsificabilidade:** a teoria deve declarar claramente onde pode falhar.
3. **Cautela científica:** resultados promissores não devem ser apresentados como prova final antes de validação independente.

---

## 18. Como citar

Enquanto o artigo formal estiver em preparação, cite este repositório como:

```text
Bauman, M. TVGD — Teoria do Vácuo Gravitacional Dispersivo:
public validation repository. GitHub repository, 2026.
```

O arquivo [`CITATION.cff`](CITATION.cff) também fornece metadados de citação para o GitHub.

---

## 19. Aviso científico

Este projeto apresenta uma teoria efetiva alternativa/em desenvolvimento.

Os códigos e resultados devem ser tratados como material científico aberto para auditoria, reprodução, crítica e tentativa de falsificação.

A TVGD não é apresentada como teoria final consolidada, mas como um programa de pesquisa quantitativo, falsificável e progressivamente covariante.

---

## 20. Contato

**Autor:** Marcelo Bauman  
**GitHub:** https://github.com/marcelobauman82-netizen  
**Repositório:** https://github.com/marcelobauman82-netizen/TVGD-Exames-Publicos
