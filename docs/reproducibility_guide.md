[reproducibility_guide.md](https://github.com/user-attachments/files/28193859/reproducibility_guide.md)
# TVGD — Guia de Reprodutibilidade

Este documento explica como um avaliador externo pode baixar o repositório, instalar dependências, verificar a estrutura e rodar os primeiros testes públicos da **Teoria do Vácuo Gravitacional Dispersivo (TVGD)**.

---

## 1. Clonar o repositório

```bash
git clone https://github.com/marcelobauman82-netizen/TVGD-Exames-Publicos.git
cd TVGD-Exames-Publicos
```

Se o nome do repositório mudar, substitua a URL acima pela URL real mostrada pelo GitHub.

---

## 2. Criar ambiente Python

No Windows:

```bash
python -m venv .venv
.venv\Scripts\activate
```

No Linux/macOS:

```bash
python -m venv .venv
source .venv/bin/activate
```

---

## 3. Instalar dependências

```bash
pip install -r requirements.txt
```

Dependências mínimas esperadas:

```text
numpy
pandas
matplotlib
scipy
astropy
requests
tqdm
```

Alguns testes avançados podem exigir pacotes adicionais, dependendo dos dados externos ou do backend cosmológico usado.

---

## 4. Rodar o quickstart

O primeiro teste recomendado é:

```bash
python cadernos/main_validation/23_reproducibility_quickstart.py
```

ou, se sua pasta estiver em inglês:

```bash
python notebooks/main_validation/23_reproducibility_quickstart.py
```

Esse script verifica:

- estrutura do repositório;
- dependências mínimas;
- presença dos scripts principais;
- sintaxe dos scripts de validação;
- possibilidade de rodar testes leves, se ativado.

---

## 5. Rodar quickstart com execução leve

```bash
python cadernos/main_validation/23_reproducibility_quickstart.py --run-basic
```

ou:

```bash
python notebooks/main_validation/23_reproducibility_quickstart.py --run-basic
```

Essa opção executa um subconjunto leve de scripts, sem baixar dados pesados.

---

## 6. Gerar painel estatístico

```bash
python cadernos/main_validation/22_model_comparison_aic_bic_summary.py
```

ou:

```bash
python notebooks/main_validation/22_model_comparison_aic_bic_summary.py
```

Saída esperada:

```text
Resultados/22_model_comparison_aic_bic_summary/
```

ou:

```text
results/22_model_comparison_aic_bic_summary/
```

Arquivos esperados:

```text
tvgd_model_comparison_summary.csv
tvgd_model_comparison_summary.json
tvgd_model_comparison_report.md
```

---

## 7. Ordem recomendada de execução

### Etapa 1 — verificação do repositório

```bash
python cadernos/main_validation/23_reproducibility_quickstart.py
python cadernos/main_validation/22_model_comparison_aic_bic_summary.py
```

### Etapa 2 — testes leves

```bash
python cadernos/main_validation/00_baseline_consistency_macro_micro.py
python cadernos/main_validation/02_black_hole_regular_core.py
python cadernos/main_validation/05_stability_causality_screening.py
python cadernos/main_validation/06_nonphantom_cosmology_background.py
python cadernos/main_validation/19_black_hole_structural_saturation.py
python cadernos/main_validation/20_falsifiability_prediction_checklist.py
```

### Etapa 3 — testes com dados externos

```bash
python cadernos/main_validation/01_sparc_rar_tvgd.py
python cadernos/main_validation/01B_sparc_ml_robustness.py
python cadernos/main_validation/01C_sparc_rar_comparison.py
python cadernos/main_validation/15_extrisk_little_things_external.py
python cadernos/main_validation/18_covrisk_rsd_fsigma8.py
```

### Etapa 4 — cosmologia e tensão de Hubble

```bash
python cadernos/main_validation/03_cosmology_background_bao.py
python cadernos/main_validation/07B_hubble_structural_phase_transition.py
python cadernos/main_validation/07E_pantheon_full_covariance_hubble.py
python cadernos/main_validation/07G_pantheon_noncalibrator_shape.py
python cadernos/main_validation/07H_hubble_structural_phase_robustness.py
```

---

## 8. Observação sobre dados externos

Alguns testes dependem de dados públicos externos, como SPARC, LITTLE THINGS, Pantheon+SH0ES, compilações de \(f\sigma_8\) e catálogos de lenteamento.

Esses dados podem não estar incluídos no repositório por tamanho, licença ou praticidade. Quando necessário, cada script deve indicar a origem dos dados ou a pasta esperada.

---

## 9. Problemas comuns

### Dependência ausente

```bash
pip install -r requirements.txt
```

### Script não encontrado

Confira se você está na raiz do repositório.

### Arquivo de dados ausente

Verifique o cabeçalho do script. Alguns testes exigem download manual dos dados.

### Git não subiu um script

Verifique o `.gitignore`. Ele não deve listar os scripts de validação pública.

---

## 10. Interpretação

Este guia não prova a TVGD. Ele apenas garante que o pacote público possa ser verificado, instalado e executado por terceiros.

A validação científica depende da reprodução dos resultados, da análise estatística e da comparação com modelos alternativos.
