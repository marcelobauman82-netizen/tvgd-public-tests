# Instruções simples — Pasta dupla_fenda

Esta pasta contém 5 testes em Python.

Eles servem para testar um protocolo de memória de fase, coerência e resposta angular usando arquivos `.s2p`.

## Aviso importante

Os testes atuais usam dados DEMO/SIMULADOS.

Eles servem para verificar se o programa está funcionando e se a ordem dos testes está correta.

Eles não são dados reais de laboratório.

## Ordem correta

Rode os testes nesta ordem:

```bash
python teste_01_gerar_arquivos_demo_memoria.py
python teste_02_static_memory.py
python teste_03_robustez_temporal_memoria_da_fase.py
python teste_04_resposta_angular_diferencial.py
python teste_05_pipeline_completo.py
