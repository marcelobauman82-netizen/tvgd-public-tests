# 1. Configure estes dados
GITHUB_USER="SEU_USUARIO_GITHUB"
REPO_NAME="tvgd-testes-bancada"
PDF_NAME="testes_bancada_TVGD_MEM_nitido(2).pdf"

# 2. Cria a pasta do projeto
mkdir "$REPO_NAME"
cd "$REPO_NAME"

# 3. Inicializa o Git
git init

# 4. Copie o PDF para esta pasta antes de continuar
# Se o PDF estiver em Downloads, por exemplo:
# cp "$HOME/Downloads/$PDF_NAME" .

# 5. Cria um README simples
cat > README.md << 'EOF'
# TVGD–MEM — Testes de Bancada

Este repositório contém o documento:

**Testes de Bancada para TVGD–MEM**

O objetivo do documento é organizar testes experimentais de bancada para avaliar a hipótese TVGD–MEM, incluindo:

- calibração de referência;
- memória de fase estática;
- robustez temporal;
- resposta angular;
- blindagem estrutural;
- dupla fenda de micro-ondas;
- comparação micro-ondas × óptica;
- controle nulo;
- teste adversarial instrumental.

Observação: os testes não provam definitivamente a TVGD, mas podem fornecer evidência experimental compatível com memória de fase, perda controlada de coerência fasorial, resposta angular e ponte entre óptica e micro-ondas.
EOF

# 6. Adiciona os arquivos ao Git
git add README.md "$PDF_NAME"

# 7. Cria o primeiro commit
git commit -m "Adicionar documento de testes de bancada TVGD-MEM"

# 8. Cria o branch principal
git branch -M main

# 9. Conecta ao repositório remoto
git remote add origin "https://github.com/$GITHUB_USER/$REPO_NAME.git"

# 10. Envia para o GitHub
git push -u origin main
