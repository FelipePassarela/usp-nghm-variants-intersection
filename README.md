# Resumo

Script para calcular frequências de genótipos a partir de um arquivo VCF, separadas por coortes. O script lê um VCF com genótipos e um ou mais arquivos CSV representando coortes (cada CSV contendo IDs de amostras), filtra as amostras de interesse e gera, para cada coorte, um CSV com as frequências de genótipos por variante.

## Requisitos

- **Python:** 3.13+

Instalação rápida (recomendado em um virtualenv):

```pwsh
python -m pip install --upgrade pip
pip install -r requirements.txt
```

Se for usar um gerenciador (ex.: `pyproject.toml`), siga seu fluxo (`poetry`, `pipx`, etc.).

## Uso

Exemplo básico (PowerShell):

```pwsh
python .\main.py .\data\COVID-SNPs3.recode.vcf \
 "data\cohorts\WES CEGH x UFES(90+).csv" \
 "data\cohorts\WES CEGH x UFES(INFECTADOS).csv" \
 -o .\frequencies\
```

- Argumentos:
  - `vcf_path` : caminho para o arquivo VCF.
  - `cohorts` : um ou mais arquivos CSV de coorte (separados por espaço).
  - `--output_dir` / `-o` : diretório de saída opcional; se não informado, as frequências são impressas no stdout.

## Formato de entrada das coortes

- O script lê coortes como arquivos CSV. Os arquivos de coorte não precisam ter cabeçalho; o script pesquisa a coluna que contém IDs de amostra.
- O script procura uma coluna cuja célula case com o padrão: `C01234-ExC123-xgenV1`. Se nenhum arquivo de coorte contiver uma coluna com esse padrão, o script lançará um `ValueError`.

Observação: Garanta que os IDs nas coortes correspondam exatamente aos nomes das amostras no cabeçalho do VCF.

## Saída

- Para cada coorte passada, o script gera um CSV no diretório de saída com o nome do arquivo de coorte. Exemplo: `entrada\coorte1.csv` -> `saida\coorte1.csv`.
- O CSV resultante contém as variantes como linhas (index no formato `chrom_pos`, por exemplo `chr1_123456`) e colunas com as frequências relativas dos genótipos: `0/0`, `0/1`, `1/1`, `./.`. Valores estão normalizados por variante (soma = 1 quando há genótipos chamados).

## Como o script funciona (resumo técnico)

- `get_cohort_ids(cohorts_path)` : lê os CSVs de coorte, detecta a coluna de IDs usando regex, alinha as colunas (preenchendo com `NaN` quando necessário) e retorna um `DataFrame` com uma coluna por coorte.
- `genotypes_matrix(vcf_path, ids_to_keep)` : abre o VCF com `pysam`, itera sobre os registros e monta um `DataFrame` com amostras nas linhas e variantes nas colunas. Genótipos são strings como `0/0`, `0/1`, `1/1` ou `./.`.
- `genotype_frequencies(df)` : conta ocorrências de cada genótipo por variante e calcula frequência relativa (contagem / total de amostras com chamada).

## Dicas e observações

- Nomes de amostras no VCF devem corresponder aos IDs dos CSVs de coorte.
- Para VCFs muito grandes, a matriz (amostras × variantes) pode consumir muita memória. Se possível, gere um VCF filtrado apenas com as variantes de interesse ou use amostras limitadas via `ids_to_keep` (passando apenas coorte(s) desejadas).
