# Gaussian Log Calculator

Extrai cargas atômicas ESP e NBO de múltiplos arquivos de saída (`.log`) do Gaussian, comparando os resultados entre várias simulações.

## Contexto
Criado como parte da minha pesquisa em química computacional: comparar cargas atômicas entre dezenas de arquivos de log do Gaussian manualmente seria ineficiente, então este script automatiza a extração, o cálculo estatístico e a geração de tabelas comparativas.

## O que faz
- Lê todos os arquivos `.log`/`.txt` de uma pasta, ordenados numericamente pelo nome (ex: `log2`, `log10`, `log20`, na ordem certa em vez de ordem alfabética)
- Extrai cargas ESP do bloco `ESP charges:` e cargas NBO do bloco `Summary of Natural Population Analysis:` de cada log
- Detecta arquivos incompletos ou com erro (crash/timeout, falta de memória, falha de convergência SCF)
- Calcula média e desvio padrão da carga de cada átomo do soluto entre todos os arquivos processados
- Gera duas tabelas HTML (`esp_cargas_table.html` e `nbo_cargas_table.html`), com erros destacados em vermelho

## Requisitos
- Python 3.6+
- Sem dependências externas — usa apenas bibliotecas padrão do Python (`os`, `re`, `sys`, `math`, `argparse`, `pathlib`)

## Como usar
```bash
python scriptLog.py <pasta> --solute <numero_de_atomos_do_soluto> --output <pasta_de_saida>
```

Ou, sem argumentos, o script pergunta interativamente:
```bash
python scriptLog.py
```

## Saída
Duas tabelas HTML na pasta de saída (padrão: `results/`), uma para cargas ESP e outra para NBO. Cada tabela traz: número do átomo, símbolo, média, desvio padrão, e uma coluna com o valor bruto de cada arquivo processado (ou o motivo do erro, se o arquivo não pôde ser lido corretamente).

## Limitações / próximos passos
- Assume um formato específico de saída do Gaussian (blocos "ESP charges:" e "Summary of Natural Population Analysis:") — não generaliza para todo tipo de cálculo
- Sem testes automatizados
- Poderia exportar também em CSV, além do HTML, para facilitar análise posterior em pandas/Excel

# Gaussian Log Calculator
 
Extracts ESP and NBO atomic charges from multiple Gaussian output (`.log`) files, comparing results across several simulations.
 
## Context
Built as part of my computational chemistry research: comparing atomic charges across dozens of Gaussian log files by hand would be unfeasible, so this script automates the extraction, statistical calculation, and generation of comparative tables.
 
## What it does
- Reads all `.log`/`.txt` files in a folder, sorted numerically by filename (e.g. `log2`, `log10`, `log20` in the correct order rather than alphabetical order)
- Extracts ESP charges from the `ESP charges:` block and NBO charges from the `Summary of Natural Population Analysis:` block of each log
- Detects and flags incomplete or errored files (crash/timeout, out of memory, SCF convergence failure)
- Calculates the mean and standard deviation of each solute atom's charge across all processed files
- Generates two HTML tables (`esp_cargas_table.html` and `nbo_cargas_table.html`), with errors highlighted in red
## Requirements
- Python 3.6+
- No external dependencies — uses only the Python standard library (`os`, `re`, `sys`, `math`, `argparse`, `pathlib`)
## Usage
```bash
python scriptLog.py <folder> --solute <number_of_solute_atoms> --output <output_folder>
```
 
Or, with no arguments, the script prompts interactively:
```bash
python scriptLog.py
```
 
## Output
Two HTML tables in the output folder (default: `results/`), one for ESP charges and one for NBO. Each table includes: atom number, symbol, mean, standard deviation, and one column with the raw value from each processed file (or the error reason, if the file couldn't be read correctly).
 
## Limitations / next steps
- Assumes a specific Gaussian output format (the "ESP charges:" and "Summary of Natural Population Analysis:" blocks) — doesn't generalize to every calculation type
- No automated tests
- Could also export to CSV, in addition to HTML, to make downstream analysis in pandas/Excel easier