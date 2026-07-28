# Gaussian Log Calculator

Extrai cargas atômicas ESP e NBO de múltiplos arquivos de saída (`.log`) do Gaussian, comparando os resultados entre várias simulações.

## Contexto
Criado como parte da minha pesquisa em química computacional: comparar cargas atômicas entre dezenas de arquivos de log do Gaussian manualmente seria ineficiente, então este script automatiza a extração, o cálculo estatístico e a geração de tabelas comparativas.

## O que faz
- Lê todos os arquivos `.log`/`.txt` de uma pasta, ordenados numericamente pelo nome (ex: `log2`, `log10`, `log20`, na ordem certa em vez de ordem alfabética)
- Extrai cargas ESP do bloco `ESP charges:` e cargas NBO do bloco `Summary of Natural Population Analysis:` de cada log
- Detecta arquivos incompletos ou com erro (crash/timeout, falta de memória, falha de convergência SCF)
- Calcula média e desvio padrão da carga de cada átomo do soluto entre todos os arquivos processados
- Gera as tabelas em HTML, em CSV, ou em ambos, à escolha do usuário
- Informa o somatório das cargas médias do soluto ao final do processamento

## Requisitos
- Python 3.9+
- Sem dependências externas — usa apenas bibliotecas padrão do Python (`os`, `csv`, `re`, `sys`, `math`, `argparse`, `pathlib`)

## Como usar
```bash
python scriptLog.py <pasta> --solute <numero_de_atomos_do_soluto> --output <pasta_de_saida>
```

Ou, sem argumentos, o script pergunta interativamente:
```bash
python scriptLog.py
```

Opções adicionais: `--formato {html,csv,ambos}` (se omitido, o script pergunta), `--decimais N` (casas decimais da carga no CSV de médias, padrão 6) e `--sep CHAR` (separador do CSV, padrão `,`).

## Saída
Na pasta de saída (padrão: `results/`), conforme o formato escolhido:

- `esp_cargas_table.html` / `nbo_cargas_table.html` — tabela para inspeção visual, com erros destacados em vermelho
- `esp_cargas_completo.csv` / `nbo_cargas_completo.csv` — mesma informação do HTML, em formato tabular
- `esp_cargas_medias.csv` / `nbo_cargas_medias.csv` — apenas os átomos do soluto, na ordem do `ljname`, com a carga arredondada na casa decimal aceita pelo input

As tabelas HTML e o CSV completo trazem: número do átomo, símbolo, média, desvio padrão, e uma coluna com o valor bruto de cada arquivo processado (ou o motivo do erro, se o arquivo não pôde ser lido corretamente). O CSV de médias traz apenas número, símbolo, carga média e desvio padrão.

## Limitações / próximos passos
- Assume um formato específico de saída do Gaussian (blocos "ESP charges:" e "Summary of Natural Population Analysis:") — não generaliza para todo tipo de cálculo
- Sem testes automatizados
- Não monta o `ljname` completo: coordenadas e parâmetros de Lennard-Jones não estão nos logs do Gaussian
- O ajuste do somatório para a carga total esperada é feito em módulo separado, ainda em desenvolvimento

---

# Gaussian Log Calculator

Extracts ESP and NBO atomic charges from multiple Gaussian output (`.log`) files, comparing results across several simulations.

## Context
Built as part of my computational chemistry research: comparing atomic charges across dozens of Gaussian log files by hand would be unfeasible, so this script automates the extraction, statistical calculation, and generation of comparative tables.

## What it does
- Reads all `.log`/`.txt` files in a folder, sorted numerically by filename (e.g. `log2`, `log10`, `log20` in the correct order rather than alphabetical order)
- Extracts ESP charges from the `ESP charges:` block and NBO charges from the `Summary of Natural Population Analysis:` block of each log
- Detects and flags incomplete or errored files (crash/timeout, out of memory, SCF convergence failure)
- Calculates the mean and standard deviation of each solute atom's charge across all processed files
- Generates tables as HTML, CSV, or both, at the user's choice
- Reports the sum of the solute's mean charges at the end of processing

## Requirements
- Python 3.9+
- No external dependencies — uses only the Python standard library (`os`, `csv`, `re`, `sys`, `math`, `argparse`, `pathlib`)

## Usage
```bash
python scriptLog.py <folder> --solute <number_of_solute_atoms> --output <output_folder>
```

Or, with no arguments, the script prompts interactively:
```bash
python scriptLog.py
```

Additional options: `--formato {html,csv,ambos}` (if omitted, the script asks), `--decimais N` (decimal places for the charge in the means CSV, default 6) and `--sep CHAR` (CSV delimiter, default `,`).

## Output
In the output folder (default: `results/`), according to the chosen format:

- `esp_cargas_table.html` / `nbo_cargas_table.html` — table for visual inspection, with errors highlighted in red
- `esp_cargas_completo.csv` / `nbo_cargas_completo.csv` — same information as the HTML, in tabular form
- `esp_cargas_medias.csv` / `nbo_cargas_medias.csv` — solute atoms only, in `ljname` order, with the charge rounded to the decimal place accepted by the input

The HTML tables and the full CSV include: atom number, symbol, mean, standard deviation, and one column with the raw value from each processed file (or the error reason, if the file couldn't be read correctly). The means CSV includes only atom number, symbol, mean charge and standard deviation.

## Limitations / next steps
- Assumes a specific Gaussian output format (the "ESP charges:" and "Summary of Natural Population Analysis:" blocks) — doesn't generalize to every calculation type
- No automated tests
- Doesn't assemble the full `ljname`: coordinates and Lennard-Jones parameters are not present in Gaussian logs
- Adjusting the sum to the expected total charge is handled by a separate module, still under development