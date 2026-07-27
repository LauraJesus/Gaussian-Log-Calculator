#!/usr/bin/env python3
import os
import csv
import re
import sys
import math
import argparse
from pathlib import Path

outliers = {
    "Pt": lambda carga:carga>0,
    #colocar outras se necessário
}


def clean_charges(cargas: list[tuple[int, str, float]],regras: dict = outliers):
    valid, invalid = [], []
    for num, symbol, charge in cargas:
        regra = regras.get(symbol)
        if regra is not None and not regra(charge):
            invalid.append((num, symbol, charge,f"{symbol} violou regra química"))
        else:
            valid.append((num, symbol, charge))

    return valid, invalid

def exit():
    sys.exit(0)

def extract_cargas_ESP(lines: list[str]) -> list[tuple[int, str, float]]:
    """
    Extrai as cargas ESP/CHELPG do bloco 'ESP charges:' do log do Gaussian.
    Formato esperado da linha: numero, simbolo, carga.
    Retorna [] se o bloco nao existir (calculo abortado antes do ajuste).
    """
    cargas = []
    inBlock = False
    for line in lines:
        if "ESP charges:" in line:
            inBlock = True
            continue
        if inBlock:
            if re.match(r'^\s+\d+\s*$', line):
                continue
            
            # ESP: "     1  Pt   0.544647" -> Num, Símbolo, Carga
            m = re.match(r'^\s+(\d+)\s+([A-Za-z]+)\s+([-\d.]+)\s*$', line)
            if m:
                cargas.append((int(m.group(1)), m.group(2), float(m.group(3))))
                continue
            
            if "Sum of ESP charges" in line:
                break
    return cargas

def extract_cargas_NBO(lines: list[str]) -> list[tuple[int, str, float]]:
    """
    Extrai as cargas NBO do bloco 'Summary of Natural Population Analysis'.
    Atencao a ordem das colunas: no NBO e simbolo, numero, carga - invertida
    em relacao ao bloco ESP. A saida e padronizada como (numero, simbolo, carga).
    """
    cargas = []
    inBlock = False
    for line in lines:
        if "Summary of Natural Population Analysis:" in line:
            inBlock = True
            continue
        if inBlock:
            if "===" in line or "Natural Population" in line or "Atom No" in line or "----" in line:
                if "===" in line and len(cargas) > 0:
                    break
                continue
            
            # NBO: "    Pt  1   -0.08039" -> Símbolo, Num, Carga
            words = line.split()
            if len(words) >= 3 and words[1].isdigit():
                try:
                    num = int(words[1])
                    symbol = words[0]
                    charge = float(words[2])
                    cargas.append((num, symbol, charge))
                except ValueError:
                    pass
    return cargas

def check_log_errors(lines: list[str]) -> str:
    """
    Distingue crash/timeout, erro de terminacao, falta de memoria e falha de
    convergencia SCF.
    """
    error_msg = "Arquivo Incompleto (Crash/Timeout)"
    last_lines = lines[-30:] if len(lines) > 30 else lines
    for line in last_lines:
        if "Error termination" in line or "Termination error" in line:
            return "Erro de terminação no Gaussian"
        elif "Out of memory" in line or "Allocation failed" in line:
            return "Falta de Memória (RAM)"
        elif "SCF Error" in line or "Convergence failure" in line:
            return "Falha na convergência (SCF)"
        elif "Normal termination" in line:
            return "Finalizado sem tags de carga"
    return error_msg

def processFolder(folderPath: str) -> tuple[dict, dict, dict, list]:
    """
    Percorre a pasta de logs e extrai as cargas de cada configuracao.
    Ordena os arquivos pelo numero no nome (log10, log20, ...), aplica a
    limpeza de outliers no ESP e registra o motivo da falha dos logs que nao
    renderam cargas.
    Retorna (esp_data, nbo_data, errors, all_files).
    """
    esp_data = {}
    nbo_data = {}
    errors = {}

    folderPath = Path(folderPath).expanduser()
    
    # Extrai todos os logs e txt ordenando pelo número no meio do nome (ex: log10, log20)
    logFilesPaths = sorted(
        list(folderPath.glob("*.log")) + list(folderPath.glob("*.txt")),
        key=lambda x: [int(c) if c.isdigit() else c for c in re.split(r'(\d+)', x.name)]
    )
    
    all_files = [f.name for f in logFilesPaths]

    if not all_files:
        print(f"[ERRO] Nenhum arquivo .log ou .txt encontrado em '{folderPath}'.")
        return esp_data, nbo_data, errors, all_files

    for fpath in logFilesPaths:
        fname = fpath.name
        print(f"Processando {fname}...")
        
        with open(fpath, 'r', encoding='utf-8', errors='replace') as f:
            lines = f.readlines()
        
        esp_cargas = extract_cargas_ESP(lines)
        esp_cargas, invalid_esp = clean_charges(esp_cargas)
        if invalid_esp:
            print(f"  AVISO: Cargas ESP inválidas em {fname}:")
            for num, symbol, charge, motivo in invalid_esp:
                print(f"    Átomo {num} ({symbol}): {charge:.6f} -> {motivo}")
        
        if esp_cargas:
            esp_data[fname] = esp_cargas
        else:
            print(f"  AVISO: Nenhuma carga ESP encontrada em {fname}")
            
        nbo_cargas = extract_cargas_NBO(lines)
        if nbo_cargas:
            nbo_data[fname] = nbo_cargas
        else:
            print(f"  AVISO: Nenhuma carga NBO encontrada em {fname}")
            
        if not esp_cargas or not nbo_cargas:
            errors[fname] = check_log_errors(lines)

    return esp_data, nbo_data, errors, all_files

def mean(values: list[float]) -> float:
    return sum(values)/len(values) if values else float('nan')

def stdev(values: list[float]) -> float:
    if len(values) < 2:
        return float('nan')
    m = mean(values)
    variance = sum((x - m) ** 2 for x in values) / (len(values) - 1)
    return math.sqrt(variance)

def nomes_curtos(all_files: list) -> list:
    """Converte 'algumacoisa_log10.log' -> 'log10' para o cabecalho das tabelas."""
    short = []
    for f in all_files:
        m = re.search(r'(\d+)\.(log|txt)$', f)
        if m:
            short.append(f"log{int(m.group(1))}")
        else:
            short.append(f)
    return short

def monta_linhas(data: dict, all_files: list, n_solute: int, errors: dict) -> tuple[list,list]:
    """Usada tanto pelo HTML quanto pelos CSVs - qualquer mudanca na regra de
    media vale para as tres saidas ao mesmo tempo.

    Retorna (short_names, linhas), onde cada linha e:
        [num_atomo, simbolo, media, desvio_padrao, carga_conf1, carga_conf2, ...]
    Cada renderizador monta o proprio cabecalho a partir de short_names."""
    if not data:
        return [], []

    max_atoms = 0
    for cargas in data.values():
        if cargas:
            max_atoms = max(max_atoms, max(c[0] for c in cargas))

    atom_index = {i: {"symbol": "", "cargas": {}} for i in range(1, max_atoms + 1)}

    for fname, cargas in data.items():
        if not cargas: continue
        for num, symbol, charge in cargas:
            atom_index[num]["symbol"] = symbol
            atom_index[num]["cargas"][fname] = charge

    short_names = nomes_curtos(all_files)
    linhas =[]
    for num in range(1, max_atoms + 1):
        symbol = atom_index[num]["symbol"]
        if not symbol:
            continue

        cargas_by_file = []
        for fname in all_files:
            val = atom_index[num]["cargas"].get(fname, None)
            if val is not None:
                cargas_by_file.append(val)
            else:
                motivo = errors.get(fname, "Dados Ausentes")
                cargas_by_file.append(motivo)

        if num <= n_solute:
            validas = [c for c in cargas_by_file if isinstance(c, float)]
            m_val= f"{mean(validas):.6f}" if validas else ""
            if len(validas) > 1:
                sd_val = f"{stdev(validas):.6f}"
            else:
                sd_val = "0.000000" if validas else ""
        else:
            m_val, sd_val = "", ""

        linhas.append([str(num), symbol, m_val, sd_val] + [
            f"{c:.6f}" if isinstance(c, float) else str(c) for c in cargas_by_file
        ])
    return short_names, linhas

def constroi_tabela(data: dict, all_files: list, n_solute: int, type_charge: str, errors: dict) -> str:
    """
    Monta a tabela formatada em HTML a partir de monta_linhas().
    """
    short_names, linhas = monta_linhas(data, all_files, n_solute, errors)
    if not linhas:
        return ""
 
    header = ["Nº Átomo", "Símbolo", "Média", "Desvio Padrão"] + short_names
 
    html = f"<html><head><meta charset='utf-8'><title>Cargas {type_charge}</title>"
    html += "<style>table {border-collapse: collapse; width: 100%; font-family: sans-serif; font-size: 14px;} th, td {border: 1px solid #ddd; padding: 8px; text-align: center;} th {background-color: #f2f2f2; position: sticky; top: 0;} .erro {color: red; font-size: 12px;}</style></head><body>"
    html += f"<h2>Tabela de Cargas: {type_charge}</h2>"
    html += "<table>"
 
    html += "<tr>" + "".join(f"<th>{th}</th>" for th in header) + "</tr>"
 
    for linha in linhas:
        html += "<tr>"
        for item in linha:
            if "Erro" in item or "Ausentes" in item:
                html += f"<td class='erro'>{item}</td>"
            else:
                html += f"<td>{item}</td>"
        html += "</tr>"
 
    html += "</table></body></html>"
    return html


def escreve_csv_completo(caminho: str, data: dict, all_files: list,
                         n_solute: int, errors: dict, sep: str = ",") -> bool:
    """CSV com todas as configuracoes, media e desvio. Espelha o HTML."""
    short_names, linhas = monta_linhas(data, all_files, n_solute, errors)
    if not linhas:
        return False
 
    header = ["num_atomo", "simbolo", "media", "desvio_padrao"] + short_names
 
    with open(caminho, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter=sep, quoting=csv.QUOTE_MINIMAL)
        w.writerow(header)
        w.writerows(linhas)
    return True
 
 
def escreve_csv_medias(caminho: str, data: dict, all_files: list,
                       n_solute: int, errors: dict, decimais: int = 6,
                       sep: str = ",", incluir_desvio: bool = True) -> tuple[bool, float]:
    """
    CSV so com os atomos do soluto, na ordem do ljname.
 
    'decimais' controla o arredondamento da carga - deve bater com a casa
    decimal que o input do DICE aceita. O desvio padrao vai numa coluna
    separada apenas como registro; ele NAO e usado no ajuste de residuo.
 
    Retorna (sucesso, soma_das_medias). A soma serve para o modulo de ajuste
    de residuo saber se ha residuo a corrigir.
    """
    _, linhas = monta_linhas(data, all_files, n_solute, errors)
    if not linhas:
        return False, 0.0
 
    header = ["num_atomo", "simbolo", "carga_media"]
    if incluir_desvio:
        header.append("desvio_padrao")
 
    saida, soma = [], 0.0
    for linha in linhas:
        if int(linha[0]) > n_solute:
            continue
        if not linha[2]:          # sem media -> atomo sem nenhuma config valida
            continue
 
        carga = round(float(linha[2]), decimais)
        soma += carga
 
        reg = [linha[0], linha[1], f"{carga:.{decimais}f}"]
        if incluir_desvio:
            reg.append(linha[3])
        saida.append(reg)
 
    with open(caminho, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter=sep, quoting=csv.QUOTE_MINIMAL)
        w.writerow(header)
        w.writerows(saida)
 
    return True, round(soma, decimais)
 
 
FORMATOS = {
    "1": ("html",  "Apenas HTML (tabela para visualizacao)"),
    "2": ("csv",   "Apenas CSV (completo + medias no padrao do input)"),
    "3": ("ambos", "HTML e CSV"),
}
 

 
def pergunta_formato(padrao: str = "3") -> str:
    """Pergunta ao usuario qual formato de saida deseja."""
    print("\nQual formato de saída você quer?")
    for k, (_, desc) in FORMATOS.items():
        marca = " (padrão)" if k == padrao else ""
        print(f"  [{k}] {desc}{marca}")
 
    while True:
        resp = input("Escolha: ").strip()
        if not resp:
            return FORMATOS[padrao][0]
        if resp.lower() == "sair":
            sys.exit(0)
        if resp in FORMATOS:
            return FORMATOS[resp][0]
        print("  Opção inválida. Digite 1, 2 ou 3.")
 
 
def pergunta_decimais(padrao: int = 6) -> int:
    """Pergunta a casa decimal aceita pelo input do DICE."""
    while True:
        resp = input(f"Quantas casas decimais o input aceita? [{padrao}]: ").strip()
        if not resp:
            return padrao
        if resp.lower() == "sair":
            sys.exit(0)
        if resp.isdigit() and 1 <= int(resp) <= 10:
            return int(resp)
        print("  Digite um inteiro entre 1 e 10.")
 


def main():
    parser = argparse.ArgumentParser(description="Extrai cargas ESP e NBO de logs Gaussian e gera TXTs.")
    parser.add_argument("pasta", nargs="?", default=None, help="Pasta contendo os arquivos .log")
    parser.add_argument("--solute", type=int, default=None, help="Número de átomos do soluto")
    parser.add_argument("--output", default="results", help="Pasta de saída para os TXTs")
    parser.add_argument("--formato", choices=["html", "csv", "ambos"], default=None,
                        help="Formato de saída. Se omitido, o programa pergunta.")
    parser.add_argument("--decimais", type=int, default=None,
                        help="Casas decimais da carga no CSV de médias (padrão 6)")
    parser.add_argument("--sep", default=",",
                        help="Separador do CSV (use ';' para Excel pt-BR)")
    args = parser.parse_args()
 
    pasta_logs = args.pasta
    if pasta_logs is None:
        while True:
            pasta_logs = input("Qual é o nome da pasta com os arquivos de log? ")
            if os.path.isdir(pasta_logs):
                break
            print(f"  [ERRO] Pasta '{pasta_logs}' não encontrada. Tente novamente.")
    else:
        if not os.path.isdir(pasta_logs):
            print(f"[ERRO] Pasta não encontrada: '{pasta_logs}'")
            sys.exit(1)
 
    n_solute = args.solute
    if n_solute is None:
        while True:
            try:
                n_solute = int(input("Quantos átomos tem o soluto? "))
                if n_solute > 0:
                    break
                print("  Digite um número inteiro positivo.")
            except ValueError:
                print("  Entrada inválida. Digite um número inteiro.")
 
    print(f"\nPasta: {pasta_logs}")
    print(f"Átomos do soluto: {n_solute}")
    print(f"Saída: {args.output}\n")
 
    print("Lendo arquivos...")
    esp_data, nbo_data, errors, all_files = processFolder(pasta_logs)
 
    if not esp_data and not nbo_data:
        print("[ERRO] Nenhum dado extraído. Verifique os arquivos na pasta.")
        sys.exit(1)
 
    os.makedirs(args.output, exist_ok=True)
 
    formato = args.formato if args.formato else pergunta_formato()
 
    decimais = args.decimais
    if formato in ("csv", "ambos") and decimais is None:
        decimais = pergunta_decimais()
 
    conjuntos = [("ESP", esp_data), ("NBO", nbo_data)]
 
    if formato in ("html", "ambos"):
        print("\nGerando tabelas HTML...")
        for tipo, data in conjuntos:
            if not data:
                continue
            html_path = os.path.join(args.output, f"{tipo.lower()}_cargas_table.html")
            with open(html_path, 'w', encoding='utf-8') as f:
                f.write(constroi_tabela(data, all_files, n_solute, tipo, errors))
            print(f"  Salvo: {html_path}")
 
    if formato in ("csv", "ambos"):
        print("\nGerando CSVs...")
        for tipo, data in conjuntos:
            if not data:
                continue
 
            csv_full = os.path.join(args.output, f"{tipo.lower()}_cargas_completo.csv")
            if escreve_csv_completo(csv_full, data, all_files, n_solute,
                                    errors, sep=args.sep):
                print(f"  Salvo: {csv_full}")
 
            csv_med = os.path.join(args.output, f"{tipo.lower()}_cargas_medias.csv")
            ok, soma = escreve_csv_medias(csv_med, data, all_files, n_solute,
                                          errors, decimais=decimais, sep=args.sep)
            if ok:
                print(f"  Salvo: {csv_med}")
                print(f"    Soma das cargas médias do soluto ({tipo}): {soma:.{decimais}f}")
                if abs(soma) > 10 ** (-decimais):
                    print(f"    ATENÇÃO: somatório != 0 -> resíduo a ser tratado "
                          f"pelo módulo de ajuste.")
 
    print("\nProcessamento concluído com sucesso!")
if __name__ == "__main__":
    main()