#!/usr/bin/env python3
import os
import re
import sys
import math
import argparse
from pathlib import Path

def exit():
    sys.exit(0)

def extract_cargas_ESP(lines: list[str]) -> list[tuple[int, str, float]]:
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

def constroi_tabela(data: dict, all_files: list, n_solute: int, type_charge: str, errors: dict) -> str:
    """
    Monta a tabela formatada em HTML
    """
    if not data:
        return ""

    short_names = []
    for f in all_files:
        m = re.search(r'(\d+)\.(log|txt)$', f)
        if m:
            short_names.append(f"log{int(m.group(1))}")
        else:
            short_names.append(f)

    max_atoms = 0
    for fname, cargas in data.items():
        if cargas:
            max_num = max(c[0] for c in cargas)
            if max_num > max_atoms:
                max_atoms = max_num

    atom_index = {i: {"symbol": "", "cargas": {}} for i in range(1, max_atoms + 1)}

    for fname, cargas in data.items():
        if not cargas: continue
        for num, symbol, charge in cargas:
            atom_index[num]["symbol"] = symbol
            atom_index[num]["cargas"][fname] = charge

    header = ["Nº Átomo", "Símbolo", "Média", "Desvio Padrão"] + short_names
    colunas = [header]
    
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
            cargas_validas = [c for c in cargas_by_file if isinstance(c, float)]
            m_val = f"{mean(cargas_validas):.6f}" if cargas_validas else ""
            sd_val = f"{stdev(cargas_validas):.6f}" if len(cargas_validas) > 1 else ("0.000000" if cargas_validas else "")
        else:
            m_val, sd_val = "", ""

        coluna = [str(num), symbol, m_val, sd_val] + [
            f"{c:.6f}" if isinstance(c, float) else str(c) for c in cargas_by_file
        ]
        colunas.append(coluna)

    html = f"<html><head><meta charset='utf-8'><title>Cargas {type_charge}</title>"
    html += "<style>table {border-collapse: collapse; width: 100%; font-family: sans-serif; font-size: 14px;} th, td {border: 1px solid #ddd; padding: 8px; text-align: center;} th {background-color: #f2f2f2; position: sticky; top: 0;} .erro {color: red; font-size: 12px;}</style></head><body>"
    html += f"<h2>Tabela de Cargas: {type_charge}</h2>"
    html += "<table>"
    
    # Cabeçalho
    html += "<tr>" + "".join(f"<th>{th}</th>" for th in colunas[0]) + "</tr>"
    
    # Linhas de dados
    for coluna in colunas[1:]:
        html += "<tr>"
        for i, item in enumerate(coluna):
            if "Erro" in item or "Ausentes" in item:
                html += f"<td class='erro'>{item}</td>"
            else:
                html += f"<td>{item}</td>"
        html += "</tr>"
        
    html += "</table></body></html>"
    return html


def main():
    parser = argparse.ArgumentParser(description="Extrai cargas ESP e NBO de logs Gaussian e gera TXTs.")
    parser.add_argument("pasta", nargs="?", default=None, help="Pasta contendo os arquivos .log")
    parser.add_argument("--solute", type=int, default=None, help="Número de átomos do soluto")
    parser.add_argument("--output", default="results", help="Pasta de saída para os TXTs")
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

    print("\nGerando tabelas HTML...")
    
    if esp_data:
        html_path = os.path.join(args.output, "esp_cargas_table.html")
        with open(html_path, 'w', encoding='utf-8') as f:
            f.write(constroi_tabela(esp_data, all_files, n_solute, "ESP", errors))
        print(f"  Salvo: {html_path}")

    if nbo_data:
        html_path = os.path.join(args.output, "nbo_cargas_table.html")
        with open(html_path, 'w', encoding='utf-8') as f:
            f.write(constroi_tabela(nbo_data, all_files, n_solute, "NBO", errors))
        print(f"  Salvo: {html_path}")

    print("\nProcessamento concluído com sucesso!")
if __name__ == "__main__":
    main()