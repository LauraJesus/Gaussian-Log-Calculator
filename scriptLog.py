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

def pergunta_regras_outliers(type_charge: str) -> dict[str, str]:
    """
    Pergunta ao usuário quais restrições de sinal aplicar.

    Retorna dicionario -> regra ('>0' ou '<0'), ex.: {'Pt': '>0', '2,3': '<0'}.
    Dict vazio significa nenhum filtro.

    O alvo pode ser um símbolo ('Pt') ou uma seleção numérica ('2,3', '4-7').
    As regras são perguntadas separadamente para ESP e NBO.
    """
    resp = input(f"\nRemover outliers por restrição de sinal em {type_charge}? "
                 f"(s/n) [n]: ").strip().lower()
    if resp != 's':
        return {}

    regras: dict[str, str] = {}
    while True:
        alvo = input("  Átomo (símbolo ou números; ENTER encerra): ").strip()
        if not alvo:
            break
        valor = input(f"  Digite o valor da restrição para '{alvo}': ").strip()
        regra = input(f"  Restrição para '{alvo}' — [1] > {valor}, [2] < {valor}: ").strip()
        if regra == '1':
            regras[alvo] = f'>{valor}'
        elif regra == '2':
            regras[alvo] = f'<{valor}'
        else:
            print("    Opção inválida. Digite 1 ou 2.")
            continue
        print(f"    Registrado: {alvo} {regras[alvo]}")

    return regras


def resolve_alvos(alvo: str, cargas: list[tuple[int, str, float]],
                  n_solute: int) -> set[int]:
    """
    Converte o alvo de uma regra no conjunto de números de átomo que ela abrange.
    """
    if alvo and all(c.isdigit() or c in ',- ' for c in alvo):
        nums = set()
        for parte in alvo.split(','):
            parte = parte.strip()
            if not parte:
                continue
            if '-' in parte:
                ini, fim = (int(p) for p in parte.split('-', 1))
                nums.update(range(min(ini, fim), max(ini, fim) + 1))
            else:
                nums.add(int(parte))
        fora = {n for n in nums if n < 1 or n > n_solute}
        if fora:
            raise ValueError(
                f"Regra '{alvo}': átomos fora do soluto (1-{n_solute}): {sorted(fora)}")
        return nums

    return {num for num, simbolo, _ in cargas
            if simbolo == alvo and num <= n_solute}


def filtra_cargas(data: dict, regras: dict[str, str],
                  n_solute: int) -> tuple[dict, dict[tuple[int, str], str]]:
    """
    Remove cargas individuais que violem as restrições de sinal.
    Retorna (data filtrado, dict (num, arquivo) -> motivo do descarte).
    """
    descartes: dict[tuple[int, str], str] = {}
    if not regras:
        return data, descartes

    filtrado: dict = {}

    for fname, cargas in data.items():
        alvos: dict[int, str] = {}
        for alvo, regra in regras.items():
            for num in resolve_alvos(alvo, cargas, n_solute):
                alvos[num] = regra

        mantidas = []
        for num, simbolo, carga in cargas:
            regra = alvos.get(num)
            if regra is None:
                mantidas.append((num, simbolo, carga))
            elif regra.startswith('>') and carga <= float(regra[1:]):
                descartes[(num, fname)] = f"Outlier ({carga:+.6f}, esperado > {regra[1:]})"
            elif regra.startswith('<') and carga >= float(regra[1:]):
                descartes[(num, fname)] = f"Outlier ({carga:+.6f}, esperado < {regra[1:]})"
            else:
                mantidas.append((num, simbolo, carga))

        filtrado[fname] = mantidas

    return filtrado, descartes


def mapa_simbolos(data: dict) -> dict[int, str]:
    """
    Monta o mapa número -> símbolo a partir dos dados ANTES do filtro.
    """
    simbolos: dict[int, str] = {}
    for cargas in data.values():
        for num, simbolo, _carga in cargas:
            simbolos.setdefault(num, simbolo)
    return simbolos


def relata_descartes(descartes: dict[tuple[int, str], str], data: dict,
                     n_solute: int, type_charge: str) -> None:
    """
    Imprime o resumo dos descartes: quantas cargas saíram e de quais átomos.

    Avisa quando um átomo perdeu todas as configurações, caso em que a média
    vira NaN.
    """
    if not descartes:
        print(f"  {type_charge}: nenhuma carga descartada.")
        return

    por_atomo: dict[int, int] = {}
    for (num, _fname) in descartes:
        por_atomo[num] = por_atomo.get(num, 0) + 1

    n_config = len(data)
    print(f"  {type_charge}: {len(descartes)} carga(s) descartada(s) "
          f"em {n_config} configuração(ões).")
    for num in sorted(por_atomo):
        restantes = n_config - por_atomo[num]
        marca = "  [ATENÇÃO] átomo sem nenhuma carga válida!" if restantes == 0 else ""
        print(f"    átomo {num:>3}: {por_atomo[num]} descartada(s), "
              f"N = {restantes}{marca}")


def processFolder(folderPath: str) -> tuple[dict, dict, dict, list]:
    esp_data = {}
    nbo_data = {}
    errors = {}

    folderPath = Path(folderPath).expanduser()
    
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

def calcula_estatisticas(data: dict, n_solute: int,
                         simbolos: dict) -> dict[int, tuple[str, int, float, float]]:
    """
    Calcula N, média e desvio padrão de cada átomo do soluto.

    Um átomo cujas cargas foram todas descartadas
    aparece com N = 0 e média NaN.

    Retorna dict num -> (símbolo, N, média, desvio padrão).
    """
    por_atomo: dict[int, list[float]] = {}
    for cargas in data.values():
        for num, _simbolo, carga in cargas:
            if num <= n_solute:
                por_atomo.setdefault(num, []).append(carga)

    nums = {n for n in simbolos if n <= n_solute} | set(por_atomo)

    estat: dict[int, tuple[str, int, float, float]] = {}
    for num in sorted(nums):
        valores = por_atomo.get(num, [])
        simbolo = simbolos.get(num, "")
        desvio = stdev(valores) if len(valores) > 1 else (0.0 if valores else float('nan'))
        estat[num] = (simbolo, len(valores), mean(valores), desvio)

    return estat


def escreve_txt_medias(estat: dict[int, tuple[str, int, float, float]],
                       caminho: str, decimais: int = 6) -> None:
    """
    Grava um TXT com uma linha por átomo do soluto, no formato 'simbolo media'.
    """
    linhas = []
    faltando = []

    for num in sorted(estat):
        simbolo, n, media, _desvio = estat[num]
        if n == 0:
            linhas.append(f"{simbolo} nan\n")
            faltando.append(num)
        else:
            linhas.append(f"{simbolo} {media:.{decimais}f}\n")

    with open(caminho, 'w', encoding='utf-8') as f:
        f.writelines(linhas)

    soma = sum(m for _s, n, m, _d in estat.values() if n > 0)
    print(f"  Salvo: {caminho}")
    print(f"    Soma das médias sobre os {len(estat)} átomos do soluto: {soma:+.6f}")
    if faltando:
        print(f"    [ATENÇÃO] átomos sem carga válida, gravados como 'nan': {faltando}")


def constroi_tabela(data: dict, all_files: list, n_solute: int, type_charge: str,
                    errors: dict, descartes: dict = None,
                    regras: dict = None, simbolos: dict = None) -> str:
    """
    Monta a tabela formatada em HTML.

    A coluna N traz quantas configurações entraram na média de cada átomo do
    soluto.
    """
    descartes = descartes or {}
    regras = regras or {}
    simbolos = simbolos or {}
    estat = calcula_estatisticas(data, n_solute, simbolos)
    if not data:
        return ""

    short_names = []
    for f in all_files:
        m = re.search(r'(\d+)\.(log|txt)$', f)
        if m:
            short_names.append(f"log{int(m.group(1))}")
        else:
            short_names.append(f)

    max_atoms = max(simbolos) if simbolos else 0
    for fname, cargas in data.items():
        if cargas:
            max_num = max(c[0] for c in cargas)
            if max_num > max_atoms:
                max_atoms = max_num

    atom_index = {i: {"symbol": simbolos.get(i, ""), "cargas": {}}
                  for i in range(1, max_atoms + 1)}

    for fname, cargas in data.items():
        if not cargas: continue
        for num, symbol, charge in cargas:
            atom_index[num]["symbol"] = symbol
            atom_index[num]["cargas"][fname] = charge

    header = ["Nº Átomo", "Símbolo", "N", "Média", "Desvio Padrão"] + short_names
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
                motivo = descartes.get((num, fname)) or errors.get(fname, "Dados Ausentes")
                cargas_by_file.append(motivo)
                
        if num <= n_solute:
            _sb, n_amostra, media, desvio = estat.get(
                num, ("", 0, float('nan'), float('nan')))
            n_val = str(n_amostra)
            m_val = f"{media:.6f}" if n_amostra else ""
            sd_val = f"{desvio:.6f}" if n_amostra else ""
        else:
            n_val, m_val, sd_val = "", "", ""

        coluna = [str(num), symbol, n_val, m_val, sd_val] + [
            f"{c:.6f}" if isinstance(c, float) else str(c) for c in cargas_by_file
        ]
        colunas.append(coluna)

    html = f"<html><head><meta charset='utf-8'><title>Cargas {type_charge}</title>"
    html += "<style>table {border-collapse: collapse; width: 100%; font-family: sans-serif; font-size: 14px;} th, td {border: 1px solid #ddd; padding: 8px; text-align: center;} th {background-color: #f2f2f2; position: sticky; top: 0;} .erro {color: red; font-size: 12px;} .outlier {color: #b35c00; background-color: #fff4e5; font-size: 12px;} .legenda {font-family: sans-serif; font-size: 13px; color: #333;}</style></head><body>"
    html += f"<h2>Tabela de Cargas: {type_charge}</h2>"

    if regras:
        itens = "".join(f"<li><b>{alvo}</b>: carga {regra}</li>"
                        for alvo, regra in regras.items())
        html += (f"<div class='legenda'><p>Restrições de sinal aplicadas ao modelo "
                 f"{type_charge}:</p><ul>{itens}</ul>"
                 f"<p>Cargas descartadas aparecem destacadas em laranja e não entram "
                 f"na média. A coluna <b>N</b> traz o número de configurações usadas "
                 f"na média de cada átomo, que varia entre átomos quando há descarte. "
                 f"O desvio padrão de um átomo filtrado é calculado sobre distribuição "
                 f"truncada e fica subestimado. Células em vermelho indicam falha de "
                 f"cálculo, não descarte por restrição.</p></div>")
    else:
        html += ("<div class='legenda'><p>Nenhuma restrição de sinal aplicada. "
                 "A coluna <b>N</b> traz o número de configurações usadas na média "
                 "de cada átomo.</p></div>")

    html += "<table>"
    
    # Cabeçalho
    html += "<tr>" + "".join(f"<th>{th}</th>" for th in colunas[0]) + "</tr>"
    
    # Linhas de dados
    for coluna in colunas[1:]:
        html += "<tr>"
        for i, item in enumerate(coluna):
            if "Outlier" in item:
                html += f"<td class='outlier'>{item}</td>"
            elif "Erro" in item or "Ausentes" in item or "Falha" in item or "Incompleto" in item:
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
    parser.add_argument("--sem-outliers", action="store_true",
                        help="Não perguntar restrições de sinal (uso não interativo)")
    parser.add_argument("--decimais", type=int, default=6,
                        help="Casas decimais no TXT de médias (padrão: 6)")
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

    regras_esp: dict[str, str] = {}
    regras_nbo: dict[str, str] = {}
    descartes_esp: dict = {}
    descartes_nbo: dict = {}

    if not args.sem_outliers:
        try:
            if esp_data:
                regras_esp = pergunta_regras_outliers("ESP")
            if nbo_data:
                regras_nbo = pergunta_regras_outliers("NBO")
        except EOFError:
            print("\n[AVISO] Entrada não interativa: nenhum filtro aplicado.")

    simbolos_esp = mapa_simbolos(esp_data)
    simbolos_nbo = mapa_simbolos(nbo_data)

    print("\nAplicando restrições...")
    try:
        if esp_data:
            esp_data, descartes_esp = filtra_cargas(esp_data, regras_esp, n_solute)
            relata_descartes(descartes_esp, esp_data, n_solute, "ESP")
        if nbo_data:
            nbo_data, descartes_nbo = filtra_cargas(nbo_data, regras_nbo, n_solute)
            relata_descartes(descartes_nbo, nbo_data, n_solute, "NBO")
    except ValueError as e:
        print(f"[ERRO] {e}")
        sys.exit(1)

    print("\nGerando saídas...")
    
    if esp_data:
        html_path = os.path.join(args.output, "esp_cargas_table.html")
        with open(html_path, 'w', encoding='utf-8') as f:
            f.write(constroi_tabela(esp_data, all_files, n_solute, "ESP", errors,
                                    descartes_esp, regras_esp, simbolos_esp))
        print(f"  Salvo: {html_path}")

        estat_esp = calcula_estatisticas(esp_data, n_solute, simbolos_esp)
        escreve_txt_medias(estat_esp,
                           os.path.join(args.output, "esp_cargas_medias.txt"),
                           args.decimais)

    if nbo_data:
        html_path = os.path.join(args.output, "nbo_cargas_table.html")
        with open(html_path, 'w', encoding='utf-8') as f:
            f.write(constroi_tabela(nbo_data, all_files, n_solute, "NBO", errors,
                                    descartes_nbo, regras_nbo, simbolos_nbo))
        print(f"  Salvo: {html_path}")

        estat_nbo = calcula_estatisticas(nbo_data, n_solute, simbolos_nbo)
        escreve_txt_medias(estat_nbo,
                           os.path.join(args.output, "nbo_cargas_medias.txt"),
                           args.decimais)

    print("\nProcessamento concluído com sucesso!")
if __name__ == "__main__":
    main()