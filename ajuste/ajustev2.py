#!/usr/bin/env python3
"""
ajuste.py — Ajuste de resíduo de carga do ciclo de polarização iterativa.

Fecha o somatório das cargas do soluto em zero, distribuindo o resíduo entre os
grupos quimicamente equivalentes proporcionalmente à largura do intervalo entre
a carga inicial e a carga média de cada grupo.

Entradas:
  - Arquivo DICE (ljname) usado na simulação inicial do ciclo, ex.:
    s0Pt001DFT.txt. As cargas iniciais do soluto (q_inicial) são lidas
    diretamente da coluna de carga do bloco da primeira molécula (o soluto,
    por convenção do DICE — Manual DICE/LAQC). As demais linhas do arquivo
    (regra de combinação, número de tipos de molécula, título, bloco do
    solvente, '$end') são preservadas tal como estão, para reconstrução do
    arquivo na saída.
  - TXT de médias gerado pelo scriptLog.py (esp_cargas_medias.txt / nbo_...),
    uma linha 'simbolo media' por átomo, na mesma ordem do bloco do soluto
    no arquivo DICE
  - Arquivo de grupos de equivalência química, declarado pelo usuário

Saídas:
  - Tabela das cargas ajustadas, impressa no terminal
  - Opcionalmente, uma cópia do arquivo DICE de entrada com as cargas do
    soluto substituídas pelas cargas ajustadas — todo o resto do arquivo
    (posições, rótulos, parâmetros de Lennard-Jones, título, bloco do
    solvente, '$end') permanece byte-a-byte idêntico ao original

Os grupos de equivalência são declarados a partir da estrutura química, e não
derivados por percepção automática de ligação. É a prática do RESP — Bayly,
Cieplak, Cornell & Kollman, J. Phys. Chem. 1993, 97, 10269 — na qual as
restrições de equivalência são impostas pelo usuário.

O mesmo arquivo de grupos serve às duas trilhas, ESP e NBO: o agrupamento
depende só da topologia do soluto, nunca das cargas. Manter a partição idêntica
nos dois modelos é o que garante que a diferença entre eles reflita apenas o
modelo de carga.
"""

import re
import sys
import argparse
from dataclasses import dataclass, field
from decimal import Decimal
from pathlib import Path


# ---------------------------------------------------------------------------
# Tabela periódica (Z -> símbolo). Fato de referência, não é escolha
# metodológica do projeto.
# ---------------------------------------------------------------------------

SIMBOLO_POR_Z = {
    1: 'H', 2: 'He', 3: 'Li', 4: 'Be', 5: 'B', 6: 'C', 7: 'N', 8: 'O', 9: 'F',
    10: 'Ne', 11: 'Na', 12: 'Mg', 13: 'Al', 14: 'Si', 15: 'P', 16: 'S',
    17: 'Cl', 18: 'Ar', 19: 'K', 20: 'Ca', 21: 'Sc', 22: 'Ti', 23: 'V',
    24: 'Cr', 25: 'Mn', 26: 'Fe', 27: 'Co', 28: 'Ni', 29: 'Cu', 30: 'Zn',
    31: 'Ga', 32: 'Ge', 33: 'As', 34: 'Se', 35: 'Br', 36: 'Kr', 37: 'Rb',
    38: 'Sr', 39: 'Y', 40: 'Zr', 41: 'Nb', 42: 'Mo', 43: 'Tc', 44: 'Ru',
    45: 'Rh', 46: 'Pd', 47: 'Ag', 48: 'Cd', 49: 'In', 50: 'Sn', 51: 'Sb',
    52: 'Te', 53: 'I', 54: 'Xe', 55: 'Cs', 56: 'Ba', 57: 'La', 58: 'Ce',
    59: 'Pr', 60: 'Nd', 61: 'Pm', 62: 'Sm', 63: 'Eu', 64: 'Gd', 65: 'Tb',
    66: 'Dy', 67: 'Ho', 68: 'Er', 69: 'Tm', 70: 'Yb', 71: 'Lu', 72: 'Hf',
    73: 'Ta', 74: 'W', 75: 'Re', 76: 'Os', 77: 'Ir', 78: 'Pt', 79: 'Au',
    80: 'Hg', 81: 'Tl', 82: 'Pb', 83: 'Bi', 84: 'Po', 85: 'At', 86: 'Rn',
    87: 'Fr', 88: 'Ra', 89: 'Ac', 90: 'Th', 91: 'Pa', 92: 'U', 93: 'Np',
    94: 'Pu', 95: 'Am', 96: 'Cm', 97: 'Bk', 98: 'Cf', 99: 'Es', 100: 'Fm',
    101: 'Md', 102: 'No', 103: 'Lr',
    # 104 é usado pelo DICE como convenção para "átomo sem massa"
    # (ponto extra de interação), não como elemento real — Manual DICE/LAQC.
    104: 'EP',
}


def simbolo_do_z(an: int) -> str:
    """Converte número atômico (coluna 'an' do arquivo DICE) em símbolo."""
    if an in SIMBOLO_POR_Z:
        return SIMBOLO_POR_Z[an]
    raise ValueError(
        f"Número atômico {an} não reconhecido. Informe o símbolo manualmente "
        f"ou atualize SIMBOLO_POR_Z."
    )


class Atomo:
    """Átomo do soluto."""

    def __init__(self, num: int, simbolo: str, carga_inicial: float):
        self.num = num
        self.simbolo = simbolo
        self.carga_inicial = carga_inicial
        self.carga_media: float = 0.0
        self.carga_final: float = 0.0
        self.grupo: str = ''


@dataclass
class TemplateDice:
    """
    Fragmentos do arquivo DICE original necessários para reescrevê-lo,
    trocando apenas a carga de cada linha do bloco do soluto (molécula 1) e
    preservando o restante do arquivo sem qualquer alteração.
    """
    cabecalho: str            # linhas 1–3 (rule, m, 'na(1) título'), verbatim
    linhas_atomo: list[str]   # uma linha crua por átomo do soluto, verbatim
    spans_carga: list[tuple[int, int]]  # (início, fim) do campo de carga em cada linha
    cauda: str                # tudo depois do bloco do soluto, verbatim


# ---------------------------------------------------------------------------
# Leitura do arquivo DICE (ljname) de entrada
# ---------------------------------------------------------------------------

def le_dice_inicial(caminho: str) -> tuple[list[Atomo], TemplateDice]:
    """
    Lê o arquivo DICE (ljname) usado na simulação inicial do ciclo e extrai
    as cargas do soluto (bloco da primeira molécula).
    Retorna (lista de Atomo, TemplateDice com os fragmentos para reconstrução
    do arquivo na saída).
    """
    with open(caminho, 'r', encoding='utf-8-sig') as f:
        linhas = f.readlines()  # mantém o '\n' de cada linha

    if len(linhas) < 3:
        raise ValueError(f"{caminho}: arquivo DICE incompleto (menos de 3 linhas).")

    linha_regra, linha_m, linha_bloco1 = linhas[0], linhas[1], linhas[2]

    partes_bloco1 = linha_bloco1.split(None, 1)
    if not partes_bloco1:
        raise ValueError(f"{caminho}:3: linha do bloco do soluto vazia.")
    try:
        na1 = int(partes_bloco1[0])
    except ValueError:
        raise ValueError(
            f"{caminho}:3: esperado 'na(1) título', não foi possível ler "
            f"na(1) em {linha_bloco1!r}"
        ) from None

    inicio_atomos = 3
    fim_atomos = 3 + na1
    if len(linhas) < fim_atomos:
        raise ValueError(
            f"{caminho}: arquivo indica {na1} átomos no soluto, mas só há "
            f"{len(linhas) - inicio_atomos} linhas após o cabeçalho."
        )

    linhas_atomo = linhas[inicio_atomos:fim_atomos]
    cauda = ''.join(linhas[fim_atomos:])

    atomos: list[Atomo] = []
    spans_carga: list[tuple[int, int]] = []

    for i, linha in enumerate(linhas_atomo, start=1):
        tokens = list(re.finditer(r'\S+', linha))
        if len(tokens) < 8:
            raise ValueError(
                f"{caminho}: linha do átomo {i} do soluto tem "
                f"{len(tokens)} campos, esperado 8 "
                f"(label an x y z q epsilon sigma): {linha!r}"
            )

        try:
            an = int(tokens[1].group())
        except ValueError:
            raise ValueError(
                f"{caminho}: átomo {i}: número atômico inválido "
                f"{tokens[1].group()!r}"
            ) from None
        simbolo = simbolo_do_z(an)

        tok_q = tokens[5]
        try:
            carga = float(tok_q.group())
        except ValueError:
            raise ValueError(
                f"{caminho}: átomo {i}: carga inválida {tok_q.group()!r}"
            ) from None

        atomos.append(Atomo(i, simbolo, carga))
        spans_carga.append((tok_q.start(), tok_q.end()))

    if not atomos:
        raise ValueError(f"{caminho}: nenhum átomo lido no bloco do soluto.")

    template = TemplateDice(
        cabecalho=linha_regra + linha_m + linha_bloco1,
        linhas_atomo=linhas_atomo,
        spans_carga=spans_carga,
        cauda=cauda,
    )
    return atomos, template


# ---------------------------------------------------------------------------
# Leitura do TXT de médias
# ---------------------------------------------------------------------------

def le_medias(caminho: str) -> tuple[dict[int, tuple[str, float]], str]:
    """
    Lê o TXT de médias gerado pelo scriptLog.py.

    Retorna (dict num -> (símbolo, carga média), modelo detectado ou '').
    """
    nome = Path(caminho).name.lower()
    if 'nbo' in nome:
        modelo = 'NBO'
    elif 'esp' in nome or 'chelpg' in nome:
        modelo = 'ESP'
    else:
        modelo = ''

    def para_float(bruto: str) -> float:
        return float('nan') if bruto.lower() == 'nan' else float(bruto)

    medias: dict[int, tuple[str, float]] = {}
    num_implicito = 0
    with open(caminho, 'r', encoding='utf-8-sig') as f:
        for n_linha, linha in enumerate(f, 1):
            linha = linha.strip()
            if not linha:
                continue
            partes = linha.split()

            if len(partes) >= 3:
                # 'num simbolo media' — se o 1º campo não for inteiro, é cabeçalho.
                try:
                    num = int(partes[0])
                except ValueError:
                    continue
                simbolo, bruto = partes[1], partes[2]
            elif len(partes) == 2:
                # 'simbolo media' — se o 2º campo não for numérico, é cabeçalho.
                simbolo, bruto = partes
                try:
                    para_float(bruto)
                except ValueError:
                    continue
                num_implicito += 1
                num = num_implicito
            else:
                raise ValueError(
                    f"{caminho}:{n_linha}: linha em formato inesperado: {linha!r}")

            try:
                valor = para_float(bruto)
            except ValueError:
                raise ValueError(
                    f"{caminho}:{n_linha}: carga inválida {bruto!r}") from None

            medias[num] = (simbolo, valor)

    return medias, modelo


# ---------------------------------------------------------------------------
# Equivalência química
# ---------------------------------------------------------------------------

def parse_selecao(texto: str, n_max: int) -> set[int]:
    """
    Converte uma seleção textual em conjunto de números de átomo.

        "1-19"          -> {1, ..., 19}
        "1,4,9"         -> {1, 4, 9}
        "1-13,15,17-19" -> {1..13, 15, 17, 18, 19}
        "todos" ou ""   -> {1, ..., n_max}
    """
    texto = texto.strip().lower()
    if texto in ('', 'todos', 'all', '*'):
        return set(range(1, n_max + 1))

    selecionados: set[int] = set()
    for parte in texto.split(','):
        parte = parte.strip()
        if not parte:
            continue
        if '-' in parte:
            ini, fim = parte.split('-', 1)
            ini, fim = int(ini.strip()), int(fim.strip())
            if ini > fim:
                ini, fim = fim, ini
            selecionados.update(range(ini, fim + 1))
        else:
            selecionados.add(int(parte))

    fora = {n for n in selecionados if n < 1 or n > n_max}
    if fora:
        raise ValueError(f"Átomos fora do intervalo 1-{n_max}: {sorted(fora)}")
    return selecionados


def le_grupos(caminho: str, n_atomos: int) -> dict[str, list[int]]:
    """
    Lê o arquivo de grupos de equivalência química declarados pelo usuário.

    Cada linha útil tem a forma 'rótulo : seleção', com a seleção na mesma
    sintaxe aceita por parse_selecao. Linhas vazias e o que vier depois de '#'
    são ignorados. Exemplo:

        Pt          : 1
        Cl          : 2,3
        O_hidroxo   : 4,6
        H_amino     : 9,11-13

    A numeração é a do bloco do soluto no arquivo DICE.

    A partição precisa ser total e exclusiva: todo átomo do soluto em exatamente
    um grupo. Falta ou repetição levanta ValueError, para que um erro de
    digitação pare o script em vez de gerar um agrupamento silenciosamente
    errado.

    Retorna dict rótulo -> lista ordenada de números de átomo.
    """
    grupos: dict[str, list[int]] = {}
    vistos: dict[int, str] = {}

    with open(caminho, 'r', encoding='utf-8-sig') as f:
        for n_linha, linha in enumerate(f, 1):
            linha = linha.split('#', 1)[0].strip()
            if not linha:
                continue
            if ':' not in linha:
                raise ValueError(f"{caminho}:{n_linha}: falta ':' em {linha!r}")

            rotulo, selecao = (p.strip() for p in linha.split(':', 1))
            if not rotulo:
                raise ValueError(f"{caminho}:{n_linha}: rótulo vazio")
            if rotulo in grupos:
                raise ValueError(f"{caminho}:{n_linha}: rótulo '{rotulo}' repetido")

            try:
                nums = sorted(parse_selecao(selecao, n_atomos))
            except ValueError as e:
                raise ValueError(f"{caminho}:{n_linha}: {e}") from None
            if not nums:
                raise ValueError(f"{caminho}:{n_linha}: grupo '{rotulo}' vazio")

            for num in nums:
                if num in vistos:
                    raise ValueError(
                        f"{caminho}:{n_linha}: átomo {num} já está em "
                        f"'{vistos[num]}' e reaparece em '{rotulo}'"
                    )
                vistos[num] = rotulo

            grupos[rotulo] = nums

    if not grupos:
        raise ValueError(f"{caminho}: nenhum grupo declarado.")

    faltando = sorted(set(range(1, n_atomos + 1)) - vistos.keys())
    if faltando:
        raise ValueError(f"{caminho}: átomos sem grupo: {faltando}")

    return grupos


def aplica_grupos(atomos: list[Atomo],
                  mapa: dict[str, list[int]]) -> dict[str, list[Atomo]]:
    """Converte o mapa de números lido do arquivo em grupos de objetos Atomo."""
    por_num = {a.num: a for a in atomos}
    grupos: dict[str, list[Atomo]] = {}

    for rotulo, nums in mapa.items():
        ausentes = [n for n in nums if n not in por_num]
        if ausentes:
            raise ValueError(
                f"Grupo '{rotulo}': átomos {ausentes} não existem no arquivo DICE."
            )
        membros = [por_num[n] for n in nums]
        for a in membros:
            a.grupo = rotulo
        grupos[rotulo] = membros

    return grupos


# ---------------------------------------------------------------------------
# Ajuste
# ---------------------------------------------------------------------------

def ajusta_residuo(grupos: dict[str, list[Atomo]], selecionados: set[int],
                   alvo: float = 0.0, criterio_grupo: str = 'media') -> dict:
    """
    Distribui o resíduo entre os grupos, proporcionalmente à largura do
    intervalo entre a carga inicial e a carga média.

    O par (q_inicial_g, q_media_g) que representa cada grupo pode ser
    definido de duas formas (parâmetro criterio_grupo):
      - 'media' (padrão)
      - 'menor_largura'
    """
    if criterio_grupo not in ('media', 'menor_largura'):
        raise ValueError(
            f"criterio_grupo deve ser 'media' ou 'menor_largura', "
            f"recebido {criterio_grupo!r}")

    linhas_relatorio = []
    soma_n_d = 0.0
    soma_n_qmed = 0.0
    soma_real_qmed = 0.0

    for rotulo, membros in grupos.items():
        n_g = len(membros)
        soma_real_qmed += sum(a.carga_media for a in membros)

        atomo_ref = None
        if criterio_grupo == 'menor_largura':
            atomo_ref = min(membros,
                            key=lambda a: abs(a.carga_media - a.carga_inicial))
            q_ini_g = atomo_ref.carga_inicial
            q_med_g = atomo_ref.carga_media
        else:
            q_ini_g = sum(a.carga_inicial for a in membros) / n_g
            q_med_g = sum(a.carga_media for a in membros) / n_g

        participa = all(a.num in selecionados for a in membros)
        if any(a.num in selecionados for a in membros) and not participa:
            raise ValueError(
                f"Seleção quebra o grupo equivalente '{rotulo}' "
                f"(átomos {[a.num for a in membros]}). Selecione o grupo inteiro."
            )

        largura = abs(q_med_g - q_ini_g) if participa else 0.0

        soma_n_d += n_g * largura
        soma_n_qmed += n_g * q_med_g

        linhas_relatorio.append({
            'grupo': rotulo, 'n': n_g, 'membros': [a.num for a in membros],
            'q_ini': q_ini_g, 'q_med': q_med_g,
            'largura': largura, 'participa': participa,
            'atomo_ref': atomo_ref.num if atomo_ref else None,
        })

    residuo = alvo - soma_n_qmed

    if abs(soma_n_d) < 1e-12:
        raise ValueError(
            "Soma das larguras é nula — não há como distribuir o resíduo. "
            "Verifique se as cargas médias diferem das iniciais."
        )

    fator = residuo / soma_n_d

    for item in linhas_relatorio:
        item['desloc'] = fator * item['largura']
        item['q_fim'] = item['q_med'] + item['desloc']
        for a in grupos[item['grupo']]:
            a.carga_final = item['q_fim']

    return {
        'linhas': linhas_relatorio,
        'residuo': residuo,
        'fator': fator,
        'soma_larguras': soma_n_d,
        'criterio_grupo': criterio_grupo,
        'soma_n_qmed': soma_n_qmed,
        'soma_real_qmed': soma_real_qmed,
        'alvo': alvo,
    }


# ---------------------------------------------------------------------------
# Saídas
# ---------------------------------------------------------------------------

def imprime_tabela(relatorio: dict, atomos: list[Atomo], modelo: str,
                   decimais: int = 6) -> None:
    """Imprime a tabela por grupo e a tabela por átomo."""
    d = decimais
    larg = d + 6

    criterio = relatorio.get('criterio_grupo', 'media')
    refs = {it['atomo_ref'] for it in relatorio['linhas'] if it.get('atomo_ref')}

    print(f"\n{'=' * 92}")
    print(f"AJUSTE DE RESÍDUO — modelo {modelo or '?'} — critério: {criterio}")
    print('=' * 92)

    if criterio == 'menor_largura':
        print("Cada grupo adota a carga do seu átomo-referência: aquele com o")
        print("menor |q_media - q_inicial| individual dentro do grupo.")

    diff_soma = relatorio['soma_n_qmed'] - relatorio['soma_real_qmed']

    print("RESÍDUO")
    print(f"  Soma das q_media, já igualadas por grupo : "
          f"{relatorio['soma_n_qmed']:+.{d}f}")
    print(f"  Carga total desejada (alvo)              : "
          f"{relatorio['alvo']:+.{d}f}")
    print(f"  Resíduo a distribuir                     : "
          f"{relatorio['residuo']:+.{d}f}")

    print("\nDIFERENÇA DAS SOMAS")
    print(f"  Soma das q_media igualadas por grupo     : "
          f"{relatorio['soma_n_qmed']:+.{d}f}")
    print(f"  Soma real das {len(atomos):<2} q_media medidas         : "
          f"{relatorio['soma_real_qmed']:+.{d}f}")
    print(f"  Diferença                                : "
          f"{diff_soma:+.{d}f}")

    print(f"\nCOMO O RESÍDUO É REPARTIDO")
    print(f"  Soma das larguras : {relatorio['soma_larguras']:.{d}f}")
    print(f"  Resíduo a distribuir: {relatorio['residuo']:+.{d}f}")
    print(f"  Fator f           : {relatorio['fator']:+.{d}f}"
          f"   (cada grupo anda {abs(relatorio['fator']) * 100:.2f}% da sua largura)")

    print(f"\n{'-' * 92}")
    cab_ref = '   átomo-ref' if criterio == 'menor_largura' else ''
    print(f"{'Grupo':<14}{'n':>3}{'q_inicial':>{larg}}{'q_media':>{larg}}"
          f"{'largura':>{larg}}{'desloc':>{larg}}{'q_final':>{larg}}{cab_ref}")
    print('-' * 92)
    for it in relatorio['linhas']:
        marca = '' if it['participa'] else '  (fora)'
        ref = f"{it['atomo_ref']:>10}" if it.get('atomo_ref') else ''
        print(f"{it['grupo']:<14}{it['n']:>3}"
              f"{it['q_ini']:>{larg}.{d}f}{it['q_med']:>{larg}.{d}f}"
              f"{it['largura']:>{larg}.{d}f}{it['desloc']:>+{larg}.{d}f}"
              f"{it['q_fim']:>{larg}.{d}f}{ref}{marca}")

    print(f"\n{'Nº':<4}{'Símb':<6}{'Grupo':<14}"
          f"{'q_inicial':>{larg}}{'q_media':>{larg}}{'q_final':>{larg}}   usada?")
    print('-' * 92)
    for a in sorted(atomos, key=lambda x: x.num):
        if criterio == 'menor_largura':
            usada = '   <-- ref' if a.num in refs else '       ---'
        else:
            usada = '     media'
        print(f"{a.num:<4}{a.simbolo:<6}{a.grupo:<14}"
              f"{a.carga_inicial:>{larg}.{d}f}{a.carga_media:>{larg}.{d}f}"
              f"{a.carga_final:>{larg}.{d}f}{usada}")

    if criterio == 'menor_largura':
        print("\n  '<-- ref': a q_media DESTE átomo foi a usada para todo o grupo.")
        print("  '---'    : a q_media deste átomo foi ignorada no cálculo")
        print("             (aparece só para conferência).")
    else:
        print("\n  'media': a q_media usada foi a média aritmética do grupo,")
        print("           não a de nenhum átomo em particular.")

    soma = sum(a.carga_final for a in atomos)
    print('-' * 92)
    print(f"Soma das q_final sobre os {len(atomos)} átomos: {soma:+.2e}  "
          f"(fecha no alvo)")
    if abs(soma) > 1e-8:
        print("[AVISO] Somatório não fechou em zero. Confira a seleção de átomos.")


def escreve_dice_final(caminho_saida: str, atomos: list[Atomo],
                       template: TemplateDice, decimais: int = 6) -> None:
    """
    Grava uma cópia do arquivo DICE de entrada com as cargas do soluto
    substituídas pelas cargas ajustadas (carga_final).

    Apenas o campo de carga (q) de cada linha do bloco do soluto é alterado;
    todos os demais campos (label, número atômico, x, y, z, epsilon, sigma),
    o espaçamento original, o título, o bloco do solvente e o '$end' são
    preservados exatamente como no arquivo de entrada.
    """
    novas_linhas = []
    for atomo, linha_original, (ini, fim) in zip(
            atomos, template.linhas_atomo, template.spans_carga):
        nova_carga = f"{atomo.carga_final:.{decimais}f}"
        nova_linha = linha_original[:ini] + nova_carga + linha_original[fim:]
        novas_linhas.append(nova_linha)

    conteudo = template.cabecalho + ''.join(novas_linhas) + template.cauda

    with open(caminho_saida, 'w', encoding='utf-8') as f:
        f.write(conteudo)

    soma = sum(Decimal(f"{a.carga_final:.{decimais}f}") for a in atomos)
    print(f"Soma das {len(atomos)} cargas gravadas: {soma}")
    if abs(soma) > Decimal('1e-5'):
        print("[AVISO] O arredondamento deixou resíduo acima de 1e-5 no arquivo. "
              "Aumente --decimais se for sensível a isso.")


def escreve_txt_cargas(caminho_saida: str, atomos: list[Atomo],
                       decimais: int = 6) -> None:
    """
    Grava um TXT com duas colunas separadas por um espaço, uma linha por
    átomo do soluto, na mesma ordem da tabela:

        q_media q_final

    Formato pensado para leitura direta (numpy.loadtxt, awk, planilha) e
    para servir de registro do ciclo: q_media é a média sobre as
    configurações Monte Carlo e q_final é a carga após a redistribuição do
    resíduo, isto é, o que vai para o arquivo DICE do ciclo seguinte.
    """
    with open(caminho_saida, 'w', encoding='utf-8') as f:
        for a in atomos:
            f.write(f"{a.carga_media:.{decimais}f} {a.carga_final:.{decimais}f}\n")

    print(f"TXT com {len(atomos)} pares (q_media q_final) gravado.")


def sugere_nome_txt(modelo: str) -> str:
    """
    Sugere um nome para o TXT de cargas, seguindo a convenção interna do
    projeto (ex.: cargas_ESP_media_final.txt). Só um ponto de partida;
    o caminho pode ser digitado livremente quando solicitado.
    """
    sufixo = modelo if modelo else 'ajustado'
    return f"cargas_{sufixo}_media_final.txt"


def sugere_nome_saida(caminho_iniciais: str, template: TemplateDice,
                      modelo: str) -> str:
    """
    Sugere um nome para o arquivo DICE atualizado, seguindo o padrão já usado
    no projeto (ex.: s1Pt001_ESP.txt / s1Pt001_NBO.txt): prefixo de ciclo
    's1', o código do complexo (extraído do título do bloco do soluto, se
    reconhecível) e o modelo de carga.

    Convenção interna do projeto, não uma exigência do DICE nem da
    literatura — é só um ponto de partida; o caminho pode ser digitado
    livremente quando solicitado.
    """
    m = re.search(r'\b([A-Za-z]+\d+)\b', template.cabecalho.splitlines()[-1])
    complexo = m.group(1) if m else Path(caminho_iniciais).stem
    sufixo = modelo if modelo else 'ajustado'
    return f"s1{complexo}_{sufixo}.txt"


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Ajusta o resíduo de carga do ciclo de polarização iterativa."
    )
    parser.add_argument('--iniciais',
                        help="Arquivo DICE (ljname) da simulação inicial do "
                             "ciclo, ex.: s0Pt001DFT.txt")
    parser.add_argument('--medias',
                        help="TXT de médias gerado pelo scriptLog.py "
                             "(esp_cargas_medias.txt / nbo_...)")
    parser.add_argument('--atomos', help="Seleção: '1-19', '1,4,9' ou 'todos'")
    modo_grupos = parser.add_mutually_exclusive_group()
    modo_grupos.add_argument('--grupos',
                             help="Arquivo de grupos de equivalência química")
    modo_grupos.add_argument('--sem-grupos', action='store_true',
                             help="Não impor equivalência (cada átomo isolado); "
                                  "use apenas como diagnóstico")
    parser.add_argument('--alvo', type=float, default=0.0,
                        help="Carga total desejada do soluto (padrão: 0)")
    parser.add_argument('--criterio-grupo', choices=['media', 'menor_largura'],
                        default=None,
                        help="Como definir (q_inicial, q_media) de cada grupo: "
                             "'media' = média aritmética dos membros; "
                             "'menor_largura' = usa o átomo do grupo com menor "
                             "|q_media - q_inicial| individual, evitando "
                             "cancelamento de sinal dentro do grupo. Se omitido, "
                             "é perguntado. Nenhum dos dois tem precedente "
                             "publicado.")
    parser.add_argument('--saida',
                        help="Caminho do arquivo DICE atualizado. Se informado, "
                             "o arquivo é gravado direto, sem perguntar "
                             "(uso em automação). Se omitido, a pergunta é "
                             "feita depois de mostrar a tabela.")
    parser.add_argument('--saida-txt',
                        help="Caminho do TXT com duas colunas 'q_media "
                             "q_final' (uma linha por átomo). Se informado, "
                             "o arquivo é gravado direto, sem perguntar "
                             "(uso em automação). Se omitido, a pergunta é "
                             "feita depois de mostrar a tabela.")
    parser.add_argument('--decimais', type=int, default=6,
                        help="Casas decimais na tabela e no arquivo (padrão: 6)")
    args = parser.parse_args()

    # --- Arquivo DICE inicial ---
    caminho_iniciais = args.iniciais
    while not caminho_iniciais or not Path(caminho_iniciais).is_file():
        if caminho_iniciais:
            print(f"  [ERRO] Arquivo não encontrado: '{caminho_iniciais}'")
        caminho_iniciais = input(
            "Arquivo DICE da simulação inicial (ljname): ").strip()

    try:
        atomos, template = le_dice_inicial(caminho_iniciais)
    except ValueError as e:
        print(f"[ERRO] {e}")
        sys.exit(1)
    print(f"\nLidos {len(atomos)} átomos do soluto em '{caminho_iniciais}'.")

    # --- TXT de médias ---
    caminho_medias = args.medias
    while not caminho_medias or not Path(caminho_medias).is_file():
        if caminho_medias:
            print(f"  [ERRO] Arquivo não encontrado: '{caminho_medias}'")
        caminho_medias = input(
            "TXT de médias (esp_cargas_medias.txt / nbo_...): ").strip()

    medias, modelo = le_medias(caminho_medias)
    print(f"Lidas {len(medias)} médias em '{caminho_medias}'"
          f"{f' — modelo {modelo}' if modelo else ''}.")

    faltando = [a.num for a in atomos if a.num not in medias]
    if faltando:
        print(f"[ERRO] Sem carga média para os átomos: {faltando}")
        sys.exit(1)

    for a in atomos:
        simbolo_medias, valor = medias[a.num]
        if simbolo_medias and simbolo_medias != a.simbolo:
            print(f"[AVISO] Átomo {a.num}: símbolo '{simbolo_medias}' no TXT de "
                  f"médias difere de '{a.simbolo}' no arquivo DICE.")
        a.carga_media = valor

    # --- Grupos ---
    if args.sem_grupos:
        caminho_grupos = None
    elif args.grupos:
        caminho_grupos = args.grupos
    else:
        # Nenhuma das duas flags foi passada: pergunta, em vez de assumir.
        caminho_grupos = input(
            "Arquivo de grupos de equivalência "
            "(ENTER para não agrupar): "
        ).strip() or None

    if caminho_grupos is None:
        grupos = {}
        for a in atomos:
            a.grupo = f"{a.simbolo}{a.num}"
            grupos[a.grupo] = [a]
        print("\n[AVISO] Sem arquivo de grupos: cada átomo é seu próprio grupo "
              "(sem equivalência química). Informe --grupos para agrupar.")
    else:
        while not Path(caminho_grupos).is_file():
            print(f"  [ERRO] Arquivo não encontrado: '{caminho_grupos}'")
            caminho_grupos = input(
                "Arquivo de grupos de equivalência "
                "(ENTER para não agrupar): "
            ).strip()
            if not caminho_grupos:
                caminho_grupos = None
                break

        if caminho_grupos is None:
            grupos = {}
            for a in atomos:
                a.grupo = f"{a.simbolo}{a.num}"
                grupos[a.grupo] = [a]
            print("\n[AVISO] Sem arquivo de grupos: cada átomo é seu próprio "
                  "grupo (sem equivalência química).")
        else:
            try:
                mapa = le_grupos(caminho_grupos, len(atomos))
                grupos = aplica_grupos(atomos, mapa)
            except ValueError as e:
                print(f"[ERRO] {e}")
                sys.exit(1)
            print(f"\nLidos {len(grupos)} grupos em '{caminho_grupos}'.")

    print(f"\nGrupos de equivalência ({len(grupos)}):")
    for rotulo, membros in grupos.items():
        print(f"  {rotulo:<18} -> átomos {[a.num for a in membros]}")

    # --- Seleção ---
    texto_sel = args.atomos
    if texto_sel is None:
        try:
            texto_sel = input(
                f"\nQuais átomos entram no ajuste? "
                f"(ex.: '1-{len(atomos)}', '1,4,9', ENTER = todos): "
            ).strip()
        except EOFError:
            texto_sel = ''

    try:
        selecionados = parse_selecao(texto_sel, len(atomos))
    except ValueError as e:
        print(f"[ERRO] {e}")
        sys.exit(1)

    # --- Critério de grupo ---
    criterio = args.criterio_grupo
    if criterio is None:
        if len(grupos) == len(atomos):
            # Cada átomo é seu próprio grupo: os dois critérios coincidem,
            # não faz sentido perguntar.
            criterio = 'media'
        else:
            print("\nCritério para definir a carga de cada grupo equivalente:")
            print("  1 - media          (média aritmética dos átomos do grupo)")
            print("  2 - menor largura  (átomo do grupo com menor "
                  "|q_media - q_inicial|)")
            escolha = ''
            while escolha not in ('1', '2'):
                try:
                    escolha = input("Escolha [1/2]: ").strip()
                except EOFError:
                    escolha = '1'
                if escolha not in ('1', '2'):
                    print("  [ERRO] Digite 1 ou 2.")
            criterio = 'media' if escolha == '1' else 'menor_largura'

    # --- Ajuste ---
    try:
        relatorio = ajusta_residuo(grupos, selecionados, args.alvo, criterio)
    except ValueError as e:
        print(f"[ERRO] {e}")
        sys.exit(1)

    imprime_tabela(relatorio, atomos, modelo, args.decimais)

    # --- Saída: arquivo DICE atualizado ---
    caminho_saida = args.saida
    if caminho_saida is None:
        try:
            resposta = input(
                "\nGostaria de gerar o arquivo DICE atualizado? [s/N]: "
            ).strip().lower()
        except EOFError:
            resposta = ''

        if resposta in ('s', 'sim', 'y', 'yes'):
            sugestao = sugere_nome_saida(caminho_iniciais, template, modelo)
            try:
                caminho_saida = input(
                    f"Caminho do arquivo DICE atualizado [{sugestao}]: "
                ).strip() or sugestao
            except EOFError:
                caminho_saida = sugestao
        else:
            caminho_saida = None

    if caminho_saida:
        escreve_dice_final(caminho_saida, atomos, template, args.decimais)
        print(f"Salvo: {caminho_saida}")
    else:
        print("\nArquivo DICE não gerado.")

    # --- Saída: TXT com q_media e q_final ---
    caminho_txt = args.saida_txt
    if caminho_txt is None:
        try:
            resposta = input(
                "\nGostaria de gerar o TXT com as cargas médias e finais? [s/N]: "
            ).strip().lower()
        except EOFError:
            resposta = ''

        if resposta in ('s', 'sim', 'y', 'yes'):
            sugestao = sugere_nome_txt(modelo)
            try:
                caminho_txt = input(
                    f"Caminho do TXT [{sugestao}]: "
                ).strip() or sugestao
            except EOFError:
                caminho_txt = sugestao
        else:
            caminho_txt = None

    if caminho_txt:
        escreve_txt_cargas(caminho_txt, atomos, args.decimais)
        print(f"Salvo: {caminho_txt}")
    else:
        print("\nTXT de cargas não gerado.")


if __name__ == "__main__":
    main()