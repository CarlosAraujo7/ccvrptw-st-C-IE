#!/usr/bin/env python3
"""
Tabela 8 do artigo, "Which workflow the results support in each situation" (17/11/2026, diário Entrada 151; plano da terceira revisão para chegar a 4,0:
tabela de decisão "quando os buffers bastam").

Tabela DERIVADA: nenhuma rodada nova e nenhum número digitado à mão. O script lê as tabelas do pacote (que o teste de regeneração já garante serem iguais
aos resultados guardados), a saber as Tabelas 4 e 7 do texto e as S10, S16, S18, S19, S20, S26, S27, S31 e S33 do suplemento, tira de cada uma as células
citadas e monta a tabela de decisão por situação. O relatório lista cada número com a célula de origem. Se uma tabela fonte mudar de forma, o script
para com uma mensagem (as buscas conferem rótulos de linha e de coluna).
Pasta das tabelas do pacote: --tabelas, ou a variável CCVRPTW_PACOTE_TABLES, ou o padrão do teste de regeneração.
Uso: python pesquisa/tabela_decisao.py [--tabelas pasta] [--tex arquivo.tex] [--out pasta]
Saída: pesquisa/resultados_decisao/relatorio_tabela_decisao.txt e a tabela LaTeX (--tex; padrão tab_decisao.tex na pasta de saída).
"""
import argparse, math, os, re

BASE = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
PACOTE_PADRAO = os.environ.get('CCVRPTW_PACOTE_TABLES', r'C:\Users\vinic\OneDrive\Documentos\Revisão Cold Chain\tables')


def corpo(tex):
    return re.search(r'\\begin\{tabular\}\{.*?\}\s*\n(.*?)\\end\{tabular\}', tex, re.S).group(1)


def linhas(tex):
    """linhas de dados da primeira tabular, cada uma como lista de células com o texto cru"""
    c = corpo(tex)
    c = re.sub(r'\\(?:toprule|midrule|bottomrule|addlinespace)(?:\[[^\]]*\])?', ' ', c)
    c = re.sub(r'\\cmidrule(?:\([a-z]+\))?\{[^}]*\}', ' ', c)
    return [[x.strip() for x in ln.split('&')] for ln in c.split('\\\\') if ln.strip()]


def nums(cel):
    """todos os números de uma célula, sem os marcadores de significância e sem o argumento de \\multicolumn"""
    c = re.sub(r'\\multicolumn\{\d+\}\{[^}]*\}', '', cel)
    c = re.sub(r'\$\^[*a]\$|\^[*a]', '', c)
    return [float(x) for x in re.findall(r'[-+]?\d+(?:\.\d+)?', c)]


def num(cel):
    v = nums(cel)
    assert v, f'célula sem número: {cel!r}'
    return v[0]


def meio(x):
    """arredonda para inteiro com meio para cima"""
    return int(math.floor(x + 0.5))


class Fontes:
    def __init__(self):
        self.itens = []

    def __call__(self, valor, origem):
        self.itens.append((valor, origem))
        return valor


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--tabelas', default=PACOTE_PADRAO)
    ap.add_argument('--tex', default=None)
    ap.add_argument('--out', default=os.path.join(BASE, 'pesquisa', 'resultados_decisao'))
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    def ler(nome):
        return linhas(open(os.path.join(a.tabelas, nome + '.tex'), encoding='utf-8').read())

    F = Fontes()

    # ---- S26: o plano nominal do PyVRP passa na certificação? (por instâncias, incerteza e família) ----
    s26 = {}
    inst = None
    for r in ler('tab_pyvrp_nominal_cert'):
        if len(r) < 7 or not re.match(r'^\d', r[2]):
            continue
        inst = r[0] or inst
        s26[(inst, r[1])] = dict(cells=num(r[2]), cert=num(r[3]), estreita=nums(r[5])[0] if nums(r[5]) else None, larga=nums(r[6])[0] if nums(r[6]) else None)
    bm, bh = s26[('benchmark', 'medium')], s26[('benchmark', 'high')]
    sp_m, sp_h = s26[(r'S\~ao Paulo', 'medium')], s26[(r'S\~ao Paulo', 'high')]
    n400h = s26[('$n=400$', 'high')]
    n200h = s26[('$n=200$', 'high')]
    assert bh['estreita'] == n200h['estreita'], 'a fração com janela estreita sob incerteza alta deveria ser a mesma no benchmark e com 200 clientes'

    # ---- S27: preço dos buffers onde o plano nominal falha (benchmark, incerteza alta) ----
    s27 = None
    inst = None
    for r in ler('tab_pyvrp_nominal_preco'):
        if len(r) < 7 or not re.match(r'^\d', r[2]):
            continue
        inst = r[0] or inst
        if inst == 'benchmark' and r[1] == 'high':
            s27 = dict(falham=num(r[2]), cal_media=nums(r[4])[0], cal_mediana=nums(r[4])[1], buf_media=nums(r[6])[0], buf_mediana=nums(r[6])[1])
    assert s27, 'S27: linha benchmark/high não encontrada'
    # quantas das células que falham (S27) têm janela estreita: a divisão do S26 em famílias estreitas e largas tem metade das células cada uma
    # (3 famílias x 3 instâncias x 5 sementes = 45); a hipótese é checada pela integralidade das contagens e pelo total de falhas do S27
    meia = bh['cells'] / 2
    cont = lambda pct: meia * pct / 100                                  # noqa: E731
    assert meia == int(meia) and all(abs(cont(p) - round(cont(p))) < 0.05 for p in (bh['estreita'], bh['larga'])), 'S26: a metade das células não dá contagens inteiras'
    falham_est, falham_larga = int(meia - round(cont(bh['estreita']))), int(meia - round(cont(bh['larga'])))
    assert falham_est + falham_larga == s27['falham'], 'S26 e S27 discordam sobre quantas células de incerteza alta falham'

    # ---- Tabela 7 (escala): cumprimento com 400 clientes ----
    esc = {}
    for r in ler('tab_escala'):
        if len(r) >= 10 and re.search(r'\d', r[4]):
            esc[r[0]] = nums(r[8])[0]
    acal400, calr400 = esc[r'\mACAL{}'], esc[r'\mCAL{}+R']
    pyv400 = esc[r'PyVRP, BUF$_{0.95}$+C+R']

    # ---- S18 (regra do limite inferior): 200 a 400 clientes ----
    s18 = None
    for r in ler('tab_regra_lcb'):
        if r and r[0].startswith('200--400'):
            s18 = dict(acalr_lcb=nums(r[3])[1], custo_acalr=num(r[7]))
    assert s18, 'S18: linha de 200 a 400 clientes não encontrada'

    # ---- S33 (comparações pareadas) e S16 (Atenas): buffer 0,95 e CAL em Atenas ----
    ath = {}
    dentro = False
    for r in ler('tab_casos_pareadas'):
        if len(r) == 1:
            dentro = r[0].startswith(r'\multicolumn{8}{@{}l}{\textit{Athens')
            continue
        if dentro and r[0] in (r'\mBUF{0.95}', r'\mCAL{}'):
            ath[r[0]] = dict(compl=num(r[1]), custo=num(r[2]), p=num(r[3]))
    assert len(ath) == 2, 'S33: bloco de Atenas incompleto'
    assert ath[r'\mCAL{}']['p'] >= 0.05, 'a diferença de custo do CAL em Atenas deixou de ser não significativa; revisar a frase'
    s16_cal = None
    for r in ler('tab_caso_atenas'):
        if r[0] == r'\mCAL{}':
            s16_cal = num(r[2])
    assert s16_cal == 100.0 and ath[r'\mCAL{}']['compl'] == 100.0, 'CAL em Atenas deixou de cumprir em todas as execuções'

    # ---- S10 (congestão ao longo do dia) ----
    s10 = {}
    plan = None
    for r in ler('tab_comparar_td'):
        if len(r) < 6 or not re.match(r'^[A-Za-z\\]', r[1]) or r[0] == 'Planned with' or r[0].startswith(' &') or r[1] == 'Method':
            continue
        plan = r[0] or plan
        s10.setdefault(plan, {})[r[1]] = dict(med=num(r[2]), alta=num(r[4]), dcm=(nums(r[3])[0] if nums(r[3]) else None),
                                              dca=(nums(r[5])[0] if nums(r[5]) else None))
    media_dia = next(v for k, v in s10.items() if k.startswith('Daily average'))
    dependente = next(v for k, v in s10.items() if k.startswith('Time-dependent'))
    td_met = [dependente[m] for m in (r'\mBUF{0.95}', r'\mCAL{}', r'\mSIM{}')]
    dif_custo = max(abs(x) for k, v in s10.items() if not k.startswith('Time-dependent') for m in v.values() for x in (m['dcm'], m['dca']) if x is not None)

    # ---- S20 (variabilidade observada): intervalos sobre todos os métodos e conjuntos ----
    taxi, truck, m3 = [], [], []
    for r in ler('tab_observado'):
        if len(r) == 7 and re.match(r'^\d', r[1]):
            taxi.append(num(r[3])); truck.append(num(r[5])); m3.append(num(r[6]))

    # ---- S31 (CV estimado) ----
    s31 = {}
    for r in ler('tab_b4'):
        if len(r) >= 7 and r[1] == 'all':
            s31['assumed' if r[0].startswith('Assumed') else 'estimated'] = num(r[3])
    assert set(s31) == {'assumed', 'estimated'}, 'S31: linhas do CV suposto e estimado não encontradas'

    # ---- S19 (arrependimento da escolha pelo valor nominal) ----
    regret = max(num(r[3]) for r in ler('tab_parteb') if len(r) == 6 and re.match(r'^\d', r[1]))

    # ---- Tabela 4 (refinamento): queda de perda de qualidade com +R ----
    queda, custo_ref = [], []
    for r in ler('tab_refino'):
        if len(r) == 6 and '+R' in r[0] and re.match(r'^\d', r[1]):
            queda.append(abs(num(r[4]))); custo_ref.append((num(r[3]), '*' in r[3]))
    assert queda and all(v <= 0 for v, sig in custo_ref if sig), 'o refinamento teria aumentado o custo de forma significativa; revisar a frase'

    # ---- tabela ----
    P = lambda x: str(meio(x))                                        # noqa: E731
    d1 = lambda x: f'{x:.1f}'                                        # noqa: E731
    assert min(s27['cal_media'], s27['cal_mediana'], s27['buf_media'], s27['buf_mediana']) > 0, 'S27: custo do buffer deixou de ser positivo; revisar a frase'
    assert s18['custo_acalr'] >= 0, 'S18: o custo médio do ACAL+R com o limite inferior deixou de ser um aumento; revisar a frase'
    assert ath[r'\mBUF{0.95}']['custo'] > 0, 'S33: o buffer 0,95 deixou de custar mais em Atenas; revisar a frase'
    v = dict(
        a1=F(P(bm['cert']), 'S26 benchmark, medium, certificado'), a2=F(P(bh['larga']), 'S26 benchmark, high, janela larga'),
        a3=F(P(sp_h['cert']), 'S26 São Paulo, high'), a4=F(P(sp_m['cert']), 'S26 São Paulo, medium'),
        b1=F(P(bh['estreita']), 'S26 benchmark, high, janela estreita'), b2=F(P(n400h['estreita']), 'S26 n=400, high, janela estreita'),
        b3=F(d1(s27['cal_media']), 'S27 benchmark, high, CAL, onde falha, média'), b4=F(d1(s27['cal_mediana']), 'S27 CAL, onde falha, mediana'),
        b5=F(d1(s27['buf_media']), 'S27 BUF0,95, onde falha, média'), b6=F(d1(s27['buf_mediana']), 'S27 BUF0,95, onde falha, mediana'),
        b7=F(str(int(s27['falham'])), 'S27 benchmark, high, células em que o plano nominal falha'),
        b8=F(str(falham_est), 'S26 (contagem): células de janela estreita que falham, de b7'),
        c1=F(P(calr400), 'Tabela 7 CAL+R, n=400, cumprimento'), c2=F(P(acal400), 'Tabela 7 ACAL, n=400, cumprimento'),
        c3=F(P(pyv400), 'Tabela 7 PyVRP BUF0,95+C+R, n=400, cumprimento'), c4=F(P(s18['acalr_lcb']), 'S18 ACAL+R, regra LCB, 200-400, cumprimento'),
        c5=F(f"{s18['custo_acalr']:.2f}", 'S18 ACAL+R, 200-400, custo médio da regra LCB'),
        d1=F(d1(ath[r'\mBUF{0.95}']['custo']), 'S33 Atenas BUF0,95, custo'), d2=F(d1(ath[r'\mCAL{}']['custo']), 'S33 Atenas CAL, custo'),
        e1=F(P(min(m['med'] for m in media_dia.values())), 'S10 média do dia, medium, mínimo'),
        e2=F(P(max(m['med'] for m in media_dia.values())), 'S10 média do dia, medium, máximo'),
        e3=F(P(min(m['alta'] for m in media_dia.values())), 'S10 média do dia, high, mínimo'),
        e4=F(P(max(m['alta'] for m in media_dia.values())), 'S10 média do dia, high, máximo'),
        e5=F(P(min(min(m['med'], m['alta']) for m in td_met)), 'S10 dependente do horário (BUF0,95, CAL, SIM+C), mínimo'),
        e6=F(P(max(max(m['med'], m['alta']) for m in td_met)), 'S10 dependente do horário (BUF0,95, CAL, SIM+C), máximo'),
        e7=F(f'{dif_custo:.2f}', 'S10 maior diferença absoluta de custo entre modelos de planejamento'),
        f1=F(P(min(taxi)), 'S20 táxi, re-certificado, mínimo'), f2=F(P(max(taxi)), 'S20 táxi, re-certificado, máximo'),
        f3=F(P(min(truck)), 'S20 caminhão, re-certificado, mínimo'), f4=F(P(max(truck)), 'S20 caminhão, re-certificado, máximo'),
        f5=F(P(min(m3)), 'S20 M3, mínimo'), f6=F(P(max(m3)), 'S20 M3, máximo'),
        f7=F(P(s31['assumed']), 'S31 CV suposto, cumprimento M3'), f8=F(P(s31['estimated']), 'S31 CV estimado, cumprimento M3'),
        g1=F(f'{regret:.3f}', 'S19 arrependimento máximo'),
        h1=F(f'{min(queda):.2f}', 'Tabela 4 linhas +R, menor queda da perda de qualidade'), h2=F(f'{max(queda):.2f}', 'Tabela 4 linhas +R, maior queda'))
    linhas_tab = [
        ('Wide windows, or medium uncertainty, with a strong solver',
         'Solve with nominal times; certify; add a buffer only where certification fails.',
         f'The nominal plan of PyVRP is certified in {v["a1"]}\\% of the benchmark cells with medium uncertainty, in {v["a2"]}\\% of those with wide windows '
         f'under high uncertainty and in {v["a3"]}--{v["a4"]}\\% in S\\~ao Paulo (Table~S26).'),
        ('Narrow windows and high uncertainty',
         'Buffer the travel times; choose the level by certified cost over a ladder rather than fixing it; certify.',
         f'The nominal plan is certified in {v["b1"]}\\% of the benchmark cells ({v["b2"]}\\% with 400 customers). Of the {v["b7"]} high-uncertainty cells where it fails, '
         f'{v["b8"]} have narrow windows; repairing it costs {v["b3"]}\\% more with the ladder (median {v["b4"]}\\%) and {v["b5"]}\\% more with a fixed 0.95 buffer '
         f'(median {v["b6"]}\\%) (Tables~S26 and~S27).'),
        ('Hundreds of customers',
         'Adapt the buffer within one search, or take the 0.95 level with an external solver; certify against the lower confidence bound.',
         f'With 400 customers the point rule leaves \\mCAL{{}}+R compliant in {v["c1"]}\\% of the runs, \\mACAL{{}} in {v["c2"]}\\% and PyVRP with the 0.95 '
         f'buffer in {v["c3"]}\\% (Table~\\ref{{tab:scale}}); under the bound \\mACAL{{}}+R complies in {v["c4"]}\\%, at a mean cost increase of '
         f'{v["c5"]}\\% over the point rule (Table~S18).'),
        ('Loose deadlines',
         'Do not fix a high buffer; calibrate it down toward the nominal level.',
         f'In Athens the 0.95 buffer costs {v["d1"]}\\% more than \\mSIM{{}}, whereas \\mCAL{{}} complies in all runs and ties it in cost '
         f'(${v["d2"]}\\%$, not significant) (Tables~S16 and~S33).'),
        ('Congestion that varies over the day',
         'Plan with time-dependent travel times; buffers and simulation protect only against random variation.',
         f'Plans built with the daily average comply in {v["e1"]}--{v["e2"]}\\% (medium uncertainty) and {v["e3"]}--{v["e4"]}\\% (high) of the runs, plans built '
         f'with time-dependent times in {v["e5"]}--{v["e6"]}\\% (\\mBUF{{0.95}}, \\mCAL{{}} and \\mSIM{{}}), and the median cost differs by at most {v["e7"]}\\% '
         f'between planning models (Table~S10).'),
        ('Uncertain variability level',
         'Estimate it from observed trips and plan with it; certify also with observed scenarios.',
         f'Re-certifying with observed scenarios restores compliance to {v["f1"]}--{v["f2"]}\\% (taxi shape) and {v["f3"]}--{v["f4"]}\\% (truck shape), but with the '
         f'observed variability no method keeps the service level ({v["f5"]}--{v["f6"]}\\%); planning with the level estimated from the trips raises it '
         f'from {v["f7"]}\\% to {v["f8"]}\\% (Tables~S20 and~S31).'),
        ('Product quality',
         'Rank plans by their nominal quality loss; keep temperature simulation out of the search.',
         f'The regret of choosing by the nominal loss is at most {v["g1"]} percentage points in every condition tested, but the proxy weakens for '
         f'quality-threshold measures (Tables~S19 and~S25).'),
        ('After any of the above',
         'Refine the certified plan; keep a move only if simulated punctuality holds.',
         f'Quality loss falls by {v["h1"]} to {v["h2"]} percentage points without raising the cost (Table~\\ref{{tab:refine}}).'),
    ]
    L = [r'\begin{table}[htbp]', r'\centering',
         r'\caption{Which workflow the results support in each situation, with the evidence behind each row. Percentages are shares of runs or cells that meet '
         r'$\alpha=0.90$ out of sample unless stated. The nominal plan is the PyVRP solution computed with nominal travel times and departing at the start of the '
         r'horizon; the ladder is the cheapest certified plan over the buffer levels $q\in\{0.80,0.90,0.95\}$. The rows summarize the evidence for the regimes tested '
         r'(synthetic benchmarks, S\~ao Paulo and Athens); they are not rules.}',
         r'\label{tab:decision}', r'\footnotesize', r'\setlength{\tabcolsep}{3.5pt}',
         r'\begin{tabular}{@{}>{\raggedright\arraybackslash}p{0.16\textwidth} >{\raggedright\arraybackslash}p{0.27\textwidth} '
         r'>{\raggedright\arraybackslash}p{0.52\textwidth}@{}}', r'\toprule',
         r'Situation & What the results support & Evidence \\', r'\midrule']
    for k, (sit, fluxo, evid) in enumerate(linhas_tab):
        L.append(f'{sit} & {fluxo} & {evid} \\\\' + (r'\addlinespace[3pt]' if k < len(linhas_tab) - 1 else ''))
    L += [r'\bottomrule', r'\end{tabular}', r'\end{table}']
    destino = a.tex or os.path.join(a.out, 'tab_decisao.tex')
    open(destino, 'w', encoding='utf-8').write('\n'.join(L) + '\n')
    rel = ['TABELA 8 (decisão): número e célula de origem', ''] + [f'  {v:>8s}  <- {o}' for v, o in F.itens]
    open(os.path.join(a.out, 'relatorio_tabela_decisao.txt'), 'w', encoding='utf-8').write('\n'.join(rel) + '\n')
    print('\n'.join(rel)); print('tabela:', destino)


if __name__ == '__main__':
    main()
