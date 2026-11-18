#!/usr/bin/env python3
"""
Tabelas S33 e S34: comparações pareadas dos casos de São Paulo e de Atenas citadas nas Seções 6.10 e 6.11 (17/11/2026, diário Entrada 149; achado 4 da
terceira revisão cega: números do texto sem tabela no PDF).

Os números já existiam em relatórios de três scripts; esta tabela os reúne a partir dos MESMOS resultados guardados, com o mesmo cálculo (copiado,
para o script rodar sozinho no teste de regeneração, com `python -I`), e confere-se contra eles:
  - `analisar_refinamento.py --dir pesquisa/resultados_caso_*/refinamento` (visão do decisor: plano escolhido de cada método contra o SIM+C;
    orçamento igual, CAL(+R) contra a melhor de quatro execuções do SIM+C(+R));
  - `caso_sp/simc100_sp.py --analisar` (SIM+C com 100 cenários de busca, só São Paulo);
  - `caso_sp/tabela_caso_sp.py` (preço da confiabilidade: plano escolhido do SIM+C contra a solução mais barata do DET; máximo do DET+C).
Diferença pareada por instância (média das sementes), em % do SIM+C para custo e CO2 e em p.p. para a perda de qualidade e a probabilidade de
cumprir; estimativa de Hodges-Lehmann (HL) e p de Wilcoxon corrigido por Holm. Famílias de Holm: painel A, por medida, sobre todos os métodos do
relatório do refinamento (9 por conjunto); painéis B e C, sobre as quatro medidas de cada comparação.
Tabela S33 = comparações pareadas do plano escolhido contra o SIM+C (e contra o SIM+C com 100 cenários); tabela S34 = orçamento igual (painel A) e preço da
confiabilidade (painel B). As duas ficam separadas porque juntas passavam de uma página.
Uso: python pesquisa/tabela_casos_sp_atenas.py [--tex S33.tex] [--tex2 S34.tex] [--out pasta]
Saída: pesquisa/resultados_casos_sp_atenas/relatorio_tabela_casos.txt, tab_casos_pareadas.tex (S33) e tab_casos_orcamento.tex (S34).
"""
import argparse, glob, json, os
from collections import defaultdict
import numpy as np
from scipy.stats import norm, wilcoxon

BASE = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
PESQ = os.path.join(BASE, 'pesquisa')
ALPHA = 0.90
ESCADA = ['C-BUF70C', 'C-BUF80C', 'C-BUF90C', 'C-BUF95C']
# conjunto -> (pasta dos resultados *_oosC.json, pasta do refinamento, rótulo)
CONJ = {'sp': ('results_c_sp', os.path.join(PESQ, 'resultados_caso_sp', 'refinamento'), r'S\~ao Paulo'),
        'td': ('results_c_sp_td', os.path.join(PESQ, 'resultados_caso_sp_td', 'refinamento'), r'S\~ao Paulo, time-dependent travel times'),
        'ath': ('results_c_atenas', os.path.join(PESQ, 'resultados_caso_atenas', 'refinamento'), 'Athens')}
LINHAS_A = [('BUF0.95+C', r'\mBUF{0.95}'), ('CAL', r'\mCAL{}'), ('CAL+R', r'\mCAL{}+R'), ('SIM4', r'\mSIM{}, best of 4')]
MEDIDAS = ((0, 'rel'), (2, 'abs'), (3, 'rel'))          # custo (%), perda (p.p.), CO2 (%); índice na tupla (f1, p_on, q_oos, f3)


def hodges_lehmann(d, nivel):
    """HL e IC bilateral (nível dado) pelos quantis das médias de Walsh (igual ao de equivalencia.py)"""
    d = np.asarray(d, float); n = len(d)
    w = np.sort((d[:, None] + d[None, :])[np.triu_indices(n)] / 2.0)
    M = len(w)
    z = norm.ppf(1 - (1 - nivel) / 2)
    C = max(int(np.floor(M / 2 - z * np.sqrt(n * (n + 1) * (2 * n + 1) / 24))), 1)
    return float(np.median(w)), float(w[C - 1]), float(w[M - C])


def holm(p):
    p = np.asarray(p, float); o = np.argsort(p); adj = np.empty(len(p)); run = 0.0
    for k, i in enumerate(o):
        run = max(run, min(1.0, (len(p) - k) * p[i])); adj[i] = run
    return adj


def pareado(planos, ks, a, b, idx, esc):
    """HL, IC 95% e p de Wilcoxon da diferença por instância (média das sementes) do método a contra o b"""
    pi, ps = defaultdict(list), defaultdict(list)
    for k in ks:
        pi[k[0]].append(planos[a][k][idx]); ps[k[0]].append(planos[b][k][idx])
    I = sorted(pi)
    x = np.array([np.mean(pi[i]) for i in I]); y = np.array([np.mean(ps[i]) for i in I])
    d = 100 * (x / y - 1) if esc == 'rel' else 100 * (x - y)
    hl, lo, hi = hodges_lehmann(d, 0.95)
    p = float(wilcoxon(d).pvalue) if np.any(d != 0) else 1.0
    return dict(n=len(I), hl=hl, lo=lo, hi=hi, p=p)


def carregar_planos(dir_ref):
    """planos escolhidos de cada método em cada execução (instância, semente), como em analisar_refinamento.py"""
    B = {(x['instance'], x['seed']): x for x in json.load(open(os.path.join(dir_ref, 'refinamento_C-BUF95C.json'), encoding='utf-8'))}
    S = {(x['instance'], x['seed']): x for x in json.load(open(os.path.join(dir_ref, 'refinamento_C-SIMC.json'), encoding='utf-8'))}
    CR = {(x['instance'], x['seed']): x for x in json.load(open(os.path.join(dir_ref, 'refinamento_CAL.json'), encoding='utf-8'))}
    planos = defaultdict(dict)
    for k, x in B.items():
        if 'C-SIMC_f1' not in x:
            continue
        g = lambda pre: (x[f'{pre}_f1'], x[f'{pre}_p_oos'], x[f'{pre}_q_oos'], x[f'{pre}_f3'])   # noqa: E731
        planos['SIM+C'][k] = g('C-SIMC'); planos['BUF0.95+C'][k] = g('base')
        planos['BUF0.95+R'][k] = g('fin'); planos['BUF0.95+R-LCB'][k] = g('finl')
        rung = [g('base')] + [g(r) for r in ('C-BUF70C', 'C-BUF80C', 'C-BUF90C') if f'{r}_f1' in x]
        planos['CAL'][k] = min(rung, key=lambda t: t[0])
        if k in S:
            planos['SIM+C+R'][k] = (S[k]['fin_f1'], S[k]['fin_p_oos'], S[k]['fin_q_oos'], S[k]['fin_f3'])
        if k in CR:
            planos['CAL+R'][k] = (CR[k]['fin_f1'], CR[k]['fin_p_oos'], CR[k]['fin_q_oos'], CR[k]['fin_f3'])
            planos['CAL+R-LCB'][k] = (CR[k]['finl_f1'], CR[k]['finl_p_oos'], CR[k]['finl_q_oos'], CR[k]['finl_f3'])
    # orçamento igual ao do CAL (4 buscas): melhor de 4 execuções do SIM+C (sementes s..s+3, mod 10), escolhida pelo custo da escolha certificada
    # antes do refinamento; SIM4+R usa o plano refinado da semente escolhida
    for (inst, sd) in list(planos['SIM+C']):
        grupo = [(inst, (sd + j) % 10) for j in range(4)]
        if not all(g in S for g in grupo):
            continue
        g = min(grupo, key=lambda g: S[g]['base_f1'])
        planos['SIM4'][(inst, sd)] = (S[g]['base_f1'], S[g]['base_p_oos'], S[g]['base_q_oos'], S[g]['base_f3'])
        planos['SIM4+R'][(inst, sd)] = (S[g]['fin_f1'], S[g]['fin_p_oos'], S[g]['fin_q_oos'], S[g]['fin_f3'])
    ks = sorted(set.intersection(*[set(v) for v in planos.values()]))
    return planos, ks


def painel_a(planos, ks):
    """cumprimento e diferença contra o SIM+C por método (Holm por medida sobre todos os métodos, como no relatório do refinamento)"""
    metodos = [m for m in planos if m != 'SIM+C']
    comp = defaultdict(dict)
    for idx, esc in MEDIDAS:
        res = [(m, pareado(planos, ks, m, 'SIM+C', idx, esc)) for m in metodos]
        for (m, r), pa in zip(res, holm([r['p'] for _, r in res])):
            r['p_holm'] = float(pa); comp[m][idx] = r
    cumpre = {m: 100 * float(np.mean([planos[m][k][1] >= ALPHA for k in ks])) for m in planos}
    return comp, cumpre


def painel_b(planos, ks):
    """orçamento igual: CAL(+R) contra SIM4(+R) nas quatro medidas, Holm sobre as quatro"""
    out = {}
    for a_, b_ in (('CAL', 'SIM4'), ('CAL+R', 'SIM4+R')):
        res = [pareado(planos, ks, a_, b_, idx, esc) for idx, esc in ((0, 'rel'), (2, 'abs'), (3, 'rel'), (1, 'abs'))]
        for r, pa in zip(res, holm([r['p'] for r in res])):
            r['p_holm'] = float(pa)
        out[a_] = res            # custo, perda, CO2, P_oos
    return out


def resumo_oos(arq):
    S = [s for s in json.load(open(arq, encoding='utf-8'))['solutions'] if s['feasible']]
    ok = [s for s in S if s['p_on'] >= ALPHA]
    cert = [s for s in S if (s.get('p_cert') or 0) >= ALPHA]
    e = min(cert, key=lambda s: s['f1']) if cert else None
    return dict(custo_ok=min(s['f1'] for s in ok) if ok else np.nan, esc=(e['f1'], e['p_on'], e['q_oos'], e['f3']) if e else None)


def painel_c():
    """SIM+C com 100 cenários de busca em São Paulo: SIM100 é a referência (como em caso_sp/simc100_sp.py)"""
    res_sp = os.path.join(BASE, 'results_c_sp')
    run = defaultdict(dict)
    for pasta, met in ((os.path.join(BASE, 'results_c_sp_simc100', 'C-SIMC100'), 'SIM100'), (os.path.join(res_sp, 'C-SIMC'), 'SIM30'),
                       (os.path.join(res_sp, 'C-BUF95C'), 'BUF0.95+C')):
        for f in glob.glob(os.path.join(pasta, '*_oosC.json')):
            nome, _, sd = os.path.basename(f)[:-len('_oosC.json')].split('__')
            run[(nome, int(sd.replace('seed', '')))][met] = resumo_oos(f)
    for (nome, sd), d in run.items():
        deg = []
        for a in ESCADA:
            f = os.path.join(res_sp, a, f'{nome}__{a}__seed{sd}_oosC.json')
            if os.path.exists(f):
                r = resumo_oos(f)
                if r['esc']:
                    deg.append(r)
        if deg:
            d['CAL'] = min(deg, key=lambda r: r['esc'][0])
    for x in json.load(open(os.path.join(CONJ['sp'][1], 'refinamento_CAL.json'), encoding='utf-8')):
        k = (x['instance'], x['seed'])
        if k in run:
            run[k]['CAL+R'] = dict(custo_ok=np.nan, esc=(x['fin_f1'], x['fin_p_oos'], x['fin_q_oos'], x['fin_f3']))
    ks = sorted(k for k, d in run.items() if 'SIM100' in d)
    cumpre, comp = {}, {}
    for m in ('SIM100', 'SIM30', 'BUF0.95+C', 'CAL', 'CAL+R'):
        v = np.array([run[k][m]['esc'] for k in ks if m in run[k] and run[k][m]['esc']])
        cumpre[m] = 100 * float(np.mean(v[:, 1] >= ALPHA))
    for m in ('SIM30', 'BUF0.95+C', 'CAL', 'CAL+R'):
        res = []
        for idx, esc in ((0, 'rel'), (2, 'abs'), (3, 'rel'), (1, 'abs')):
            pi, ps = defaultdict(list), defaultdict(list)
            for k in ks:
                if m in run[k] and run[k][m]['esc'] and run[k]['SIM100']['esc']:
                    pi[k[0]].append(run[k][m]['esc'][idx]); ps[k[0]].append(run[k]['SIM100']['esc'][idx])
            I = sorted(pi)
            x = np.array([np.mean(pi[i]) for i in I]); y = np.array([np.mean(ps[i]) for i in I])
            dd = 100 * (x / y - 1) if esc == 'rel' else 100 * (x - y)
            hl, lo, hi = hodges_lehmann(dd, 0.95)
            res.append(dict(n=len(I), hl=hl, lo=lo, hi=hi, p=float(wilcoxon(dd).pvalue) if np.any(dd != 0) else 1.0))
        for r, pa in zip(res, holm([r['p'] for r in res])):
            r['p_holm'] = float(pa)
        comp[m] = res            # custo, CO2? ordem: custo, perda, CO2, P_oos
    return cumpre, comp, len(ks)


def painel_d(pasta):
    """preço da confiabilidade e custo do DET+C (como em caso_sp/tabela_caso_sp.py)"""
    R = {}
    for f in glob.glob(os.path.join(BASE, pasta, 'C-*', '*_oosC.json')):
        d = json.load(open(f, encoding='utf-8')); m = d['meta']
        R[(m['algorithm'], m['instance'], m['seed'])] = d['solutions']
    run = defaultdict(dict)
    for (a, i, s), sols in R.items():
        S = [x for x in sols if x['feasible']]
        ok = [x for x in S if x['p_on'] >= ALPHA]
        cert = [x for x in S if (x.get('p_cert') or 0) >= ALPHA]
        run[(i, s)][a] = dict(custo_ok=min(x['f1'] for x in ok) if ok else np.nan, esc=min(cert, key=lambda x: x['f1']) if cert else None,
                              mais_barata=min(S, key=lambda x: x['f1']) if S else None)
    out = {}
    for nivel in ('MEDIUM', 'HIGH', 'DATA'):
        ks = [k for k in run if k[0].endswith(nivel)]
        if not ks:
            continue
        pdet, psim, pc, pq, pe, dc = [], [], [], [], [], []
        for k in sorted(ks):
            d = run[k]
            if 'C-DET' in d and 'C-SIMC' in d and d['C-SIMC']['esc'] and d['C-DET']['mais_barata']:
                a_, b_ = d['C-SIMC']['esc'], d['C-DET']['mais_barata']
                pc.append(100 * (a_['f1'] / b_['f1'] - 1)); pq.append(100 * (a_['q_oos'] - b_['q_oos']))
                pe.append(100 * (a_['f3'] / b_['f3'] - 1)); pdet.append(b_['p_on']); psim.append(a_['p_on'])
            if 'C-DETC' in d and 'C-SIMC' in d and np.isfinite(d['C-DETC']['custo_ok']) and np.isfinite(d['C-SIMC']['custo_ok']):
                dc.append(100 * (d['C-DETC']['custo_ok'] / d['C-SIMC']['custo_ok'] - 1))
        out[nivel] = dict(n=len(ks), p_det=float(np.median(pdet)), p_sim=float(np.median(psim)), custo=float(np.median(pc)), perda=float(np.median(pq)),
                          co2=float(np.median(pe)), detc_med=float(np.median(dc)), detc_max=float(np.max(dc)))
    return out


def fh(x):
    return '0.00' if abs(x) < 0.005 else f'{x:+.2f}'


def fp(p):
    return r'$<$0.001' if p < 0.001 else (f'{p:.3f}' if p < 0.1 else f'{p:.2f}')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--tex', default=None, help='Tabela S33 (comparações pareadas contra o SIM+C)')
    ap.add_argument('--tex2', default=None, help='Tabela S34 (orçamento igual e preço da confiabilidade)')
    ap.add_argument('--out', default=os.path.join(PESQ, 'resultados_casos_sp_atenas'))
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    rel = []

    def log(t=''):
        print(t, flush=True); rel.append(t)

    A, Bq, nks = {}, {}, {}
    for c in ('sp', 'td', 'ath'):
        planos, ks = carregar_planos(CONJ[c][1])
        nks[c] = len({k[0] for k in ks})
        A[c] = painel_a(planos, ks); Bq[c] = painel_b(planos, ks)
        log(f'=== {c}: {len(ks)} execuções em comum, {nks[c]} instâncias ===')
        for m, _ in LINHAS_A:
            r = A[c][0][m]
            log(f"  {m:10s} cumpre {A[c][1][m]:5.1f}%  custo {r[0]['hl']:+.3f} (p_Holm {r[0]['p_holm']:.3g})  perda {r[2]['hl']:+.3f} (p_Holm {r[2]['p_holm']:.3g})  "
                f"CO2 {r[3]['hl']:+.3f} (p_Holm {r[3]['p_holm']:.3g})")
        for a_, res in Bq[c].items():
            log(f"  orçamento igual {a_}: " + '  '.join(f"{n} {r['hl']:+.3f} (p_Holm {r['p_holm']:.3g})" for n, r in zip(('custo', 'perda', 'CO2', 'P_oos'), res)))
    cumpre_c, comp_c, n_c = painel_c()
    log(f'=== SIM100 em São Paulo ({n_c} execuções) ===  cumpre ' + '  '.join(f'{m} {v:.1f}%' for m, v in cumpre_c.items()))
    for m, res in comp_c.items():
        log(f"  {m:10s} contra SIM100: " + '  '.join(f"{n} {r['hl']:+.3f} (p_Holm {r['p_holm']:.3g})" for n, r in zip(('custo', 'perda', 'CO2', 'P_oos'), res)))
    D = {'sp': painel_d('results_c_sp'), 'ath': painel_d('results_c_atenas')}
    for c, dd in D.items():
        for nivel, v in dd.items():
            log(f"=== preço ({c}, {nivel}, {v['n']} execuções) ===  P {v['p_det']:.3f} → {v['p_sim']:.3f}; custo {v['custo']:+.2f}%  perda {v['perda']:+.3f} p.p.  "
                f"CO2 {v['co2']:+.2f}%; DET+C contra SIM+C: mediana {v['detc_med']:+.2f}%, máximo {v['detc_max']:+.2f}%")
    open(os.path.join(a.out, 'relatorio_tabela_casos.txt'), 'w', encoding='utf-8').write('\n'.join(rel) + '\n')

    # ---- tabelas LaTeX: S33 (comparações pareadas contra o SIM+C) e S34 (orçamento igual e preço da confiabilidade) ----
    def linha_a(rot, cumpre, r):
        cel = [f'{cumpre:.1f}']
        if r is None:
            cel += ['--'] * 6
        else:
            for idx in (0, 2, 3):
                cel += [fh(r[idx]['hl']), fp(r[idx]['p_holm'])]
        return rot + ' & ' + ' & '.join(cel) + r' \\'

    def cab(c):
        return r'\multicolumn{8}{@{}l}{\textit{' + CONJ[c][2] + f' ({nks[c]} ' + ('days' if c == 'ath' else 'instances') + r')}} \\'
    T1 = [r' & & \multicolumn{2}{c}{Cost (\%)} & \multicolumn{2}{c}{Quality loss (p.p.)} & \multicolumn{2}{c}{CO$_2$ (\%)} \\',
          r'\cmidrule(lr){3-4}\cmidrule(lr){5-6}\cmidrule(l){7-8}',
          r'Method & Complies (\%) & HL & $p$ & HL & $p$ & HL & $p$ \\', r'\midrule']
    for c in ('sp', 'td', 'ath'):
        T1.append(cab(c))
        T1.append(linha_a(r'\mSIM{}', A[c][1]['SIM+C'], None))
        for m, rot in LINHAS_A:
            T1.append(linha_a(rot, A[c][1][m], A[c][0][m]))
        T1.append(r'\midrule')
    T1.append(r'\multicolumn{8}{@{}l}{\textit{S\~ao Paulo, against \mSIM{} with 100 search scenarios (' + f'{nks["sp"]} instances)' + r'}} \\')
    T1.append(linha_a(r'\mSIM{}, 100 scenarios', cumpre_c['SIM100'], None))
    for m, rot in (('SIM30', r'\mSIM{}, 30 scenarios'), ('BUF0.95+C', r'\mBUF{0.95}'), ('CAL', r'\mCAL{}'), ('CAL+R', r'\mCAL{}+R')):
        r = comp_c[m]
        T1.append(linha_a(rot, cumpre_c[m], {0: r[0], 2: r[1], 3: r[2]}))
    T2 = [r'\multicolumn{9}{@{}l}{\textit{A. Equal budget (four searches)}} \\',
          r' & \multicolumn{2}{c}{Cost (\%)} & \multicolumn{2}{c}{Quality loss (p.p.)} & \multicolumn{2}{c}{CO$_2$ (\%)} & \multicolumn{2}{c}{Punctuality (p.p.)} \\',
          r'\cmidrule(lr){2-3}\cmidrule(lr){4-5}\cmidrule(lr){6-7}\cmidrule(l){8-9}',
          r'Comparison & HL & $p$ & HL & $p$ & HL & $p$ & HL & $p$ \\', r'\midrule']
    for c in ('sp', 'ath'):
        T2.append(r'\multicolumn{9}{@{}l}{\textit{' + CONJ[c][2] + f' ({nks[c]} ' + ('days' if c == 'ath' else 'instances') + r')}} \\')
        for a_, rot in (('CAL', r'\mCAL{} against \mSIM{}, best of 4'), ('CAL+R', r'\mCAL{}+R against \mSIM{}+R, best of 4')):
            res = Bq[c][a_]
            T2.append(rot + ' & ' + ' & '.join(f"{fh(r['hl'])} & {fp(r['p_holm'])}" for r in res) + r' \\')
        if c == 'sp':
            T2.append(r'\midrule')
    T3 = [r'\multicolumn{7}{@{}l}{\textit{B. Price of reliability (medians over runs)}} \\',
          r' & & \multicolumn{2}{c}{On-time probability} & \multicolumn{3}{c}{Change from \mDET{}} \\',
          r'\cmidrule(lr){3-4}\cmidrule(l){5-7}',
          r'Set & Runs & \mDET{} & \mSIM{} & Cost (\%) & Loss (p.p.) & CO$_2$ (\%) \\', r'\midrule']
    for c, nivel, rot in (('sp', 'MEDIUM', r'S\~ao Paulo, medium'), ('sp', 'HIGH', r'S\~ao Paulo, high'), ('ath', 'DATA', 'Athens')):
        v = D[c][nivel]
        T3.append(f"{rot} & {v['n']} & {v['p_det']:.3f} & {v['p_sim']:.3f} & {v['custo']:+.2f} & {v['perda']:+.2f} & {v['co2']:+.2f}" + r' \\')
    T3 += [r'\midrule', r'\multicolumn{7}{@{}l}{\textit{\mDETC{} against \mSIM{}: cost of the cheapest compliant solution (\%)}} \\',
           r' & & \multicolumn{2}{c}{Median} & \multicolumn{3}{c}{Maximum} \\', r'\midrule']
    for c, nivel, rot in (('sp', 'MEDIUM', r'S\~ao Paulo, medium'), ('sp', 'HIGH', r'S\~ao Paulo, high'), ('ath', 'DATA', 'Athens')):
        v = D[c][nivel]
        T3.append(f"{rot} & {v['n']} & \\multicolumn{{2}}{{c}}{{{v['detc_med']:+.2f}}} & \\multicolumn{{3}}{{c}}{{{v['detc_max']:+.2f}}}" + r' \\')
    tex1 = (r"""\begin{table}[htbp]
\centering
\caption{S\~ao Paulo and Athens: the plan a decision-maker would choose (the cheapest certified one; the refined plan for +R) with each method against \mSIM{} (paired comparisons quoted in Sections~6.10 and~6.11; 10 seeds per instance; $\alpha=0.90$; evaluation on 1{,}000 common scenarios). Complies: \% of runs whose plan meets $\alpha$ out of sample. Differences to \mSIM{}: Hodges--Lehmann (HL) estimate of the paired difference per instance (mean over seeds) in \% for cost and CO$_2$ and in percentage points for quality loss, with the Holm-adjusted $p$ of the Wilcoxon signed-rank test (Holm within each measure over the nine methods of the refinement analysis, of which four are shown). Best of 4: four runs of \mSIM{}, the budget of \mCAL{}. The last block takes \mSIM{} with 100 search scenarios as the reference (Holm over the four measures of each method).}
\label{tab:sp-ath}
\footnotesize\setlength{\tabcolsep}{3.5pt}
\begin{tabular}{@{}l r rr rr rr@{}}
\toprule
""" + '\n'.join(T1) + r"""
\bottomrule
\end{tabular}
\end{table}
""")
    tex2 = (r"""\begin{table}[htbp]
\centering
\caption{S\~ao Paulo and Athens: equal budget and price of reliability (quoted in Sections~6.10 and~6.11; 10 seeds per instance; $\alpha=0.90$; evaluation on 1{,}000 common scenarios). A: the chosen plan of \mCAL{} against the best of four runs of \mSIM{}, and the refined plans of both; Hodges--Lehmann (HL) estimate of the paired difference per instance (mean over seeds), with the Holm-adjusted $p$ of the Wilcoxon signed-rank test over the four measures of each comparison; Punctuality is the difference in the out-of-sample probability that no customer is late (the joint on-time probability) of the chosen plan, in percentage points. B: the plan chosen from \mSIM{} against the cheapest plan of \mDET{}, which ignores punctuality (medians over runs, 80 per uncertainty level in S\~ao Paulo and 90 in Athens), and the cost of the cheapest compliant solution of \mDETC{} relative to \mSIM{}.}
\label{tab:sp-ath-b}
\footnotesize\setlength{\tabcolsep}{3.5pt}
\begin{tabular}{@{}l rr rr rr rr@{}}
\toprule
""" + '\n'.join(T2) + r"""
\bottomrule
\end{tabular}

\medskip
\begin{tabular}{@{}l r rr rrr@{}}
\toprule
""" + '\n'.join(T3) + r"""
\bottomrule
\end{tabular}
\end{table}
""")
    destino = a.tex or os.path.join(a.out, 'tab_casos_pareadas.tex')
    open(destino, 'w', encoding='utf-8').write(tex1)
    destino2 = a.tex2 or os.path.join(a.out, 'tab_casos_orcamento.tex')
    open(destino2, 'w', encoding='utf-8').write(tex2)
    print('tabelas:', destino, destino2)


if __name__ == '__main__':
    main()
