#!/usr/bin/env python3
"""
Análise do B3 (16/11/2026, diário Entrada 130): frota dedicada por categoria contra o multicompartimento.
Plano fixado antes de rodar (Entrada 130); este script só executa o plano.

Algoritmo C-BUF90C; n = 50; 18 instâncias (seis famílias x LOW/MEDIUM/HIGH); sementes 0-4.
Variantes: MIX (Z = 1, compartimento único, carga mista, set-point dominante), Z = 2, Z = 3 (results_c_b2), D2 e D4 (results_c_b3).
Por (instância, semente) em que a variante e a base MIX têm solução certificada fora da amostra (P >= α): diferença da melhor perda de
qualidade (p.p.), do melhor custo (%) e do melhor CO2 (%) contra a base MIX; mediana sobre os pares (igual à Tabela S15).
Também: % de execuções com solução que cumpre; % do valor entregue acima do limite legal da própria categoria (conformidade_legal_b2).
Comparação direta de D2 e D4 contra Z = 2 (por instância; HL, Wilcoxon, Holm nas três medidas).
Uso: python pesquisa/analisar_frota_dedicada.py [--sem-legal]
Saída: pesquisa/resultados_c_b3/relatorio_frota_dedicada.txt, resumo_frota_dedicada.json, tab_frota_dedicada.tex
"""
import argparse, glob, json, os, sys
from collections import defaultdict
from multiprocessing import Pool
import numpy as np
from scipy.stats import wilcoxon

BASE = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
sys.path.insert(0, os.path.join(BASE, 'src')); sys.path.insert(0, os.path.join(BASE, 'pesquisa'))
from equivalencia import hodges_lehmann    # noqa: E402
from analisar_refinamento import holm      # noqa: E402

ALPHA = 0.90
OUT = os.path.join(BASE, 'pesquisa', 'resultados_c_b3')
VAR = [('MIX', os.path.join('results_c_b2', 'C-BUF90C-MIX')), ('Z=2', os.path.join('results_c_b2', 'C-BUF90C_Z2')),
       ('Z=3', os.path.join('results_c_b2', 'C-BUF90C_Z3')), ('D2', os.path.join('results_c_b3', 'C-BUF90C-D2')),
       ('D4', os.path.join('results_c_b3', 'C-BUF90C-D4'))]
SEM = (0, 1, 2, 3, 4)


def carregar(pasta):
    """{(instância, semente): dict(ok, q, c, e)} só para n = 50 e sementes 0-4."""
    d = {}
    for f in glob.glob(os.path.join(BASE, pasta, '*_oosC.json')):
        r = json.load(open(f, encoding='utf-8')); m = r['meta']
        if m['seed'] not in SEM or '_050_' not in m['instance']:
            continue
        B = [s for s in r['solutions'] if s['feasible'] and s['p_on'] >= ALPHA]
        d[(m['instance'], m['seed'])] = dict(ok=bool(B), q=min((s['q_oos'] for s in B), default=np.nan),
                                             c=min((s['f1'] for s in B), default=np.nan), e=min((s['f3'] for s in B), default=np.nan))
    return d


def legal(pasta, workers):
    import conformidade_legal_b2 as cl
    arqs = [f for f in glob.glob(os.path.join(BASE, pasta, '*_solutions.json'))
            if '_050_' in os.path.basename(f) and int(f.split('seed')[-1].split('_')[0]) in SEM]
    with Pool(workers) as pool:
        res = pool.map(cl._job, sorted(arqs), chunksize=4)
    v = [100 * np.mean(r['frac']) for r in res if r['frac']]
    return float(np.mean(v)) if v else float('nan'), len(v)


def pares(A, B):
    """pares (instância, semente) com solução certificada nas duas variantes"""
    return [k for k in B if k in A and A[k]['ok'] and B[k]['ok']]


def difs(A, B, campo, rel):
    ks = pares(A, B)
    v = defaultdict(list)
    for k in ks:
        v[k[0]].append(100 * (B[k][campo] / A[k][campo] - 1) if rel else 100 * (B[k][campo] - A[k][campo]))
    return {i: float(np.mean(x)) for i, x in v.items()}, ks


def main(sem_legal, workers):
    D = {n: carregar(p) for n, p in VAR}
    os.makedirs(OUT, exist_ok=True)
    out = ['B3 — frota dedicada por categoria contra o multicompartimento; C-BUF90C, n = 50, 18 instâncias, sementes 0-4; α = %.2f' % ALPHA,
           'execuções por variante: ' + ', '.join(f'{n} {len(d)}' for n, d in D.items())]
    resumo = {'linhas': {}, 'direto': {}}
    LEG = {}
    if not sem_legal:
        for n, p in VAR:
            LEG[n] = legal(p, workers)
    for n, _ in VAR[1:]:
        ks = pares(D['MIX'], D[n])
        if not ks:
            continue
        dq = np.array([100 * (D[n][k]['q'] - D['MIX'][k]['q']) for k in ks])
        dc = np.array([100 * (D[n][k]['c'] / D['MIX'][k]['c'] - 1) for k in ks])
        de = np.array([100 * (D[n][k]['e'] / D['MIX'][k]['e'] - 1) for k in ks])
        ok1 = 100 * np.mean([D['MIX'][k]['ok'] for k in D[n] if k in D['MIX']]); okz = 100 * np.mean([D[n][k]['ok'] for k in D[n]])
        resumo['linhas'][n] = dict(pares=len(ks), ok_mix=ok1, ok=okz, dq=float(np.median(dq)), dc=float(np.median(dc)), de=float(np.median(de)),
                                   dq_q1=float(np.percentile(dq, 25)), dq_q3=float(np.percentile(dq, 75)),
                                   dc_q1=float(np.percentile(dc, 25)), dc_q3=float(np.percentile(dc, 75)),
                                   legal=LEG.get(n, (float('nan'), 0))[0], legal_mix=LEG.get('MIX', (float('nan'), 0))[0])
        out.append(f'  {n:4s} pares {len(ks):3d}  cumpre MIX {ok1:5.1f}% / {n} {okz:5.1f}%  Δ perda {np.median(dq):+6.2f} p.p. [{np.percentile(dq, 25):+.2f}, {np.percentile(dq, 75):+.2f}]  '
                   f'Δ custo {np.median(dc):+6.2f}% [{np.percentile(dc, 25):+.2f}, {np.percentile(dc, 75):+.2f}]  Δ CO2 {np.median(de):+6.2f}%  '
                   f'legal {LEG.get(n, (float("nan"), 0))[0]:.1f}% (MIX {LEG.get("MIX", (float("nan"), 0))[0]:.1f}%)')
    out.append('\nComparação direta contra Z = 2 (por instância; HL [IC 95%]; Wilcoxon; Holm nas 3 medidas):')
    for n in ('D2', 'D4', 'Z=3'):
        res = []
        for nm, campo, rel in (('Δ perda (p.p.)', 'q', False), ('Δ custo (%)', 'c', True), ('Δ CO2 (%)', 'e', True)):
            dd, ks = difs(D['Z=2'], D[n], campo, rel)
            I = sorted(dd)
            if len(I) < 3:
                continue
            v = np.array([dd[i] for i in I])
            hl, lo, hi = hodges_lehmann(v, 0.95)
            p = wilcoxon(v).pvalue if np.any(np.abs(v) > 1e-12) else 1.0
            res.append([nm, len(I), hl, lo, hi, p, int(np.sum(v < -1e-12)), int(np.sum(v > 1e-12))])
        for r, pa in zip(res, holm([r[5] for r in res])):
            r.append(pa)
        out.append(f'  {n} contra Z=2:')
        for nm, ni, hl, lo, hi, p, mn, mx, pa in res:
            out.append(f'      {nm:15s} n={ni:2d}  HL {hl:+7.3f} [{lo:+7.3f}, {hi:+7.3f}]  p={p:.3g}  p_Holm={pa:.3g}  ({mn} menores, {mx} maiores)')
            resumo['direto'][f'{n}|{nm}'] = dict(n=ni, hl=hl, lo=lo, hi=hi, p=p, p_holm=pa)
    txt = '\n'.join(out); print(txt)
    open(os.path.join(OUT, 'relatorio_frota_dedicada.txt'), 'w', encoding='utf-8').write(txt + '\n')
    json.dump(resumo, open(os.path.join(OUT, 'resumo_frota_dedicada.json'), 'w', encoding='utf-8'), indent=1)
    tabela(resumo)


def tabela(resumo):
    if not resumo['linhas']:
        return
    NOME = {'Z=2': r'$Z=2$ zones', 'Z=3': r'$Z=3$ zones', 'D2': r'D2: frozen / chilled', 'D4': r'D4: one per category'}
    L = [r'\begin{table}[htbp]', r'\centering',
         r'\caption{Dedicated vehicles per product class versus multi-compartment vehicles, \mBUF{0.90}, $n=50$, 18 instances $\times$ 5 seeds. '
         r'Reference: a single compartment carrying the same mixed load at the dominant set-point (Z=1, as in Table~S15). Dedicated fleets: each route carries one class '
         r'of product and runs one zone at the set-point of the coldest category present (D2: frozen / chilled; D4: one class per category). Differences are medians over '
         r'the pairs in which both variants have a plan certified out of sample, for the best value of each objective among compliant solutions (interquartile range of the '
         r'cost difference in brackets). Compliance: \% of runs with a certified plan. Legal: mean share of delivered value above the legal limit of its own category.}',
         r'\label{tab:b3}', r'\footnotesize\setlength{\tabcolsep}{3pt}', r'\begin{tabular}{@{}l r r r r r r@{}}', r'\toprule',
         r'Variant & Pairs & Compl. (\%) & $\Delta$ loss (p.p.) & $\Delta f_1$ (\%) & $\Delta f_3$ (\%) & Legal (\%) \\', r'\midrule']
    for n in ('Z=2', 'Z=3', 'D2', 'D4'):
        r = resumo['linhas'].get(n)
        if r is None:
            continue
        L.append(f'{NOME[n]} & {r["pares"]} & {r["ok"]:.0f} & {r["dq"]:+.2f} & {r["dc"]:+.2f} [{r["dc_q1"]:+.1f}, {r["dc_q3"]:+.1f}] & {r["de"]:+.2f} & {r["legal"]:.1f} \\\\')
    r0 = next(iter(resumo['linhas'].values()))
    L += [r'\midrule', f'$Z=1$, mixed load & -- & {r0["ok_mix"]:.0f} & -- & -- & -- & {r0["legal_mix"]:.1f} \\\\', r'\bottomrule', r'\end{tabular}', r'\end{table}']
    open(os.path.join(OUT, 'tab_frota_dedicada.tex'), 'w', encoding='utf-8').write('\n'.join(L) + '\n')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--sem-legal', action='store_true')
    ap.add_argument('--workers', type=int, default=2)
    a = ap.parse_args()
    main(a.sem_legal, a.workers)
