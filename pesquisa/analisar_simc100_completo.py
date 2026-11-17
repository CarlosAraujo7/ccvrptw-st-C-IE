#!/usr/bin/env python3
"""
B2 (16/11/2026, diário Entrada 128): SIM+C com 100 cenários de busca em todo o benchmark (n = 25, 50 e 100).
Plano fixado antes de rodar (Entrada 128); este script só executa o plano.

Execuções: results_c_extra/C-SIMC100 (n = 50: Entrada 81; n = 25 e 100: B2), sementes 0-4, MEDIUM e HIGH (36 instâncias).
Plano escolhido = solução certificada mais barata (P_cert >= α) da frente; medido nos 1.000 cenários comuns
(avaliação fora da amostra de analisar_c.py; cria os *_oosC.json que faltam).
Comparações por instância (média das sementes 0-4): SIM+C com 100 cenários contra SIM+C com 30, contra BUF0.95+C e
contra o CAL (4 buscas). Hodges-Lehmann com IC de 95%, Wilcoxon e Holm nas quatro medidas.
Também: % que cumpre α fora da amostra, execuções sem plano certificado, otimismo da amostra de busca
(média sobre a frente de P_on dentro da amostra − P_on fora) e razão de tempo de execução.
Uso: python pesquisa/analisar_simc100_completo.py [--sem-oos]
Saída: pesquisa/resultados_robustez/relatorio_simc100_completo.txt, resumo_simc100_completo.json, tab_simc100_completo.tex
"""
import argparse, glob, json, os, sys
from collections import defaultdict
from multiprocessing import Pool
import numpy as np
from scipy.stats import wilcoxon

BASE = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
sys.path.insert(0, os.path.join(BASE, 'pesquisa'))
import calibracao_preditiva_base as b      # noqa: E402
import analisar_c as ac                    # noqa: E402
from equivalencia import hodges_lehmann    # noqa: E402
from analisar_refinamento import holm      # noqa: E402

ALPHA = b.ALPHA
OUT = os.path.join(BASE, 'pesquisa', 'resultados_robustez')
D100 = os.path.join(BASE, 'results_c_extra', 'C-SIMC100')
D30 = os.path.join(BASE, 'results_c', 'C-SIMC')
SEM = (0, 1, 2, 3, 4)
MEDIDAS = (('custo (%)', 0, 'rel'), ('perda (p.p.)', 2, 'abs'), ('CO2 (%)', 3, 'rel'), ('P_oos (p.p.)', 1, 'abs'))
GRUPOS = (('n = 25', lambda i: '_025_' in i), ('n = 50', lambda i: '_050_' in i), ('n = 100', lambda i: '_100_' in i),
          ('todas', lambda i: True))


def tup(s):
    return None if s is None else (s['f1'], s['p_on'], s['q_oos'], s['f3'])


def carregar100():
    P = {}; otim = {}; tempo = {}
    for f in glob.glob(os.path.join(D100, '*_oosC.json')):
        d = json.load(open(f, encoding='utf-8')); m = d['meta']
        if m['seed'] not in SEM or not m['instance'].endswith(('_MEDIUM', '_HIGH')):
            continue
        S = [s for s in d['solutions'] if s['feasible']]
        cert = [s for s in S if (s.get('p_cert') or 0) >= ALPHA]
        k = (m['instance'], m['seed'])
        P[k] = tup(min(cert, key=lambda s: s['f1']) if cert else None)
        o = [s['p_on_insample'] - s['p_on'] for s in S if s.get('p_on_insample') is not None]
        otim[k] = float(np.mean(o)) if o else float('nan')
        tempo[k] = json.load(open(f.replace('_oosC.json', '_metrics.json'), encoding='utf-8'))['elapsed_s']
    return P, otim, tempo


def carregar30():
    """Planos escolhidos e tempos dos métodos do motor com 30 cenários / buffers, mesmas instâncias e sementes."""
    R = b.carregar('results_c')
    P = defaultdict(dict); otim = {}; tempo = {}
    for k, d in R.items():
        if k[1] not in SEM or not k[0].endswith(('_MEDIUM', '_HIGH')):
            continue
        P['SIM+C'][k] = tup(d.get('C-SIMC'))
        P['BUF0.95+C'][k] = tup(d.get('C-BUF95C'))
        P['CAL'][k] = tup(b.melhor([d.get(a) for a in b.DEG]))
    for f in glob.glob(os.path.join(D30, '*_oosC.json')):
        d = json.load(open(f, encoding='utf-8')); m = d['meta']
        if m['seed'] not in SEM or not m['instance'].endswith(('_MEDIUM', '_HIGH')):
            continue
        S = [s for s in d['solutions'] if s['feasible']]
        o = [s['p_on_insample'] - s['p_on'] for s in S if s.get('p_on_insample') is not None]
        k = (m['instance'], m['seed'])
        otim[k] = float(np.mean(o)) if o else float('nan')
        tempo[k] = json.load(open(f.replace('_oosC.json', '_metrics.json'), encoding='utf-8'))['elapsed_s']
    return P, otim, tempo


def por_inst(D, filtro, idx):
    pi = defaultdict(list)
    for (i, _), t in D.items():
        if t is not None and filtro(i):
            pi[i].append(t[idx])
    return {i: float(np.mean(v)) for i, v in pi.items()}


def comparar(Pm, Pr, filtro):
    res = []
    for nm, idx, esc in MEDIDAS:
        x, y = por_inst(Pm, filtro, idx), por_inst(Pr, filtro, idx)
        I = sorted(set(x) & set(y))
        if len(I) < 3:
            return None
        a = np.array([x[i] for i in I]); c = np.array([y[i] for i in I])
        dd = 100 * (a / c - 1) if esc == 'rel' else 100 * (a - c)
        hl, lo, hi = hodges_lehmann(dd, 0.95)
        p = wilcoxon(dd).pvalue if np.any(np.abs(dd) > 1e-12) else 1.0
        res.append([nm, len(I), hl, lo, hi, p, int(np.sum(dd < -1e-12)), int(np.sum(dd > 1e-12))])
    for r, pa in zip(res, holm([r[5] for r in res])):
        r.append(pa)
    return res


def cumpre(D, filtro):
    v = [t[1] for k, t in D.items() if t is not None and filtro(k[0])]
    sem = sum(1 for k, t in D.items() if t is None and filtro(k[0]))
    return (100 * float(np.mean(np.array(v) >= ALPHA)) if v else float('nan')), len(v), sem


def main(sem_oos):
    if not sem_oos:
        files = sorted(glob.glob(os.path.join(D100, '*_solutions.json')))
        with Pool(3) as pool:
            pool.map(ac._job, files, chunksize=2)
    P100, o100, t100 = carregar100()
    P, o30, t30 = carregar30()
    P['SIM+C 100'] = P100
    os.makedirs(OUT, exist_ok=True)
    out = [f'B2 — SIM+C com 100 cenários de busca (SIM+C 100) contra SIM+C com 30, BUF0.95+C e CAL; sementes 0-4, MEDIUM e HIGH; α = {ALPHA}',
           f'execuções de SIM+C 100: {len(P100)}; instâncias: {len({k[0] for k in P100})}']
    resumo = {'comp': {}, 'cumpre': {}, 'tempo': {}, 'otim': {}}
    for g, filtro in GRUPOS:
        ks = [k for k in P100 if filtro(k[0])]
        if not ks:
            continue
        out.append(f'\n== {g} ({len({k[0] for k in ks})} instâncias, {len(ks)} execuções) ==')
        for m in ('SIM+C', 'SIM+C 100', 'BUF0.95+C', 'CAL'):
            Dm = {k: v for k, v in P[m].items() if k[0] in {x[0] for x in ks}}     # mesmas instâncias (e sementes 0-4) do SIM+C 100
            c, nc, sp = cumpre(Dm, filtro)
            resumo['cumpre'][f'{g}|{m}'] = c
            v = np.array([t for t in Dm.values() if t is not None])
            out.append(f'  {m:10s} execuções {len(Dm):3d} sem plano {sp:2d}  cumpre {c:5.1f}%  custo {v[:, 0].mean():9.1f}  '
                       f'perda {100 * v[:, 2].mean():.3f}%  CO2 {v[:, 3].mean():.1f}')
        com = [k for k in ks if k in t30]
        if com:
            r_t = np.mean([t100[k] / t30[k] for k in com]); r_t2 = np.mean([t100[k] for k in com]) / np.mean([t30[k] for k in com])
            oo100 = np.nanmean([o100[k] for k in com]); oo30 = np.nanmean([o30[k] for k in com])
            out.append(f'  tempo médio por execução: 100 cenários {np.mean([t100[k] for k in com]):.0f} s, 30 cenários {np.mean([t30[k] for k in com]):.0f} s '
                       f'(razão das médias {r_t2:.2f}); otimismo da amostra de busca: 100 cen. {oo100:+.4f}, 30 cen. {oo30:+.4f}')
            resumo['tempo'][g] = dict(t100=float(np.mean([t100[k] for k in com])), t30=float(np.mean([t30[k] for k in com])), razao=float(r_t2))
            resumo['otim'][g] = dict(o100=float(oo100), o30=float(oo30))
        for ref in ('SIM+C', 'BUF0.95+C', 'CAL'):
            r = comparar(P['SIM+C 100'], P[ref], filtro)
            if r is None:
                continue
            out.append(f'  SIM+C 100 contra {ref} (por instância; HL [IC 95%]; p; p_Holm; menores/maiores)')
            for nm, ni, hl, lo, hi, p, mn, mx, pa in r:
                out.append(f'      {nm:14s} n={ni:2d}  HL {hl:+7.3f} [{lo:+7.3f}, {hi:+7.3f}]  p={p:.3g}  p_Holm={pa:.3g}  ({mn}/{mx})')
                resumo['comp'][f'{g}|{ref}|{nm}'] = dict(n=ni, hl=hl, lo=lo, hi=hi, p=p, p_holm=pa, menores=mn, maiores=mx)
    txt = '\n'.join(out); print(txt)
    open(os.path.join(OUT, 'relatorio_simc100_completo.txt'), 'w', encoding='utf-8').write(txt + '\n')
    json.dump(resumo, open(os.path.join(OUT, 'resumo_simc100_completo.json'), 'w', encoding='utf-8'), indent=1)
    tabela(resumo)


def tabela(resumo):
    """Tabela do suplemento: SIM+C 100 contra SIM+C 30 e contra o CAL, por tamanho (HL; * = p_Holm < 0,05)."""
    gs = [g for g, _ in GRUPOS if f'{g}|SIM+C|custo (%)' in resumo['comp']]
    if not gs:
        return

    def cel(g, ref, nm):
        d = resumo['comp'].get(f'{g}|{ref}|{nm}')
        if d is None:
            return '--'
        return ('0.00' if abs(d['hl']) < 0.005 else f'{d["hl"]:+.2f}') + ('$^*$' if d['p_holm'] < 0.05 else '')
    ROT = {'n = 25': r'$n=25$', 'n = 50': r'$n=50$', 'n = 100': r'$n=100$', 'todas': 'All'}
    L = [r'\begin{table}[htbp]', r'\centering',
         r'\caption{Search sample of the simheuristic: \mSIM{} with 100 travel-time scenarios inside the search (SIM+C 100) against \mSIM{} with 30 scenarios and against \mCAL{} (four searches). '
         r'Medium and high uncertainty, seeds 0--4, by instance size (number of instances in brackets). Hodges--Lehmann estimate of the paired difference per instance '
         r'(mean over seeds); $^*$: $p_{\mathrm{Holm}}<0.05$ (Wilcoxon signed-rank test, Holm correction over the four measures including punctuality). '
         r'Complies: \% of the chosen plans meeting $\alpha=0.90$ out of sample, SIM+C 100 / reference. The runs with 50 customers are those of the first check of the search sample.}',
         r'\label{tab:simc100-full}', r'\small\setlength{\tabcolsep}{4pt}', r'\begin{tabular}{@{}l l r r r r@{}}', r'\toprule',
         r'Instances & Against & Cost (\%) & Loss (p.p.) & CO$_2$ (\%) & Compl. (\%) \\', r'\midrule']
    for g in gs:
        ni = resumo['comp'][f'{g}|SIM+C|custo (%)']['n']
        for j, (ref, rot) in enumerate((('SIM+C', r'\mSIM{} (30 scenarios)'), ('CAL', r'\mCAL{}'))):
            c1, c0 = resumo['cumpre'].get(f'{g}|SIM+C 100'), resumo['cumpre'].get(f'{g}|{ref}')
            L.append(f'{(ROT[g] + f" ({ni})") if j == 0 else ""} & {rot} & {cel(g, ref, "custo (%)")} & {cel(g, ref, "perda (p.p.)")} & {cel(g, ref, "CO2 (%)")} & '
                     f'{c1:.1f} / {c0:.1f} \\\\')
        L.append(r'\midrule')
    L[-1] = r'\bottomrule'
    L += [r'\end{tabular}', r'\end{table}']
    open(os.path.join(OUT, 'tab_simc100_completo.tex'), 'w', encoding='utf-8').write('\n'.join(L) + '\n')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--sem-oos', action='store_true', help='não recalcula a avaliação fora da amostra')
    a = ap.parse_args()
    main(a.sem_oos)
