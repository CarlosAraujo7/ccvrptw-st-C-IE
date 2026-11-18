#!/usr/bin/env python3
"""
B4 (17/11/2026, diário Entrada 139): variabilidade de viagem ESTIMADA dos dados observados e usada no planejamento.

Pergunta (P2 e P10 da revisão): se o CV dos tempos de viagem é estimado dos dados (táxis de Nova York) e usado no
planejamento (buffers e certificação), os planos cumprem α em cenários observados, e a que custo?
  κ_obs = 0,27: média do CV observado por faixa de duração (0,29 / 0,27 / 0,25 / 0,20) sobre os trechos dos planos
  escolhidos do CAL nas 6 instâncias MEDIUM com n = 50 (média 0,272).
Braços (6 instâncias MEDIUM, n = 50, sementes 0-4):
  A  escada C-BUF70C/80C/90C/95C planejada e certificada com κ = 0,10 (execuções existentes, results_c)
  B  o mesmo com cópias das instâncias com κ = 0,27 (instances/b4_cv_observado; saída em results_c_b4)
Plano escolhido = CAL (certificado mais barato entre os quatro degraus). Avaliação nos 1.000 cenários de
reavaliar_observado.py (mesma semente por instância), iguais para os dois braços:
  M3        CV observado pela duração nominal de cada trecho do plano + fator do dia (medida principal)
  M2 (27)   forma observada + fator do dia, CV 0,27 (verificação)
Uso: python pesquisa/tempos_observados/b4_cv_estimado.py --preparar
     python pesquisa/tempos_observados/b4_cv_estimado.py --rodar --workers 3
     python pesquisa/tempos_observados/b4_cv_estimado.py --analisar
Saída: pesquisa/tempos_observados/resultados_b4/ (relatorio_b4.txt, resumo_b4.json, tab_b4.tex)
"""
import argparse, glob, json, math, os, sys, time
from collections import defaultdict
from multiprocessing import Pool
import numpy as np

AQUI = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.normpath(os.path.join(AQUI, '..', '..'))
sys.path.insert(0, os.path.join(BASE, 'src')); sys.path.insert(0, os.path.join(BASE, 'pesquisa')); sys.path.insert(0, AQUI)
KAPPA_OBS = 0.27
INST_ORIG = os.path.join(BASE, 'instances', 'generated')
INST_B4 = os.path.join(BASE, 'instances', 'b4_cv_observado')
RES_A = os.path.join(BASE, 'results_c')
RES_B = os.path.join(BASE, 'results_c_b4')
OUTD = os.path.join(AQUI, 'resultados_b4')
TLC = os.path.join(AQUI, 'resultados')
ESC = ['C-BUF70C', 'C-BUF80C', 'C-BUF90C', 'C-BUF95C']
NOMES = [f'CCVRPTW-ST_{f}_050_MEDIUM' for f in ('C1', 'C2', 'R1', 'R2', 'RC1', 'RC2')]
SEMENTES = (0, 1, 2, 3, 4)
ALPHA = 0.90; K = 1000
CV_DUR = [(10, 0.29), (20, 0.27), (40, 0.25), (math.inf, 0.20)]


def preparar():
    os.makedirs(INST_B4, exist_ok=True)
    for nm in NOMES:
        d = json.load(open(os.path.join(INST_ORIG, nm + '.json'), encoding='utf-8'))
        antes = d['thermal_profile']['sigma_t_pct']
        d['thermal_profile']['sigma_t_pct'] = KAPPA_OBS
        json.dump(d, open(os.path.join(INST_B4, nm + '.json'), 'w', encoding='utf-8'), ensure_ascii=False)
        print(f'{nm}: κ {antes} -> {KAPPA_OBS}')


def _rodar(t):
    import experiment_runner as er
    return er._worker(t)


def rodar(workers):
    import experiment_runner as er
    T = []
    for alg in ESC:
        cfg = er.ALGO_CONFIG[alg]; out_dir = os.path.join(RES_B, alg)
        for nm in NOMES:
            for s in SEMENTES:
                key = f'{nm}__{alg}__seed{s}'
                if os.path.exists(os.path.join(out_dir, key + '_metrics.json')):
                    continue
                T.append({'inst_path': os.path.join(INST_B4, nm + '.json'), 'algo': alg, 'algo_label': alg, 'seed': s, 'pop_size': 100,
                          'n_gen': 200, 'mc_samples': cfg['mc_samples'], 'runner': cfg['runner'], 'out_dir': out_dir, 'result_key': key,
                          'use_multizone': False, 'n_zones': 1, 'dedup': cfg.get('dedup', False), 'c_params': cfg.get('c_params')})
    print(f'{len(T)} execuções a rodar', flush=True)
    t0 = time.time()
    with Pool(workers) as pool:
        for i, r in enumerate(pool.imap_unordered(_rodar, T), 1):
            print(f"  [{i}/{len(T)}] {r.get('key', '')} {r.get('status', '')} {r.get('elapsed', '')} [{(time.time() - t0) / 60:.1f} min]", flush=True)
    print('FIM B4 (execuções)', flush=True)


def escolhida(arq):
    if not os.path.exists(arq):
        return None
    S = json.load(open(arq, encoding='utf-8'))
    c = [s for s in S if s['feasible'] and (s.get('p_cert') or 0) >= ALPHA]
    return min(c, key=lambda s: s['f1_cost_BRL']) if c else None


def multiplicadores(sig, e):
    x = np.exp(sig * e)
    return x / x.mean(axis=0, keepdims=True)


def analisar():
    import analisar_c as ac
    import ccvrp_c as cc
    A = np.load(os.path.join(TLC, 'amostras_tlc.npz'))
    rho = json.load(open(os.path.join(TLC, 'variabilidade_tlc.json'), encoding='utf-8'))['fator_dia']['fracao_variancia']
    os.makedirs(OUTD, exist_ok=True)
    sk27 = math.sqrt(math.log1p(KAPPA_OBS ** 2))
    ac.INST_DIR = INST_ORIG; ac._cache.clear()
    cel = {}      # (braço, instância, semente) -> dict(custo, P3, P2)
    sem_plano = {'A': 0, 'B': 0}
    for nm in NOMES:
        inst, T, sc, _ = ac._ctx(nm, 1)
        n = len(inst['customers'])
        rng = np.random.default_rng(99991 + sum(map(ord, nm)))
        e = rng.choice(A['z_res'], size=(K, n + 1)); ed = rng.choice(A['z_dia'], size=(K, 1))
        mix = math.sqrt(rho) * ed + math.sqrt(1 - rho) * e
        tt2 = multiplicadores(sk27, mix)
        for braco, res in (('A', RES_A), ('B', RES_B)):
            for sd in SEMENTES:
                lad = [escolhida(os.path.join(res, a, f'{nm}__{a}__seed{sd}_solutions.json')) for a in ESC]
                lad = [x for x in lad if x]
                if not lad:
                    sem_plano[braco] += 1
                    continue
                x = min(lad, key=lambda s: s['f1_cost_BRL'])
                sig = np.full(n + 1, sk27)
                for r in x['routes']:
                    prev = 0
                    for c in r:
                        tau = T[prev][c]
                        cv = next(v for lim, v in CV_DUR if tau < lim)
                        sig[c] = math.sqrt(math.log1p(cv * cv)); prev = c
                tt3 = multiplicadores(sig[None, :], mix)
                P = {}
                for nome, tt in (('P3', tt3), ('P2', tt2)):
                    on = np.ones(K, dtype=bool)
                    for r, d0 in zip(x['routes'], x['partida_min']):
                        on &= cc.pontual_rota(r, inst, T, tt, d0)
                    P[nome] = float(on.mean())
                cel[(braco, nm, sd)] = dict(custo=x['f1_cost_BRL'], **P)
    out = [f'B4 — CV estimado dos dados (κ = {KAPPA_OBS}) no planejamento; 6 instâncias MEDIUM com n = 50, sementes 0-4; α = {ALPHA}; K = {K}',
           f'execuções sem plano certificado: braço A {sem_plano["A"]} de 30, braço B {sem_plano["B"]} de 30']
    resumo = {'sem_plano': sem_plano}
    for braco, rot in (('A', 'A: planejado e certificado com κ = 0,10 (original)'), ('B', f'B: planejado e certificado com κ = {KAPPA_OBS} (estimado dos dados)')):
        ks = [k for k in cel if k[0] == braco]
        p3 = np.array([cel[k]['P3'] for k in ks]); p2 = np.array([cel[k]['P2'] for k in ks]); c_ = np.array([cel[k]['custo'] for k in ks])
        d = dict(n=len(ks), cumpre_M3=100 * float(np.mean(p3 >= ALPHA)), P_M3=float(p3.mean()), cumpre_M2=100 * float(np.mean(p2 >= ALPHA)),
                 P_M2=float(p2.mean()), custo=float(c_.mean()))
        resumo[braco] = d
        out.append(f'  {rot}: {len(ks)} planos; cumpre sob M3 {d["cumpre_M3"]:5.1f}% (P média {d["P_M3"]:.3f}); sob M2(0,27) {d["cumpre_M2"]:5.1f}% (P média {d["P_M2"]:.3f}); custo médio {d["custo"]:.1f}')
    pares = [(k[1], k[2]) for k in cel if k[0] == 'A' and ('B', k[1], k[2]) in cel]
    dc = np.array([100 * (cel[('B', nm, sd)]['custo'] / cel[('A', nm, sd)]['custo'] - 1) for nm, sd in pares])
    dp3 = np.array([100 * (cel[('B', nm, sd)]['P3'] - cel[('A', nm, sd)]['P3']) for nm, sd in pares])
    por_inst = defaultdict(list)
    for (nm, sd), v in zip(pares, dc):
        por_inst[nm].append(v)
    mi = np.array([np.mean(v) for v in por_inst.values()])
    resumo['custo'] = dict(pares=len(pares), mediana=float(np.median(dc)), q1=float(np.percentile(dc, 25)), q3=float(np.percentile(dc, 75)),
                           media=float(dc.mean()), por_instancia=[float(x) for x in mi], dP3_mediana=float(np.median(dp3)))
    out.append(f'  custo do B contra o A, por (instância, semente) ({len(pares)} pares): mediana {np.median(dc):+.2f}% [quartis {np.percentile(dc, 25):+.2f}; {np.percentile(dc, 75):+.2f}], '
               f'média {dc.mean():+.2f}%; por instância (média das sementes): ' + ', '.join(f'{x:+.2f}' for x in mi))
    out.append(f'  ganho de P sob M3 (B − A), mediana dos pares: {np.median(dp3):+.1f} p.p.')
    for nm in NOMES:
        ka = [cel[('A', nm, s)]['P3'] for s in SEMENTES if ('A', nm, s) in cel]; kb = [cel[('B', nm, s)]['P3'] for s in SEMENTES if ('B', nm, s) in cel]
        ka2 = [cel[('A', nm, s)]['P2'] for s in SEMENTES if ('A', nm, s) in cel]; kb2 = [cel[('B', nm, s)]['P2'] for s in SEMENTES if ('B', nm, s) in cel]
        out.append(f'    {nm.replace("CCVRPTW-ST_", ""):16s} P média sob M3: A {np.mean(ka):.3f}  B {np.mean(kb):.3f}   |  sob M2(0,27): A {np.mean(ka2):.3f}  B {np.mean(kb2):.3f}')
    # acrescentado em 17/11 depois de ver o resultado (diário Entrada 142): distribuição de P por braço, para interpretar o M2
    for braco in ('A', 'B'):
        for nome in ('P3', 'P2'):
            p = np.array([cel[k][nome] for k in cel if k[0] == braco])
            marg = int(np.sum((p >= 0.85) & (p < ALPHA)))
            out.append(f'  braço {braco}, {nome}: P mínimo {p.min():.3f}, quartis {np.percentile(p, 25):.3f} / {np.median(p):.3f} / {np.percentile(p, 75):.3f}; '
                       f'{int(np.sum(p < ALPHA))} planos abaixo de {ALPHA} (dos quais {marg} entre 0,85 e {ALPHA})')
    resumo['por_plano'] = [dict(braco=k[0], instancia=k[1], semente=k[2], custo=v['custo'], P3=v['P3'], P2=v['P2']) for k, v in sorted(cel.items())]
    # divisão por largura de janela, acrescentada em 17/11 depois de ver os totais (exploratória; diário Entrada 142)
    def janela(nm):
        return 'narrow' if nm.split('_')[1].endswith('1') else 'wide'

    def estrato(braco, g):
        ks = [k for k in cel if k[0] == braco and (g is None or janela(k[1]) == g)]
        p3 = np.array([cel[k]['P3'] for k in ks]); p2 = np.array([cel[k]['P2'] for k in ks])
        return dict(n=len(ks), c3=100 * float(np.mean(p3 >= ALPHA)), P3=float(p3.mean()), c2=100 * float(np.mean(p2 >= ALPHA)), P2=float(p2.mean()))

    def custo_estrato(g):
        v = np.array([100 * (cel[('B', nm, sd)]['custo'] / cel[('A', nm, sd)]['custo'] - 1) for nm, sd in pares if g is None or janela(nm) == g])
        return float(np.median(v)), float(np.percentile(v, 25)), float(np.percentile(v, 75)), float(v.mean())

    resumo['por_janela'] = {}
    for braco in ('A', 'B'):
        for g in (None, 'narrow', 'wide'):
            e = estrato(braco, g)
            if braco == 'B':
                e['custo'] = dict(zip(('mediana', 'q1', 'q3', 'media'), custo_estrato(g)))
            resumo['por_janela'][f'{braco}|{g or "all"}'] = e
            out.append(f'  janela {g or "todas":6s} braço {braco}: {e["n"]} planos; cumpre M3 {e["c3"]:5.1f}% (P média {e["P3"]:.3f}); cumpre M2 {e["c2"]:5.1f}% (P média {e["P2"]:.3f})'
                       + (f'; custo do B contra o A: mediana {e["custo"]["mediana"]:+.2f}% [{e["custo"]["q1"]:+.2f}; {e["custo"]["q3"]:+.2f}], média {e["custo"]["media"]:+.2f}%' if braco == 'B' else ''))
    txt = '\n'.join(out); print(txt)
    open(os.path.join(OUTD, 'relatorio_b4.txt'), 'w', encoding='utf-8').write(txt + '\n')
    json.dump(resumo, open(os.path.join(OUTD, 'resumo_b4.json'), 'w', encoding='utf-8'), indent=1)
    linhas = []
    for braco, rot in (('A', r'Assumed, $0.10$'), ('B', r'Estimated from data, $0.27$')):
        for g, glab in ((None, 'all'), ('narrow', 'narrow'), ('wide', 'wide')):
            e = resumo['por_janela'][f'{braco}|{g or "all"}']
            if braco == 'A':
                cs = '--'
            else:
                cc_ = e['custo']; cs = '%+.2f\\%% (%+.2f; %+.2f)' % (cc_['mediana'], cc_['q1'], cc_['q3'])
            linhas.append('%s & %s & %d & %.1f & %.3f & %.1f & %s \\\\' % (rot if g is None else '', glab, e['n'], e['c3'], e['P3'], e['c2'], cs))
        if braco == 'A':
            linhas.append(r'\midrule')
    L = [r'\begin{table}[htbp]', r'\centering',
         r'\caption{Planning with the travel-time variability estimated from observed trips. The coefficient of variation is the mean observed CV by trip-duration band of the legs of the chosen plans (New York taxis, $0.27$), '
         r'used for the buffers and the certification in the planned-with-data plans, against the assumed $0.10$. Six medium-uncertainty instances with $n=50$, five seeds; chosen plan: cheapest certified rung of '
         r'$q\in\{0.70,0.80,0.90,0.95\}$ (\mCAL{}); 1{,}000 common scenarios per instance drawn from the observed residuals and day factor, with the observed CV by duration band (M3) or with CV $0.27$ (M2). '
         r'Complies: \% of plans with $\hat P\ge0.90$; $\bar P$: mean probability. Cost: change of the planned-with-data plan against the assumed-CV plan, median (quartiles) over the pairs. '
         r'Windows: narrow, families C1, R1 and RC1; wide, families C2, R2 and RC2 (a split added after seeing the totals).}',
         r'\label{tab:b4}', r'\footnotesize\setlength{\tabcolsep}{3pt}', r'\begin{tabular}{@{}l l r r r r r@{}}', r'\toprule',
         r'Planned with CV & Windows & Plans & Compl. M3 (\%) & $\bar P$ M3 & Compl. M2 (\%) & Cost vs assumed \\', r'\midrule'] + linhas + \
        [r'\bottomrule', r'\end{tabular}', r'\end{table}']
    open(os.path.join(OUTD, 'tab_b4.tex'), 'w', encoding='utf-8').write('\n'.join(L) + '\n')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--preparar', action='store_true'); ap.add_argument('--rodar', action='store_true')
    ap.add_argument('--analisar', action='store_true'); ap.add_argument('--workers', type=int, default=3)
    a = ap.parse_args()
    if a.preparar:
        preparar()
    if a.rodar:
        rodar(a.workers)
    if a.analisar:
        analisar()
