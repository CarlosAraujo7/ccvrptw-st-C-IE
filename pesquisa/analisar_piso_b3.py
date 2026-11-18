#!/usr/bin/env python3
"""
B3 repetido com limite inferior de temperatura (17/11/2026, diário Entrada 153). PÓS-HOC.

Contexto. Na Entrada 138 o compartimento único no set-point mais frio (Z1c) bateu as duas zonas nas três medidas, e a Entrada 141 mostrou
que isso vinha de o modelo não penalizar o excesso de frio. Aqui o modelo ganha um limite inferior (src/thermal_route.py, variável
CCVRP_PISO_FRACAO = 1): um produto refrigerado cuja temperatura cai abaixo do ponto de congelamento antes da entrega perde o valor.
As buscas Z1c, Z = 2 e Z = 3 (C-BUF90C, n = 50, 18 instâncias x 5 sementes) foram refeitas com o limite (results_c_b3_lb). D2, D4 e o
compartimento único no set-point dominante (MIX) não mudam: nenhuma das suas zonas fica abaixo de -1 °C, então nenhum refrigerado
chega ao ponto de congelamento (conferido aqui: congelado = 0 em todas as soluções gravadas, e no piloto de 17/11 a frente do D2 saiu idêntica).
O controle do Z = 2 sem o limite, com o código atual, reproduziu o B2 gravado (frente idêntica no piloto).

Painel A: o modelo sem limite inferior (resultados gravados em 10, 16 e 17/11). Painel B: com o limite inferior, buscas refeitas.
Medidas, por (instância, semente) com plano certificado fora da amostra (P >= alfa) na variante e na base MIX: diferença da melhor perda de
qualidade (p.p.; a do painel B já conta o valor perdido por congelamento), do melhor custo (%) e do melhor CO2 (%) contra MIX; mediana
sobre os pares. Também: % de execuções com plano certificado; % do valor entregue acima do limite legal da própria categoria (Legal);
% do valor entregue dos refrigerados cuja temperatura cai abaixo do ponto de congelamento antes da entrega (Congelado; cenário nominal,
tau_p = 10 h, partida no set-point da própria categoria).
Sensibilidade (só avaliação, sem rodar de novo): as frentes gravadas do Z1c e do Z = 2 reavaliadas com fração perdida de 0,1 a 1 (--sens).

Uso: python pesquisa/analisar_piso_b3.py [--tex arquivo.tex] [--out pasta] [--workers 2] [--sens]
Saída: pesquisa/resultados_c_b3_lb/relatorio_piso_b3.txt, resumo_piso_b3.json, tab_frota_dedicada.tex
"""
import argparse, glob, json, os, sys
from multiprocessing import Pool
import numpy as np
from scipy.stats import wilcoxon

BASE = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
sys.path.insert(0, os.path.join(BASE, 'src')); sys.path.insert(0, os.path.join(BASE, 'pesquisa'))
from equivalencia import hodges_lehmann       # noqa: E402
from analisar_refinamento import holm         # noqa: E402
import analisar_frota_dedicada as af          # noqa: E402

ALPHA, SEM = af.ALPHA, af.SEM
OUT = os.path.join(BASE, 'pesquisa', 'resultados_c_b3_lb')
LB = os.environ.get('PISO_B3_LB', 'results_c_b3_lb')
A = [('MIX', os.path.join('results_c_b2', 'C-BUF90C-MIX')), ('Z=2', os.path.join('results_c_b2', 'C-BUF90C_Z2')),
     ('Z=3', os.path.join('results_c_b2', 'C-BUF90C_Z3')), ('D2', os.path.join('results_c_b3', 'C-BUF90C-D2')),
     ('D4', os.path.join('results_c_b3', 'C-BUF90C-D4')), ('Z1c', os.path.join('results_c_b3', 'C-BUF90C-Z1c'))]
B = [('Z=2', os.path.join(LB, 'C-BUF90C-Z2')), ('Z=3', os.path.join(LB, 'C-BUF90C-Z3')), ('Z1c', os.path.join(LB, 'C-BUF90C-Z1c'))]


def _arquivos(pasta):
    return sorted(f for f in glob.glob(os.path.join(BASE, pasta, '*_solutions.json'))
                  if '_050_' in os.path.basename(f) and int(f.split('seed')[-1].split('_')[0]) in SEM)


def _job_cong(path):
    """Fração do valor entregue dos refrigerados que cai abaixo do ponto de congelamento antes da entrega (nominal, tau_p = 10 h)."""
    os.environ['CCVRP_PISO_FRACAO'] = '1'
    import nsga2_simheuristic as ns
    import thermal_route as tr
    m = json.load(open(path_metrics(path), encoding='utf-8'))
    inst = ns.load_instance(os.path.join(BASE, 'instances', 'generated', m['instance'] + '.json'))
    _, ttime = ns.build_matrices(inst); cust = inst['customers']
    mz = bool(m.get('use_multizone')); nz = m['n_zones'] if mz else 1
    nom = tr.cenario_nominal(inst, n_zones=nz)
    ci = ns._cat_index_map(inst) if mz else None
    fr = []
    for s in json.load(open(path, encoding='utf-8')):
        if not s['feasible']:
            continue
        zones = ns.assign_zones(s['routes'], s['zone_block'], inst, ci) if mz else None
        cong = tot = 0.0
        for i, (r, d0) in enumerate(zip(s['routes'], s['partida_min'])):
            if mz:
                zi = zones[i]
                zoc = [zi['zone_of_category'][cust[c - 1]['product_category']] for c in r]
                d = tr.route_quality(r, inst, ttime, nom, zoc, zi['zone_temp'], per_customer_kinetics=True, detalhe=True, t0_min=d0)
            else:
                d = tr.route_quality(r, inst, ttime, nom, per_customer_kinetics=True, detalhe=True, t0_min=d0)
            v = d['values']; esf = np.array([c in tr.PISO_C for c in d['cat']])
            cong += float((d['cong'][0] * v)[esf].sum()); tot += float(v[esf].sum())
        if tot > 0:
            fr.append(cong / tot)
    return fr


def path_metrics(path):
    return path.replace('_solutions.json', '_metrics.json')


def congelado(pasta, workers):
    with Pool(workers) as pool:
        res = pool.map(_job_cong, _arquivos(pasta), chunksize=4)
    v = [100 * np.mean(r) for r in res if r]
    return (float(np.mean(v)) if v else float('nan')), len(v)


def linha(D, n, base, leg, cong):
    ks = af.pares(base, D[n])
    if not ks:
        return None
    dq = np.array([100 * (D[n][k]['q'] - base[k]['q']) for k in ks])
    dc = np.array([100 * (D[n][k]['c'] / base[k]['c'] - 1) for k in ks])
    de = np.array([100 * (D[n][k]['e'] / base[k]['e'] - 1) for k in ks])
    return dict(pares=len(ks), ok_mix=100 * np.mean([base[k]['ok'] for k in D[n] if k in base]), ok=100 * np.mean([D[n][k]['ok'] for k in D[n]]),
                dq=float(np.median(dq)), dc=float(np.median(dc)), de=float(np.median(de)),
                dc_q1=float(np.percentile(dc, 25)), dc_q3=float(np.percentile(dc, 75)), legal=leg[0], cong=cong[0])


def direto(X, Y):
    """Comparação por instância (HL, Wilcoxon, Holm nas três medidas) da variante Y contra X."""
    res = []
    for nm, campo, rel in (('Δ perda (p.p.)', 'q', False), ('Δ custo (%)', 'c', True), ('Δ CO2 (%)', 'e', True)):
        dd, ks = af.difs(X, Y, campo, rel)
        I = sorted(dd)
        if len(I) < 3:
            continue
        v = np.array([dd[i] for i in I])
        hl, lo, hi = hodges_lehmann(v, 0.95)
        p = wilcoxon(v).pvalue if np.any(np.abs(v) > 1e-12) else 1.0
        res.append([nm, len(I), hl, lo, hi, p, int(np.sum(v < -1e-12)), int(np.sum(v > 1e-12))])
    for r, pa in zip(res, holm([r[5] for r in res])):
        r.append(pa)
    return res


def _nom_fracao(args):
    """Reavalia (sem gravar) as soluções certificadas de um arquivo no cenário NOMINAL com fração perdida phi (a certificação vem do arquivo
    *_oosC.json gravado; ela não depende de phi): devolve (instância, semente, melhor q, melhor custo, melhor CO2, há certificada)."""
    path, phi = args[0], args[1]
    tau = args[2] if len(args) > 2 else None             # constante de tempo do produto (h); None = a do estudo, 10 h
    os.environ['CCVRP_PISO_FRACAO'] = repr(float(phi))
    import nsga2_simheuristic as ns
    import thermal_route as tr
    tr.TAU_P_PISO_H = float(tau) if tau else 10.0        # o worker é reaproveitado: sempre refazer o valor
    m = json.load(open(path_metrics(path), encoding='utf-8'))
    oos = {s['sol_id']: s for s in json.load(open(path.replace('_solutions.json', '_oosC.json'), encoding='utf-8'))['solutions']}
    inst = ns.load_instance(os.path.join(BASE, 'instances', 'generated', m['instance'] + '.json'))
    _, ttime = ns.build_matrices(inst); cust = inst['customers']
    mz = bool(m.get('use_multizone')); nz = m['n_zones'] if mz else 1
    nom = tr.cenario_nominal(inst, n_zones=nz)
    ci = ns._cat_index_map(inst) if mz else None
    B_ = []
    for s in json.load(open(path, encoding='utf-8')):
        if not s['feasible'] or oos[s['sol_id']]['p_on'] < ALPHA:
            continue
        zones = ns.assign_zones(s['routes'], s['zone_block'], inst, ci) if mz else None
        acc = tot = 0.0
        for i, (r, d0) in enumerate(zip(s['routes'], s['partida_min'])):
            zi = zones[i]
            zoc = [zi['zone_of_category'][cust[c - 1]['product_category']] for c in r]
            q = float(tr.route_quality(r, inst, ttime, nom, zoc, zi['zone_temp'], per_customer_kinetics=True, t0_min=d0)[0])
            v = sum(cust[c - 1].get('load_value_BRL', 1.0) for c in r)
            acc += q * v; tot += v
        B_.append((float(1 - acc / tot), s['f1_cost_BRL'], s['f3_co2_kg']))
    return (m['instance'], m['seed'], min((b[0] for b in B_), default=np.nan), min((b[1] for b in B_), default=np.nan),
            min((b[2] for b in B_), default=np.nan), bool(B_))


def sensibilidade(workers, fracoes=(0.0, 0.1, 0.25, 0.5, 1.0)):
    """Frentes gravadas do Z1c e do Z = 2 (buscas sem limite inferior) reavaliadas no cenário nominal com fração perdida phi: o phi em que o Z1c deixa de vencer."""
    out = []
    for phi in fracoes:
        ref = {}
        for nome, pasta in (('Z1c', A[5][1]), ('Z=2', A[1][1])):
            tarefas = [(f, phi) for f in _arquivos(pasta)]
            with Pool(workers) as pool:
                r = pool.map(_nom_fracao, tarefas, chunksize=2)
            ref[nome] = {(i, s): dict(ok=ok, q=q, c=c, e=e) for i, s, q, c, e, ok in r}
        dd, ks = af.difs(ref['Z=2'], ref['Z1c'], 'q', False)
        v = np.array([dd[i] for i in sorted(dd)])
        hl, lo, hi = hodges_lehmann(v, 0.95)
        out.append(dict(phi=phi, n=len(v), hl=hl, lo=lo, hi=hi, mais_zero=int(np.sum(v > 1e-12)), menos_zero=int(np.sum(v < -1e-12))))
    return out


def sens_tau(workers, taus=(10.0, 4.0, 2.0, 1.0)):
    """Planos REFEITOS com o limite (Z1c e Z = 2) reavaliados no cenário nominal, com perda total, a constantes de tempo do produto de 10 h (a do estudo)
    a 1 h: quanto da vantagem restante do Z1c depende de o produto refrigerado levar horas para congelar numa zona a -18 °C (a cinética de degradação
    segue a temperatura do ar do compartimento; o congelamento segue a do produto). Não refaz a busca: é um limite superior do efeito de tau."""
    out = []
    for tau in taus:
        ref = {}
        for nome, pasta in (('Z1c', B[2][1]), ('Z=2', B[0][1])):
            tarefas = [(f, 1.0, tau) for f in _arquivos(pasta)]
            with Pool(workers) as pool:
                r = pool.map(_nom_fracao, tarefas, chunksize=2)
            ref[nome] = {(i, s): dict(ok=ok, q=q, c=c, e=e) for i, s, q, c, e, ok in r}
        dd, ks = af.difs(ref['Z=2'], ref['Z1c'], 'q', False)
        v = np.array([dd[i] for i in sorted(dd)])
        hl, lo, hi = hodges_lehmann(v, 0.95)
        out.append(dict(tau=tau, n=len(v), hl=hl, lo=lo, hi=hi, mais_zero=int(np.sum(v > 1e-12)), menos_zero=int(np.sum(v < -1e-12)),
                        q_z1c=float(100 * np.nanmedian([x['q'] for x in ref['Z1c'].values()])), q_z2=float(100 * np.nanmedian([x['q'] for x in ref['Z=2'].values()]))))
    return out


def main(tex, out_dir, workers, sens):
    DA = {n: af.carregar(p) for n, p in A}
    DB = {n: af.carregar(p) for n, p in B}
    base = DA['MIX']
    os.makedirs(out_dir, exist_ok=True)
    out = ['B3 repetido com limite inferior de temperatura (pós-hoc, Entrada 153); C-BUF90C, n = 50, 18 instâncias, sementes 0-4; α = %.2f' % ALPHA,
           'execuções por variante (painel A): ' + ', '.join(f'{n} {len(d)}' for n, d in DA.items()),
           'execuções por variante (painel B, limite inferior): ' + ', '.join(f'{n} {len(d)}' for n, d in DB.items())]
    LEG, CON = {}, {}
    for n, p in A:
        LEG[('A', n)] = af.legal(p, workers); CON[('A', n)] = congelado(p, workers)
    for n, p in B:
        LEG[('B', n)] = af.legal(p, workers); CON[('B', n)] = congelado(p, workers)
    resumo = {'A': {}, 'B': {}, 'direto': {}, 'sens': []}
    for pn, D, vs in (('A', DA, [x[0] for x in A[1:]]), ('B', DB, [x[0] for x in B])):
        out.append(f'\nPainel {pn}:')
        for n in vs:
            r = linha(D, n, base, LEG[(pn, n)], CON[(pn, n)])
            if r is None:
                continue
            resumo[pn][n] = r
            out.append(f'  {n:4s} pares {r["pares"]:3d}  cumpre MIX {r["ok_mix"]:5.1f}% / {n} {r["ok"]:5.1f}%  Δ perda {r["dq"]:+6.2f} p.p.  Δ custo {r["dc"]:+6.2f}% '
                       f'[{r["dc_q1"]:+.1f}, {r["dc_q3"]:+.1f}]  Δ CO2 {r["de"]:+6.2f}%  legal {r["legal"]:.1f}%  congelado {r["cong"]:.1f}%')
    resumo['mix'] = dict(ok=100 * np.mean([v['ok'] for v in base.values()]), legal=LEG[('A', 'MIX')][0], cong=CON[('A', 'MIX')][0])
    out.append(f'  MIX  cumpre {resumo["mix"]["ok"]:.1f}%  legal {resumo["mix"]["legal"]:.1f}%  congelado {resumo["mix"]["cong"]:.1f}%')
    out.append('\nCongelado nas variantes que não mudaram (deve ser 0): ' + ', '.join(f'{n} {CON[("A", n)][0]:.2f}%' for n in ('D2', 'D4')))
    out.append('\nComparação direta, por instância (HL [IC 95%]; Wilcoxon; Holm nas 3 medidas):')
    pares_dir = [('A', 'Z1c contra Z=2 (painel A, sem limite inferior)', DA['Z=2'], DA['Z1c']),
                 ('B', 'Z1c contra Z=2 (painel B, com limite inferior)', DB['Z=2'], DB['Z1c']),
                 ('B', 'Z=3 contra Z=2 (painel B)', DB['Z=2'], DB['Z=3']),
                 ('AB', 'Z=2 com limite inferior contra Z=2 sem ele (o preço de não esfriar demais)', DA['Z=2'], DB['Z=2']),
                 ('AB', 'Z1c com limite inferior contra Z1c sem ele', DA['Z1c'], DB['Z1c'])]
    for pn, rot, X, Y in pares_dir:
        res = direto(X, Y)
        out.append(f'  {rot}:')
        for nm, ni, hl, lo, hi, p, mn, mx, pa in res:
            out.append(f'      {nm:15s} n={ni:2d}  HL {hl:+7.3f} [{lo:+7.3f}, {hi:+7.3f}]  p={p:.3g}  p_Holm={pa:.3g}  ({mn} menores, {mx} maiores)')
            resumo['direto'][f'{rot}|{nm}'] = dict(n=ni, hl=hl, lo=lo, hi=hi, p=p, p_holm=pa)
    if sens:
        out.append('\nSensibilidade (só avaliação das frentes gravadas, sem limite inferior na busca): Δ perda do Z1c contra o Z=2 (HL, p.p.), com fração perdida phi:')
        for s in sensibilidade(workers):
            resumo['sens'].append(s)
            out.append(f'  phi = {s["phi"]:.2f}: HL {s["hl"]:+.3f} [{s["lo"]:+.3f}, {s["hi"]:+.3f}] em {s["n"]} instâncias ({s["mais_zero"]} com Z1c pior, {s["menos_zero"]} com Z1c melhor)')
        out.append('\nSensibilidade à constante de tempo do produto (planos refeitos com o limite, só reavaliados, perda total, cenário nominal): Δ perda do Z1c contra o Z=2 (HL, p.p.):')
        resumo['sens_tau'] = []
        for s in sens_tau(workers):
            resumo['sens_tau'].append(s)
            out.append(f'  tau = {s["tau"]:4.1f} h: HL {s["hl"]:+.3f} [{s["lo"]:+.3f}, {s["hi"]:+.3f}] em {s["n"]} instâncias ({s["mais_zero"]} com Z1c pior, {s["menos_zero"]} com Z1c melhor); '
                       f'perda mediana do melhor plano certificado: Z1c {s["q_z1c"]:.2f}%, Z=2 {s["q_z2"]:.2f}%')
    txt = '\n'.join(out); print(txt)
    open(os.path.join(out_dir, 'relatorio_piso_b3.txt'), 'w', encoding='utf-8').write(txt + '\n')
    json.dump(resumo, open(os.path.join(out_dir, 'resumo_piso_b3.json'), 'w', encoding='utf-8'), indent=1, default=float)
    tabela(resumo, tex or os.path.join(out_dir, 'tab_frota_dedicada.tex'))


NOME = {'Z1c': r'$Z=1$, coldest set-point', 'Z=2': r'$Z=2$ zones', 'Z=3': r'$Z=3$ zones', 'D2': r'D2: frozen / chilled', 'D4': r'D4: one per category'}


def tabela(resumo, destino):
    def lin(n, r):
        return (f'{NOME[n]} & {r["pares"]} & {r["ok"]:.0f} & {r["dq"]:+.2f} & {r["dc"]:+.2f} [{r["dc_q1"]:+.1f}, {r["dc_q3"]:+.1f}] & {r["de"]:+.2f} & '
                f'{r["legal"]:.1f} & {r["cong"]:.1f} \\\\')
    L = [r'\begin{table}[htbp]', r'\centering',
         r'\caption{Dedicated vehicles per product class versus multi-compartment vehicles, \mBUF{0.90}, $n=50$, 18 instances $\times$ 5 seeds. '
         r'Reference: a single compartment carrying the same mixed load at the dominant set-point (Z=1, as in Table~S15). Z=1, coldest set-point: the same mixed routes, '
         r'with the single compartment run at the coldest set-point of any category on board. Dedicated fleets: each route carries one class of product and runs one zone '
         r'at the set-point of the coldest category present (D2: frozen / chilled; D4: one class per category). Panel A: the model of the first version, which does not '
         r'penalize over-cooling. Panel B: the same searches repeated with a lower temperature limit, added after the result of panel A: a refrigerated product whose temperature '
         r'falls below its freezing point ($-1.5$\,\textdegree C for meat, $-0.5$\,\textdegree C for dairy and produce) before delivery loses its value; the dedicated fleets and the '
         r'single compartment at the dominant set-point never reach those limits, so their rows are the same in both panels. Differences are medians over the pairs in which both '
         r'variants have a plan certified out of sample, for the best value of each objective among compliant solutions (in panel B the loss includes the value lost by freezing; '
         r'interquartile range of the cost difference in brackets). Compliance: \% of runs with a certified plan. Legal: mean share of delivered value above the legal limit of its own '
         r'category. Frozen: mean share of the delivered value of refrigerated products that falls below its freezing point before delivery (nominal scenario).}',
         r'\label{tab:b3}', r'\footnotesize\setlength{\tabcolsep}{3pt}', r'\begin{tabular}{@{}l r r r r r r r@{}}', r'\toprule',
         r'Variant & Pairs & Compl. (\%) & $\Delta$ loss (p.p.) & $\Delta f_1$ (\%) & $\Delta f_3$ (\%) & Legal (\%) & Frozen (\%) \\', r'\midrule',
         r'\multicolumn{8}{@{}l}{\textit{A. Over-cooling not penalized}} \\']
    for n in ('Z1c', 'Z=2', 'Z=3', 'D2', 'D4'):
        if n in resumo['A']:
            L.append(lin(n, resumo['A'][n]))
    L += [r'\midrule', r'\multicolumn{8}{@{}l}{\textit{B. Lower temperature limit (freezing point), searches repeated}} \\']
    for n in ('Z1c', 'Z=2', 'Z=3'):
        if n in resumo['B']:
            L.append(lin(n, resumo['B'][n]))
    m = resumo['mix']
    L += [r'\midrule', f'$Z=1$, mixed load & -- & {m["ok"]:.0f} & -- & -- & -- & {m["legal"]:.1f} & {m["cong"]:.1f} \\\\', r'\bottomrule', r'\end{tabular}', r'\end{table}']
    open(destino, 'w', encoding='utf-8').write('\n'.join(L) + '\n')
    print('tabela:', destino)


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--tex', default=None)
    ap.add_argument('--out', default=OUT)
    ap.add_argument('--workers', type=int, default=2)
    ap.add_argument('--sens', action='store_true')
    a = ap.parse_args()
    main(a.tex, a.out, a.workers, a.sens)
