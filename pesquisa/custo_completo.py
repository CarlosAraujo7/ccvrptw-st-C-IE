#!/usr/bin/env python3
"""
A1 do ciclo C1 (diário, Entrada 121, 16/11/2026): custo completo da confiabilidade.

Pergunta de revisor: "confiabilidade custa < 0,5%" vale só para f1 = custo fixo + R$/km. O que muda se as horas de rota
(motorista) e a energia de refrigeração forem precificadas? Sem reotimizar: usa as soluções guardadas.

Para cada execução (mesma seleção da Tabela S13): a solução mais barata que cumpre α fora da amostra em
results_c/<alg>/*_oosC.json; o número de rotas R vem de *_solutions.json (mesmo sol_id). Com f1 e f3 guardados:
    L = (f1 − c_fix·R − c_extra·max(0, R − K)) / c_km          (km totais)
    D = (f3 − e_km·L) / (e_el·P_ref)                            (horas totais de rota, cenário nominal)
Custo completo: f1 + c·D, com c = w + p·P_ref (w: R$/h do motorista, premissa; p: R$/kWh, premissa).
Verificação de consistência: D/R deve ficar entre 0,3 h e 16 h (horizonte) em todas as soluções.

Saídas: resultados_custo_completo/{custo_completo.json, relatorio_custo_completo.txt} e a tabela LaTeX (--tex).
"""
import argparse, glob, json, os, re
from collections import defaultdict
import numpy as np

BASE = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
OUT = os.path.join(BASE, 'pesquisa', 'resultados_custo_completo')
TEX_PADRAO = r'C:\Users\vinic\OneDrive\Documentos\Revisão Cold Chain\tables\tab_preco_completo.tex'
NIV = ['LOW', 'MEDIUM', 'HIGH']
W = [15.0, 30.0, 60.0]      # R$/h do motorista (premissa; ×0,5, ×1, ×2)
P_KWH = 1.0                 # R$/kWh (premissa)
src = open(os.path.join(BASE, 'src', 'nsga2_simheuristic.py'), encoding='utf-8').read()
EXTRA = float(re.search(r'EXTRA_VEHICLE_COST\s*=\s*([0-9.]+)', src).group(1))
_cache = {}


def params(nome):
    if nome not in _cache:
        d = json.load(open(os.path.join(BASE, 'instances', 'generated', nome + '.json'), encoding='utf-8'))
        fl, co = d['fleet'], d['co2_params']
        _cache[nome] = dict(K=fl['num_vehicles'], c_fix=fl['fixed_cost_BRL'], c_km=fl['variable_cost_BRL_km'],
                            e_km=co['alpha_kg_per_km'], e_el=co['beta_kg_per_kWh'], P=co['P_refrigeration_kW'])
    return _cache[nome]


def derivar(f1, f3, R, p):
    L = (f1 - p['c_fix'] * R - EXTRA * max(0, R - p['K'])) / p['c_km']
    D = (f3 - p['e_km'] * L) / (p['e_el'] * p['P'])
    return L, D


def melhores(pattern, alpha, checks):
    out = defaultdict(list)
    for f in glob.glob(pattern):
        d = json.load(open(f, encoding='utf-8')); m = d['meta']
        S = [s for s in d['solutions'] if s['feasible'] and (alpha is None or s['p_on'] >= alpha)]
        if not S:
            continue
        s = min(S, key=lambda x: x['f1'])
        sols = json.load(open(f.replace('_oosC.json', '_solutions.json'), encoding='utf-8'))
        R = len(next(x for x in sols if x['sol_id'] == s['sol_id'])['routes'])
        p = params(m['instance'])
        L, D = derivar(s['f1'], s['f3'], R, p)
        checks.append((D / R, L, R, m['instance']))
        out[m['instance']].append((s['f1'], D, R))
    return {k: np.mean(v, axis=0) for k, v in out.items()}


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--tex', default=TEX_PADRAO); a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    checks = []
    det = melhores(os.path.join(BASE, 'results_c', 'C-DET', '*_oosC.json'), None, checks)
    linhas, res, rel = [], {}, []
    for rot, alg, alpha in (('\\mSIM', 'C-SIMC', 0.90), ('\\mBUF{0.95}', 'C-BUF95C', 0.90), ('\\mSIM', 'C-SIMC95', 0.95)):
        m = melhores(os.path.join(BASE, 'results_c', alg, '*_oosC.json'), alpha, checks)
        for u in NIV:
            ks = [k for k in m if k.endswith(u) and k in det]
            dD = np.median([m[k][1] - det[k][1] for k in ks])
            d1 = np.median([100 * (m[k][0] / det[k][0] - 1) for k in ks])
            dr = np.median([m[k][2] - det[k][2] for k in ks])
            cel = []
            for w in W:
                c = w + P_KWH * 3.5
                cel.append(np.median([100 * ((m[k][0] + c * m[k][1]) - (det[k][0] + c * det[k][1])) / det[k][0] for k in ks]))
            res[f'{alg}_{alpha}_{u}'] = dict(n=len(ks), dD_h=dD, dR=dr, df1_pct=d1, dfull_pct=dict(zip(map(str, W), cel)))
            rel.append(f'{alg} alpha={alpha} {u:6s} n={len(ks)} ΔD={dD:+.2f} h  ΔR={dr:+.2f}  Δf1={d1:+.2f}%  Δ(f1+c·D)/f1_DET = ' +
                       '  '.join(f'w={w:.0f}: {x:+.2f}%' for w, x in zip(W, cel)))
            linhas.append(f"{alpha:.2f} & {rot}{{}} & {u.capitalize()} & {dD:+.2f} & {d1:+.2f} & " + ' & '.join(f'{x:+.2f}' for x in cel) + r' \\')
    r = np.array([c[0] for c in checks])
    bad = int(np.sum((r < 0.3) | (r > 16)))
    rel.insert(0, f'verificação: {len(r)} soluções; D/R min {r.min():.2f} h, mediana {np.median(r):.2f} h, máx {r.max():.2f} h; fora de [0,3; 16] h: {bad}')
    open(os.path.join(OUT, 'relatorio_custo_completo.txt'), 'w', encoding='utf-8').write('\n'.join(rel) + '\n')
    json.dump(res, open(os.path.join(OUT, 'custo_completo.json'), 'w', encoding='utf-8'), indent=1)
    cab = ' & '.join(f'$w={w:.0f}$' for w in W)
    tex = r"""\begin{table}[htbp]
\centering
\caption{Price of reliability with driver time and refrigeration energy priced (extends Table~S13). Cheapest out-of-sample-compliant solution versus the cheapest solution of \mDET{}, medians over the 18 instances of each uncertainty level. $\Delta D$: change in total route time (h per plan, nominal scenario, derived from the stored $f_1$ and $f_3$). $\Delta f_1$: change in the distance-based cost of the paper. Last three columns: change in $f_1+c\,D$ as a percentage of the \mDET{} value of $f_1$, with $c=w+1.0\times3.5$\,BRL/h, where $w$ is the driver cost in BRL/h and 3.5\,kW is the refrigeration power priced at 1.0\,BRL/kWh; $w$ and the energy price are assumptions of this study.}
\label{tab:price-full}
\small
\begin{tabular}{@{}c l l r r r r r@{}}
\toprule
 & & & & & \multicolumn{3}{c}{$\Delta(f_1+cD)$ (\% of $f_1^{\mathrm{DET}}$)} \\
\cmidrule(l){6-8}
$\alpha$ & Method & Level & $\Delta D$ (h) & $\Delta f_1$ (\%) & """ + cab + r""" \\
\midrule
""" + '\n'.join(linhas) + r"""
\bottomrule
\end{tabular}
\end{table}
"""
    open(a.tex, 'w', encoding='utf-8').write(tex)
    print('\n'.join(rel))


if __name__ == '__main__':
    main()
