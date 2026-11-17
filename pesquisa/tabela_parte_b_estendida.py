#!/usr/bin/env python3
"""Tabela S25 (A3, diário Entradas 123-125): sensibilidade térmica estendida, a partir de resultados_parte_b_estendida/parte_b.csv
e (se existir) resultados_parte_b_estendida_poshoc/parte_b.csv. Agregação como na Tabela S19: E[perda] = média das instâncias;
ρ_min = mínimo; arrependimento = máximo; violação = média; P̂ = média. Condições sem variação (ρ indefinido) saem com '--'."""
import argparse, csv, json, os
from collections import OrderedDict
import numpy as np

AQUI = os.path.dirname(os.path.abspath(__file__))
TEX_PADRAO = r'C:\Users\vinic\OneDrive\Documentos\Revisão Cold Chain\tables\tab_parteb_estendida.tex'
CSVS = [os.path.join(AQUI, 'resultados_parte_b_estendida', 'parte_b.csv'),
        os.path.join(AQUI, 'resultados_parte_b_estendida_poshoc', 'parte_b.csv')]
ROT = {'base': 'Base case (reproduces Table S19)',
       'Ea x0.8 (rate at set-point fixed)': r'Activation energy $\times0.8$',
       'Ea x1.2 (rate at set-point fixed)': r'Activation energy $\times1.2$',
       'zero-order kinetics': 'Zero-order kinetics',
       'persistent thermostat bias (sigma x2)': r'Persistent thermostat bias, $\sigma_T\times2$',
       'quality threshold Q<0.80': r'Threshold: quality $<0.80$ (never reached)',
       'combined stress': r'Combined stress, threshold $0.80$',
       'quality threshold Q<0.95 (post hoc)': r'Threshold: quality $<0.95$ (post hoc)',
       'quality threshold Q<0.90 (post hoc)': r'Threshold: quality $<0.90$ (post hoc)',
       'combined stress, threshold 0.90 (post hoc)': r'Combined stress, threshold $0.90$ (post hoc)'}


def fmt(x, f):
    return '--' if x is None or (isinstance(x, float) and np.isnan(x)) else format(x, f)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--tex', default=TEX_PADRAO); a = ap.parse_args()
    por = OrderedDict()
    for path in CSVS:
        if os.path.exists(path):
            for r in csv.DictReader(open(path, encoding='utf-8')):
                por.setdefault(r['condicao'], []).append(r)
    linhas, res = [], {}
    for c, L in por.items():
        f = lambda k: np.array([float(x[k]) for x in L])
        rho = f('rho_q'); arr = f('arrep_q_abs')
        v = dict(n=len(L), q=100 * f('q_esp').mean(),
                 rho=None if np.all(np.isnan(rho)) else float(np.nanmin(rho)),
                 arrep=None if np.all(np.isnan(arr)) or np.all(np.isnan(rho)) else 100 * float(np.nanmax(arr)),
                 quebra=100 * f('quebra_esp').mean(), p=float(f('p_pontual').mean()))
        res[c] = v
        linhas.append(f"{ROT.get(c, c)} & {v['q']:.3f} & {fmt(v['rho'], '.4f')} & {fmt(v['arrep'], '.4f')} & {v['quebra']:.2f} & {v['p']:.3f} " + r'\\')
    json.dump(res, open(os.path.join(AQUI, 'resultados_parte_b_estendida', 'resumo.json'), 'w', encoding='utf-8'), indent=1)
    tex = r"""\begin{table}[htbp]
\centering
\caption{Extended sensitivity of the thermal conclusion (extends Table~S19): re-evaluation of the final fronts of \mDET{}, \mBUF{0.90} and \mSIM{} (18 instances with $n=50$; up to 150 distinct solutions each; same solutions, scenarios and measures as Table~S19; the base case reproduces its first row). Activation energy: scaled in every category, with the pre-exponential factor adjusted to keep the rate at the set-point. Zero-order kinetics: loss $=\min(1,\Lambda)$. Persistent bias: the thermostat deviation of a scenario is the same along the whole route. Threshold: the ``loss'' is the expected share of delivered value whose retained quality is below the threshold, a non-linear measure; the threshold of $0.80$ is never reached (intra-day losses are about 3\%), so its correlation is undefined (--), and the thresholds $0.95$ and $0.90$ were added post hoc after that result. $E[\text{loss}]$: mean expected loss; $\rho_{\min}$: smallest Spearman correlation, over instances, between nominal and expected loss; Regret: largest loss, over instances, of choosing by the nominal instead of the expected value within cost budgets; Breach: expected share of delivered value whose product temperature exceeded its modeled category threshold; $\hat P$: mean probability that no customer is served more than 5\,min late.}
\label{tab:partB-ext}
\small
\begin{tabular}{@{}p{0.40\textwidth} r r r r r@{}}
\toprule
Condition & $E[\text{loss}]$ (\%) & $\rho_{\min}$ & Regret (p.p.) & Breach (\%) & $\hat P$ \\
\midrule
""" + '\n'.join(linhas) + r"""
\bottomrule
\end{tabular}
\end{table}
"""
    open(a.tex, 'w', encoding='utf-8').write(tex)
    for c, v in res.items():
        print(f"{c:46s} n={v['n']:3d} E[loss]={v['q']:.3f}%  rho_min={fmt(v['rho'], '.4f')}  regret={fmt(v['arrep'], '.4f')} p.p.  breach={v['quebra']:.3f}%  P={v['p']:.3f}")


if __name__ == '__main__':
    main()
