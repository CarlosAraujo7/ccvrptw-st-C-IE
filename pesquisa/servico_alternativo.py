#!/usr/bin/env python3
"""
A2 do ciclo C1 (diário, Entradas 121-123, 16/11/2026): métricas de nível de serviço alternativas à probabilidade CONJUNTA do plano.

Pergunta de revisor: "por que o nível de serviço conjunto do plano inteiro é a métrica certa? O que muda com nível por cliente,
fração de entregas no prazo ou atraso esperado?" Sem reotimizar: avalia o plano escolhido (o mais barato certificado, como em
validar_teoria_buffer.py) de cada método nos MESMOS 1.000 cenários comuns, no benchmark e em São Paulo, incerteza MEDIUM e HIGH.

Desenho fixado antes de calcular (métricas por plano; depois, média ou % sobre os planos de cada método):
  P_conj     = Pr[todos os clientes no prazo]            (a métrica do artigo)      → % de planos com P_conj ≥ α
  mín π_j    = pior cliente                                                           → % de planos com mín π_j ≥ α (nível POR CLIENTE)
  média π_j  = fração esperada de entregas no prazo
  atraso     = E[max(0, chegada − prazo − δ)] por cliente (min)
  nº atrasos = Σ_j (1 − π_j) clientes atrasados esperados por plano
Métodos: DET+C, BUF0.90+C, BUF0.95+C, CAL e SIM+C. α = 0,90; δ = 5 min.
Saídas: resultados_servico_alternativo/{servico_alternativo.json, relatorio_servico_alternativo.txt} e a tabela LaTeX (--tex).
"""
import argparse, glob, json, os, sys
from collections import defaultdict
import numpy as np

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, AQUI)
import validar_teoria_buffer as vtb   # noqa: E402  (mesmo contexto, cenários e escolha do plano)
ac, cc = vtb.ac, vtb.cc
BASE = vtb.BASE
OUT = os.path.join(BASE, 'pesquisa', 'resultados_servico_alternativo')
TEX_PADRAO = r'C:\Users\vinic\OneDrive\Documentos\Revisão Cold Chain\tables\tab_servico_alt.tex'
ALPHA, DELTA = vtb.ALPHA, vtb.DELTA
NOME = {'C-DETC': r'\mDETC{}', 'C-BUF90C': r'\mBUF{0.90}', 'C-BUF95C': r'\mBUF{0.95}', 'CAL': r'\mCAL{}', 'C-SIMC': r'\mSIM{}'}
ORDEM = ['C-DETC', 'C-BUF90C', 'C-BUF95C', 'CAL', 'C-SIMC']


def avalia(R, D, inst, T, tt):
    cust = inst['customers']; K = tt.shape[0]
    ont, late = [], []
    for r, d in zip(R, D):
        t = np.full(K, float(d)); prev = 0
        for c in r:
            t = t + T[prev][c] * tt[:, c]
            ont.append((t - cust[c - 1]['due_date_min']) <= DELTA)
            late.append(np.maximum(0.0, t - cust[c - 1]['due_date_min'] - DELTA))
            t = np.maximum(t, cust[c - 1]['ready_time_min']) + cust[c - 1]['service_time_min']; prev = c
    ont = np.array(ont).T; late = np.array(late).T        # K × clientes
    pi = ont.mean(axis=0)
    return dict(P=float(ont.all(axis=1).mean()), mn=float(pi.min()), mean=float(pi.mean()),
                atraso=float(late.mean()), nlate=float((1 - pi).sum()))


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--tex', default=TEX_PADRAO); a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    L = defaultdict(lambda: defaultdict(list))
    for grupo, res, idir in vtb.GRUPOS:
        ac.INST_DIR = idir; ac._cache.clear()
        nomes = sorted({os.path.basename(f).split('__')[0] for f in glob.glob(os.path.join(res, 'C-SIMC', '*_solutions.json'))})
        nomes = [n for n in nomes if n.upper().endswith(('MEDIUM', 'HIGH'))]
        for nome in nomes:
            inst, T, sc, _ = ac._ctx(nome, 1); tt = sc['tt']
            for sd in range(10):
                esc = {}
                for al in ['C-DETC', 'C-BUF70C', 'C-BUF80C', 'C-BUF90C', 'C-BUF95C', 'C-SIMC']:
                    f = os.path.join(res, al, f'{nome}__{al}__seed{sd}_solutions.json')
                    if os.path.exists(f):
                        e = vtb.escolha(json.load(open(f, encoding='utf-8')))
                        if e:
                            esc[al] = e
                lad = [esc[al] for al, _ in vtb.ESCADA if al in esc]
                if lad:
                    esc['CAL'] = min(lad, key=lambda x: x['f1_cost_BRL'])
                for al in ORDEM:
                    if al in esc:
                        L[grupo][al].append(avalia(esc[al]['routes'], esc[al]['partida_min'], inst, T, tt))
            print(grupo, nome, flush=True)
    res_json, linhas, rel = {}, [], []
    for grupo, titulo in (('benchmark', r'Benchmark (medium and high uncertainty)'), ('SP', r'S\~ao Paulo (medium and high uncertainty)')):
        linhas.append(r'\midrule' if linhas else '')
        linhas.append(r'\multicolumn{8}{@{}l}{\textit{' + titulo + r'}} \\')
        for al in ORDEM:
            X = L[grupo][al]
            if not X:
                continue
            v = dict(planos=len(X), conj=100 * np.mean([x['P'] >= ALPHA for x in X]), cli=100 * np.mean([x['mn'] >= ALPHA for x in X]),
                     media=100 * np.mean([x['mean'] for x in X]), pior=100 * np.median([x['mn'] for x in X]),
                     atraso=np.mean([x['atraso'] for x in X]), nlate=np.mean([x['nlate'] for x in X]))
            res_json[f'{grupo}_{al}'] = v
            rel.append(f"{grupo:9s} {al:9s} planos {v['planos']:4d}  P_conj>=α {v['conj']:5.1f}%  mín π>=α {v['cli']:5.1f}%  média π {v['media']:5.2f}%  "
                       f"pior cliente (mediana) {v['pior']:5.1f}%  atraso {v['atraso']:.3f} min/cliente  atrasos esperados {v['nlate']:.2f}/plano")
            linhas.append(f"{NOME[al]} & {v['planos']} & {v['conj']:.1f} & {v['cli']:.1f} & {v['media']:.2f} & {v['pior']:.1f} & {v['atraso']:.3f} & {v['nlate']:.2f} " + r'\\')
    open(os.path.join(OUT, 'relatorio_servico_alternativo.txt'), 'w', encoding='utf-8').write('\n'.join(rel) + '\n')
    json.dump(res_json, open(os.path.join(OUT, 'servico_alternativo.json'), 'w', encoding='utf-8'), indent=1)
    tex = r"""\begin{table}[htbp]
\centering
\caption{Alternative service-level measures for the chosen plan of each method (the cheapest certified plan), on the 1{,}000 common out-of-sample scenarios. Joint: \% of plans with $\Pr[\text{all customers on time}]\ge0.90$ (the measure of the paper). Per customer: \% of plans in which every customer is on time with probability at least $0.90$. Mean $\pi_j$: expected \% of deliveries on time. Worst $\pi_j$: median over plans of the least punctual customer (\%). Tardiness: expected minutes beyond the due date plus $\delta=5$\,min, per customer. Late: expected number of late customers per plan. Plans: chosen plans, one per run with a certified plan (360 runs on the benchmark and 160 in S~ao Paulo; the remaining runs have none).}
\label{tab:service-alt}
\small
\setlength{\tabcolsep}{3.5pt}
\begin{tabular}{@{}l r rr r r r r@{}}
\toprule
Method & Plans & Joint & Per customer & Mean $\pi_j$ & Worst $\pi_j$ & Tardiness & Late \\
\midrule
""" + '\n'.join(x for x in linhas if x) + r"""
\bottomrule
\end{tabular}
\end{table}
"""
    open(a.tex, 'w', encoding='utf-8').write(tex)
    print('\n'.join(rel))


if __name__ == '__main__':
    main()
