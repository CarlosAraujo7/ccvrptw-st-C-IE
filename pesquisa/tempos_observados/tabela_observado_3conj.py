#!/usr/bin/env python3
"""
Tabela S20 com os TRÊS conjuntos (diário, Entrada 121, N2): benchmark, São Paulo e São Paulo com tempos dependentes do horário.
O texto da Seção 6.13 já usava as faixas dos três; tabela_observado.py imprimia só dois. Mesmos JSON de entrada
(resultados/resumo_reavaliacao_observada.json, resumo_recertificacao_observada.json, resumo_caminhoes.json).
Saída: --tex (padrão: tables/tab_observado.tex do pacote C&IE).
"""
import argparse, json, os

AQUI = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(AQUI, 'resultados')
TEX_PADRAO = r'C:\Users\vinic\OneDrive\Documentos\Revisão Cold Chain\tables\tab_observado.tex'
NOME = {'C-DETC': r'\mDETC{}', 'C-BUF90C': r'\mBUF{0.90}', 'C-BUF95C': r'\mBUF{0.95}', 'CAL': r'\mCAL{}',
        'C-SIMC': r'\mSIM{}', 'ACAL': r'\mACAL{}'}
CONJ = [('benchmark (média e alta)', 'Benchmark (medium and high uncertainty; 360 runs)'),
        ('SP', r'S\~ao Paulo (160 runs)'),
        ('SP com horário', r'S\~ao Paulo, time-dependent travel times (160 runs)')]


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--tex', default=TEX_PADRAO); a = ap.parse_args()
    R = json.load(open(os.path.join(RES, 'resumo_reavaliacao_observada.json'), encoding='utf-8'))
    C = json.load(open(os.path.join(RES, 'resumo_recertificacao_observada.json'), encoding='utf-8'))
    T = json.load(open(os.path.join(RES, 'resumo_caminhoes.json'), encoding='utf-8'))
    L = [r'\begin{table}[htbp]', r'\centering',
         r'\caption{Chosen plans under observed travel-time variability: percentage of runs whose chosen plan meets '
         r'$\alpha=0.90$ on 1{,}000 scenarios. M0: lognormal model of the paper. Taxi: instance CV with the shape and day '
         r'factor observed in New York taxi trips (M2). Truck: instance CV with the shape observed in the DT-CARGO truck trips. '
         r'Re-cert.: the stored front re-certified on 500 independent scenarios of the same observed model (runs without a '
         r're-certified plan are excluded). M3: CV observed in the taxi trips for the arc duration (20\% to 29\%).}',
         r'\label{tab:observed}', r'\small', r'\setlength{\tabcolsep}{4pt}',
         r'\begin{tabular}{@{}l r rr rr r@{}}', r'\toprule',
         r' & & \multicolumn{2}{c}{Taxi} & \multicolumn{2}{c}{Truck} & \\',
         r'\cmidrule(lr){3-4}\cmidrule(lr){5-6}',
         r'Method & M0 & orig. & re-cert. & orig. & re-cert. & M3 \\']
    for k, titulo in CONJ:
        L += [r'\midrule', r'\multicolumn{7}{@{}l}{\textit{' + titulo + r'}} \\']
        for alg in ('C-DETC', 'C-BUF90C', 'C-BUF95C', 'CAL', 'C-SIMC', 'ACAL'):
            r = R[k][alg]; c = C[k][alg]['M2']; t = T[k][alg]
            cel = [f"{r['M0']['cumpre']:.0f}", f"{r['M2']['cumpre']:.0f}", f"{c['cumpre']:.0f}",
                   f"{t['original']:.0f}", f"{t['recertificada']:.0f}", f"{r['M3']['cumpre']:.0f}"]
            L.append(NOME[alg] + ' & ' + ' & '.join(cel) + r' \\')
    L += [r'\bottomrule', r'\end{tabular}', r'\end{table}']
    open(a.tex, 'w', encoding='utf-8').write('\n'.join(L) + '\n')
    print('tabela:', a.tex)


if __name__ == '__main__':
    main()
