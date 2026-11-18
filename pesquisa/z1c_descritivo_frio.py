#!/usr/bin/env python3
"""
Descritivo do Z1c (17/11/2026, diário Entrada 141): quanto da carga NÃO congelada viaja na zona única a -18 °C.

No Z1c (compartimento único, rotas mistas, set-point da categoria mais fria a bordo) toda rota com algum cliente de
congelados roda a -18 °C. O modelo não mede o excesso de frio (a violação do limite inferior de produtos refrigerados);
este script só conta, nas soluções guardadas, a fração das rotas com congelados e a fração da massa não congelada
transportada nessas rotas. Descritivo, pós-hoc (não faz parte do desenho fixado na Entrada 138).
Uso: python pesquisa/z1c_descritivo_frio.py
Saída: pesquisa/resultados_c_b3/descritivo_frio_z1c.txt
Desde a Entrada 153 também roda sobre o Z1c refeito com o limite inferior de temperatura:
  python pesquisa/z1c_descritivo_frio.py --pasta results_c_b3_lb/C-BUF90C-Z1c --saida pesquisa/resultados_c_b3_lb/descritivo_frio_z1c_limite.txt
"""
import argparse, glob, json, os, sys
import numpy as np

BASE = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
sys.path.insert(0, os.path.join(BASE, 'src'))
import nsga2_simheuristic as ns   # noqa: E402

PASTA = os.path.join(BASE, 'results_c_b3', 'C-BUF90C-Z1c')
SAIDA = os.path.join(BASE, 'pesquisa', 'resultados_c_b3', 'descritivo_frio_z1c.txt')


def main(pasta=PASTA, saida=SAIDA, rotulo='Z1c'):
    arqs = sorted(glob.glob(os.path.join(pasta, '*_solutions.json')))
    cache, por_exec = {}, []
    for a in arqs:
        nm = os.path.basename(a).split('__')[0]
        if nm not in cache:
            cache[nm] = ns.load_instance(os.path.join(BASE, 'instances', 'generated', nm + '.json'))['customers']
        cust = cache[nm]
        fr_r, fr_m = [], []
        for s in json.load(open(a, encoding='utf-8')):
            if not s['feasible']:
                continue
            com, m_nao, m_nao_com = 0, 0.0, 0.0
            for r in s['routes']:
                cats = [cust[c - 1]['product_category'] for c in r]
                tem = 'frozen' in cats
                com += tem
                for c, k in zip(r, cats):
                    if k != 'frozen':
                        q = cust[c - 1]['demand_kg']
                        m_nao += q
                        m_nao_com += q if tem else 0.0
            fr_r.append(com / len(s['routes']))
            fr_m.append(m_nao_com / m_nao if m_nao else 0.0)
        if fr_r:
            por_exec.append((np.mean(fr_r), np.mean(fr_m)))
    v = np.array(por_exec)
    linhas = [f'{rotulo} (C-BUF90C, n = 50): {len(arqs)} arquivos, {len(v)} execuções com solução viável',
              'média sobre as execuções da média sobre as soluções viáveis de cada frente final',
              'rotas com algum cliente de congelados (zona a -18 °C): média %.1f%% (mín %.1f%%, máx %.1f%%)' % (100 * v[:, 0].mean(), 100 * v[:, 0].min(), 100 * v[:, 0].max()),
              'massa não congelada transportada nessas rotas: média %.1f%% (mín %.1f%%, máx %.1f%%)' % (100 * v[:, 1].mean(), 100 * v[:, 1].min(), 100 * v[:, 1].max())]
    txt = '\n'.join(linhas)
    print(txt)
    os.makedirs(os.path.dirname(saida), exist_ok=True)
    open(saida, 'w', encoding='utf-8').write(txt + '\n')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--pasta', default=None, help='pasta das frentes (relativa à raiz do projeto)')
    ap.add_argument('--saida', default=None, help='arquivo de saída (relativo à raiz do projeto)')
    ap.add_argument('--rotulo', default=None)
    a = ap.parse_args()
    main(os.path.join(BASE, a.pasta) if a.pasta else PASTA, os.path.join(BASE, a.saida) if a.saida else SAIDA,
         a.rotulo or ('Z1c com limite inferior' if a.pasta else 'Z1c'))
