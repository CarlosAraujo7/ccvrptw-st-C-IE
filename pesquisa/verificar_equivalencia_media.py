#!/usr/bin/env python3
"""
O veredito de equivalência de custo resiste a estimadores diferentes do Hodges-Lehmann? (17/11/2026, diário Entrada 148;
terceira revisão cega.)

Só resultados já guardados (nenhuma rodada nova): `pesquisa/resultados_c/por_execucao.csv` (alfa 0,9; 54 instâncias; 10 sementes), o mesmo arquivo
de `equivalencia.py` e `margem_equivalencia.py`. Para cada método contra o SIM+C, d_i é a diferença relativa pareada (%) entre as médias das
sementes de cada instância (como nos dois scripts citados). O artigo declara equivalência quando o IC de 90% de Hodges-Lehmann (HL) de d fica
dentro de ±0,5%; o HL descreve a instância típica e ignora a cauda. Aqui o mesmo critério é aplicado a outros estimadores de posição:
  - HL (o do artigo; repetido para conferir a Tabela S22/S32);
  - média, com IC de 90% por bootstrap de instâncias (percentis 5 e 95 de 20.000 reamostragens);
  - mediana e média aparada de 10%, também por bootstrap;
  - média, com IC de 90% de t de Student pareado.
"Margem mínima" é o maior limite absoluto do IC de 90%: o método é equivalente a ±m% quando a margem mínima é menor que m.
Também: p bruto de Wilcoxon (sinais e postos, bilateral, zeros descartados) e de uma permutação de sinais da média (100.000 sorteios), a lista das
execuções cujo custo passa de 1,5 vez a mediana do SIM+C na mesma instância e o efeito de retirar a instância com a maior diferença.
Saída: pesquisa/resultados_robustez/relatorio_equivalencia_media.txt. Uso: python pesquisa/verificar_equivalencia_media.py [--csv caminho] [--out pasta]
"""
import argparse, csv, os
from collections import defaultdict
import numpy as np
from scipy import stats
from scipy.stats import norm

BASE = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
REF = 'C-SIMC'
METODOS = [('C-BUF95C', 'BUF0.95+C'), ('C-BUF90C', 'BUF0.90+C'), ('C-CAL-BEST', 'CAL'), ('C-MOEADC', 'MOEA/D+C'), ('C-BUF70C', 'BUF0.70+C')]
LIMITE = 0.5
B_BOOT = 20000
N_PERM = 100000
SEMENTE = 20261117


def hodges_lehmann(d, nivel):
    """HL e IC bilateral (nível dado) pelos quantis das médias de Walsh (igual ao de equivalencia.py e margem_equivalencia.py)."""
    d = np.asarray(d, float); n = len(d)
    w = np.sort((d[:, None] + d[None, :])[np.triu_indices(n)] / 2.0)
    M = len(w)
    z = norm.ppf(1 - (1 - nivel) / 2)
    C = max(int(np.floor(M / 2 - z * np.sqrt(n * (n + 1) * (2 * n + 1) / 24))), 1)
    return float(np.median(w)), float(w[C - 1]), float(w[M - C])


def carregar(caminho):
    """{algoritmo: {instância: {semente: custo da solução mais barata que cumpre}}}"""
    t = defaultdict(lambda: defaultdict(dict))
    for r in csv.DictReader(open(caminho, encoding='utf-8')):
        v = r['custo_min_cumpre']
        if r['alpha'] == '0.9' and v not in ('', 'nan') and np.isfinite(float(v)):
            t[r['algorithm']][r['instance']][int(r['seed'])] = float(v)
    return t


def diferencas(x, base):
    ks = sorted(k for k in x if k in base)
    return ks, np.array([100 * (np.mean(list(x[k].values())) / np.mean(list(base[k].values())) - 1) for k in ks])


def ic90_media(d, k):
    """IC de 90% da média por bootstrap de instâncias, com o mesmo esquema de margem_equivalencia.py (um fluxo [SEMENTE, k] por método)"""
    idx = np.random.default_rng([SEMENTE, k]).integers(0, len(d), size=(B_BOOT, len(d)))
    mb = d[idx].mean(axis=1)
    return float(np.percentile(mb, 5)), float(np.percentile(mb, 95))


def ic90_boot(d, est, rng):
    idx = rng.integers(0, len(d), size=(B_BOOT, len(d)))
    v = np.array([est(d[i]) for i in idx])
    return float(np.percentile(v, 5)), float(np.percentile(v, 95))


def margem(lo, hi):
    return max(abs(lo), abs(hi))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--csv', default=os.path.join(BASE, 'pesquisa', 'resultados_c', 'por_execucao.csv'))
    ap.add_argument('--out', default=os.path.join(BASE, 'pesquisa', 'resultados_robustez'))
    a = ap.parse_args()
    dados = carregar(a.csv)
    base = dados[REF]
    med_sim = {k: float(np.median(list(v.values()))) for k, v in base.items()}
    rng = np.random.default_rng(SEMENTE)
    os.makedirs(a.out, exist_ok=True)
    rel = open(os.path.join(a.out, 'relatorio_equivalencia_media.txt'), 'w', encoding='utf-8')

    def log(t=''):
        print(t, flush=True); rel.write(t + '\n')

    log(f'EQUIVALÊNCIA DE CUSTO (±{LIMITE}%) POR OUTROS ESTIMADORES. Custo da solução mais barata que cumpre, alfa 0,90, contra o SIM+C; '
        f'{B_BOOT} reamostragens bootstrap de instâncias, semente {SEMENTE}.')
    log('"marg" = maior limite absoluto do IC de 90% (equivalente a ±m% se marg < m). Diferenças em %.\n')
    log(f'{"método":11s} {"n":>3s} | {"HL":>6s} {"marg":>5s} | {"média":>6s} {"IC90":>16s} {"marg":>5s} | {"mediana":>7s} {"marg":>5s} | '
        f'{"aparada":>7s} {"marg":>5s} | {"t:marg":>6s} | {"p Wilcoxon":>10s} {"p perm.média":>12s} | {"d>+0,5":>6s} {"d<-0,5":>6s} {"máx d":>6s}')
    res = {}
    for k_met, (alg, nome) in enumerate(METODOS):
        ks, d = diferencas(dados[alg], base)
        hl, _, _ = hodges_lehmann(d, 0.95)
        _, lo90, hi90 = hodges_lehmann(d, 0.90)
        lo_m, hi_m = ic90_media(d, k_met)
        lo_md, hi_md = ic90_boot(d, np.median, rng)
        lo_t, hi_t = ic90_boot(d, lambda v: stats.trim_mean(v, 0.10), rng)
        se = d.std(ddof=1) / np.sqrt(len(d)); tq = stats.t.ppf(0.95, len(d) - 1)
        marg_t = margem(d.mean() - tq * se, d.mean() + tq * se)
        try:
            p_w = float(stats.wilcoxon(d).pvalue)
        except ValueError:
            p_w = float('nan')
        sinais = rng.choice([-1.0, 1.0], size=(N_PERM, len(d)))
        p_perm = float((1 + np.sum(np.abs((sinais * d).mean(axis=1)) >= abs(d.mean()) - 1e-12)) / (N_PERM + 1))
        res[alg] = dict(d=d, ks=ks)
        log(f'{nome:11s} {len(d):3d} | {hl:+6.2f} {margem(lo90, hi90):5.2f} | {d.mean():+6.2f} [{lo_m:+6.2f},{hi_m:+6.2f}] {margem(lo_m, hi_m):5.2f} | '
            f'{np.median(d):+7.2f} {margem(lo_md, hi_md):5.2f} | {stats.trim_mean(d, 0.10):+7.2f} {margem(lo_t, hi_t):5.2f} | {marg_t:6.2f} | '
            f'{p_w:10.2g} {p_perm:12.3g} | {int((d > 0.5).sum()):6d} {int((d < -0.5).sum()):6d} {d.max():+6.1f}')

    log('\nSem a instância de maior diferença (um a um):')
    for k_met, (alg, nome) in enumerate(METODOS):
        d = res[alg]['d']; ks = res[alg]['ks']
        j = int(np.argmax(d)); d2 = np.delete(d, j)
        hl, _, _ = hodges_lehmann(d2, 0.95); _, lo90, hi90 = hodges_lehmann(d2, 0.90)
        lo_m, hi_m = ic90_media(d2, 100 + k_met)
        log(f'  {nome:11s} sem {ks[j]} (d = {d[j]:+.1f}%): média {d2.mean():+.2f}, IC90 da média [{lo_m:+.2f}, {hi_m:+.2f}] (marg {margem(lo_m, hi_m):.2f}); '
            f'HL {hl:+.2f} (marg {margem(lo90, hi90):.2f})')

    log('\nExecuções (instância, semente) com custo maior que 1,5 vez a mediana do SIM+C na instância:')
    for alg, nome in [('C-SIMC', 'SIM+C')] + METODOS:
        x = dados[alg]
        lista = [(k, s, c, c / med_sim[k]) for k, sd in x.items() if k in med_sim for s, c in sd.items() if c / med_sim[k] > 1.5]
        n_runs = sum(len(sd) for k, sd in x.items() if k in med_sim)
        top = '; '.join(f'{k.replace("CCVRPTW-ST_", "")} s{s}: {c:,.0f} ({r:.2f}x)' for k, s, c, r in sorted(lista, key=lambda t: -t[3])[:3])
        log(f'  {nome:11s} {len(lista):3d} de {n_runs} execuções com custo' + (f'; as maiores: {top}' if lista else ''))
    rel.close()


if __name__ == '__main__':
    main()
