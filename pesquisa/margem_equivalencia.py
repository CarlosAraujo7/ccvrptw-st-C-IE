#!/usr/bin/env python3
"""
Justificativa da margem de equivalência de ±0,5% pelo custo (17/11/2026, diário Entrada 147).

Crítica da revisão cega de 17/11: a margem de ±0,5% da Tabela S22 aparecia no artigo sem justificativa. Este script usa SÓ resultados já
guardados (nenhuma rodada nova): o custo da solução mais barata que cumpre, por instância e semente, do bloco 1 (alfa = 0,90; 54 instâncias;
10 sementes; `pesquisa/resultados_c/por_execucao.csv`, o mesmo arquivo de `equivalencia.py`). Calcula, para o SIM+C, o BUF 0,95+C, o BUF 0,90+C,
o CAL e o MOEA/D+C:
  A. a variação do custo entre execuções do mesmo método e da mesma instância com sementes diferentes: coeficiente de variação (CV)
     sobre as 10 sementes e diferença relativa absoluta |d| entre duas sementes (todos os pares ordenados);
  B. a calibração A/A do TOST: cada método contra ele mesmo, com as 10 sementes divididas em dois grupos de 5 de todas as 126 maneiras
     (as médias por instância de cada grupo são comparadas como na Tabela S22). A "margem mínima" é a menor margem que o TOST aceita, isto é,
     o maior limite absoluto do IC de 90% de Hodges–Lehmann. Com 5 contra 5 sementes o ruído é MAIOR que o da comparação principal (10 contra 10);
  C. cada método contra o SIM+C com as 10 sementes, como na Tabela S22: estimativa HL, diferença média, margem mínima e veredito do TOST com margens
     de 0,25%, 0,5% e 1%. O BUF 0,70 entra de propósito: sua diferença típica (HL) é pequena, mas a média é dominada por poucas instâncias, o que
     marca o limite do estimador (descreve a instância típica, não a cauda).
Origem da margem (diário, Entrada 84): fixada em 08/10/2026 em dados arquivados da versão anterior do experimento, porque com 1% até CAL e MOEA/D,
significativamente diferentes, seriam declarados equivalentes.
Ampliação de 17/11/2026 (diário, Entrada 149; achados 1 e 3 da terceira revisão cega): o painel C também aplica o critério à MÉDIA das 54 diferenças
relativas, com IC de 90% por bootstrap de instâncias (percentis 5 e 95 de B_BOOT reamostragens; semente fixa por método), porque o HL descreve a instância
típica e ignora a cauda: pelo HL o CAL é equivalente a ±0,5% (margem 0,39%), pela média não (uma execução em 540 custa 2,3 vezes a mediana do SIM+C na
instância). A legenda traz essa execução extrema e a âncora em reais da margem (0,5% do custo mediano do plano; custos da Tabela S12: 150 BRL por rota,
2,5 BRL/km; valores em reais do modelo, só diferenças relativas têm interpretação).
Saída: pesquisa/resultados_robustez/relatorio_margem_equivalencia.txt e margem_equivalencia.json; tabela LaTeX com --tex (padrão: tab_margem.tex na
mesma pasta). Uso: python pesquisa/margem_equivalencia.py [--csv caminho] [--out pasta] [--tex arquivo.tex]
"""
import argparse, csv, itertools, json, os
from collections import defaultdict
import numpy as np
from scipy.stats import norm

BASE = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
LIMITE = 0.5                      # margem do artigo (%)
MARGENS = (0.25, 0.5, 1.0)        # margens do painel C (%)
MARGENS_AA = (0.5, 0.25, 0.1)     # margens do painel B (%)
REF = 'C-SIMC'
B_BOOT = 20000                    # reamostragens do bootstrap de instâncias (média das diferenças)
SEMENTE = 20261117                # semente do bootstrap; um fluxo por método, [SEMENTE, índice do método]
C_FIX, C_KM = 150.0, 2.5          # custo fixo por rota (BRL) e custo por km (BRL/km): Tabela S12, geradores de instâncias
# (código no CSV, macro do LaTeX)
METODOS = [('C-SIMC', r'\mSIM{}'), ('C-BUF95C', r'\mBUF{0.95}'), ('C-BUF90C', r'\mBUF{0.90}'),
           ('C-CAL-BEST', r'\mCAL{}'), ('C-MOEADC', r'\mMOEAD{}')]
CONTRA_SIM = [('C-BUF95C', r'\mBUF{0.95}'), ('C-BUF90C', r'\mBUF{0.90}'), ('C-CAL-BEST', r'\mCAL{}'),
              ('C-MOEADC', r'\mMOEAD{}'), ('C-BUF70C', r'\mBUF{0.70}')]


def hodges_lehmann(d, nivel):
    """Estimador HL e IC bilateral (nível dado) pelos quantis das médias de Walsh (igual ao de equivalencia.py)."""
    d = np.asarray(d, float); n = len(d)
    w = np.sort((d[:, None] + d[None, :])[np.triu_indices(n)] / 2.0)
    M = len(w)
    z = norm.ppf(1 - (1 - nivel) / 2)
    C = max(int(np.floor(M / 2 - z * np.sqrt(n * (n + 1) * (2 * n + 1) / 24))), 1)
    return float(np.median(w)), float(w[C - 1]), float(w[M - C])


def carregar(caminho):
    """{algoritmo: {instância: {semente: custo da solução mais barata que cumpre}}}; execuções sem solução que cumpre ficam de fora"""
    t = defaultdict(lambda: defaultdict(dict))
    for r in csv.DictReader(open(caminho, encoding='utf-8')):
        v = r['custo_min_cumpre']
        if r['alpha'] == '0.9' and v not in ('', 'nan') and np.isfinite(float(v)):
            t[r['algorithm']][r['instance']][int(r['seed'])] = float(v)
    return t


def matriz(por_inst, insts, sementes):
    """matriz instâncias x sementes; exige que todas as execuções tenham custo (vale para os métodos dos painéis A e B)"""
    faltam = [(i, s) for i in insts for s in sementes if s not in por_inst.get(i, {})]
    assert not faltam, f'execuções sem custo: {faltam[:3]}'
    return np.array([[por_inst[i][s] for s in sementes] for i in insts])


def variacao(m):
    """A. CV entre sementes por instância e |d| entre duas execuções (pares ordenados de sementes diferentes)"""
    cv = 100 * np.std(m, axis=1, ddof=1) / np.mean(m, axis=1)
    d = np.concatenate([np.abs(100 * (m[:, a] / m[:, b] - 1)) for a, b in itertools.permutations(range(m.shape[1]), 2)])
    return dict(cv_med=float(np.median(cv)), cv_p90=float(np.percentile(cv, 90)), d_med=float(np.median(d)),
                d_p75=float(np.percentile(d, 75)), d_p90=float(np.percentile(d, 90)), d_gt=float(100 * np.mean(d > LIMITE)),
                n_pares=int(len(d)))


def aa(m):
    """B. o método contra ele mesmo: todas as divisões das sementes em dois grupos iguais (a semente 0 sempre no grupo A)"""
    ns = m.shape[1]
    hl_abs, marg = [], []
    for comb in itertools.combinations(range(ns), ns // 2):
        if 0 not in comb:
            continue
        outras = [s for s in range(ns) if s not in comb]
        d = 100 * (m[:, outras].mean(axis=1) / m[:, list(comb)].mean(axis=1) - 1)
        h, lo, hi = hodges_lehmann(d, 0.90)
        hl_abs.append(abs(h)); marg.append(max(abs(lo), abs(hi)))
    marg = np.array(marg)
    return dict(divisoes=len(marg), hl_max=float(max(hl_abs)), marg_med=float(np.median(marg)), marg_max=float(marg.max()),
                equiv={str(g): int(np.sum(marg < g)) for g in MARGENS_AA})


def extremo(x, base):
    """execução (instância, semente, custo, razão) do método com o maior custo relativo à mediana do SIM+C na instância"""
    med = {k: float(np.median(list(v.values()))) for k, v in base.items()}
    pior = None
    for k, sd in x.items():
        if k not in med:
            continue
        for s, c in sd.items():
            r = c / med[k]
            if pior is None or r > pior[3]:
                pior = (k, s, c, r)
    return pior


def custo_plano(base):
    """âncora em reais da margem: custo do plano (SIM+C, solução mais barata que cumpre, média das sementes por instância)"""
    m = np.array([np.mean(list(v.values())) for v in base.values()])
    q1, med, q3 = (float(v) for v in np.percentile(m, [25, 50, 75]))
    por_n = {str(n): float(np.median([np.mean(list(v.values())) for k, v in base.items() if f'_{n:03d}_' in k])) for n in (25, 50, 100)}
    return dict(mediana=med, q1=q1, q3=q3, minimo=float(m.min()), maximo=float(m.max()), por_n=por_n,
                meia_pct=0.005 * med, meia_pct_q1=0.005 * q1, meia_pct_q3=0.005 * q3,
                km=0.005 * med / C_KM, fracao_fixo=0.005 * med / C_FIX)


def contra_sim(x, base, ns, k_met=0):
    """C. método contra o SIM+C: média das sementes com custo por instância (como em equivalencia.py), diferença relativa pareada (%);
    HL (instância típica) e média com IC de 90% por bootstrap de instâncias (a margem pela média usa o mesmo critério: maior limite absoluto)"""
    ks = sorted(k for k in x if k in base)
    d = np.array([100 * (np.mean(list(x[k].values())) / np.mean(list(base[k].values())) - 1) for k in ks])
    hl, lo95, hi95 = hodges_lehmann(d, 0.95)
    _, lo90, hi90 = hodges_lehmann(d, 0.90)
    marg = max(abs(lo90), abs(hi90))
    idx = np.random.default_rng([SEMENTE, k_met]).integers(0, len(d), size=(B_BOOT, len(d)))
    mb = d[idx].mean(axis=1)
    lo_m, hi_m = float(np.percentile(mb, 5)), float(np.percentile(mb, 95))
    marg_m = max(abs(lo_m), abs(hi_m))
    sem_custo = sum(ns - len(x[k]) for k in ks)
    # concentração: instâncias com diferença > 2% e a parte da soma dos valores absolutos que as 5 maiores diferenças carregam
    grandes = [k for k, v in zip(ks, d) if v > 2]
    ab = np.sort(np.abs(d))[::-1]
    top5 = float(100 * ab[:5].sum() / ab.sum()) if ab.sum() > 0 else float('nan')
    estreitas_altas = sum(1 for k in grandes if k.split('_')[-1] == 'HIGH' and k.split('_')[1] in ('C1', 'R1', 'RC1'))
    return dict(n=len(d), hl=hl, ic95=(lo95, hi95), ic90=(lo90, hi90), media=float(d.mean()), margem_min=marg,
                ic90_media=(lo_m, hi_m), margem_media=marg_m,
                execucoes_sem_custo=int(sem_custo), n_gt2=len(grandes), n_gt2_alta_estreita=int(estreitas_altas), top5=top5,
                equiv={str(g): bool(marg < g) for g in MARGENS}, equiv_media={str(g): bool(marg_m < g) for g in MARGENS})


def sim_nao(v, sim='yes', nao='no'):
    return sim if v else nao


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--csv', default=os.path.join(BASE, 'pesquisa', 'resultados_c', 'por_execucao.csv'))
    ap.add_argument('--out', default=os.path.join(BASE, 'pesquisa', 'resultados_robustez'))
    ap.add_argument('--tex', default=None)
    a = ap.parse_args()
    dados = carregar(a.csv)
    insts = sorted(dados[REF])
    sementes = sorted(set.intersection(*[set(dados[REF][i]) for i in insts]))
    ni, ns = len(insts), len(sementes)
    mats = {alg: matriz(dados[alg], insts, sementes) for alg, _ in METODOS}      # painéis A e B: todas as execuções têm custo
    os.makedirs(a.out, exist_ok=True)
    rel = open(os.path.join(a.out, 'relatorio_margem_equivalencia.txt'), 'w', encoding='utf-8')

    def log(t=''):
        print(t, flush=True); rel.write(t + '\n')

    res = dict(instancias=ni, sementes=ns, limite=LIMITE, A={}, B={}, C={})
    log(f'MARGEM DE EQUIVALÊNCIA DO CUSTO (±{LIMITE}%): variação entre sementes. Custo da solução mais barata que cumpre, alfa 0,90, '
        f'{ni} instâncias, {ns} sementes.')
    log('\nA. Variação entre execuções do mesmo método e da mesma instância (sementes diferentes), %')
    for alg, nome in METODOS:
        v = variacao(mats[alg]); res['A'][alg] = v
        log(f"  {alg:12s} CV mediano {v['cv_med']:.2f}  CV p90 {v['cv_p90']:.2f} | |d| mediano {v['d_med']:.2f}  p75 {v['d_p75']:.2f}  "
            f"p90 {v['d_p90']:.2f}  fração > {LIMITE}%: {v['d_gt']:.0f}% de {v['n_pares']} pares")
    log(f'\nB. Calibração A/A: o método contra ele mesmo, {ns} sementes em dois grupos de {ns // 2} (todas as divisões), %')
    for alg, nome in METODOS:
        v = aa(mats[alg]); res['B'][alg] = v
        eq = '  '.join(f"±{g}%: {v['equiv'][str(g)]}/{v['divisoes']}" for g in MARGENS_AA)
        log(f"  {alg:12s} |HL| máx {v['hl_max']:.3f}  margem mínima mediana {v['marg_med']:.3f}  máxima {v['marg_max']:.3f}  equivalente em {eq}")
    log(f'\nC. Cada método contra o SIM+C, {ns} sementes (como na Tabela S22), %. HL = instância típica; média = média das {ni} diferenças com '
        f'IC de 90% por bootstrap de instâncias ({B_BOOT} reamostragens, semente {SEMENTE}); margem mínima = maior limite absoluto do IC de 90%')
    for k_met, (alg, nome) in enumerate(CONTRA_SIM):
        v = contra_sim(dados[alg], dados[REF], ns, k_met); res['C'][alg] = v
        eq = '  '.join(f"±{g}%: {sim_nao(v['equiv'][str(g)])}/{sim_nao(v['equiv_media'][str(g)])}" for g in MARGENS)
        log(f"  {alg:12s} n={v['n']}  HL {v['hl']:+.2f}  IC95 [{v['ic95'][0]:+.2f}, {v['ic95'][1]:+.2f}]  margem {v['margem_min']:.2f}  | "
            f"média {v['media']:+.2f}  IC90 [{v['ic90_media'][0]:+.2f}, {v['ic90_media'][1]:+.2f}]  margem {v['margem_media']:.2f}  | "
            f"equivalente (HL/média) {eq}"
            + (f"  ({v['execucoes_sem_custo']} execuções sem solução que cumpre ficam fora da média da instância)" if v['execucoes_sem_custo'] else ''))
        log(f"{'':14s}concentração: {v['n_gt2']} instâncias com diferença > 2% ({v['n_gt2_alta_estreita']} com incerteza alta e janelas estreitas); "
            f"as 5 maiores diferenças em valor absoluto carregam {v['top5']:.0f}% da soma dos valores absolutos")
        v['extremo'] = extremo(dados[alg], dados[REF])
        e = v['extremo']
        log(f"{'':14s}execução mais cara em relação à mediana do SIM+C na instância: {e[0].replace('CCVRPTW-ST_', '')} semente {e[1]}, custo {e[2]:,.0f} "
            f"({e[3]:.2f} vez a mediana)")
    res['custo_plano'] = custo_plano(dados[REF])
    cp = res['custo_plano']
    log(f"\nÂncora em reais (valores do modelo; custos da Tabela S12: {C_FIX:.0f} BRL por rota, {C_KM} BRL/km): custo do plano do SIM+C (solução mais barata que cumpre, "
        f"média das sementes por instância): mediana {cp['mediana']:,.0f}, quartis {cp['q1']:,.0f} e {cp['q3']:,.0f}, mínimo {cp['minimo']:,.0f}, máximo "
        f"{cp['maximo']:,.0f}; por tamanho (n = 25, 50, 100): {cp['por_n']['25']:,.0f}, {cp['por_n']['50']:,.0f}, {cp['por_n']['100']:,.0f}. "
        f"0,5% da mediana = {cp['meia_pct']:.1f} BRL (quartis: {cp['meia_pct_q1']:.1f} a {cp['meia_pct_q3']:.1f}); corresponde a {cp['km']:.1f} km "
        f"e a {100 * cp['fracao_fixo']:.0f}% do custo fixo de uma rota")
    rel.close()
    json.dump(res, open(os.path.join(a.out, 'margem_equivalencia.json'), 'w', encoding='utf-8'), indent=1)

    # tabela LaTeX (suplemento): painéis A e B (7 colunas) e painel C (9 colunas) em duas tabulars dentro da mesma tabela flutuante
    L = []
    L.append(r'\multicolumn{7}{@{}l}{\textit{A. Variation between runs of the same method on the same instance (different seeds)}} \\')
    L.append(r'Method & CV med. & CV p90 & $|d|$ med. & $|d|$ p75 & $|d|$ p90 & $|d|>0.5\%$ \\')
    L.append(r'\midrule')
    for alg, nome in METODOS:
        v = res['A'][alg]
        L.append(f"{nome} & {v['cv_med']:.2f} & {v['cv_p90']:.2f} & {v['d_med']:.2f} & {v['d_p75']:.2f} & {v['d_p90']:.2f} & {v['d_gt']:.0f} \\\\")
    L.append(r'\midrule')
    L.append(r'\multicolumn{7}{@{}l}{\textit{B. Each method against itself, seeds split 5 + 5 in the 126 possible ways (splits declared equivalent)}} \\')
    L.append(r'Method & $|\mathrm{HL}|$ max & Margin med. & Margin max & $\pm0.5\%$ & $\pm0.25\%$ & $\pm0.1\%$ \\')
    L.append(r'\midrule')
    for alg, nome in METODOS:
        v = res['B'][alg]
        L.append(f"{nome} & {v['hl_max']:.2f} & {v['marg_med']:.2f} & {v['marg_max']:.2f} & "
                 + ' & '.join(str(v['equiv'][str(g)]) for g in MARGENS_AA) + r' \\')
    L2 = []
    L2.append(r'\multicolumn{9}{@{}l}{\textit{C. Each method against \mSIM{}, ten seeds, as in Table~\ref{tab:equiv}}} \\')
    L2.append(r' & \multicolumn{2}{c}{Hodges--Lehmann} & \multicolumn{3}{c}{Mean (bootstrap)} & \multicolumn{3}{c}{Equivalent within (HL / mean)} \\')
    L2.append(r'\cmidrule(lr){2-3}\cmidrule(lr){4-6}\cmidrule(l){7-9}')
    L2.append(r'Method & Est. & Margin & Est. & 90\% interval & Margin & $\pm0.25\%$ & $\pm0.5\%$ & $\pm1\%$ \\')
    L2.append(r'\midrule')
    for alg, nome in CONTRA_SIM:
        v = res['C'][alg]
        ver = ' & '.join(f"{sim_nao(v['equiv'][str(g)])} / {sim_nao(v['equiv_media'][str(g)])}" for g in MARGENS)
        L2.append(f"{nome} & {v['hl']:+.2f} & {v['margem_min']:.2f} & {v['media']:+.2f} & {v['ic90_media'][0]:+.2f}, {v['ic90_media'][1]:+.2f} & "
                  f"{v['margem_media']:.2f} & {ver} \\\\")
    cp = res['custo_plano']
    milhar = lambda x: f'{x:,.0f}'.replace(',', '{,}')           # noqa: E731
    e = res['C']['C-CAL-BEST']['extremo']
    _, fam, n_ext, nivel_ext = e[0].split('_')
    legenda_c = (r'C: each method against \mSIM{} with all seeds. HL is the Hodges--Lehmann estimate of the typical instance and its margin the larger '
                 r'absolute limit of the 90\% Hodges--Lehmann interval; the mean of the 54 relative differences has a 90\% percentile-bootstrap interval '
                 f'over instances ({milhar(B_BOOT)} resamples) and its margin is the larger absolute limit of that interval. The last three columns give the '
                 r'verdict for each margin as Hodges--Lehmann / mean. The 0.70 buffer shows the limit of the Hodges--Lehmann estimator: its typical '
                 r'difference is small, but a few instances dominate its mean. \mCAL{} is equivalent within $\pm0.5\%$ by the typical instance but not by '
                 f'the mean, because one run in {ni * ns} of \\mCAL{{}} ({fam}, $n={int(n_ext)}$, {nivel_ext.lower()} uncertainty, seed {e[1]}) costs '
                 f'{e[3]:.1f} times the median of \\mSIM{{}} on that instance. The median plan of \\mSIM{{}} costs BRL {milhar(cp["mediana"])} '
                 f'(quartiles {milhar(cp["q1"])} and {milhar(cp["q3"])}), so 0.5\\% is BRL {cp["meia_pct"]:.0f} ({cp["meia_pct_q1"]:.0f} to '
                 f'{cp["meia_pct_q3"]:.0f}), about {cp["km"]:.0f}\\,km of driving at {C_KM}\\,BRL/km or {100 * cp["fracao_fixo"]:.0f}\\% of the fixed cost '
                 r'of a route; monetary values are in the units of the model, and the margin is a statistical resolution, not a threshold of '
                 r'operational relevance.')
    tex = (r"""\begin{table}[htbp]
\centering
\caption{Calibration of the $\pm0.5\%$ equivalence margin against the variation of the cost between seeds (cost of the cheapest compliant solution, $\alpha=0.90$, 54 instances, 10 seeds; all values in \% except counts). A: coefficient of variation (CV) of the cost over the ten seeds of each instance (median and 90th percentile over instances) and absolute relative difference $|d|$ between two runs with different seeds (median, 75th and 90th percentiles); the last column is the share of pairs with $|d|>0.5\%$. B: each method against itself, with the seeds divided in two groups of five in all 126 ways and the per-instance means compared as in Table~\ref{tab:equiv}; the margin is the smallest for which the two one-sided tests declare equivalence (the larger absolute limit of the 90\% Hodges--Lehmann interval), and the last three columns count the splits declared equivalent within each margin. Five against five seeds is noisier than the ten against ten of the main analysis. """ + legenda_c + r"""}
\label{tab:margem}
\small
\begin{tabular}{@{}l r r r r r r@{}}
\toprule
""" + '\n'.join(L) + r"""
\bottomrule
\end{tabular}

\medskip
\begin{tabular}{@{}l rr rrr ccc@{}}
\toprule
""" + '\n'.join(L2) + r"""
\bottomrule
\end{tabular}
\end{table}
""")
    destino = a.tex or os.path.join(a.out, 'tab_margem.tex')
    open(destino, 'w', encoding='utf-8').write(tex)
    print('tabela:', destino)


if __name__ == '__main__':
    main()
