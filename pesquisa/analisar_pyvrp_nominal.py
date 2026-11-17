#!/usr/bin/env python3
"""
Análise do B1 (diário Entradas 128-133): PyVRP SEM buffer, isto é, o plano nominal do resolvedor externo.

Fonte do plano nominal: iteração 0 do PyVRP-SIM (resultados_pyvrp_sim/pyvrp_sim.jsonl): custo de distância resolvido com os
tempos NOMINAIS, mesmo tempo por resolução, rotas saindo no instante 0, certificado em 500 cenários (a mesma amostra do pipeline).
Cobre as 332 células (instância x semente) do B1: benchmark (36 instâncias MEDIUM/HIGH) e São Paulo (16), 5 sementes; n = 200 e 400
(24 instâncias), 3 sementes. As resoluções novas de pipeline_nominal.jsonl são a amostra de replicação (Entrada 133).

Blocos (plano fixado nas Entradas 128 e 133):
  A  certificação do plano nominal, por conjunto e nível; separação pelo tipo de família (exploratório, post hoc)
  B  preço da confiabilidade: custo do plano certificado do PyVRP contra o plano nominal da mesma célula (HL por instância;
     descritivo separado pelo resultado da certificação do plano nominal)
  C  post hoc: CAL com um degrau nominal
  D  decomposição aditiva (log da razão de custos) da diferença SIM+C contra PyVRP CAL+R: simulação na busca contra buffer
     (mesmo motor), motor (mesmo buffer 0,95) e escada com refinamento do PyVRP
  E  replicação e carga
Uso: python pesquisa/analisar_pyvrp_nominal.py
Saídas em pesquisa/resultados_pyvrp_nominal/: relatorio_pyvrp_nominal.txt, resumo_pyvrp_nominal.json,
  tab_pyvrp_nominal_cert.tex, tab_pyvrp_nominal_preco.tex, tab_pyvrp_nominal_decomp.tex
"""
import csv, json, math, os, sys
from collections import defaultdict
import numpy as np
from scipy import stats
from scipy.stats import wilcoxon

BASE = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
sys.path.insert(0, os.path.join(BASE, 'pesquisa'))
import calibracao_preditiva_base as b      # noqa: E402
from equivalencia import hodges_lehmann    # noqa: E402
from analisar_refinamento import holm      # noqa: E402

ALPHA = b.ALPHA
OUT = os.path.join(BASE, 'pesquisa', 'resultados_pyvrp_nominal')
SIM = os.path.join(BASE, 'pesquisa', 'resultados_pyvrp_sim', 'pyvrp_sim.jsonl')
NOM = os.path.join(OUT, 'pipeline_nominal.jsonl')
PIPE_BENCH = os.path.join(BASE, 'pesquisa', 'resultados_pyvrp_pipeline', 'pipeline.jsonl')
PIPE_ESC = os.path.join(BASE, 'pesquisa', 'resultados_escala', 'pipeline_escala.jsonl')
# rótulo (LaTeX), filtro de instância, pasta do motor, pasta do refinamento, arquivo do PyVRP com buffer
CONJ = [
    ('benchmark', lambda i: '_SP' not in i and '_200_' not in i and '_400_' not in i and i.endswith(('_MEDIUM', '_HIGH')),
     'results_c', 'resultados_refinamento', PIPE_BENCH),
    ('S\\~ao Paulo', lambda i: '_SP' in i, 'results_c_sp', os.path.join('resultados_caso_sp', 'refinamento'), PIPE_BENCH),
    ('$n=200$', lambda i: '_200_' in i, 'results_c_escala', os.path.join('resultados_escala', 'refinamento'), PIPE_ESC),
    ('$n=400$', lambda i: '_400_' in i, 'results_c_escala', os.path.join('resultados_escala', 'refinamento'), PIPE_ESC),
]
_cache = {}


def ler(p):
    return [json.loads(x) for x in open(p, encoding='utf-8')] if os.path.exists(p) else []


def nivel(i):
    return 'MEDIUM' if i.endswith('_MEDIUM') else 'HIGH'


def tipo(i):
    """famílias 1 (C1, R1, RC1): horizonte curto e janelas estreitas; famílias 2: horizonte longo e janelas largas"""
    if '_SP' in i:
        return None
    return 'estreita' if i.replace('CCVRPTW-ST_', '').split('_')[0].endswith('1') else 'larga'


def tup(s):
    return None if s is None else (s['f1'], s['p_on'], s['q_oos'], s['f3'])


def planos_motor(res, refin):
    if (res, refin) in _cache:
        return _cache[(res, refin)]
    R = b.carregar(res)
    P = defaultdict(dict)
    for k, d in R.items():
        P['SIM+C'][k] = tup(d.get('C-SIMC'))
        P['BUF0.95+C'][k] = tup(d.get('C-BUF95C'))
        P['CAL'][k] = tup(b.melhor([d.get(a) for a in b.DEG]))
    p = os.path.join(BASE, 'pesquisa', refin, 'refinamento_CAL.json')
    if os.path.exists(p):
        for x in json.load(open(p, encoding='utf-8')):
            P['CAL+R'][(x['instance'], x['seed'])] = (x['fin_f1'], x['fin_p_oos'], x['fin_q_oos'], x['fin_f3'])
    _cache[(res, refin)] = P
    return P


def buffer_pyvrp(arq, filtro):
    """{(inst, semente): {q: linha}} só com resoluções viáveis; planos PY-CAL(+R) e PY-BUF0.95(+R) (None sem plano certificado)"""
    por = defaultdict(dict)
    for x in ler(arq):
        if filtro(x['instance']) and x['viavel']:
            por[(x['instance'], x['seed'])][x['q']] = x
    P = defaultdict(dict); cal = {}
    for k, d in por.items():
        cert = [x for x in d.values() if x['base_p_cert'] >= ALPHA]
        e = min(cert, key=lambda x: x['base_f1']) if cert else None
        cal[k] = e
        P['PY-CAL'][k] = None if e is None else (e['base_f1'], e['base_p_oos'], e['base_q_oos'], e['base_f3'])
        P['PY-CAL+R'][k] = None if e is None else (e['fin_f1'], e['fin_p_oos'], e['fin_q_oos'], e['fin_f3'])
        x = d.get(0.95)
        ok = x is not None and x['base_p_cert'] >= ALPHA
        P['PY-BUF0.95'][k] = (x['base_f1'], x['base_p_oos'], x['base_q_oos'], x['base_f3']) if ok else None
        P['PY-BUF0.95+R'][k] = (x['fin_f1'], x['fin_p_oos'], x['fin_q_oos'], x['fin_f3']) if ok else None
    return por, cal, P


def nominal_k0(filtro):
    """iteração 0 do PyVRP-SIM: {(inst, semente): dict(f1, f3, p_cert, viavel)}"""
    d = {}
    for x in ler(SIM):
        if not filtro(x['instance']):
            continue
        it = [i for i in x['its'] if i['k'] == 0][0]
        d[(x['instance'], x['seed'])] = it
    return d


def por_inst(valores, filtro_chave=None):
    """{instância: média sobre as células}"""
    pi = defaultdict(list)
    for (i, s), v in valores.items():
        if v is not None and (filtro_chave is None or filtro_chave((i, s))):
            pi[i].append(v)
    return {i: float(np.mean(v)) for i, v in pi.items()}


def hl_inst(x, y):
    """HL por instância da diferença percentual 100*(x/y - 1); x e y são {instância: custo médio}"""
    I = sorted(set(x) & set(y))
    if len(I) < 3:
        return None
    dd = np.array([100 * (x[i] / y[i] - 1) for i in I])
    hl, lo, hi = hodges_lehmann(dd, 0.95)
    p = wilcoxon(dd).pvalue if np.any(np.abs(dd) > 1e-12) else 1.0
    return dict(n=len(I), hl=hl, lo=lo, hi=hi, p=p, menores=int(np.sum(dd < -1e-12)), maiores=int(np.sum(dd > 1e-12)))


def resumo_estat(v):
    v = np.asarray(v, float)
    return (len(v), float(v.mean()), float(np.median(v))) if len(v) else (0, float('nan'), float('nan'))


def main():
    out = ['B1 — PyVRP sem buffer (plano nominal = iteração 0 do PyVRP-SIM); α = %.2f; certificação em 500 cenários' % ALPHA]
    R = {'cert': {}, 'preco': {}, 'hl': {}, 'cal0': {}, 'decomp': {}, 'replicacao': {}}
    for rot, filtro, res, refin, pipe in CONJ:
        k0 = nominal_k0(filtro)
        if not k0:
            continue
        mot = planos_motor(res, refin)
        por, cal, Pp = buffer_pyvrp(pipe, filtro)
        out.append(f'\n==== {rot}: {len(k0)} células (instância x semente), sementes {sorted({k[1] for k in k0})} ====')
        # ---- A. certificação do plano nominal
        out.append('A. Certificação do plano nominal (P_cert >= α em 500 cenários)')
        for lv in ('MEDIUM', 'HIGH'):
            ks = [k for k in k0 if nivel(k[0]) == lv]
            if not ks:
                continue
            v = [k for k in ks if k0[k]['viavel']]
            pc = np.array([k0[k]['p_cert'] for k in v]); cert = pc >= ALPHA
            linha = (f'  {lv:6s} células {len(ks):3d} (viáveis {len(v):3d})  certificado {100 * cert.mean():5.1f}% ({int(cert.sum())})  '
                     f'P_cert média {pc.mean():.3f}  mediana {np.median(pc):.3f}  mín {pc.min():.3f}')
            d = dict(n=len(ks), viaveis=len(v), cert=100 * float(cert.mean()), n_cert=int(cert.sum()), pbar=float(pc.mean()), pmed=float(np.median(pc)))
            if tipo(ks[0][0]) is not None:
                for t_ in ('estreita', 'larga'):
                    kt = [k for k in v if tipo(k[0]) == t_]
                    ct = np.array([k0[k]['p_cert'] >= ALPHA for k in kt])
                    linha += f' | {t_}: {100 * ct.mean():5.1f}% de {len(kt)}'
                    d[t_] = dict(n=len(kt), cert=100 * float(ct.mean()), n_cert=int(ct.sum()))
                for fam in ('C1', 'R1', 'RC1', 'C2', 'R2', 'RC2'):
                    kt = [k for k in v if k[0].replace('CCVRPTW-ST_', '').split('_')[0] == fam]
                    if kt:
                        d[f'fam_{fam}'] = dict(n=len(kt), cert=100 * float(np.mean([k0[k]['p_cert'] >= ALPHA for k in kt])))
            out.append(linha)
            R['cert'][f'{rot}|{lv}'] = d
        # ---- B. preço da confiabilidade
        out.append('B. Preço da confiabilidade: custo do plano certificado do PyVRP contra o plano nominal da mesma célula (%)')
        nom_f1 = {k: (v['f1'] if v['viavel'] else None) for k, v in k0.items()}
        cel = defaultdict(lambda: defaultdict(list))
        for k, it in k0.items():
            if not it['viavel'] or k not in cal:
                continue
            e = cal[k]
            cn = it['p_cert'] >= ALPHA
            lv = nivel(k[0])
            if e is not None:
                cel[lv]['cal_todas'].append(100 * (e['base_f1'] / it['f1'] - 1))
                cel[lv]['calR_todas'].append(100 * (e['fin_f1'] / it['f1'] - 1))
                cel[lv]['cal_passa' if cn else 'cal_falha'].append(100 * (e['base_f1'] / it['f1'] - 1))
            x = por[k].get(0.95)
            if x is not None and x['base_p_cert'] >= ALPHA:
                cel[lv]['b95_todas'].append(100 * (x['base_f1'] / it['f1'] - 1))
                cel[lv]['b95_passa' if cn else 'b95_falha'].append(100 * (x['base_f1'] / it['f1'] - 1))
            # C. post hoc: CAL com um degrau nominal (a) e política "nominal primeiro" (b): para no plano nominal se ele passa
            if e is not None:
                c0 = min(it['f1'] if cn else math.inf, e['base_f1'])
                cel[lv]['cal0_vs_cal'].append(100 * (c0 / e['base_f1'] - 1))
                cel[lv]['cal0_vs_nom'].append(100 * (c0 / it['f1'] - 1))
                cnf = it['f1'] if cn else e['base_f1']
                cel[lv]['nf_vs_cal'].append(100 * (cnf / e['base_f1'] - 1))
                cel[lv]['nf_solves'].append(1 if cn else 4)
        for lv in ('MEDIUM', 'HIGH'):
            c = cel.get(lv)
            if not c:
                continue
            lin = [f'  {lv:6s}']
            d = {}
            for nome, rotulo in (('cal_todas', 'CAL, todas'), ('cal_passa', 'CAL, nominal passa'), ('cal_falha', 'CAL, nominal falha'),
                                 ('calR_todas', 'CAL+R, todas'), ('b95_todas', 'BUF0.95, todas'), ('b95_passa', 'BUF0.95, passa'), ('b95_falha', 'BUF0.95, falha')):
                n, m, md = resumo_estat(c.get(nome, []))
                d[nome] = dict(n=n, media=m, mediana=md)
                lin.append(f'{rotulo}: n={n} média {m:+.2f} mediana {md:+.2f}')
            out.append(' | '.join(lin))
            R['preco'][f'{rot}|{lv}'] = d
            n, m, md = resumo_estat(c.get('cal0_vs_cal', [])); n2, m2, md2 = resumo_estat(c.get('cal0_vs_nom', []))
            n3, m3, md3 = resumo_estat(c.get('nf_vs_cal', [])); sol = float(np.mean(c.get('nf_solves', [np.nan])))
            out.append(f'      C. (post hoc) CAL com degrau nominal: contra CAL média {m:+.3f} mediana {md:+.3f}; '
                       f'contra o plano nominal média {m2:+.3f} mediana {md2:+.3f} (n={n})')
            out.append(f'         política "nominal primeiro" (para no plano nominal se ele passa; senão a escada completa): custo contra o CAL média {m3:+.3f} '
                       f'mediana {md3:+.3f}; resoluções por célula {sol:.2f} contra 3 do CAL')
            R['cal0'][f'{rot}|{lv}'] = dict(n=n, vs_cal_media=m, vs_cal_mediana=md, vs_nom_media=m2, vs_nom_mediana=md2,
                                            nf_vs_cal_media=m3, nf_vs_cal_mediana=md3, nf_solves=sol)
        out.append('   HL por instância do custo (%), plano certificado do PyVRP contra o nominal (Holm entre as 3 comparações):')
        for lv in ('MEDIUM', 'HIGH', 'todos'):
            filt = (lambda k, lv=lv: lv == 'todos' or nivel(k[0]) == lv)
            base = por_inst({k: v for k, v in nom_f1.items() if v is not None}, filt)
            comps = []
            for nome, P_ in (('PY-CAL', Pp['PY-CAL']), ('PY-CAL+R', Pp['PY-CAL+R']), ('PY-BUF0.95', Pp['PY-BUF0.95'])):
                x = por_inst({k: (t[0] if t is not None else None) for k, t in P_.items() if k in nom_f1}, filt)
                r_ = hl_inst(x, base)
                if r_:
                    comps.append((nome, r_))
            for (nome, r_), pa in zip(comps, holm([c_[1]['p'] for c_ in comps])):
                r_['p_holm'] = pa
                out.append(f'      {lv:6s} {nome:11s} n={r_["n"]:2d}  HL {r_["hl"]:+7.3f} [{r_["lo"]:+7.3f}, {r_["hi"]:+7.3f}]  '
                           f'p={r_["p"]:.3g}  p_Holm={pa:.3g}  ({r_["menores"]} menores, {r_["maiores"]} maiores)')
                R['hl'][f'{rot}|{lv}|{nome}'] = r_
        # ---- D. decomposição aditiva (por célula com os quatro planos)
        cels = [k for k in mot['SIM+C'] if mot['SIM+C'][k] is not None and mot['BUF0.95+C'].get(k) is not None
                and Pp['PY-BUF0.95'].get(k) is not None and Pp['PY-CAL+R'].get(k) is not None]
        for lv_nome, filt in (('todos', lambda k: True), ('MEDIUM', lambda k: nivel(k[0]) == 'MEDIUM'), ('HIGH', lambda k: nivel(k[0]) == 'HIGH')):
            ck = [k for k in cels if filt(k)]
            if len({k[0] for k in ck}) < 3:
                continue
            term = defaultdict(lambda: defaultdict(list))
            for k in ck:
                a, b95c, pb95, pcalr = mot['SIM+C'][k][0], mot['BUF0.95+C'][k][0], Pp['PY-BUF0.95'][k][0], Pp['PY-CAL+R'][k][0]
                term[k[0]]['total'].append(math.log(a / pcalr)); term[k[0]]['sim'].append(math.log(a / b95c))
                term[k[0]]['motor'].append(math.log(b95c / pb95)); term[k[0]]['escada'].append(math.log(pb95 / pcalr))
            d = {'n_inst': len(term), 'n_cel': len(ck)}
            lin = (f'  D. decomposição ({lv_nome}; {len(term)} instâncias, {len(ck)} células; '
                   f'100·média do log da razão de custos [IC 95% sobre as instâncias]):')
            for nome in ('total', 'sim', 'motor', 'escada'):
                v = np.array([np.mean(t[nome]) for t in term.values()]) * 100
                m = float(v.mean()); se = float(v.std(ddof=1) / math.sqrt(len(v))) if len(v) > 1 else float('nan')
                tcrit = stats.t.ppf(0.975, len(v) - 1) if len(v) > 1 else float('nan')
                d[nome] = dict(media=m, lo=m - tcrit * se, hi=m + tcrit * se)
                lin += f' {nome} {m:+.2f} [{m - tcrit * se:+.2f}, {m + tcrit * se:+.2f}];'
            out.append(lin)
            R['decomp'][f'{rot}|{lv_nome}'] = d
        # ---- D'. mesma decomposição em Hodges-Lehmann por instância (consistente com o resto do artigo; as partes não somam exatamente)
        sem = {k[1] for k in k0}

        def custo_inst(P_, sementes=None):
            return por_inst({k: (t_[0] if t_ is not None else None) for k, t_ in P_.items()
                             if filtro(k[0]) and (sementes is None or k[1] in sementes)})
        x_sim = custo_inst(mot['SIM+C']); x_b95c = custo_inst(mot['BUF0.95+C'])
        x_pb95 = custo_inst(Pp['PY-BUF0.95'], sem); x_pcalr = custo_inst(Pp['PY-CAL+R'], sem)
        pares_d = (('total', x_pcalr, x_sim, 'PyVRP CAL+R / SIM+C'), ('sim', x_b95c, x_sim, 'BUF0.95+C / SIM+C (mesmo motor)'),
                   ('motor', x_pb95, x_b95c, 'PyVRP BUF0.95 / BUF0.95+C (mesmo buffer)'), ('escada', x_pcalr, x_pb95, 'PyVRP CAL+R / PyVRP BUF0.95'))
        rs = [(nome, hl_inst(a_, b_), rot_) for nome, a_, b_, rot_ in pares_d]
        rs = [(n_, r_, ro_) for n_, r_, ro_ in rs if r_]
        for (n_, r_, ro_), pa in zip(rs, holm([r_['p'] for _, r_, _ in rs])):
            r_['p_holm'] = pa
            out.append(f"  D'. HL por instância, {ro_:42s} n={r_['n']:2d}  HL {r_['hl']:+7.3f} [{r_['lo']:+7.3f}, {r_['hi']:+7.3f}]  p={r_['p']:.3g}  p_Holm={pa:.3g}")
            R.setdefault('decomp_hl', {})[f'{rot}|{n_}'] = r_
    # ---- E. replicação e carga
    out.append('\n==== E. Replicação das resoluções novas (pipeline_nominal.jsonl) contra a iteração 0 do PyVRP-SIM, e carga ====')
    novo = [x for x in ler(NOM) if x['tag'] == 'nom']
    k0_all = {}
    for rot, filtro, *_ in CONJ:
        k0_all.update(nominal_k0(filtro))
    pares = [(x, k0_all[(x['instance'], x['seed'])]) for x in novo if (x['instance'], x['seed']) in k0_all]
    if pares:
        for nn in sorted({x['n'] for x, _ in pares}):
            xs = [(x, it) for x, it in pares if x['n'] == nn]
            dc = np.array([100 * (x['base_f1'] / it['f1'] - 1) for x, it in xs]); dp = np.array([x['base_p_cert'] - it['p_cert'] for x, it in xs])
            ig = int(np.sum(np.abs(dc) < 0.005))
            mesma = sum((x['base_p_cert'] >= ALPHA) == (it['p_cert'] >= ALPHA) for x, it in xs)
            out.append(f'  n = {nn:3d}: {len(xs):2d} pares; custo idêntico em {ig}; |Δcusto| máx. {np.max(np.abs(dc)):.3f}%; '
                       f'|ΔP_cert| máx. {np.max(np.abs(dp)):.3f}; mesma decisão de certificação em {mesma} de {len(xs)}')
            R['replicacao'][str(nn)] = dict(pares=len(xs), identicos=ig, dcusto_max=float(np.max(np.abs(dc))), mesma_decisao=int(mesma))
        itn = [x for x in novo if x.get('n_iter')]
        if itn:
            out.append(f'  esforço e carga nas {len(itn)} resoluções novas: iterações do PyVRP mediana {np.median([x["n_iter"] for x in itn]):.0f}; '
                       f'CPU do sistema mediana {np.median([x["sys_cpu_pct"] for x in itn]):.0f}% (mín {min(x["sys_cpu_pct"] for x in itn):.0f}%)')
        tol = [x for x in novo if 'viavel_tol' in x]
        out.append(f'  viabilidade: {sum(1 for x in tol if x["viavel_tol"])} de {len(tol)} resoluções com o campo de tolerância são viáveis; '
                   f'estritas: {sum(1 for x in tol if x["viavel"])} de {len(tol)}; sem o campo: {len(novo) - len(tol)}')
    ctrl = [x for x in ler(NOM) if x['tag'] == 'ctrl']
    antigo = {(x['instance'], x['q'], x['seed']): x for x in ler(PIPE_BENCH) + ler(PIPE_ESC)}
    cc = [(x, antigo[(x['instance'], x['q'], x['seed'])]) for x in ctrl if (x['instance'], x['q'], x['seed']) in antigo]
    if cc:
        dc = np.array([100 * (x['base_f1'] / o['base_f1'] - 1) for x, o in cc])
        out.append(f'  controle q = 0,90 (semente 0): {len(cc)} reexecuções; custo idêntico em {int(np.sum(np.abs(dc) < 1e-9))}; |Δcusto| máx. {np.max(np.abs(dc)):.3f}%')
        R['replicacao']['controle'] = dict(n=len(cc), identicos=int(np.sum(np.abs(dc) < 1e-9)), max=float(np.max(np.abs(dc))))
    txt = '\n'.join(out); print(txt)
    open(os.path.join(OUT, 'relatorio_pyvrp_nominal.txt'), 'w', encoding='utf-8').write(txt + '\n')
    json.dump(R, open(os.path.join(OUT, 'resumo_pyvrp_nominal.json'), 'w', encoding='utf-8'), indent=1)
    tabelas(R)


def fmt(v, casas=2, sinal=True):
    if v is None or (isinstance(v, float) and math.isnan(v)):
        return '--'
    if sinal and abs(v) < 0.5 * 10 ** (-casas):
        return '0.00'
    return f'{v:+.{casas}f}' if sinal else f'{v:.{casas}f}'


def tabelas(R):
    rots = [c[0] for c in CONJ]
    # S26: certificação do plano nominal
    L = [r'\begin{table}[htbp]', r'\centering',
         r'\caption{PyVRP without buffers: the plan solved with nominal travel times (no buffer, departures at the start of the horizon) is certified on 500 scenarios. '
         r'Cells: instance $\times$ seed (5 seeds; 3 with $n=200$ and 400). Certified: \% of the cells whose plan reaches $\hat P\ge\alpha=0.90$; $\bar P$: mean of $\hat P$. '
         r'Narrow / wide: the same percentage for the families with narrow time windows and short horizon (C1, R1, RC1) and with wide windows and long horizon (C2, R2, RC2); '
         r'this split was added after seeing the totals (exploratory). No cell is dropped.}',
         r'\label{tab:pyvrp-nominal}', r'\small\setlength{\tabcolsep}{4pt}', r'\begin{tabular}{@{}l l r r r r r@{}}', r'\toprule',
         r'Instances & Uncertainty & Cells & Certified (\%) & $\bar P$ & Narrow (\%) & Wide (\%) \\', r'\midrule']
    prim = True
    for rot in rots:
        linhas = [(lv, R['cert'].get(f'{rot}|{lv}')) for lv in ('MEDIUM', 'HIGH')]
        linhas = [(lv, d) for lv, d in linhas if d]
        if not linhas:
            continue
        if not prim:
            L.append(r'\midrule')
        prim = False
        for j, (lv, d) in enumerate(linhas):
            es = fmt(d['estreita']['cert'], 1, False) if 'estreita' in d else '--'
            la = fmt(d['larga']['cert'], 1, False) if 'larga' in d else '--'
            L.append(f'{rot if j == 0 else ""} & {lv.lower()} & {d["n"]} & {d["cert"]:.1f} & {d["pbar"]:.3f} & {es} & {la} \\\\')
    L += [r'\bottomrule', r'\end{tabular}', r'\end{table}']
    open(os.path.join(OUT, 'tab_pyvrp_nominal_cert.tex'), 'w', encoding='utf-8').write('\n'.join(L) + '\n')
    # S27: preço da confiabilidade
    L = [r'\begin{table}[htbp]', r'\centering',
         r'\caption{Price of certifying the plan of the external solver: cost of the certified plan with buffers, relative to the nominal plan of the same cell (\%). '
         r'\mCAL{}: cheapest certified rung of $q\in\{0.80,0.90,0.95\}$; \mBUF{0.95}: the fixed rung. Mean (median) over cells; Fails: only the cells whose nominal plan '
         r'misses the certification (number of cells in the third column). Negative values are search noise, because the nominal plan comes from one solve and \mCAL{} from '
         r'the best of three. The split by the result of the certification of the nominal plan is descriptive.}',
         r'\label{tab:pyvrp-price}', r'\small\setlength{\tabcolsep}{3pt}', r'\begin{tabular}{@{}l l r r r r r@{}}', r'\toprule',
         r' & & & \multicolumn{2}{c}{PyVRP, \mCAL{}} & \multicolumn{2}{c}{PyVRP, \mBUF{0.95}} \\', r'\cmidrule(lr){4-5}\cmidrule(lr){6-7}',
         r'Instances & Uncertainty & Fails & All cells & Fails & All cells & Fails \\', r'\midrule']
    prim = True
    for rot in rots:
        linhas = [(lv, R['preco'].get(f'{rot}|{lv}')) for lv in ('MEDIUM', 'HIGH')]
        linhas = [(lv, d) for lv, d in linhas if d]
        if not linhas:
            continue
        if not prim:
            L.append(r'\midrule')
        prim = False
        for j, (lv, d) in enumerate(linhas):
            def cel(nome, d=d):
                x = d.get(nome)
                if not x or x['n'] == 0:
                    return '--'
                return f'{fmt(x["media"])} ({fmt(x["mediana"])})'
            nf = d['cal_falha']['n'] if 'cal_falha' in d else 0
            L.append(f'{rot if j == 0 else ""} & {lv.lower()} & {nf} & {cel("cal_todas")} & {cel("cal_falha")} & {cel("b95_todas")} & {cel("b95_falha")} \\\\')
    L += [r'\bottomrule', r'\end{tabular}', r'\end{table}']
    open(os.path.join(OUT, 'tab_pyvrp_nominal_preco.tex'), 'w', encoding='utf-8').write('\n'.join(L) + '\n')
    # S28: decomposição (HL por instância; as partes não somam exatamente)
    L = [r'\begin{table}[htbp]', r'\centering',
         r'\caption{Where the cost gap between \mSIM{} and the external solver comes from. Hodges--Lehmann estimate of the paired difference per instance in the cost of the chosen plan (\%, mean over seeds; '
         r'negative: the first plan is cheaper); $^*$: $p_{\mathrm{Holm}}<0.05$ (Wilcoxon signed-rank test, Holm correction over the four comparisons). '
         r'Total: PyVRP, \mCAL{}+R against \mSIM{}. The other columns separate it into three steps: a 0.95 buffer in our solver against simulation inside its search; '
         r'PyVRP against our solver at the same 0.95 buffer (a cost-only solver against a three-objective search); the ladder with refinement against the fixed buffer in PyVRP. '
         r'The Hodges--Lehmann estimates of the steps do not add up exactly to the total.}',
         r'\label{tab:pyvrp-decomp}', r'\small\setlength{\tabcolsep}{4pt}', r'\begin{tabular}{@{}l r r r r r@{}}', r'\toprule',
         r'Instances & Inst. & Total & Buffer vs simulation & Solver & Ladder, +R \\', r'\midrule']
    for rot in rots:
        if f'{rot}|total' not in R.get('decomp_hl', {}):
            continue
        d = R['decomp_hl']

        def cel(nome, rot=rot, d=d):
            x = d.get(f'{rot}|{nome}')
            return '--' if x is None else fmt(x['hl']) + ('$^*$' if x['p_holm'] < 0.05 else '')
        L.append(f'{rot} & {d[f"{rot}|total"]["n"]} & {cel("total")} & {cel("sim")} & {cel("motor")} & {cel("escada")} \\\\')
    L += [r'\bottomrule', r'\end{tabular}', r'\end{table}']
    open(os.path.join(OUT, 'tab_pyvrp_nominal_decomp.tex'), 'w', encoding='utf-8').write('\n'.join(L) + '\n')


if __name__ == '__main__':
    main()
