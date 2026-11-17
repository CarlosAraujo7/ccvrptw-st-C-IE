#!/usr/bin/env python3
"""
B1 (16/11/2026, diário Entrada 128): pipeline do PyVRP SEM buffer (nível nominal) + certificação + refinamento.

Pergunta: do ganho de custo do PyVRP com buffers certificados sobre o SIM+C, quanto é o preço da confiabilidade
(buffer e certificação) e quanto é o motor? O nível nominal do MESMO pipeline (tempos médios, sem folga; rotas saem no
instante 0) dá o custo que o PyVRP teria se ignorasse a incerteza.

Usa pipeline_pyvrp_buffer.tarefa SEM nenhuma alteração, com q = 0,0 (buffer_q falso => tempos de planejamento nominais).
Mesma resolução (PyVRP 0.13.3; tempo por resolução de referencia_pyvrp.TEMPO; 240 s e 480 s na escala), mesma amostra
de refinamento e de certificação, mesma porteira (P_cert >= alfa = 0,90) e mesma avaliação nos 1.000 cenários comuns.
Conjuntos e sementes iguais aos das rodadas já feitas:
  bench   36 instâncias MEDIUM/HIGH, sementes 0-4 (180 resoluções)
  sp      16 instâncias do caso SP, sementes 0-4 (80)
  escala  24 instâncias (n = 200 e 400), sementes 0-2 (72)
  --controle  reexecuta o degrau q = 0,90, semente 0, nas 36 instâncias do benchmark e nas 24 da escala (controle de carga
              da máquina; o da escala foi acrescentado na Entrada 129, antes de rodar a escala)
Campos extras gravados por resolução: n_iter (iterações do PyVRP), sys_cpu_pct (CPU média do sistema durante a
resolução), conjunto, tag ('nom' ou 'ctrl').

Uso: python pesquisa/pyvrp_nominal.py --piloto
     python pesquisa/pyvrp_nominal.py --rodar --conjuntos bench sp --controle --workers 3
     python pesquisa/pyvrp_nominal.py --rodar --conjuntos escala --workers 3
     python pesquisa/pyvrp_nominal.py --analisar
Saída: pesquisa/resultados_pyvrp_nominal/ (jsonl retomável)
"""
import argparse, glob, json, os, random, sys, time
from multiprocessing import Pool

BASE = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
sys.path.insert(0, os.path.join(BASE, 'src')); sys.path.insert(0, os.path.join(BASE, 'pesquisa'))
sys.path.insert(0, os.path.join(BASE, 'pesquisa', 'caso_sp'))
OUT = os.path.join(BASE, 'pesquisa', 'resultados_pyvrp_nominal')
ARQ = os.path.join(OUT, 'pipeline_nominal.jsonl')
ARQ_PILOTO = os.path.join(OUT, 'piloto.jsonl')
SEM = {'bench': (0, 1, 2, 3, 4), 'sp': (0, 1, 2, 3, 4), 'escala': (0, 1, 2)}
Q_CTRL = 0.90


def lista_tarefas(conjuntos, controle):
    T = []
    if 'bench' in conjuntos:
        for p in sorted(glob.glob(os.path.join(BASE, 'instances', 'generated', 'CCVRPTW-ST_*.json'))):
            if p.endswith(('_MEDIUM.json', '_HIGH.json')):
                T += [('bench', p, 0.0, s, 'nom') for s in SEM['bench']]
                if controle:
                    T.append(('bench', p, Q_CTRL, 0, 'ctrl'))
    if 'sp' in conjuntos:
        for p in sorted(glob.glob(os.path.join(BASE, 'instances', 'caso_sp', 'CCVRPTW-ST_SP*.json'))):
            T += [('sp', p, 0.0, s, 'nom') for s in SEM['sp']]
    if 'escala' in conjuntos:
        for p in sorted(glob.glob(os.path.join(BASE, 'instances', 'escala', 'CCVRPTW-ST_*.json'))):
            T += [('escala', p, 0.0, s, 'nom') for s in SEM['escala']]
            if controle:                      # acréscimo registrado na Entrada 129 (antes de rodar a escala)
                T.append(('escala', p, Q_CTRL, 0, 'ctrl'))
    return T


def _chave(t):
    return (os.path.basename(t[1])[:-5], t[2], t[3], t[4])


def _tarefa(t):
    """Uma resolução (conjunto, caminho, q, semente, rótulo) com os campos extras de carga."""
    conj, path, q, sd, tag = t
    import psutil
    import referencia_pyvrp as rp
    rp.TEMPO.update({200: 240, 400: 480})
    if conj == 'escala':                  # mesmo contexto de rodar_escala._init
        import rodar_caso_sp
        rodar_caso_sp.aplicar_matrizes()
    import pipeline_pyvrp_buffer as pp
    from pyvrp import Model
    orig = Model.solve
    cap = {}

    def solve(self, *a, **k):
        res = orig(self, *a, **k)
        cap['n_iter'] = int(res.num_iterations)
        cap['feasivel_pyvrp'] = bool(res.is_feasible())
        cap['excesso_carga'] = int(sum(res.best.excess_load()))
        cap['time_warp'] = int(res.best.time_warp())
        cap['rotas'] = [[int(v) for v in r.visits()] for r in res.best.routes()]
        return res
    Model.solve = solve
    psutil.cpu_percent(interval=None)     # inicia a janela de medição
    t0 = time.time()
    try:
        r = pp.tarefa((path, q, sd, conj == 'sp'))
    finally:
        Model.solve = orig
    r.update(conjunto=conj, tag=tag, n_iter=cap.get('n_iter'), sys_cpu_pct=round(psutil.cpu_percent(interval=None), 1),
             wall_total_s=round(time.time() - t0, 1), feasivel_pyvrp=cap.get('feasivel_pyvrp'),
             excesso_carga=cap.get('excesso_carga'), time_warp=cap.get('time_warp'))
    r['viavel_tol'] = _viavel_tol(path, q, cap)
    return r


def _viavel_tol(path, q, cap, tol=1e-6):
    """Viabilidade com tolerância de ponto flutuante (Entrada 131): a regra estrita de pipeline_pyvrp_buffer.tarefa
    (`load > cap`, `t > due`) rejeita rota exatamente cheia por erro de soma (619,0000000000001 > 619,0)."""
    import nsga2_simheuristic as ns
    import ccvrp_c as cc
    inst = ns.load_instance(path); n = len(inst['customers'])
    ctx = cc.Contexto(inst, 0, buffer_q=q, decoder='split', p_ls=0.0)
    R = cap.get('rotas', [])
    if not cap.get('feasivel_pyvrp') or sorted(c for r in R for c in r) != list(range(1, n + 1)):
        return False
    tt = ctx.ttime_plan
    for r in R:
        load = 0.0; t = 0.0; prev = 0
        for c in r:
            load += ctx.dem[c]
            if load > ctx.cap + tol:
                return False
            t += tt[prev][c]
            if t > ctx.due[c] + tol:
                return False
            t = max(t, ctx.rdy[c]) + ctx.svc[c]; prev = c
    return True


def _rodar_lista(T, arq, workers):
    os.makedirs(OUT, exist_ok=True)
    print(f'{len(T)} resoluções a fazer; {workers} processos', flush=True)
    t0 = time.time()
    with Pool(workers) as pool, open(arq, 'a', encoding='utf-8') as fh:
        for i, r in enumerate(pool.imap_unordered(_tarefa, T), 1):
            fh.write(json.dumps(r) + '\n'); fh.flush()
            print(f"  [{i}/{len(T)}] {r['conjunto']} {r['instance']} q={r['q']} s{r['seed']} {r['tag']} viável={r['viavel']} "
                  f"iter={r['n_iter']} cpu_sys={r['sys_cpu_pct']}% custo={r['base_f1']:.1f} "
                  f"P_cert {r['base_p_cert']:.3f} -> {r['ref_p_cert']:.3f} P_oos {r['base_p_oos']:.3f} "
                  f"[{(time.time() - t0) / 60:.1f} min]", flush=True)
    print('FIM', flush=True)


def rodar(conjuntos, workers, controle):
    feitas = set()
    if os.path.exists(ARQ):
        for L in open(ARQ, encoding='utf-8'):
            x = json.loads(L); feitas.add((x['instance'], x['q'], x['seed'], x['tag']))
    T = [t for t in lista_tarefas(conjuntos, controle) if _chave(t) not in feitas]
    random.Random(2026).shuffle(T)
    print(f'{len(feitas)} já feitas', flush=True)
    _rodar_lista(T, ARQ, workers)


def piloto(workers):
    """Seis resoluções para conferir o caminho de código e medir tempos; NÃO entram na análise."""
    g = os.path.join(BASE, 'instances', 'generated'); e = os.path.join(BASE, 'instances', 'escala')
    s = os.path.join(BASE, 'instances', 'caso_sp')
    T = [('bench', os.path.join(g, 'CCVRPTW-ST_C1_025_HIGH.json'), 0.0, 0, 'piloto'),
         ('bench', os.path.join(g, 'CCVRPTW-ST_C1_050_HIGH.json'), 0.0, 0, 'piloto'),
         ('bench', os.path.join(g, 'CCVRPTW-ST_C1_100_HIGH.json'), 0.0, 0, 'piloto'),
         ('bench', os.path.join(g, 'CCVRPTW-ST_C1_050_HIGH.json'), Q_CTRL, 0, 'piloto'),
         ('sp', os.path.join(s, 'CCVRPTW-ST_SP1_050_HIGH.json'), 0.0, 0, 'piloto'),
         ('escala', os.path.join(e, 'CCVRPTW-ST_C1_200_HIGH.json'), 0.0, 0, 'piloto')]
    _rodar_lista(T, ARQ_PILOTO, workers)


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--piloto', action='store_true')
    ap.add_argument('--rodar', action='store_true')
    ap.add_argument('--analisar', action='store_true')
    ap.add_argument('--conjuntos', nargs='+', default=['bench', 'sp'], choices=['bench', 'sp', 'escala'])
    ap.add_argument('--controle', action='store_true')
    ap.add_argument('--workers', type=int, default=3)
    a = ap.parse_args()
    if a.piloto:
        piloto(a.workers)
    if a.rodar:
        rodar(a.conjuntos, a.workers, a.controle)
    if a.analisar:
        raise SystemExit('análise: ver analisar_pyvrp_nominal.py (escrito depois dos dados, com o plano da Entrada 128)')
