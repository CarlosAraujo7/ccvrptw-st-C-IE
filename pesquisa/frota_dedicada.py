#!/usr/bin/env python3
"""
B3 (16/11/2026, diário Entrada 130): frota dedicada por categoria como referência do multicompartimento.

Variantes (cada rota carrega uma única CLASSE de produto):
  D4  uma classe por categoria (frozen, meat, dairy, produce); cada veículo tem uma zona no set-point da sua categoria
  D2  duas classes: frozen e refrigerados (meat, dairy, produce podem dividir rota); a zona roda no set-point da
      categoria mais fria presente na rota (regra de zona do B2)
A avaliação é a do modelo multizona do B2 com n_zones = 1 (uma zona por veículo): custo, qualidade, CO2 e pontualidade
saem das mesmas funções das linhas Z = 1, 2 e 3 da Tabela S15.

Pureza SEM alterar src/: dentro do processo de cada execução, o módulo ccvrp_c é remendado
  (a) split_decode aplicado separadamente à subsequência de cada classe do cromossomo;
  (b) vizinhança granular da busca local: os k vizinhos mais próximos DENTRO da classe;
  (c) _eliminar_rota só reinsere em rotas da mesma classe.
Conferência (aborta se falhar): todas as rotas de todas as soluções das frentes finais são puras.

Uso:  python pesquisa/frota_dedicada.py --piloto
      python pesquisa/frota_dedicada.py --rodar --variantes D4 D2 --tamanhos 50 --workers 3
      python pesquisa/frota_dedicada.py --conferir          (pureza das frentes já gravadas)
Saída: results_c_b3/C-BUF90C-D4 e -D2 (mesmos arquivos de qualquer execução do motor).
"""
import argparse, glob, json, os, sys, time
from multiprocessing import Pool

BASE = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
sys.path.insert(0, os.path.join(BASE, 'src')); sys.path.insert(0, os.path.join(BASE, 'pesquisa'))
RES = os.path.join(BASE, 'results_c_b3')
ALGO = 'C-BUF90C'
CLASSE = {'D4': {'frozen': 'frozen', 'meat': 'meat', 'dairy': 'dairy', 'produce': 'produce'},
          'D2': {'frozen': 'frozen', 'meat': 'chilled', 'dairy': 'chilled', 'produce': 'chilled'}}
SEMENTES = (0, 1, 2, 3, 4)
_ORIG = {}


def aplicar(variante):
    """Remenda ccvrp_c neste processo (idempotente) para impor rotas puras por classe."""
    import ccvrp_c as cc
    from collections import defaultdict
    if not _ORIG:
        _ORIG.update(split=cc.split_decode, elim=cc._eliminar_rota, init=cc.Contexto.__init__)
    mapa = CLASSE[variante]

    def init(self, inst, *a, **k):
        _ORIG['init'](self, inst, *a, **k)
        cs = inst['customers']; n = len(cs)
        self.classe = [''] + [mapa[c['product_category']] for c in cs]
        self.classes = sorted(set(self.classe[1:]))
        kv = len(self.viz[1]) if n else 10
        self.viz = [[]] + [sorted((j for j in range(1, n + 1) if j != i and self.classe[j] == self.classe[i]),
                                  key=lambda j: self.dist[i][j])[:kv] for i in range(1, n + 1)]

    def split(chrom, ctx):
        rotas = []; viavel = True
        for k in ctx.classes:
            sub = [c for c in chrom if ctx.classe[c] == k]
            if sub:
                r, f = _ORIG['split'](sub, ctx)
                rotas += r; viavel = viavel and f
        return rotas, viavel

    def elim(R, ctx):
        grupos = defaultdict(list)
        for i, r in enumerate(R):
            grupos[ctx.classe[r[0]]].append(i)
        for k in sorted(grupos):
            idx = grupos[k]
            if len(idx) < 2:
                continue
            novo = _ORIG['elim']([R[i] for i in idx], ctx)
            if novo is not None:
                fora = set(idx)
                return [r for i, r in enumerate(R) if i not in fora] + novo
        return None
    cc.Contexto.__init__ = init; cc.split_decode = split; cc._eliminar_rota = elim
    cc._frota_dedicada = variante


def _init(variante):
    aplicar(variante)


def _rodar(t):
    import experiment_runner as er
    return er._worker(t)


def tarefas(variantes, tamanhos, niveis):
    import experiment_runner as er
    cfg = er.ALGO_CONFIG[ALGO]
    T = []
    for v in variantes:
        out_dir = os.path.join(RES, f'{ALGO}-{v}')
        for f in sorted(glob.glob(os.path.join(BASE, 'instances', 'generated', 'CCVRPTW-ST_*.json'))):
            nm = os.path.basename(f)[:-5]; partes = nm.split('_')
            if int(partes[2]) not in tamanhos or (niveis and partes[3] not in niveis):
                continue
            for s in SEMENTES:
                key = f'{nm}__{ALGO}-{v}__seed{s}'
                if os.path.exists(os.path.join(out_dir, key + '_metrics.json')):
                    continue
                T.append({'inst_path': f, 'algo': ALGO, 'algo_label': f'{ALGO}-{v}', 'seed': s, 'pop_size': 100, 'n_gen': 200,
                          'mc_samples': cfg['mc_samples'], 'runner': cfg['runner'], 'out_dir': out_dir, 'result_key': key,
                          'use_multizone': True, 'n_zones': 1, 'dedup': cfg.get('dedup', False), 'c_params': cfg.get('c_params'),
                          'variante': v})
    return T


def _tarefa(t):
    aplicar(t['variante'])
    return _rodar(t)


def conferir(variantes, silencioso=False):
    """Pureza das rotas de todas as soluções das frentes gravadas."""
    import nsga2_simheuristic as ns
    ok = True; n_sol = n_rotas = 0
    for v in variantes:
        mapa = CLASSE[v]
        for f in sorted(glob.glob(os.path.join(RES, f'{ALGO}-{v}', '*_solutions.json'))):
            m = json.load(open(f.replace('_solutions.json', '_metrics.json'), encoding='utf-8'))
            inst = ns.load_instance(os.path.join(BASE, 'instances', 'generated', m['instance'] + '.json'))
            cs = inst['customers']
            for s in json.load(open(f, encoding='utf-8')):
                n_sol += 1
                for r in s['routes']:
                    n_rotas += 1
                    if len({mapa[cs[c - 1]['product_category']] for c in r}) > 1:
                        ok = False
                        print(f'ROTA IMPURA em {os.path.basename(f)}: {r}')
    if not silencioso:
        print(f'pureza: {n_sol} soluções, {n_rotas} rotas, {"todas puras" if ok else "HÁ ROTAS IMPURAS"}')
    return ok


def piloto():
    """Uma instância (n = 50, MEDIUM), semente 0, D4 e D2; confere a pureza e mede o tempo."""
    global RES
    RES = os.path.join(BASE, 'results_c_b3_piloto')
    T = [t for t in tarefas(['D4', 'D2'], [50], ['MEDIUM']) if t['seed'] == 0 and '_C1_' in t['inst_path']]
    print(f'{len(T)} tarefas de piloto', flush=True)
    t0 = time.time()
    with Pool(2) as pool:
        for r in pool.imap_unordered(_tarefa, T):
            print(f"  {r.get('key', '')} {r.get('status', '')} {r.get('elapsed', '')} [{(time.time() - t0) / 60:.1f} min]", flush=True)
    conferir(['D4', 'D2'])


def avaliar_oos(variantes, workers):
    """Avaliação fora da amostra (os mesmos 1.000 cenários comuns) das frentes gravadas: cria os *_oosC.json."""
    import analisar_c as ac
    arqs = sorted(f for v in variantes for f in glob.glob(os.path.join(RES, f'{ALGO}-{v}', '*_solutions.json')))
    with Pool(workers) as pool:
        pool.map(ac._job, arqs, chunksize=2)
    print(f'avaliação fora da amostra: {len(arqs)} arquivos', flush=True)


def rodar(variantes, tamanhos, niveis, workers):
    os.makedirs(RES, exist_ok=True)
    T = tarefas(variantes, tamanhos, niveis)
    T.sort(key=lambda t: -int(os.path.basename(t['inst_path']).split('_')[2]))
    print(f'{len(T)} execuções a rodar', flush=True)
    t0 = time.time()
    with Pool(workers) as pool:
        for i, r in enumerate(pool.imap_unordered(_tarefa, T), 1):
            print(f"  [{i}/{len(T)}] {r.get('key', '')} {r.get('status', '')} {r.get('elapsed', '')} [{(time.time() - t0) / 60:.1f} min]", flush=True)
    print('execuções concluídas; conferindo a pureza das frentes...', flush=True)
    if not conferir(variantes):
        raise SystemExit('FALHA: há rotas impuras; resultados não devem ser usados')
    avaliar_oos(variantes, workers)
    print('FIM B3 (execuções e avaliação fora da amostra)', flush=True)


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--piloto', action='store_true')
    ap.add_argument('--rodar', action='store_true')
    ap.add_argument('--conferir', action='store_true')
    ap.add_argument('--variantes', nargs='+', default=['D4', 'D2'], choices=['D4', 'D2'])
    ap.add_argument('--tamanhos', nargs='+', type=int, default=[50])
    ap.add_argument('--niveis', nargs='+', default=None)
    ap.add_argument('--workers', type=int, default=3)
    a = ap.parse_args()
    if a.piloto:
        piloto()
    if a.rodar:
        rodar(a.variantes, a.tamanhos, a.niveis, a.workers)
    if a.conferir:
        conferir(a.variantes)
