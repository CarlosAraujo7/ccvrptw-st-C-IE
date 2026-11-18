#!/usr/bin/env python3
"""
Teste do limite inferior de temperatura (17/11/2026, diário Entrada 153): src/thermal_route.py com CCVRP_PISO_FRACAO.
  1. Desligado (variável ausente ou 0), route_quality devolve exatamente o mesmo que a versão anterior (backup de 16/10/2026),
     em rotas aleatórias, com cenário nominal e com 200 cenários simulados, com e sem zonas.
  2. Ligado (fração 1): produto refrigerado numa zona a -18 °C congela (Q = 0) e o congelado não; na zona a 2 °C nenhum congela.
  3. A fração é linear: com 0,5 a qualidade dos produtos que congelam é metade da original.
  4. Congelamento só depende da exposição: numa rota curta (menos de 2 h) o produto refrigerado ainda não congela.
  5. (Entrada 154) A constante de tempo do produto é configurável por CCVRP_PISO_TAU_H: ausente ou inválida = 10 h, idêntico à versão
     anterior (backup de 18/11/2026); 2 h dá exatamente o mesmo que trocar a constante na versão anterior; um produto que responde mais
     depressa nunca congela menos.
Uso: python pesquisa/test_piso_temperatura.py
"""
import importlib.util, os, sys
import numpy as np

BASE = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
sys.path.insert(0, os.path.join(BASE, 'src'))
import nsga2_simheuristic as ns   # noqa: E402
import thermal_route as tr        # noqa: E402


def carregar_antigo():
    from importlib.machinery import SourceFileLoader
    carregador = SourceFileLoader('thermal_route_antigo', os.path.join(BASE, 'src', 'thermal_route.py.bak_20261117_piso'))
    spec = importlib.util.spec_from_loader('thermal_route_antigo', carregador)
    m = importlib.util.module_from_spec(spec); carregador.exec_module(m)
    return m


def zonas(inst, rota, T0):
    """uma zona única na temperatura T0 (como o Z1c) para toda a rota"""
    zt = {0: {'T_set_C': T0, 'T_min_C': T0 - 3.0, 'T_max_C': T0 + 5.0, 'sigma_T_C': 2.0}}
    return [0] * len(rota), zt


def main():
    ant = carregar_antigo()
    inst = ns.load_instance(os.path.join(BASE, 'instances', 'generated', 'CCVRPTW-ST_C1_050_HIGH.json'))
    _, ttime = ns.build_matrices(inst)
    cust = inst['customers']
    cat = {c['id']: c['product_category'] for c in cust}
    rng = np.random.default_rng(5)
    nom = tr.cenario_nominal(inst, n_zones=1)
    sc = tr.sortear_cenarios(inst, 200, np.random.default_rng(7), n_zones=1, tempo_viagem=True)
    ok = True

    def conf(nome, cond):
        nonlocal ok
        print(('PASSA ' if cond else 'FALHA ') + nome)
        ok = ok and cond

    os.environ.pop('CCVRP_PISO_FRACAO', None)
    iguais = True
    for _ in range(25):
        rota = [int(x) for x in rng.choice(np.arange(1, 51), size=int(rng.integers(3, 9)), replace=False)]
        zc, zt = zonas(inst, rota, float(rng.choice([-18.0, 2.0, 4.0])))
        for cenario in (nom, sc):
            a = tr.route_quality(rota, inst, ttime, cenario, zc, zt, per_customer_kinetics=True)
            b = ant.route_quality(rota, inst, ttime, cenario, zc, zt, per_customer_kinetics=True)
            c1 = tr.route_quality(rota, inst, ttime, cenario, per_customer_kinetics=True)
            c2 = ant.route_quality(rota, inst, ttime, cenario, per_customer_kinetics=True)
            iguais = iguais and np.array_equal(a, b) and np.array_equal(c1, c2)
    conf('1. desligado: idêntico à versão anterior (25 rotas, nominal e 200 cenários, com e sem zonas)', iguais)

    # rota longa (tempo de sobra) com os quatro tipos de produto
    por_cat = {k: [c for c, v in cat.items() if v == k] for k in ('frozen', 'meat', 'dairy', 'produce')}
    rota = [por_cat['frozen'][0], por_cat['meat'][0], por_cat['dairy'][0], por_cat['produce'][0], por_cat['meat'][1], por_cat['dairy'][1]]
    # força um cronograma longo: janelas de serviço da instância; mede a duração total da rota pelo detalhe
    zc, zt = zonas(inst, rota, -18.0)
    base = ant.route_quality(rota, inst, ttime, nom, zc, zt, per_customer_kinetics=True, detalhe=True)
    os.environ['CCVRP_PISO_FRACAO'] = '1'
    d1 = tr.route_quality(rota, inst, ttime, nom, zc, zt, per_customer_kinetics=True, detalhe=True)
    cats = d1['cat']
    gelou = d1['cong'][0]
    conf('2a. zona a -18 °C: o congelado nunca é penalizado', not any(g for g, k in zip(gelou, cats) if k == 'frozen'))
    refrig = [g for g, k in zip(gelou, cats) if k != 'frozen']
    print('     produtos refrigerados que congelaram na zona a -18 °C (ordem da rota):', [bool(g) for g in gelou], cats)
    conf('2b. zona a -18 °C: ao menos o último cliente refrigerado da rota congelou', bool(gelou[-1]) or bool(gelou[-2]) or bool(gelou[-3]))
    conf('2c. onde congelou, Q = 0; onde não, Q igual ao original', all((d1['Q'][0, j] == 0.0) if gelou[j] else (d1['Q'][0, j] == base['Q'][0, j]) for j in range(len(rota))))

    zc2, zt2 = zonas(inst, rota, 2.0)
    base2 = ant.route_quality(rota, inst, ttime, nom, zc2, zt2, per_customer_kinetics=True, detalhe=True)
    d2 = tr.route_quality(rota, inst, ttime, nom, zc2, zt2, per_customer_kinetics=True, detalhe=True)
    conf('2d. zona a 2 °C: nenhum produto congela e Q não muda', (not d2['cong'].any()) and np.array_equal(d2['Q'], base2['Q']))
    ds = tr.route_quality(rota, inst, ttime, sc, zc2, zt2, per_customer_kinetics=True, detalhe=True)
    conf('2e. zona a 2 °C, 200 cenários: ainda nenhum congela', not ds['cong'].any())

    os.environ['CCVRP_PISO_FRACAO'] = '0.5'
    dm = tr.route_quality(rota, inst, ttime, nom, zc, zt, per_customer_kinetics=True, detalhe=True)
    conf('3. fração 0,5: Q dos que congelam é metade do original', all(np.isclose(dm['Q'][0, j], 0.5 * base['Q'][0, j]) for j in range(len(rota)) if gelou[j]))

    # 4. rota curta: um só cliente refrigerado numa zona a -18 °C não chega a 2 h de exposição
    os.environ['CCVRP_PISO_FRACAO'] = '1'
    uma = [por_cat['meat'][0]]
    zc3, zt3 = zonas(inst, uma, -18.0)
    viagem_min = ttime[0][uma[0]]
    pronto = cust[uma[0] - 1]['ready_time_min']
    t0_justo = max(0.0, pronto - viagem_min)                # sai de modo a chegar na abertura da janela: sem espera
    d3 = tr.route_quality(uma, inst, ttime, nom, zc3, zt3, per_customer_kinetics=True, detalhe=True, t0_min=t0_justo)
    d3b = tr.route_quality(uma, inst, ttime, nom, zc3, zt3, per_customer_kinetics=True, detalhe=True, t0_min=0.0)
    espera_h = max(0.0, pronto - viagem_min) / 60.0
    print(f'     rota de um cliente: viagem {viagem_min / 60:.2f} h, janela abre em {pronto / 60:.2f} h; '
          f'saindo justo: Tp_min = {d3["Tp_min"][0, 0]:.2f} °C (congelou = {bool(d3["cong"][0, 0])}); '
          f'saindo em 0 (espera de {espera_h:.2f} h): Tp_min = {d3b["Tp_min"][0, 0]:.2f} °C (congelou = {bool(d3b["cong"][0, 0])})')
    conf('4. a exposição manda: viagem curta sem espera não congela; com a espera de mais de 2 h congela',
         (not d3['cong'][0, 0]) and (bool(d3b['cong'][0, 0]) == (espera_h + viagem_min / 60.0 > 2.0)))
    # 5. constante de tempo do produto configurável (Entrada 154): comparação com a versão de 18/11/2026, que tem a constante fixa
    from importlib.machinery import SourceFileLoader
    carregador = SourceFileLoader('thermal_route_bak_tau', os.path.join(BASE, 'src', 'thermal_route.py.bak_20261118_tau'))
    spec = importlib.util.spec_from_loader('thermal_route_bak_tau', carregador)
    bak = importlib.util.module_from_spec(spec); carregador.exec_module(bak)
    rng5 = np.random.default_rng(11)
    casos = []
    for _ in range(25):
        rota5 = [int(x) for x in rng5.choice(np.arange(1, 51), size=int(rng5.integers(3, 9)), replace=False)]
        zc5, zt5 = zonas(inst, rota5, float(rng5.choice([-18.0, 2.0, 4.0])))
        casos.append((rota5, zc5, zt5))

    def varrer(mod):
        saida = []
        for rota5, zc5, zt5 in casos:
            for cenario in (nom, sc):
                d = mod.route_quality(rota5, inst, ttime, cenario, zc5, zt5, per_customer_kinetics=True, detalhe=True)
                saida.append((d['Q'], d['Tp_min'], d['cong']))
        return saida

    def iguais_a(a, b):
        return all(np.array_equal(x[0], y[0]) and np.array_equal(x[1], y[1], equal_nan=True) and np.array_equal(x[2], y[2]) for x, y in zip(a, b))

    os.environ['CCVRP_PISO_FRACAO'] = '1'
    os.environ.pop('CCVRP_PISO_TAU_H', None)
    bak.TAU_P_PISO_H = 10.0
    padrao = varrer(tr)
    conf('5a. variável ausente: idêntico à versão com a constante de 10 h (25 rotas, nominal e 200 cenários)', iguais_a(padrao, varrer(bak)))
    invalidos = True
    for v in ('abc', '0', '-3', ''):
        os.environ['CCVRP_PISO_TAU_H'] = v
        invalidos = invalidos and iguais_a(padrao, varrer(tr))
    conf('5b. valor inválido, zero ou negativo: volta às 10 h', invalidos)
    os.environ['CCVRP_PISO_TAU_H'] = '2'
    rapido = varrer(tr)
    bak.TAU_P_PISO_H = 2.0
    conf('5c. CCVRP_PISO_TAU_H = 2: idêntico a trocar a constante da versão anterior por 2 h', iguais_a(rapido, varrer(bak)))
    mais = all(np.all(r[2] | ~p[2]) for r, p in zip(rapido, padrao))                      # quem congela com 10 h congela com 2 h
    frio = all(np.all((r[1] <= p[1] + 1e-12) | np.isnan(p[1])) for r, p in zip(rapido, padrao))
    n_mais = sum(int(r[2].sum()) for r in rapido); n_pad = sum(int(p[2].sum()) for p in padrao)
    print(f'     produtos que congelam nas 25 rotas (cenários somados): {n_pad} com 10 h, {n_mais} com 2 h')
    conf('5d. um produto que responde mais depressa nunca congela menos (e, aqui, congela mais)', mais and frio and n_mais > n_pad)
    os.environ.pop('CCVRP_PISO_TAU_H', None)
    os.environ.pop('CCVRP_PISO_FRACAO', None)
    print('\nRESULTADO:', 'TODOS PASSAM' if ok else 'HÁ FALHAS')
    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()
