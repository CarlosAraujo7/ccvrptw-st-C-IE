#!/usr/bin/env python3
"""
A3 do ciclo C1 (diário, Entradas 121-123, 16/11/2026): sensibilidade térmica ESTENDIDA da Tabela S19.

Pergunta de revisor: "'temperatura não muda a decisão' vale para outra cinética, outro ruído, ou uma medida com limiar?"
Reaproveita parte_b_caracterizacao.py (mesmas 18 instâncias com n = 50, mesmas frentes de DET, BUF0.90+C e SIM+C, até 150 soluções
distintas por instância, K = 200 cenários, números aleatórios comuns) e troca só as condições. Desenho fixado antes de calcular:
  base                 → deve REPRODUZIR a linha "base" da Tabela S19 (E[perda] 2,78%, ρ_min 0,9915, arrependimento 0,0044 p.p.)
  Ea ×0,8 e Ea ×1,2    → energia de ativação escalada em TODAS as categorias, com A0 ajustado para manter a taxa no ponto de ajuste
  ordem zero           → perda = min(1, Λ) em vez de 1 − exp(−Λ) (mesmas taxas)
  viés persistente     → o desvio do termostato de cada cenário é o MESMO em toda a rota (correlação 1), com σ_T × 2
  limiar de qualidade  → "perda" = fração do valor entregue com qualidade retida abaixo de 0,80 (medida não linear)
  estresse combinado   → Ea ×1,2, viés persistente, σ_T × 2 e limiar
Medidas como na Tabela S19: E[perda]; ρ_min (Spearman mínimo, nominal × esperado, entre instâncias); arrependimento máximo
(escolher pelo nominal em vez do esperado, dentro de orçamentos de custo); violação do limiar de temperatura; P̂ de pontualidade.
Uso: python pesquisa/parte_b_estendida.py   Saída: pesquisa/resultados_parte_b_estendida/parte_b.csv
"""
import math, os, sys
import numpy as np

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, AQUI)
import parte_b_caracterizacao as pb   # noqa: E402
tr, cc = pb.tr, pb.cc
_preparar_orig = pb.preparar

pb.BASE_COND.update(ea=1.0, ordem0=False, vies=False, cliff=None)


def condicoes():
    return [('base', {}),
            ('Ea x0.8 (rate at set-point fixed)', {'ea': 0.8}),
            ('Ea x1.2 (rate at set-point fixed)', {'ea': 1.2}),
            ('zero-order kinetics', {'ordem0': True}),
            ('persistent thermostat bias (sigma x2)', {'vies': True, 'escala_sigmaT': 2.0}),
            ('quality threshold Q<0.80', {'cliff': 0.80}),
            ('combined stress', {'ea': 1.2, 'vies': True, 'escala_sigmaT': 2.0, 'cliff': 0.80})]


def preparar(inst, cond):
    inst = _preparar_orig(inst, cond)
    f = cond.get('ea', 1.0)
    if f != 1.0:
        tp = inst['thermal_profile']

        def esc(d, tset):
            ea0 = d['Ea_J_mol']; ea1 = ea0 * f
            d['A0'] = d['A0'] * math.exp((ea1 - ea0) / (tr.R_GAS * (tset + 273.15)))   # mantém k(T_set)
            d['Ea_J_mol'] = ea1
        if 'Ea_J_mol' in tp:
            esc(tp, tp['T_set_C'])
        for c in tp.get('per_category', {}).values():
            esc(c, c.get('T_set_C', tp['T_set_C']))
    return inst


def avaliar(sol, inst, ttime, sc, cond, lim):
    acc = 0.0; tot = 0.0; viol = 0.0; on = None
    sc = pb.ruido_dependente(sol, inst, ttime, sc, cond)
    if cond.get('vies') and np.any(sc['eps']):
        eps = sc['eps'].copy()
        for r in sol['routes']:
            for c in r:
                eps[:, c, :] = eps[:, r[0], :]
        sc = {**sc, 'eps': eps}
    for r, d0 in zip(sol['routes'], sol['partida_min']):
        d = tr.route_quality(r, inst, ttime, sc, t_start_h=cond['t_start'], tau_p_h=cond['tau_p'],
                             per_customer_kinetics=cc.CINETICA_POR_CLIENTE, detalhe=True, t0_min=d0)
        v = d['values']; Q = d['Q']
        if cond.get('ordem0'):
            Q = np.clip(1.0 + np.log(np.maximum(Q, 1e-12)), 0.0, 1.0)
        if cond.get('cliff') is not None:
            acc = acc + ((Q >= cond['cliff']) * v[None, :]).sum(1)
        else:
            acc = acc + (Q * v[None, :]).sum(1)
        tot += v.sum()
        viol = viol + ((d['Tp_max'] > lim) * v[None, :]).sum(1)
        if 'tt' in sc:
            o = (d['late'] <= cc.DELTA_MIN).all(1)
            on = o if on is None else (on & o)
    return 1 - acc / tot, viol / tot, (on.astype(float) if on is not None else None)


pb.condicoes, pb.preparar, pb.avaliar = condicoes, preparar, avaliar

def condicoes_poshoc():
    """PÓS-HOC (diário, Entrada 125): limiares mais próximos da cauda, porque 0,80 nunca é atingido."""
    return [('quality threshold Q<0.95 (post hoc)', {'cliff': 0.95}),
            ('quality threshold Q<0.90 (post hoc)', {'cliff': 0.90}),
            ('combined stress, threshold 0.90 (post hoc)', {'ea': 1.2, 'vies': True, 'escala_sigmaT': 2.0, 'cliff': 0.90})]


if __name__ == '__main__':
    poshoc = '--poshoc' in sys.argv
    if poshoc:
        pb.condicoes = condicoes_poshoc
    out = os.path.join(AQUI, 'resultados_parte_b_estendida_poshoc' if poshoc else 'resultados_parte_b_estendida')
    sys.argv = [sys.argv[0], '--dirs', 'results_c', '--algs', 'C-DET', 'C-BUF90C', 'C-SIMC', '--tamanho', '050',
                '--K', '200', '--out', out]
    pb.main()
