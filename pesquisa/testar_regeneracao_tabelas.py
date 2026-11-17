#!/usr/bin/env python3
"""
A7 (diário Entrada 137; Lei Geral I.4): teste de regeneração das tabelas do artigo.

Roda o gerador de cada tabela a partir dos resultados guardados, lê a tabela gerada e compara os NÚMEROS do corpo da tabela
(o que está entre \\begin{tabular} e \\end{tabular}) com os da tabela do pacote de submissão (tables/). Legendas, rótulos e
numeração não entram na comparação (foram ajustados à mão no pacote para a C&IE). PASSA = mesma sequência de números.
Os geradores antigos escrevem em pesquisa/latex_ejor_c/tables/; esse diretório é copiado antes e restaurado depois, de modo que
o teste não deixa rastro ali. Os resultados derivados em pesquisa/resultados_* (relatórios) podem ser reescritos com o mesmo conteúdo.

Uso: python pesquisa/testar_regeneracao_tabelas.py [--pacote "<pasta tables do pacote>"] [--rapido | --tudo] [--so tab_x tab_y ...]
  --rapido  (padrão) só os geradores que levam alguns minutos ou menos; --tudo inclui os lentos.
Saída: pesquisa/resultados_teste_regeneracao/relatorio_teste_regeneracao.txt; código de saída 1 se alguma tabela falhar.
"""
import argparse, os, re, shutil, subprocess, sys, tempfile, time

BASE = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
PACOTE_PADRAO = os.environ.get('CCVRPTW_PACOTE_TABLES', r'C:\Users\vinic\OneDrive\Documentos\Revisão Cold Chain\tables')
LTX = os.path.join('pesquisa', 'latex_ejor_c', 'tables')
RES = 'pesquisa'

# (tabelas geradas, comando após `python -X utf8 -I`, ambiente extra, onde achar cada tabela, rápido?, tempo máximo em s)
# 'onde': '{tmp}' = pasta temporária do teste (opção --tex), ou um caminho relativo ao projeto
REGISTRO = [
    # S22-S25 e observado (geradores com --tex)
    (['tab_preco_completo'], ['pesquisa/custo_completo.py', '--tex', '{tmp}/tab_preco_completo.tex'], {}, {'tab_preco_completo': '{tmp}'}, True, 900),
    (['tab_servico_alt'], ['pesquisa/servico_alternativo.py', '--tex', '{tmp}/tab_servico_alt.tex'], {}, {'tab_servico_alt': '{tmp}'}, True, 900),
    (['tab_parteb_estendida'], ['pesquisa/tabela_parte_b_estendida.py', '--tex', '{tmp}/tab_parteb_estendida.tex'], {}, {'tab_parteb_estendida': '{tmp}'}, True, 300),
    (['tab_observado'], ['pesquisa/tempos_observados/tabela_observado_3conj.py', '--tex', '{tmp}/tab_observado.tex'], {}, {'tab_observado': '{tmp}'}, True, 300),
    # S26-S30 (Bloco B)
    (['tab_pyvrp_nominal_cert', 'tab_pyvrp_nominal_preco', 'tab_pyvrp_nominal_decomp'], ['pesquisa/analisar_pyvrp_nominal.py'], {},
     {t: 'pesquisa/resultados_pyvrp_nominal' for t in ('tab_pyvrp_nominal_cert', 'tab_pyvrp_nominal_preco', 'tab_pyvrp_nominal_decomp')}, True, 600),
    (['tab_simc100_completo'], ['pesquisa/analisar_simc100_completo.py', '--sem-oos'], {}, {'tab_simc100_completo': 'pesquisa/resultados_robustez'}, True, 600),
    (['tab_frota_dedicada'], ['pesquisa/analisar_frota_dedicada.py', '--workers', '2'], {}, {'tab_frota_dedicada': 'pesquisa/resultados_c_b3'}, False, 1800),
    # geradores antigos (escrevem em latex_ejor_c/tables)
    (['tab_equivalencia'], ['pesquisa/equivalencia.py'], {}, {'tab_equivalencia': LTX}, True, 600),
    (['tab_refino'], ['pesquisa/analisar_refinamento.py', '--dir', 'pesquisa/resultados_refinamento', '--tabela', '{tmp}/tab_refino.tex'], {}, {'tab_refino': '{tmp}'}, True, 600),
    (['tab_robustez'], ['pesquisa/gerar_tab_robustez.py'], {}, {'tab_robustez': LTX}, True, 600),
    (['tab_acal'], ['pesquisa/resumo_acal.py'], {}, {'tab_acal': LTX}, True, 900),
    (['tab_escala'], ['pesquisa/escala/analisar_escala.py'], {}, {'tab_escala': LTX}, True, 1200),
    (['tab_degrau_nominal'], ['pesquisa/analisar_degrau_nominal.py'], {}, {'tab_degrau_nominal': LTX}, True, 1200),
    (['tab_pyvrp_pipeline'], ['pesquisa/pipeline_pyvrp_buffer.py', '--analisar'], {}, {'tab_pyvrp_pipeline': LTX}, True, 900),
    (['tab_pyvrp_sim'], ['pesquisa/pyvrp_sim.py', '--analisar'], {}, {'tab_pyvrp_sim': LTX}, True, 900),
    (['tab_caso_sp'], ['pesquisa/caso_sp/tabela_caso_sp.py'], {}, {'tab_caso_sp': LTX}, False, 1800),
    (['tab_caso_sp_pico'], ['pesquisa/caso_sp/tabela_caso_sp.py'], {'CCVRP_SP_VARIANTE': 'pico'}, {'tab_caso_sp_pico': LTX}, False, 1800),
    (['tab_caso_sp_td'], ['pesquisa/caso_sp/tabela_caso_sp.py'], {'CCVRP_SP_VARIANTE': 'td'}, {'tab_caso_sp_td': LTX}, False, 1800),
    (['tab_caso_atenas'], ['pesquisa/caso_sp/tabela_caso_sp.py'], {'CCVRP_SP_VARIANTE': 'atenas'}, {'tab_caso_atenas': LTX}, False, 1800),
    (['tab_comparar_td'], ['pesquisa/caso_sp/comparar_td.py', '--tabela'], {}, {'tab_comparar_td': LTX}, False, 1800),
    (['tab_regra_lcb'], ['pesquisa/regra_lcb.py'], {}, {'tab_regra_lcb': LTX}, False, 3600),
    (['tab_teoria'], ['pesquisa/validar_teoria_buffer.py'], {}, {'tab_teoria': LTX}, False, 3600),
    (['tab_conformidade', 'tab_custo', 'tab_hv', 'tab_preco', 'tab_b2', 'tab_clima', 'tab_exato', 'tab_parteb'],
     ['pesquisa/gerar_tabelas_c.py'], {}, {t: LTX for t in ('tab_conformidade', 'tab_custo', 'tab_hv', 'tab_preco', 'tab_b2', 'tab_clima',
                                                              'tab_exato', 'tab_parteb')}, False, 3600),
    (['tab_tempo'], ['pesquisa/tempo_microbenchmark.py', '--tabela'], {}, {'tab_tempo': LTX}, True, 600),
]
NAO_COBERTAS = ['tab_params', 'tab_sp_dados']


def numeros(tex):
    """sequência de números dos corpos das tabelas (entre \\begin{tabular} e \\end{tabular}), sem comandos LaTeX"""
    corpos = re.findall(r'\\begin\{tabular\}.*?\\end\{tabular\}', tex, flags=re.S)
    out = []
    for c in corpos:
        c = re.sub(r'\\begin\{tabular\}\{[^}]*\}', ' ', c)
        c = re.sub(r'\\(mBUFnc|mBUF|mCAL|mSIM|mACAL|mDETC|mDET|mMOEAD)\{[^}]*\}', ' ', c)     # macros com argumento numérico
        c = re.sub(r'\\[a-zA-Z]+', ' ', c)
        out += re.findall(r'[-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?', c)
    return out


def comparar(gerada, pacote):
    a, b = numeros(open(gerada, encoding='utf-8').read()), numeros(open(pacote, encoding='utf-8').read())
    if a == b:
        return True, f'{len(a)} números idênticos'
    dif = [(i, x, y) for i, (x, y) in enumerate(zip(a, b)) if x != y][:4]
    return False, f'gerada {len(a)} números, pacote {len(b)}; primeiras diferenças (posição, gerada, pacote): {dif}'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--pacote', default=PACOTE_PADRAO)
    ap.add_argument('--tudo', action='store_true')
    ap.add_argument('--so', nargs='+', default=None)
    a = ap.parse_args()
    tmp = tempfile.mkdtemp(prefix='teste_tabelas_')
    snap = os.path.join(tmp, 'snapshot_latex'); shutil.copytree(os.path.join(BASE, LTX), snap)
    linhas = [f'TESTE DE REGENERAÇÃO DAS TABELAS — pacote: {a.pacote}', time.strftime('%d/%m/%Y %H:%M'), '']
    resultados = {}
    try:
        for tabelas, cmd, env, onde, rapido, limite in REGISTRO:
            if a.so and not set(tabelas) & set(a.so):
                continue
            if not a.tudo and not a.so and not rapido:
                for t in tabelas:
                    resultados[t] = ('PULADA', 'gerador lento (use --tudo)')
                continue
            args = [x.replace('{tmp}', tmp) for x in cmd]
            t0 = time.time()
            try:
                r = subprocess.run([sys.executable, '-X', 'utf8', '-I'] + args, cwd=BASE, capture_output=True, text=True, encoding='utf-8',
                                   errors='replace', timeout=limite, env={**os.environ, **env})
                ok_exec, msg_exec = r.returncode == 0, (r.stderr or '').strip().split('\n')[-1][:160]
            except subprocess.TimeoutExpired:
                ok_exec, msg_exec = False, f'tempo máximo de {limite} s excedido'
            dt = time.time() - t0
            for t in tabelas:
                base_onde = onde[t].replace('{tmp}', tmp)
                gerada = os.path.join(base_onde if os.path.isabs(base_onde) else os.path.join(BASE, base_onde), t + '.tex')
                pacote = os.path.join(a.pacote, t + '.tex')
                if not ok_exec:
                    resultados[t] = ('ERRO', f'o gerador falhou ({msg_exec}); {dt:.0f} s')
                elif not os.path.exists(gerada):
                    resultados[t] = ('ERRO', f'tabela gerada não encontrada: {gerada}')
                elif not os.path.exists(pacote):
                    resultados[t] = ('SEM PACOTE', f'a tabela não existe no pacote ({pacote})')
                else:
                    ok, msg = comparar(gerada, pacote)
                    resultados[t] = ('PASSA' if ok else 'FALHA', f'{msg}; {dt:.0f} s')
            print(f'{tabelas}: ' + ', '.join(f'{resultados[t][0]}' for t in tabelas), flush=True)
    finally:
        # restaura o diretório do LaTeX antigo exatamente como estava
        for nome in os.listdir(os.path.join(BASE, LTX)):
            if nome not in os.listdir(snap):
                os.remove(os.path.join(BASE, LTX, nome))
        for nome in os.listdir(snap):
            shutil.copy2(os.path.join(snap, nome), os.path.join(BASE, LTX, nome))
    for t, (st, msg) in sorted(resultados.items()):
        linhas.append(f'{st:10s} {t:26s} {msg}')
    c = {k: sum(1 for s, _ in resultados.values() if s == k) for k in ('PASSA', 'FALHA', 'ERRO', 'SEM PACOTE', 'PULADA')}
    linhas += ['', 'Resumo: ' + ', '.join(f'{k} {v}' for k, v in c.items()), 'Tabelas sem gerador testado: ' + '; '.join(NAO_COBERTAS)]
    txt = '\n'.join(linhas); print('\n' + txt)
    saida = os.path.join(BASE, 'pesquisa', 'resultados_teste_regeneracao'); os.makedirs(saida, exist_ok=True)
    open(os.path.join(saida, 'relatorio_teste_regeneracao.txt'), 'w', encoding='utf-8').write(txt + '\n')
    shutil.rmtree(tmp, ignore_errors=True)
    sys.exit(1 if (c['FALHA'] or c['ERRO']) else 0)


if __name__ == '__main__':
    main()
