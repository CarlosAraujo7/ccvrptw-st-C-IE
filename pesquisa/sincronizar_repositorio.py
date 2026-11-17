#!/usr/bin/env python3
"""
Copia para repositorio_github/ccvrptw-st os arquivos do projeto acrescentados ou alterados nas melhorias de
16-17/10/2026 (tempos dependentes do horário, caso de Atenas, teoria do buffer, degrau nominal; diário
Entradas 101-107). Lista explícita: nada de dados brutos (TomTom, CERTH, NOAA), backups ou logs.
Uso: python pesquisa/sincronizar_repositorio.py
"""
import os, shutil

BASE = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
REPO = os.path.join(BASE, 'repositorio_github', 'ccvrptw-st')
ARQUIVOS = [
    'src/ccvrp_c.py', 'src/thermal_route.py', 'src/tempo_dependente.py',
    'pesquisa/teste_tempo_dependente.py', 'pesquisa/validar_teoria_buffer.py', 'pesquisa/analisar_degrau_nominal.py',
    'pesquisa/pipeline_pyvrp_buffer.py', 'pesquisa/sincronizar_repositorio.py',
    'pesquisa/caso_sp/rodar_caso_sp.py', 'pesquisa/caso_sp/tabela_caso_sp.py', 'pesquisa/caso_sp/gerar_caso_sp_td.py',
    'pesquisa/caso_sp/comparar_td.py', 'pesquisa/caso_sp/extrair_tomtom.py', 'pesquisa/caso_sp/rodar_caso_sp_td.ps1',
    'pesquisa/caso_certh/gerar_caso_atenas.py', 'pesquisa/caso_certh/rodar_caso_atenas.ps1',
    'pesquisa/clima/processar_isd.py', 'pesquisa/clima/perfis_16716199999.json', 'pesquisa/clima/resumo_clima_atenas.txt',
    # calibração (17/10/2026, Entradas 108-110)
    'src/experiment_runner.py', 'pesquisa/rodar_acal.py', 'pesquisa/avaliar_acal.py', 'pesquisa/resumo_acal.py',
    'pesquisa/rodar_acal_teste.ps1', 'pesquisa/calibracao_preditiva_base.py', 'pesquisa/calibracao_preditiva.py',
    'pesquisa/calibracao_preditiva_dev.py', 'pesquisa/refinamento_simulacao.py',
    # tempos de viagem observados (17/10/2026, Entrada 113)
    'pesquisa/tempos_observados/analisar_tlc.py', 'pesquisa/tempos_observados/reavaliar_observado.py',
    'pesquisa/tempos_observados/recertificar_observado.py', 'pesquisa/tempos_observados/tabela_observado.py',
    'pesquisa/tempos_observados/analisar_dtcargo.py', 'pesquisa/tempos_observados/reavaliar_caminhoes.py',
    # escala (Entrada 114)
    'pesquisa/escala/gerar_escala.py', 'pesquisa/escala/rodar_escala.py', 'pesquisa/escala/rodar_escala.ps1',
    'pesquisa/escala/analisar_escala.py',      # análise (18/10/2026, Entrada 116)
    # regra LCB e figura-síntese (18/10/2026, Entrada 117)
    'pesquisa/regra_lcb.py', 'pesquisa/rodar_regra_lcb.py', 'pesquisa/figura_sintese.py',
    # B3/B4 (18/10/2026, Entrada 118)
    'pesquisa/pyvrp_sim.py', 'pesquisa/rodar_ciclo_b3b4.ps1', 'pesquisa/rodar_ciclo_b3b4.cmd',
    # ciclo C1 da revisão para a C&IE: bloco A (16/11/2026, Entradas 121-127) e bloco B (16-17/11/2026, Entradas 128-136)
    'pesquisa/custo_completo.py', 'pesquisa/servico_alternativo.py', 'pesquisa/parte_b_estendida.py',
    'pesquisa/tabela_parte_b_estendida.py', 'pesquisa/tempos_observados/tabela_observado_3conj.py',
    'pesquisa/equivalencia.py', 'pesquisa/COMO_REGERAR_AS_TABELAS.md',
    'pesquisa/pyvrp_nominal.py', 'pesquisa/analisar_pyvrp_nominal.py', 'pesquisa/frota_dedicada.py',
    'pesquisa/analisar_frota_dedicada.py', 'pesquisa/analisar_simc100_completo.py', 'pesquisa/orquestrar_bloco_b.sh',
    'pesquisa/testar_regeneracao_tabelas.py',     # teste de regeneração das tabelas (A7, Entrada 137)
    'pesquisa/gerar_tabelas_c.py',                # saída heterogênea de tempos renomeada (Entrada 137)
    # lançadores alterados desde a última sincronização
    'pesquisa/caso_sp/rodar_caso_sp_v2.ps1', 'pesquisa/rodar_ciclo_20261015.ps1', 'pesquisa/rodar_execucao_c_v2.ps1',
    'pesquisa/rodar_medicoes_maquina_livre.ps1', 'pesquisa/rodar_pos_execucao.ps1',
]
PASTAS = ['instances/caso_sp_td', 'instances/caso_sp_media_operacao', 'instances/caso_atenas',
          'pesquisa/resultados_caso_sp_td', 'pesquisa/resultados_caso_sp_media', 'pesquisa/resultados_caso_atenas',
          'pesquisa/resultados_teoria_buffer', 'pesquisa/resultados_degrau_nominal', 'pesquisa/resultados_acal',
          'pesquisa/resultados_calibracao_preditiva', 'pesquisa/tempos_observados/resultados', 'instances/escala',
          'pesquisa/resultados_escala', 'pesquisa/resultados_regra_lcb', 'pesquisa/resultados_refinamento',
          'pesquisa/resultados_caso_sp/refinamento', 'pesquisa/resultados_pyvrp_sim',
          # ciclo C1 (resumos pequenos; os resultados brutos por execução vão para um arquivo permanente, não para o git)
          'pesquisa/resultados_custo_completo', 'pesquisa/resultados_servico_alternativo', 'pesquisa/resultados_parte_b_estendida',
          'pesquisa/resultados_parte_b_estendida_poshoc', 'pesquisa/resultados_pyvrp_nominal', 'pesquisa/resultados_c_b3',
          'pesquisa/resultados_robustez', 'pesquisa/resultados_teste_regeneracao']   # o LaTeX gerado fica fora do repositório (.gitignore)
IGNORAR = ('.bak', '__pycache__', 'log_', '_backup_')


def main():
    n = 0
    for a in ARQUIVOS:
        o = os.path.join(BASE, a); d = os.path.join(REPO, a)
        os.makedirs(os.path.dirname(d), exist_ok=True); shutil.copy2(o, d); n += 1
    for p in PASTAS:
        for raiz, _, fs in os.walk(os.path.join(BASE, p)):
            if any(x in raiz for x in IGNORAR):
                continue
            for f in fs:
                if any(x in f for x in IGNORAR):
                    continue
                o = os.path.join(raiz, f); d = os.path.join(REPO, os.path.relpath(o, BASE))
                os.makedirs(os.path.dirname(d), exist_ok=True); shutil.copy2(o, d); n += 1
    print(n, 'arquivos copiados para', REPO)


if __name__ == '__main__':
    main()
