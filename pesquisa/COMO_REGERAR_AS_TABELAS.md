# Como regerar as tabelas do artigo (mapa tabela → script)

Atualizado em 16/11/2026 (diário, Entradas 121-124). Cada tabela do artigo e do suplemento é gerada por um script a partir de resultados brutos guardados; nenhuma é digitada à mão. Rodar na raiz do projeto (`C:\Pesquisas\CCVRPTW-ST`), com Python 3.12 e `pip install -r repositorio_github/ccvrptw-st/requirements.txt`.

| Tabela (arquivo) | Script | Entrada principal |
|---|---|---|
| `tab_conformidade`, `tab_custo`, `tab_hv`, `tab_preco`, `tab_b2`, `tab_clima`, `tab_exato`, `tab_parteb` | `pesquisa/gerar_tabelas_c.py` | `pesquisa/resultados_c/por_execucao.csv`, `results_c/*/*_oosC.json`, `results_c_b2` |
| `tab_tempo` | `pesquisa/tempo_microbenchmark.py --tabela` | `results_c_tempo/` (microbenchmark controlado em máquina livre; NÃO usar os `elapsed_s` de `results_c`, medidos em sessões com cargas diferentes: o teste de regeneração apontou a divergência e `gerar_tabelas_c.py` agora grava essa versão como `tab_tempo_execucao_principal.tex`) |
| `tab_equivalencia` | `pesquisa/equivalencia.py` | `pesquisa/resultados_robustez/` |
| `tab_robustez` | `pesquisa/gerar_tab_robustez.py` | `pesquisa/resultados_robustez/decisor.json` |
| `tab_refino` | `pesquisa/analisar_refinamento.py` | `pesquisa/resultados_refinamento/refinamento_*.json` |
| `tab_acal` | `pesquisa/resumo_acal.py` | `results_c_acal_*` |
| `tab_regra_lcb` | `pesquisa/regra_lcb.py` | `results_c*`, `results_c_escala`, `pesquisa/resultados_refinamento` |
| `tab_caso_sp`, `tab_caso_sp_pico` | `pesquisa/caso_sp/tabela_caso_sp.py` | `results_c_sp`, `results_c_sp_pico` |
| `tab_comparar_td` | `pesquisa/caso_sp/comparar_td.py` | `results_c_sp_td` e variantes |
| `tab_escala` | `pesquisa/escala/analisar_escala.py` | `results_c_escala` |
| `tab_degrau_nominal` | `pesquisa/analisar_degrau_nominal.py` | resultados existentes |
| `tab_pyvrp_pipeline` | `pesquisa/pipeline_pyvrp_buffer.py` | resultados do PyVRP |
| `tab_pyvrp_sim` | `pesquisa/pyvrp_sim.py` | `pyvrp_sim.jsonl` |
| `tab_teoria` | `pesquisa/validar_teoria_buffer.py` | `results_c`, `results_c_sp` |
| `tab_observado` (três conjuntos) | `pesquisa/tempos_observados/tabela_observado_3conj.py` | `pesquisa/tempos_observados/resultados/*.json` |

Novos neste ciclo (ciclo C1, diário Entradas 121-124; saem direto no pacote C&IE `tables/`):

| Tabela | Script | Observação |
|---|---|---|
| `tab_preco_completo` (S23) | `pesquisa/custo_completo.py` | Deriva as horas de rota de f1, f3 e do número de rotas; verifica D/R em [0,3; 16] h |
| `tab_servico_alt` (S24) | `pesquisa/servico_alternativo.py` | Reavalia o plano escolhido nos 1.000 cenários comuns |
| `tab_parteb_estendida` (S25) | `pesquisa/parte_b_estendida.py`, depois `pesquisa/tabela_parte_b_estendida.py` | A linha "base" deve reproduzir a Tabela S19 |

Bloco B (diário Entradas 128-130; saem para o pacote C&IE `tables/` depois de conferidas):

| Tabela | Script | Observação |
|---|---|---|
| `tab_pyvrp_nominal_cert`, `tab_pyvrp_nominal_escada` (B1) | `pesquisa/pyvrp_nominal.py` (execução; retomável) e `pesquisa/analisar_pyvrp_nominal.py` | PyVRP sem buffer; o controle de reexecução (q = 0,90, semente 0) e o número de iterações ficam no relatório |
| `tab_simc100_completo` (B2) | `src/experiment_runner.py --algorithms C-SIMC100 --sizes 25 100 ...` e `pesquisa/analisar_simc100_completo.py` | SIM+C com 100 cenários de busca; n = 50 vem da Entrada 81 |
| `tab_frota_dedicada` (B3) | `pesquisa/frota_dedicada.py` e `pesquisa/analisar_frota_dedicada.py` | pureza das rotas conferida pelo script; legal via `conformidade_legal_b2.py` |
| (orquestração) | `pesquisa/orquestrar_bloco_b.sh` | encadeia B3, B2 e a fase 2 do B1 depois da fase 1 |

## Avisos de reprodutibilidade

- **Dentro de um mesmo ambiente**, a mesma semente reproduz o mesmo resultado (`pesquisa/test_determinismo_det.py`, 4/4 PASS). **Entre ambientes** (versão de numpy/BLAS, hardware, sistema operacional) as frentes de Pareto podem diferir (Entradas 18, 20 e 21 do diário); por isso o repositório e o Zenodo guardam os **resultados brutos**, e as tabelas se regeneram deles sem reotimizar.
- Laços limitados por tempo de parede (protocolo oficial de comparação) tornam a execução do zero sensível à carga da máquina; replay bit a bit só com o modo de iterações fixas (Lei Geral, II.9).
- Dados de terceiros: OpenStreetMap e DT-CARGO (ODbL), Atenas (CC BY 4.0), TomTom Traffic Index e CNPJ da Receita Federal (termos próprios). Ver `DATA_LICENSE.md` do repositório.
- **Teste de regeneração (A7; Entrada 137):** `python pesquisa/testar_regeneracao_tabelas.py` roda os geradores, lê as tabelas geradas e compara os NÚMEROS do corpo de cada uma com os da tabela do pacote de submissão (`--pacote` indica a pasta `tables/`; legendas e numeração não entram, porque foram ajustadas à mão para a C&IE). Por padrão roda só os geradores rápidos; `--tudo` inclui os lentos (`gerar_tabelas_c.py`, `regra_lcb.py`, `validar_teoria_buffer.py`, as variantes do caso SP e `comparar_td.py`) e `--so tab_x` roda uma tabela. O diretório `pesquisa/latex_ejor_c/tables/` é copiado antes e restaurado depois. Relatório em `pesquisa/resultados_teste_regeneracao/relatorio_teste_regeneracao.txt`. Sem gerador testado: `tab_params` e `tab_sp_dados` (tabelas de parâmetros e de dados, sem script gerador).
