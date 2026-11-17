#!/bin/bash
# Encadeia as etapas do Bloco B depois da fase 1 do B1 (diário Entradas 128-132), no máximo 3 processos por vez:
#   (espera a fase 1)  ->  B3 (frota dedicada, n = 50)  ->  B2 (SIM+C com 100 cenários, n = 25 e 100)
# A fase 2 do B1 (nominal na escala, 72 resoluções + 24 de controle, ≈ 9,6 h de CPU) NÃO entra: as resoluções nominais da escala
# repetem a iteração 0 do PyVRP-SIM (mesmas instâncias, sementes e tempos), já gravada em resultados_pyvrp_sim/pyvrp_sim.jsonl
# (Entrada 132). Para rodá-la mesmo assim: python -X utf8 -I pesquisa/pyvrp_nominal.py --rodar --conjuntos escala --controle --workers 3
# Todas as etapas são retomáveis; se uma falhar a cadeia para (&&) e o erro fica no log da etapa.
# Uso (na raiz do projeto, em segundo plano): bash pesquisa/orquestrar_bloco_b.sh
cd /c/Pesquisas/CCVRPTW-ST || exit 1
export PYTHONIOENCODING=utf-8
PID=$(powershell.exe -NoProfile -Command "(Get-CimInstance Win32_Process -Filter \"Name='python.exe'\" | Where-Object { \$_.CommandLine -like '*pyvrp_nominal.py*--conjuntos bench sp*' -and \$_.CommandLine -notlike '*spawn_main*' } | Select-Object -First 1).ProcessId" | tr -d '\r\n ')
echo "$(date) esperando a fase 1 do B1 (PID '$PID')"
if [ -n "$PID" ]; then
  powershell.exe -NoProfile -Command "Wait-Process -Id $PID -ErrorAction SilentlyContinue"
fi
echo "$(date) fase 1 terminou; iniciando B3"
python -X utf8 -I pesquisa/frota_dedicada.py --rodar --variantes D4 D2 --tamanhos 50 --workers 3 > pesquisa/log_b3_rodada.txt 2>&1 \
 && echo "$(date) B3 concluído; iniciando B2" \
 && python -X utf8 -I src/experiment_runner.py --algorithms C-SIMC100 --sizes 25 100 --uncertainty MEDIUM HIGH --seeds 5 --workers 3 --resume --results-dir results_c_extra > pesquisa/log_b2_rodada.txt 2>&1 \
 && echo "$(date) B2 concluído; FIM DA ORQUESTRAÇÃO (a escala do B1 não foi rodada: ver Entrada 132)"
echo "$(date) orquestrador terminou (código $?)"
