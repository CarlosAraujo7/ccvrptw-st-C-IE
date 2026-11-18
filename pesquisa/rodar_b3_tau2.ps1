# B3 com constante de tempo do produto de 2 h (18/11/2026, diário Entrada 154): Z1c e Z = 2 (180 execuções) com o limite inferior de temperatura
# (CCVRP_PISO_FRACAO = 1) e o produto respondendo em 2 h (CCVRP_PISO_TAU_H = 2); saída em results_c_b3_lb_tau2. O plano da análise está no cabeçalho de
# pesquisa/analisar_piso_b3.py (painel C), gravado antes de a primeira execução existir.
# Como rodar: powershell -ExecutionPolicy Bypass -File .\pesquisa\rodar_b3_tau2.ps1
# Retomável: o script de busca pula o que já tem _metrics.json; ao fim confere a pureza das frentes e faz a avaliação fora da amostra.
# Mantém o computador acordado enquanto roda. Lançado fora da árvore de processos do app (WMI), para sobreviver ao fechamento da sessão.
$ErrorActionPreference = "Continue"
$projeto = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path
Set-Location $projeto
$env:PYTHONIOENCODING = "utf-8"
$env:CCVRP_PISO_FRACAO = "1"
$env:CCVRP_PISO_TAU_H = "2"
$res = "results_c_b3_lb_tau2"
New-Item -ItemType Directory -Force -Path $res | Out-Null
$log = Join-Path $projeto "$res\log_rodada.txt"
$py = (Get-Command python -ErrorAction SilentlyContinue).Source
if (-not $py) { $py = Join-Path $env:LOCALAPPDATA "Programs\Python\Python312\python.exe" }
Add-Type -Name Energia -Namespace CCVRPTau2 -MemberDefinition '[DllImport("kernel32.dll")] public static extern uint SetThreadExecutionState(uint esFlags);'
[CCVRPTau2.Energia]::SetThreadExecutionState([uint32]2147483649) | Out-Null

"INICIO B3 tau_p = 2 h $(Get-Date)  (python: $py)" | Tee-Object -FilePath $log -Append
& $py -X utf8 pesquisa\frota_dedicada.py --rodar --variantes Z1c Z2 --tamanhos 50 --workers 2 --res $res 2>&1 | Tee-Object -FilePath $log -Append
"FIM B3 tau_p = 2 h $(Get-Date)" | Tee-Object -FilePath $log -Append
