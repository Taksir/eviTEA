# Run the clean incumbent (R9E) engine end to end on one development sample and
# compare every stage with the frozen Route-B artifacts. Resumable: rerun the
# same command after an interruption and completed stages are reused from the
# work directory.
#
#   ISO1   (frozen R9E 192 rows, b3e48681...; stage oracles TE_locator_route_b_v2\reports)
#   AKA017 (frozen R9E 792 rows, 1d1110f9...; stage oracles AKA017\route_b_r9e_validation_full55x)
#
#   powershell -ExecutionPolicy Bypass -File .\tools\windows\VALIDATE_INCUMBENT.ps1 -Sample ISO1
#   powershell -ExecutionPolicy Bypass -File "...\VALIDATE_INCUMBENT.ps1" -Sample AKA017 -Workers 6
#
# Expected verdict: INCUMBENT_EQUIVALENCE_PASS in <Work>\INCUMBENT_VALIDATION_REPORT.json
[CmdletBinding(PositionalBinding = $false)]
param(
    [string]$Python = "python",
    [ValidateSet("ISO1", "AKA017")]
    [string]$Sample = "ISO1",
    [string]$Work = "",
    [int]$Workers = 0,
    [Parameter(ValueFromRemainingArguments = $true)]
    [string[]]$Extra = @()
)

$ErrorActionPreference = "Stop"
$env:PYTHONDONTWRITEBYTECODE = "1"
if (-not $Work) { throw "Pass -Work <directory> (a writable work directory)." }
$EviTEARoot = (Resolve-Path (Join-Path $PSScriptRoot "..\..")).Path
$Script = Join-Path $EviTEARoot "validation\incumbent_stage_checks\run_incumbent_validation.py"
$ArgsList = @($Script, "--sample", $Sample, "--work", $Work)
if ($Workers -gt 0) { $ArgsList += @("--workers", "$Workers") }
$ArgsList += $Extra
& $Python @ArgsList
exit $LASTEXITCODE
