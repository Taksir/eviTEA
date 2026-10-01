# End-to-end validation of the public `evitea call` pipeline (clean incumbent +
# v3 extension) on ISO1 and AKA017.  Prepares reference/TE assets once (adopts
# the frozen dm6 Bowtie2 index; BUILDS both v3 codebooks with the EviTEA
# language stages and checks them against the frozen codebooks), runs
# `evitea call --profile <sample>` from each prepared msBWT, and compares every
# layer hash (R9E, bilateral, asymmetric, full memory, INTERNAL21) and
# calls.native.tsv with the frozen callsets.
#
# Resumable: re-run the same command after an interruption.
# Expected verdict: EVITEA_CALL_END_TO_END_EQUIVALENCE_PASS in
#   <Work>\EVITEA_CALL_END_TO_END_EQUIVALENCE.json
#
#   powershell -ExecutionPolicy Bypass -File .\tools\windows\VALIDATE_EVITEA_CALL.ps1
#   ... -Samples ISO1          (ISO1 only; also proves the ISO1 R9E incumbent)
[CmdletBinding(PositionalBinding = $false)]
param(
    [string]$Python = "python",
    [string]$EdlibPython = "",
    [Parameter(Mandatory = $true)][string]$Work,
    [Parameter(Mandatory = $true)][string]$RouteBV3,
    [int]$Threads = 8,
    [int]$Workers = 0,
    [string[]]$Samples = @("ISO1", "AKA017"),
    [ValidateSet("built", "supplied")]
    [string]$Languages = "built"
)

$ErrorActionPreference = "Stop"
$env:PYTHONDONTWRITEBYTECODE = "1"
$EviTEARoot = (Resolve-Path (Join-Path $PSScriptRoot "..\..")).Path
$Script = Join-Path $EviTEARoot "validation\end_to_end\run_call_equivalence.py"
$ArgsList = @($Script, "--python", $Python, "--work", $Work, "--route-b-v3", $RouteBV3,
              "--threads", "$Threads", "--languages", $Languages,
              "--samples") + $Samples
if ($EdlibPython) { $ArgsList += @("--edlib-python", $EdlibPython) }
if ($Workers -gt 0) { $ArgsList += @("--workers", "$Workers") }
& $Python @ArgsList
exit $LASTEXITCODE
