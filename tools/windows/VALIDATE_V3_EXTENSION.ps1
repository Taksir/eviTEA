# Run the extension layers (paired-end anchor support, full memory, internal 21-mer) on
# AKA017 and ISO1 from the prepared msBWTs; compare every stage and the
# bilateral / asymmetric / full-memory / internal 21-mer layers with the fixed
# reference results (-RouteBV2, -RouteBV3).  Resumable (reuses finished stages in -Work).
#
#   powershell -ExecutionPolicy Bypass -File .\tools\windows\VALIDATE_V3_EXTENSION.ps1
[CmdletBinding(PositionalBinding = $false)]
param(
    [string]$Python = "python",
    [string]$EdlibPython = "",
    [Parameter(Mandatory = $true)][string]$Work,
    [Parameter(Mandatory = $true)][string]$RouteBV2,
    [Parameter(Mandatory = $true)][string]$RouteBV3,
    [int]$Threads = 8
)

$ErrorActionPreference = "Stop"
$env:PYTHONDONTWRITEBYTECODE = "1"
$EviTEARoot = (Resolve-Path (Join-Path $PSScriptRoot "..\..")).Path
$Script = Join-Path $EviTEARoot "validation\v3_layer_checks\run_v3_extension_equivalence.py"
$EdlibArgs = @()
if ($EdlibPython) { $EdlibArgs = @("--edlib-python", $EdlibPython) }
& $Python $Script --route-b-v2 $RouteBV2 --route-b-v3 $RouteBV3 --work $Work --threads $Threads @EdlibArgs
exit $LASTEXITCODE
