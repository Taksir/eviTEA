# Validation-only probe: run the ISO-1 paired-end anchor support transfer procedure
# unchanged on AKA-017 and compare with the fixed paired-end anchor support candidates of AKA-017.
# Resumable: rerun the same command after an interruption.
#
#   powershell -ExecutionPolicy Bypass -File .\tools\windows\PROBE_ISO_PROCEDURE_AKA017.ps1
[CmdletBinding(PositionalBinding = $false)]
param(
    [string]$Python = "python",
    [Parameter(Mandatory = $true)][string]$EdlibPython,
    [Parameter(Mandatory = $true)][string]$Work,
    [Parameter(Mandatory = $true)][string]$RouteBV3,
    [int]$Threads = 8,
    [ValidateSet("iso", "aka-historical")]
    [string]$HostReadNames = "iso"
)

$ErrorActionPreference = "Stop"
$env:PYTHONDONTWRITEBYTECODE = "1"
$EviTEARoot = (Resolve-Path (Join-Path $PSScriptRoot "..\..")).Path
$Script = Join-Path $EviTEARoot "validation\v3_layer_checks\run_iso_procedure_on_aka017.py"
& $Python $Script --route-b-v3 $RouteBV3 --work $Work --edlib-python $EdlibPython --threads $Threads --host-read-names $HostReadNames
exit $LASTEXITCODE
