# Language builders (probe sets / codebooks) against the fixed reference
# tables (-RouteBV3): every stage in isolation, then the whole chain from the
# regenerated G1A catalogue.  Expected result LANGUAGE_BUILDER_EQUIVALENCE_PASS
# in <Cache>\LANGUAGE_EQUIVALENCE.json (final terminal d631f16a..., internal-only
# 44b98960...).  About 10 minutes; no sample data are read.
#
#   powershell -ExecutionPolicy Bypass -File .\tools\windows\VALIDATE_LANGUAGES.ps1
[CmdletBinding(PositionalBinding = $false)]
param(
    [string]$Python = "python",
    [Parameter(Mandatory = $true)][string]$Cache,
    [Parameter(Mandatory = $true)][string]$RouteBV3
)

$ErrorActionPreference = "Stop"
$env:PYTHONDONTWRITEBYTECODE = "1"
$EviTEARoot = (Resolve-Path (Join-Path $PSScriptRoot "..\..")).Path
$Script = Join-Path $EviTEARoot "validation\language_checks\run_language_equivalence.py"
& $Python $Script --route-b-v3 $RouteBV3 --cache $Cache
exit $LASTEXITCODE
