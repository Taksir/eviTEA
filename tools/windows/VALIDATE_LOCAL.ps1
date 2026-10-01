[CmdletBinding(PositionalBinding = $false)]
param(
    [string]$Python = "python",
    [switch]$RunLiveMsbwt
)

$ErrorActionPreference = "Stop"
$EviTEARoot = (Resolve-Path (Join-Path $PSScriptRoot "..\..")).Path
Push-Location $EviTEARoot
try {
    if ($RunLiveMsbwt) {
        $env:EVITEA_RUN_MSBWT_LIVE = "1"
    }
    & $Python -m pytest
    if ($LASTEXITCODE -ne 0) {
        exit $LASTEXITCODE
    }
    & $Python -m build
    if ($LASTEXITCODE -ne 0) {
        exit $LASTEXITCODE
    }
}
finally {
    Pop-Location
}
