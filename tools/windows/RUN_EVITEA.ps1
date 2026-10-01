[CmdletBinding(PositionalBinding = $false)]
param(
    [string]$Python = "python",
    [Parameter(ValueFromRemainingArguments = $true)]
    [string[]]$EviTEAArgs
)

$ErrorActionPreference = "Stop"
& $Python -m evitea @EviTEAArgs
if ($LASTEXITCODE -ne 0) {
    exit $LASTEXITCODE
}
