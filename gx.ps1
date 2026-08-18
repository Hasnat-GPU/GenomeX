# GenomeX launcher for Windows PowerShell: the cockpit stays here, the science
# runs in WSL. Override with $env:GENOMEX_WSL_DISTRO / $env:GENOMEX_ENV.
param([Parameter(ValueFromRemainingArguments = $true)][string[]]$Args)

$distro  = if ($env:GENOMEX_WSL_DISTRO) { $env:GENOMEX_WSL_DISTRO } else { "Ubuntu-24.04" }
$envName = if ($env:GENOMEX_ENV) { $env:GENOMEX_ENV } else { "gx" }

$repo = (& wsl.exe -d $distro wslpath -u "$PSScriptRoot").Trim()
$inner = "export MAMBA_ROOT_PREFIX=`$HOME/micromamba; cd '$repo' && " +
         "`$HOME/.local/bin/micromamba run -n '$envName' python -m genomex " + ($Args -join " ")
& wsl.exe -d $distro -- bash -lc $inner
