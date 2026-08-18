#!/usr/bin/env bash
# GenomeX launcher. Runs the pipeline in the Linux environment that holds the
# bioinformatics tools, whether you invoke it from Linux, WSL, or Git Bash on
# Windows. Override with GENOMEX_WSL_DISTRO / GENOMEX_ENV if yours differ.
set -euo pipefail

ENV_NAME="${GENOMEX_ENV:-gx}"
DISTRO="${GENOMEX_WSL_DISTRO:-Ubuntu-24.04}"
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

micromamba_bin() {
  command -v micromamba || echo "$HOME/.local/bin/micromamba"
}

if [ -r /proc/version ] && ! grep -qi microsoft /proc/version 2>/dev/null || [ "$(uname -s)" = "Linux" ]; then
  # Already on Linux (native or inside WSL): run here.
  cd "$here"
  exec "$(micromamba_bin)" run -n "$ENV_NAME" python -m genomex "$@"
fi

# Windows shell: translate this directory into a WSL path and hop over.
win_path="$(cygpath -w "$here" 2>/dev/null || echo "$here")"
wsl_path="$(wsl.exe -d "$DISTRO" wslpath -u "$win_path" | tr -d '\r')"
exec wsl.exe -d "$DISTRO" -- bash -lc \
  "export MAMBA_ROOT_PREFIX=\$HOME/micromamba; cd '$wsl_path' && \
   \$HOME/.local/bin/micromamba run -n '$ENV_NAME' python -m genomex $*"
