#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
mkdir -p "${HOME}/.local/bin"
ln -sfn "${ROOT}/bin/marv-gen" "${HOME}/.local/bin/marv-gen"
chmod +x "${ROOT}/bin/marv-gen"
if ! echo "$PATH" | grep -q "${HOME}/.local/bin"; then
  echo 'export PATH="$HOME/.local/bin:$PATH"' >> "${HOME}/.bashrc"
  export PATH="$HOME/.local/bin:$PATH"
fi
echo "Installed: $(command -v marv-gen)"
echo "Pipeline:  ${ROOT}"
echo ""
echo "Activate Java 17 for Nextflow:"
echo "  conda activate nf-java17"
echo "Then run:"
echo "  marv-gen --help"
