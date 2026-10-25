#!/usr/bin/env bash
# Applies the Scoary 1.6.16 SciPy compatibility patch inside a conda environment.
set -euo pipefail
ENV="${1:-vv_assoc}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PATCH="$HERE/scoary-1.6.16_scipy-binomtest.patch"

PKGDIR="$(conda run -n "$ENV" python -c 'import scoary, os; print(os.path.dirname(os.path.dirname(scoary.__file__)))')"
VERSION="$(conda run -n "$ENV" python -c 'from importlib.metadata import version; print(version("scoary"))')"

if [[ "$VERSION" != "1.6.16" ]]; then
    echo "Scoary $VERSION found in '$ENV'; this patch is for 1.6.16 only." >&2
    exit 1
fi
if grep -q "Compatibility shim (mics-vvulnificus)" "$PKGDIR/scoary/methods.py"; then
    echo "Scoary in '$ENV' is already patched."
    exit 0
fi
patch -p1 -d "$PKGDIR" < "$PATCH"
echo "Scoary in '$ENV' patched: $PKGDIR/scoary/methods.py"
