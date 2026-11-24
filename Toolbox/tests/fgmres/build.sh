#!/usr/bin/env bash
# Build and run the small standalone FGMRES test (does not build suanPan).
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
root="$(cd "$here/../../.." && pwd)"

if [[ -z "${MKL_LIB_DIR:-}" ]]; then
    for candidate in /opt/intel/oneapi/mkl/latest/lib/intel64 "${MKLROOT:-/nonexistent}/lib/intel64" /usr/lib/x86_64-linux-gnu; do
        if compgen -G "$candidate/libmkl_rt.so*" > /dev/null; then
            MKL_LIB_DIR="$candidate"
            break
        fi
    done
fi

# oneMKL can also be installed with `pip install mkl` in a virtual environment.
if [[ -z "${MKL_LIB_DIR:-}" && -n "${VIRTUAL_ENV:-}" ]]; then
    MKL_LIB_DIR="$VIRTUAL_ENV/lib"
fi

if [[ -z "${MKL_LIB_DIR:-}" ]] || ! compgen -G "$MKL_LIB_DIR/libmkl_rt.so*" > /dev/null; then
    echo "oneMKL runtime not found. Set MKL_LIB_DIR to the directory containing libmkl_rt.so." >&2
    exit 1
fi

mkl_lib="$(basename "$(find "$MKL_LIB_DIR" -maxdepth 1 -name 'libmkl_rt.so*' | head -n 1)")"
out="${TMPDIR:-/tmp}/fgmres_test"

g++ -std=c++20 -O1 -Wall -Wextra \
    -I"$here/shims" -I"$root" -I"$root/Include" \
    "$here/test_fgmres.cpp" \
    -L"$MKL_LIB_DIR" "-l:$mkl_lib" -Wl,-rpath,"$MKL_LIB_DIR" \
    -o "$out"

"$out"
