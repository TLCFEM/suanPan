#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# linux-patcher.sh
#
# Usage:
#   linux-patcher.sh <target_binary>
#
# Description:
#   - Copies all non-system shared library dependencies of <target_binary> into ../lib
#   - Patches the rpath of <target_binary> to use '$ORIGIN/../lib'
#   - Patches the rpath of each copied library in ../lib to use '$ORIGIN'
#
# Requirements:
#   - ldd, awk, realpath, patchelf must be installed and available in PATH
#
# Example:
#   ./linux-patcher.sh ./bin/suanPan
#
# -----------------------------------------------------------------------------

set -e

if [ $# -ne 1 ]; then
  echo "Usage: $0 <target_binary>"
  exit 1
fi

for TOOL in ldd awk realpath patchelf; do
  if ! command -v "$TOOL" >/dev/null 2>&1; then
    echo "Error: Required tool '$TOOL' is not installed or not in PATH."
    exit 1
  fi
done

TARGET="$1"
LIBDIR="$(realpath "$(dirname "$TARGET")/../lib")"

mkdir -p "$LIBDIR"

echo "Scanning $TARGET ..."

ldd "$TARGET" | awk '/=>/ { print $3 }' | while read -r DEPENDENCY; do
  if [ -z "$DEPENDENCY" ] || [ ! -f "$DEPENDENCY" ]; then continue; fi
  FILENAME=$(basename "$DEPENDENCY")
  case "$FILENAME" in
  libgcc_s* | libgfortran* | libgomp* | libm* | libquadmath* | libstdc++*) ;;
  *) continue ;;
  esac
  DEST="$LIBDIR/$FILENAME"
  if [ ! -f "$DEST" ]; then
    echo "Copying $FILENAME to $LIBDIR"
    cp "$DEPENDENCY" "$DEST"
  fi
done

echo "Patching $TARGET ..."

patchelf --set-rpath "$ORIGIN/../lib" "$TARGET"

for DEPENDENCY in "$LIBDIR"/*.so*; do
  if [ -f "$DEPENDENCY" ]; then
    echo "Patching $DEPENDENCY ..."
    patchelf --set-rpath "$ORIGIN" "$DEPENDENCY"
  fi
done

echo "Done."
