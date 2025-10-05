#!/usr/bin/env bash

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
  TARGET="$LIBDIR/$FILENAME"
  if [ ! -f "$TARGET" ]; then
    echo "Copying $FILENAME to $LIBDIR"
    cp "$DEPENDENCY" "$TARGET"
  fi
done

echo "Patching $TARGET ..."

patchelf --set-rpath '$ORIGIN/../lib' "$TARGET"

echo "Done."
