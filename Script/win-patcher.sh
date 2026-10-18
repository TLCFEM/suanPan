#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# win-patcher.sh - Windows DLL cleanup utility
#
# Prunes unnecessary DLLs from a Windows executable's directory by analyzing
# its import table and removing any DLLs not actually needed.
#
# Usage: ./win-patcher.sh <path/to/app.exe>
#
# Example:
#   ./win-patcher.sh "C:/path/to/app/bin/suanPan.exe"
# -----------------------------------------------------------------------------

set -euo pipefail
shopt -s nullglob

usage() {
  echo "Usage: $0 <path/to/app.exe>"
  exit 2
}

[[ $# -eq 1 ]] || usage
TARGET=$1
[[ -f "$TARGET" ]] || {
  echo "Not a file: $TARGET"
  exit 1
}

DIR=$(cd "$(dirname "$TARGET")" && pwd)
BASE=$(basename "$TARGET")

have() { command -v "$1" >/dev/null 2>&1; }

get_deps_with_objdump() {
  objdump -p "$TARGET" 2>/dev/null | awk '/DLL Name:/ {print tolower($3)}'
}

get_deps_with_dumpbin() {
  dumpbin /nologo /dependents "$TARGET" 2>/dev/null |
    awk '
      /Image has the following dependencies/ {insec=1; next}
      insec && NF {print tolower($1)}
      insec && /^$/ {insec=0}
    '
}

get_deps() {
  if have objdump; then
    get_deps_with_objdump
  elif have dumpbin; then
    get_deps_with_dumpbin
  else
    echo "Error: need objdump (MinGW/binutils) or dumpbin (VS tools) in PATH." >&2
    exit 1
  fi
}

declare -A NEED=()
while IFS= read -r d; do
  [[ -n "${d:-}" ]] && NEED["$d"]=1
done < <(get_deps)

echo "Target: $BASE"
echo "Directory: $DIR"
echo "Dependencies found:"
for k in "${!NEED[@]}"; do echo "  - $k"; done | sort

for dll in "$DIR"/*.dll "$DIR"/*.DLL; do
  [[ -f "$dll" ]] || continue
  name=$(basename "$dll")
  lower=$(printf '%s' "$name" | tr '[:upper:]' '[:lower:]')
  if [[ -n "${NEED[$lower]+y}" ]]; then
    continue
  fi
  echo "Remove: $name"
  rm -f -- "$dll"
done

echo "Done."
