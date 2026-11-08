#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# RemoveBinary.sh - Remove compiled binaries from build folders
#
# Removes suanPan binaries from cmake-build* and build* directories
# relative to the script location.
#
# Usage: ./RemoveBinary.sh
# -----------------------------------------------------------------------------

set -e

shopt -s nullglob

SCRIPT_DIR=$(realpath "$0")
SCRIPT_DIR=$(dirname "$SCRIPT_DIR")
SCRIPT_DIR=$(dirname "$SCRIPT_DIR")
cd "$SCRIPT_DIR" || exit

path_array=(cmake-build*)
path_array+=(build*)

for path in "${path_array[@]}"; do
  path_file="${path}/suanPan"
  if [ -f "$path_file" ]; then
    rm "$path_file"
  fi
done

echo "Successfully removed all compiled binaries in following folders."
echo "${path_array[@]}"
