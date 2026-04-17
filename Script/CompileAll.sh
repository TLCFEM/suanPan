#!/usr/bin/env bash

set -e

if [ -f "/opt/intel/oneapi/setvars.sh" ]; then
  source /opt/intel/oneapi/setvars.sh
fi

SCRIPT_DIR=$(realpath "$0")
SCRIPT_DIR=$(dirname "$SCRIPT_DIR")
SCRIPT_DIR=$(dirname "$SCRIPT_DIR")
cd "$SCRIPT_DIR" || exit

folders=()

while IFS= read -r -d '' folder; do
  folders+=("$folder")
done < <(find . -mindepth 1 -maxdepth 1 -type d -name "cmake-build*" -print0)

while IFS= read -r -d '' folder; do
  folders+=("$folder")
done < <([ -d ./build ] && find ./build -mindepth 1 -maxdepth 1 -type d -print0)

for folder in "${folders[@]}"; do
  (
    echo "Compiling $folder"
    cd "$folder" || exit
    cmake --build . --target all -j "$(nproc)"
  )
done
