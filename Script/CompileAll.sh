#!/usr/bin/env bash

set -e

if [ -f "/opt/intel/oneapi/setvars.sh" ]; then
  source /opt/intel/oneapi/setvars.sh
fi

SCRIPT_DIR=$(realpath "$0")
SCRIPT_DIR=$(dirname "$SCRIPT_DIR")
SCRIPT_DIR=$(dirname "$SCRIPT_DIR")
cd "$SCRIPT_DIR" || exit

find . \
  \( -maxdepth 1 -mindepth 1 -type d -name "cmake-build*" \) -o \
  \( -maxdepth 2 -mindepth 2 -type d -path "./build/*" \) \
  -print0 |
while IFS= read -r -d '' folder; do
  (
    echo "Compiling $folder"
    cd "$folder" || exit
    cmake --build . --target all -j "$(nproc)"
  )
done
