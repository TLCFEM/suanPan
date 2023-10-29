#!/bin/bash

if [ -f "/opt/intel/oneapi/setvars.sh" ]; then
  source /opt/intel/oneapi/setvars.sh
fi

SCRIPT_DIR=$(realpath "$0")
SCRIPT_DIR=$(dirname "$SCRIPT_DIR")
SCRIPT_DIR=$(dirname "$SCRIPT_DIR")
cd "$SCRIPT_DIR" || exit

# shellcheck disable=SC2044
for folder in $(find . -maxdepth 1 -type d -name "cmake-build*"); do
  (
    echo "Compiling $folder"
    cd "$folder" || exit
    cmake --build . --target all -j "$(nproc)"
  )
done
