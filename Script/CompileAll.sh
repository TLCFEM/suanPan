#!/bin/bash

if [ -f "/opt/intel/oneapi/setvars.sh" ]; then
  source /opt/intel/oneapi/setvars.sh
fi

# shellcheck disable=SC2044
for folder in $(find . -maxdepth 1 -type d -name "cmake-build*"); do
  (
    echo "Compiling $folder"
    cd "$folder" || exit
    make -j"$(nproc)"
  )
done
