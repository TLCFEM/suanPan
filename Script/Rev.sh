#!/bin/bash

if [ -f "./suanPan.cpp" ]; then
  file_path="./Toolbox/revision.h"
elif [ -f "../suanPan.cpp" ]; then
  file_path="../Toolbox/revision.h"
else
  echo "Error: Cannot locate thr root directory."
  exit
fi

git_rev=$(git rev-parse --short=8 HEAD)

echo "constexpr auto SUANPAN_REVISION = \"$git_rev\";" > $file_path

echo "Revision tag set to $git_rev"
