#!/bin/bash

git --version >/dev/null 2>&1

if [ ! $? ]; then
  echo "Git is not installed, not setting revision."
  exit
fi

if [ ! -d ".git" ]; then
  echo "Not a git repository, can't set revision."
  exit
fi

if [ -f "./suanPan.cpp" ]; then
  file_path="./Toolbox/revision.h"
elif [ -f "../suanPan.cpp" ]; then
  file_path="../Toolbox/revision.h"
else
  echo "Error: Cannot locate the root directory."
  exit
fi

git_rev=$(git rev-parse --short=8 HEAD)

echo "constexpr auto SUANPAN_REVISION = \"$git_rev\";" >$file_path

echo "Revision tag set to $git_rev"
