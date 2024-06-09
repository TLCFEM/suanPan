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

SCRIPT_DIR=$(realpath "$0")
SCRIPT_DIR=$(dirname "$SCRIPT_DIR")
SCRIPT_DIR=$(dirname "$SCRIPT_DIR")
cd "$SCRIPT_DIR" || exit

# sleep random time to avoid git conflict in parallel execution
sleep $((RANDOM % 4)).$(((RANDOM % 10) + 1))s

git_rev=$(git rev-parse --short=8 HEAD)

file_path="./Toolbox/revision.h"

echo "constexpr auto SUANPAN_REVISION = \"$git_rev\";" >$file_path

echo "Revision tag set to $git_rev"
