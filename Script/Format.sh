#!/usr/bin/env bash

ROOT_DIR=$(dirname "$(dirname "$(realpath "$0")")")

EXCLUDE_DIRS=("Include" "Toolbox/superlu*" "cmake*" "build*")
EXTENSIONS=("cpp" "c" "h" "hpp")

PRUNE_EXPR=""
for DIR in "${EXCLUDE_DIRS[@]}"; do
  PRUNE_EXPR="$PRUNE_EXPR -path '$ROOT_DIR/$DIR' -prune -o"
done

NAME_EXPR=""
for EXT in "${EXTENSIONS[@]}"; do
  if [ -z "$NAME_EXPR" ]; then
    NAME_EXPR="-name '*.$EXT'"
  else
    NAME_EXPR="$NAME_EXPR -o -name '*.$EXT'"
  fi
done

eval "find $ROOT_DIR $PRUNE_EXPR \( $NAME_EXPR \) -type f -exec cf -i {} +"
