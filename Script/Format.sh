#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# Format.sh - Format source code with clang-format or cf
#
# Formats C/C++/Header files across the project using clang-format or cf.
# Excludes third-party code and build directories.
#
# Usage: ./Format.sh
# -----------------------------------------------------------------------------
ROOT_DIR=$(dirname "$(dirname "$(realpath "$0")")")

EXCLUDE_PATTERNS=("Include" "Toolbox/superlu*" "cmake*" "build*")
EXTENSIONS=("cpp" "c" "h" "hpp")

if command -v cf >/dev/null 2>&1; then
  FORMATTER=(cf)
elif command -v clang-format >/dev/null 2>&1; then
  FORMATTER=(clang-format)
else
  echo "Error: neither 'cf' nor 'clang-format' is available." >&2
  exit 1
fi

PRUNE_ARGS=()
for PATTERN in "${EXCLUDE_PATTERNS[@]}"; do
  for DIR in "$ROOT_DIR"/$PATTERN; do
    [ -d "$DIR" ] && PRUNE_ARGS+=(-path "$DIR" -prune -o)
  done
done

NAME_ARGS=()
for EXT in "${EXTENSIONS[@]}"; do
  NAME_ARGS+=(-name "*.$EXT" -o)
done
unset 'NAME_ARGS[${#NAME_ARGS[@]}-1]'

find "$ROOT_DIR" "${PRUNE_ARGS[@]}" \( "${NAME_ARGS[@]}" \) -type f -print0 | xargs -0 -r -n1 -P "$(nproc)" "${FORMATTER[@]}" -i
