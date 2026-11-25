#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# Coverage.sh - Run examples and collect timing data
#
# Executes example files with the compiled suanPan binary, records execution
# times, and generates a summary report. Creates a Results.tar.gz archive.
# Example, smoke-command, and test-suite failures are all reported and make
# the script exit with a nonzero status.
#
# Usage: ./Coverage.sh <build_folder>
#   build_folder: Directory containing the compiled suanPan binary
#
# Example: ./Coverage.sh ./build
# -----------------------------------------------------------------------------
set -uo pipefail

if [ $# -eq 0 ]; then
  echo "Usage: $0 <build_folder>"
  exit 1
fi

cd "$1" || exit 1

if [ ! -f "suanPan" ]; then
  echo "This script must be executed in the folder contains suanPan."
  exit 1
fi

if [ ! -d "../Example" ]; then
  echo "This script must be executed in the sibling folder of Example folder."
  exit 1
fi

# GNU time is not installed everywhere. Prefer it when available, otherwise
# measure wall-clock time with nanosecond date output (or whole seconds).
find_time_binary() {
  local candidate
  for candidate in "$(type -P time 2>/dev/null)" /usr/bin/time /usr/local/bin/gtime /opt/homebrew/bin/gtime; do
    if [ -n "$candidate" ] && [ -x "$candidate" ] && "$candidate" -f '%e' /bin/true >/dev/null 2>&1; then
      echo "$candidate"
      return 0
    fi
  done
  return 1
}

TIME_BIN=$(find_time_binary) || TIME_BIN=""

mapfile -t files < <(find ../Example -name "*.supan" | grep -v "Misc")

declare -A timings

ROOTDIR=$(pwd)
TIME_TMP="$ROOTDIR/.coverage_time"
mkdir -p Results
trap 'rm -f "$TIME_TMP"' EXIT

LAST_STATUS=0
LAST_ELAPSED="n/a"

# Runs "$@" with its output sent to $1. Sets LAST_STATUS and LAST_ELAPSED.
run_timed() {
  local log_file=$1
  shift
  local start_ns start_sec end_ns end_sec

  if [ -n "$TIME_BIN" ]; then
    "$TIME_BIN" -f '%e' "$@" >"$log_file" 2>"$TIME_TMP"
    LAST_STATUS=$?
    LAST_ELAPSED=$(grep -E '^[0-9]+(\.[0-9]+)?$' "$TIME_TMP" | tail -n 1)
    [ -n "$LAST_ELAPSED" ] || LAST_ELAPSED=$(tail -n 1 "$TIME_TMP")
  else
    start_sec=$SECONDS
    start_ns=$(date +%s%N 2>/dev/null || echo 0)
    "$@" >"$log_file" 2>&1
    LAST_STATUS=$?
    end_ns=$(date +%s%N 2>/dev/null || echo 0)
    end_sec=$SECONDS
    if [[ "$start_ns" =~ ^[0-9]+$ && "$end_ns" =~ ^[0-9]+$ ]] && [ "$start_ns" -gt 0 ] && [ "$end_ns" -ge "$start_ns" ]; then
      LAST_ELAPSED=$(awk -v s="$start_ns" -v e="$end_ns" 'BEGIN { printf "%.2f", (e - s) / 1000000000 }')
    else
      LAST_ELAPSED=$((end_sec - start_sec))
    fi
  fi
}

example_failures=0
for file in "${files[@]}"; do
  echo "Processing $file ..."

  FILEPATH=$(realpath "$file")
  FILEBASE=$(basename "$FILEPATH")

  cd "$(dirname "$FILEPATH")" || exit 1

  # Truncate rather than append so repeated runs do not duplicate output.
  run_timed "$ROOTDIR/Results/$FILEBASE" "$ROOTDIR/suanPan" -nc -f "$FILEBASE"

  if [ "$LAST_STATUS" -ne 0 ]; then
    echo "Error processing $file (exit code $LAST_STATUS)." >&2
    example_failures=$((example_failures + 1))
  fi

  cd "$ROOTDIR" || exit 1

  timings["$file"]="$LAST_ELAPSED"
done

echo -e "\n================================= Processing Time Summary ================================="
printf "%-80s %10s\n" "File" "Time (s)"
printf "%-80s %10s\n" "--------------------------------------------------------------------------------" "----------"

for file in "${files[@]}"; do
  printf "%-80s %10s\n" "$file" "${timings["$file"]}"
done

# These -t invocations exercise CLI parsing only (test_mode() is currently a
# no-op); the actual Catch2 suite is run via -ct below.
smoke_failures=0
for args in "-t" "-vb -t" "-nc -vb -t" "-nc -vb -np -t" "-v"; do
  # shellcheck disable=SC2086
  if ! ./suanPan $args >/dev/null 2>&1; then
    echo "Smoke command failed: suanPan $args" >&2
    smoke_failures=$((smoke_failures + 1))
  fi
done

test_failures=0
echo -e "\nRunning Catch2 test suite (output saved to Results/tests.log) ..."
if ./suanPan -ct ~"Large Mixed Precision" ~"Large Sparse Solve Type" >"$ROOTDIR/Results/tests.log" 2>&1; then
  tail -n 3 "$ROOTDIR/Results/tests.log"
else
  test_failures=1
  echo "Catch2 test suite failed; tail of Results/tests.log:" >&2
  tail -n 40 "$ROOTDIR/Results/tests.log" >&2
fi

if ! tar -czf Results.tar.gz Results; then
  echo "Failed to create Results.tar.gz." >&2
  exit 1
fi

if [ "$example_failures" -ne 0 ] || [ "$smoke_failures" -ne 0 ] || [ "$test_failures" -ne 0 ]; then
  echo "Coverage run finished with $example_failures failed example(s), $smoke_failures failed smoke command(s), and $test_failures failed test suite(s)." >&2
  exit 1
fi

echo "Coverage run completed successfully."
