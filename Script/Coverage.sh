#!/bin/bash

if [ $# -eq 0 ]; then
  echo "Usage: $0 <build_folder> <log_file>"
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

mapfile -t files < <(find ../Example -name "*.supan")

if [ $# -eq 2 ]; then
  log_file=$(realpath "$2")
else
  log_file="/dev/null"
fi

declare -A timings

true >"$log_file"

current_dir=$(pwd)

for file in "${files[@]}"; do
  echo "Processing $file ..." | tee -a "$log_file"

  file_path=$(realpath "$file")
  file_base=$(basename "$file_path")

  cd $(dirname "$file_path") || exit 1

  exec 3>&1 4>&2
  TIME_OUTPUT=$({ /usr/bin/time -f "%e" "$current_dir/suanPan" -f "$file_base" >>"$log_file"; } 2>&1 1>&3)
  exec 3>&- 4>&-

  cd "$current_dir" || exit 1

  # shellcheck disable=SC2181
  if [ $? -ne 0 ]; then
    echo "Error processing $file." >&2
    exit 1
  fi

  timings["$file"]="$TIME_OUTPUT"
done

echo -e "\n================================= Processing Time Summary ================================="
printf "%-80s %10s\n" "File" "Time (s)"
printf "%-80s %10s\n" "--------------------------------------------------------------------------------" "----------"

for file in "${files[@]}"; do
  printf "%-80s %10s\n" "$file" "${timings["$file"]}"
done

{
  ./suanPan -t
  ./suanPan -vb -t
  ./suanPan -nc -vb -t
  ./suanPan -nc -vb -np -t
  ./suanPan -v
  ./suanPan -ct ~"Large Mixed Precision" ~"Large Sparse Solve Type"
  ./suanPan -c -f ../Example/Misc/Converter/TEST.inp
  ./suanPan -s -f ../Example/Misc/Converter/TEST.inp
} >>"$log_file"

rm ../Example/Misc/Converter/TEST_out.inp
rm ../Example/Misc/Converter/TEST_out.supan
