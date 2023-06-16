#!/bin/bash

if [ $# -eq 0 ]; then
  echo "Usage: Coverage.sh <build_folder> <log_file>"
  exit
fi

cd "$1" || exit

if [ ! -f "suanPan" ]; then
  echo "This script must be executed in the folder contains suanPan."
  exit
fi

if [ ! -d "../Example" ]; then
  echo "This script must be executed in the sibling folder of Example folder."
  exit
fi

files=$(find ../Example -name "*.supan")

if [ $# -eq 2 ]; then
  log_file=$2
else
  log_file="/dev/null"
fi

cp ../Example/Solver/EZ .
cp ../Example/Material/Table .
cp ../Example/Material/C.txt .
cp ../Example/Material/T.txt .
cp ../Example/Material/CYCLE.txt .
cp ../Example/Material/EHIST .
cp ../Example/Material/example .
cp ../Example/Section/HIST .

for file in $files; do
  echo "Processing $file ..."
  time ./suanPan -f "$file" >>"$log_file"
done

{
  ./suanPan -t
  ./suanPan -vb -t
  ./suanPan -nc -vb -t
  ./suanPan -nc -vb -np -t
  ./suanPan -v
  ./suanPan -ctest ~"Large Mixed Precision" ~"Large Sparse Solve Type"
  ./suanPan -c -f ../Example/Misc/Converter/TEST.inp
  ./suanPan -s -f ../Example/Misc/Converter/TEST.inp
} >>"$log_file"

rm ../Example/Misc/Converter/TEST_out.inp
rm ../Example/Misc/Converter/TEST_out.supan
