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

for file in $files; do
  echo "Processing $file ..."
  ./suanPan -f "$file" >>$log_file
done

{
  ./suanPan -t
  ./suanPan -v
  ./suanPan -ctest
  ./suanPan -c -f ../Example/Misc/Converter/TEST.inp
  ./suanPan -s -f ../Example/Misc/Converter/TEST.inp
} >>$log_file

rm ../Example/Misc/Converter/TEST_out.inp
rm ../Example/Misc/Converter/TEST_out.supan
