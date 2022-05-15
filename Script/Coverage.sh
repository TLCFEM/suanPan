#!/bin/bash

if [ $# -eq 0 ]; then
  echo "Usage: Coverage.sh <build_folder>"
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

cp ../Example/Solver/EZ .
cp ../Example/Material/Table .
cp ../Example/Material/C.txt .
cp ../Example/Material/T.txt .
cp ../Example/Material/CYCLE.txt .
./suanPan -t >/dev/null
./suanPan -v >/dev/null
./suanPan -ctest >/dev/null
./suanPan -c -f ../Example/Misc/Converter/TEST.inp >/dev/null
./suanPan -s -f ../Example/Misc/Converter/TEST.inp >/dev/null

for file in $files; do
  echo "Processing $file ..."
  ./suanPan -f "$file" >/dev/null
done

rm ../Example/Misc/Converter/TEST_out.inp
rm ../Example/Misc/Converter/TEST_out.supan
