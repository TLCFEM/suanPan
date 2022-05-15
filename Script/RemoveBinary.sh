#!/bin/bash

shopt -s nullglob

if [ -f "suanPan.cpp" ]; then
  cd .
elif [ -f "../suanPan.cpp" ]; then
  cd ..
else
  echo "Please invoke this script in the root folder."
  exit 1
fi

path_array=(cmake-build*)
path_array+=(build*)

for path in "${path_array[@]}"; do
  path_file="${path}/suanPan"
  if [ -f "$path_file" ]; then
    rm "$path_file"
  fi
done

echo "Successfully removed all compiled binaries in following folders."
echo "${path_array[@]}"
