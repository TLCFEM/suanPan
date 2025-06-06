#!/bin/bash

if [ $# -ne 3 ]; then
  echo "This script updates all version tags in the source tree."
  echo "Please invoke under the root folder of the project."
  echo "    Usage: $0 major minor patch"
  exit 0
fi

SCRIPT_DIR=$(realpath "$0")
SCRIPT_DIR=$(dirname "$SCRIPT_DIR")
SCRIPT_DIR=$(dirname "$SCRIPT_DIR")
cd "$SCRIPT_DIR" || exit

# ci/cd
yml_file=$(find ".github/workflows" -name "*.yml")
for file in $yml_file; do
  if [ -f "$file" ]; then
    sed -i "s/suanPan-[0-9]\.[0-9]\.[0-9]/suanPan-$1\.$2\.$3/g" "$file"
    sed -i "s/suanpan.[0-9]\.[0-9]\.[0-9]/suanpan.$1\.$2\.$3/g" "$file"
  fi
done

# chocolatey
file_name="Enhancement/suanPan.nuspec"
if [ -f $file_name ]; then
  sed -i "s/<version>[0-9]\.[0-9]\.[0-9]<\/version>/<version>$1\.$2\.$3<\/version>/g" $file_name
fi

# msvc
file_name="Resource/suanPan.rc"
if [ -f $file_name ]; then
  sed -i "s/[0-9],[0-9],[0-9],0/$1,$2,$3,0/g" $file_name
  sed -i "s/[0-9]\.[0-9]\.[0-9]\.0/$1\.$2\.$3\.0/g" $file_name
fi

# application
file_name="Toolbox/argument.cpp"
if [ -f $file_name ]; then
  sed -i "s/constexpr auto SUANPAN_MAJOR = [0-9];/constexpr auto SUANPAN_MAJOR = $1;/g" $file_name
  sed -i "s/constexpr auto SUANPAN_MINOR = [0-9];/constexpr auto SUANPAN_MINOR = $2;/g" $file_name
  sed -i "s/constexpr auto SUANPAN_PATCH = [0-9];/constexpr auto SUANPAN_PATCH = $3;/g" $file_name
fi

# cpack
file_name="CMakeLists.txt"
if [ -f $file_name ]; then
  sed -i "s/set(CPACK_PACKAGE_VERSION \"[0-9]\.[0-9]\.[0-9]\")/set(CPACK_PACKAGE_VERSION \"$1\.$2\.$3\")/g" $file_name
fi

# snapcraft
file_name="snapcraft.yaml"
if [ -f $file_name ]; then
  sed -i "s/version: \"[0-9]\.[0-9]\.[0-9]\"/version: \"$1\.$2\.$3\"/g" $file_name
fi

# citation
file_name="CITATION.cff"
if [ -f $file_name ]; then
  sed -i "s/version: \"[0-9]\.[0-9]\.[0-9]\"/version: \"$1\.$2\.$3\"/g" $file_name
fi
