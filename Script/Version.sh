#!/bin/bash

if [ $# -ne 3 ]; then
  echo "This script updates all version tags in the source tree."
  echo "Please invoke under the root folder of the project."
  echo "    Usage: $0 major minor patch"
  exit 0
fi

folder_name="."

if [ ! -f "suanPan.cpp" ]; then
  if [ -f "Version.sh" ]; then
    folder_name=".."
  else
    echo "Please invoke under the root folder of the project."
    exit 0
  fi
fi

# ci/cd
yml_file=$(find "$folder_name/.github/workflows" -name "*.yml")
for file in $yml_file; do
  if [ -f "$file" ]; then
    sed -i "s/suanPan-[0-9]\.[0-9]\.[0-9]/suanPan-$1\.$2\.$3/g" "$file"
  fi
done

# inno setup
file_name="$folder_name/Enhancement/suanPan.iss"
if [ -f $file_name ]; then
  sed -i "s/#define MyAppVersion \"[0-9]\.[0-9]\"/#define MyAppVersion \"$1\.$2\"/g" $file_name
fi

# chocolatey
file_name="$folder_name/Enhancement/suanpan.nuspec"
if [ -f $file_name ]; then
  sed -i "s/<version>[0-9]\.[0-9]<\/version>/<version>$1\.$2<\/version>/g" $file_name
fi

# msvc
file_name="$folder_name/Resource/suanPan.rc"
if [ -f $file_name ]; then
  sed -i "s/[0-9],[0-9],[0-9],0/$1,$2,$3,0/g" $file_name
  sed -i "s/[0-9]\.[0-9]\.[0-9]\.0/$1\.$2\.$3\.0/g" $file_name
fi

# application
file_name="$folder_name/Toolbox/argument.cpp"
if [ -f $file_name ]; then
  sed -i "s/constexpr auto SUANPAN_MAJOR = [0-9];/constexpr auto SUANPAN_MAJOR = $1;/g" $file_name
  sed -i "s/constexpr auto SUANPAN_MINOR = [0-9];/constexpr auto SUANPAN_MINOR = $2;/g" $file_name
  sed -i "s/constexpr auto SUANPAN_PATCH = [0-9];/constexpr auto SUANPAN_PATCH = $3;/g" $file_name
fi

# cpack
file_name="$folder_name/CMakeLists.txt"
if [ -f $file_name ]; then
  sed -i "s/set(CPACK_PACKAGE_VERSION \"[0-9]\.[0-9]\.[0-9]\")/set(CPACK_PACKAGE_VERSION \"$1\.$2\.$3\")/g" $file_name
fi

# snapcraft
file_name="$folder_name/snapcraft.yaml"
if [ -f $file_name ]; then
  sed -i "s/version: \"[0-9]\.[0-9]\"/version: \"$1\.$2\"/g" $file_name
fi
