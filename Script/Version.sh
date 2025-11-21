#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# Version.sh
#
# Usage:
#   Version.sh <major> <minor> <patch>
#
# Description:
#   - Updates all version tags in the source tree to the specified major, minor, and patch numbers.
#   - Modifies CI/CD workflow files, Chocolatey spec, MSVC resource, application source, CPack config,
#     Snapcraft config, and citation files.
#   - Should be run from the root folder of the project.
#
# Requirements:
#   - bash, realpath, sed, find must be available in PATH.
#
# Example:
#   ./Version.sh 2 1 0
#
# -----------------------------------------------------------------------------

set -e

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
    sed -i "s/suanPan-\d\.\d\.\d/suanPan-$1\.$2\.$3/g" "$file"
    sed -i "s/suanpan.\d\.\d\.\d/suanpan.$1\.$2\.$3/g" "$file"
  fi
done

# chocolatey
file_name="Enhancement/suanPan.nuspec"
if [ -f $file_name ]; then
  sed -i "s/<version>\d\.\d\.\d<\/version>/<version>$1\.$2\.$3<\/version>/g" $file_name
fi

# msvc
file_name="Resource/suanPan.rc"
if [ -f $file_name ]; then
  sed -i "s/\d,\d,\d,0/$1,$2,$3,0/g" $file_name
  sed -i "s/\d\.\d\.\d\.0/$1\.$2\.$3\.0/g" $file_name
fi

# application
file_name="Toolbox/argument.cpp"
if [ -f $file_name ]; then
  sed -i "s/constexpr auto SUANPAN_MAJOR = \d;/constexpr auto SUANPAN_MAJOR = $1;/g" $file_name
  sed -i "s/constexpr auto SUANPAN_MINOR = \d;/constexpr auto SUANPAN_MINOR = $2;/g" $file_name
  sed -i "s/constexpr auto SUANPAN_PATCH = \d;/constexpr auto SUANPAN_PATCH = $3;/g" $file_name
fi

# cpack
file_name="CMakeLists.txt"
if [ -f $file_name ]; then
  sed -i "s/set(CPACK_PACKAGE_VERSION \"\d\.\d\.\d\")/set(CPACK_PACKAGE_VERSION \"$1\.$2\.$3\")/g" $file_name
fi

# snapcraft
file_name="snapcraft.yaml"
if [ -f $file_name ]; then
  sed -i "s/version: \"\d\.\d\.\d\"/version: \"$1\.$2\.$3\"/g" $file_name
fi

# citation
file_name="CITATION.cff"
if [ -f $file_name ]; then
  sed -i "s/version: \"\d\.\d\.\d\"/version: \"$1\.$2\.$3\"/g" $file_name
fi

# rpm spec
file_name="suanPan.spec.rpkg"
if [ -f $file_name ]; then
  sed -i "s/\d\.\d\.\d/$1\.$2\.$3/g" $file_name
fi
