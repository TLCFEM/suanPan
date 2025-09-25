#!/bin/bash

set -e

if [ $# -gt 0 ]; then
  pushd "$1" >/dev/null
fi

BIN_DIR="bin"
LIB_DIR="lib"

echo "Bundling all required dylibs for binaries and existing dylibs..."

copy_deps() {
  local target="$1"
  local hostdir=$(dirname "$target")
  otool -L "$target" | tail -n +2 | awk '{print $1}' | while read -r dep; do
    [[ "$dep" == /usr/lib/* || "$dep" == /System/* || "$dep" == "$target" ]] && continue
    depname=$(basename "$dep")
    if [[ "$dep" == @rpath/* ]]; then
      abspath="$hostdir/${dep#@rpath/}"
    else
      abspath="$dep"
    fi
    if [[ -f "$abspath" && ! -f "$LIB_DIR/$depname" ]]; then
      echo "Copying $abspath to $LIB_DIR/"
      cp "$abspath" "$LIB_DIR/"
      copy_deps "$abspath"
      touch "$LIB_DIR/.new_dylib_copied"
    fi
  done
}

while :; do
  rm -f "$LIB_DIR/.new_dylib_copied"

  for target in "$BIN_DIR"/* "$LIB_DIR"/*.dylib; do
    [[ -f "$target" ]] && file "$target" | grep -q 'Mach-O' && copy_deps "$target"
  done

  [ -f "$LIB_DIR/.new_dylib_copied" ] || break
done

for target in "$BIN_DIR"/*; do
  [[ -f "$target" ]] && file "$target" | grep -q 'Mach-O' || continue
  echo "Patching $target"
  otool -L "$target" | tail -n +2 | awk '{print $1}' | while read -r dep; do
    [[ "$dep" == /usr/lib/* || "$dep" == /System/* ]] && continue
    install_name_tool -change "$dep" "@executable_path/../lib/$(basename "$dep")" "$target"
  done
done

for target in "$LIB_DIR"/*.dylib; do
  [[ -f "$target" ]] && file "$target" | grep -q 'Mach-O' || continue
  echo "Patching $target"
  install_name_tool -id "@loader_path/$(basename "$target")" "$target"
  otool -L "$target" | tail -n +2 | awk '{print $1}' | while read -r dep; do
    [[ "$dep" == /usr/lib/* || "$dep" == /System/* ]] && continue
    depname=$(basename "$dep")
    if [[ -f "$LIB_DIR/$depname" ]]; then
      install_name_tool -change "$dep" "@loader_path/$depname" "$target"
    fi
  done
done

echo "Verification:"

for target in "$BIN_DIR"/* "$LIB_DIR"/*.dylib; do
  [[ -f "$target" ]] && file "$target" | grep -q 'Mach-O' || continue
  echo "== $target =="
  otool -L "$target"
done

if [ $# -gt 0 ]; then
  popd >/dev/null
fi

echo "Done!"
