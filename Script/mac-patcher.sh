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
	otool -L "$target" | tail -n +2 | awk '{print $1}' | while read -r dep; do
		if [[ "$dep" == /usr/lib/* || "$dep" == /System/* ]]; then
			continue
		fi
		depname=$(basename "$dep")
		if [[ -f "$dep" && ! -f "$LIB_DIR/$depname" ]]; then
			echo "Copying $dep to $LIB_DIR/"
			cp "$dep" "$LIB_DIR/"
			touch "$LIB_DIR/.new_dylib_copied"
		fi
	done
}

while :; do
	rm -f "$LIB_DIR/.new_dylib_copied"

	for binfile in "$BIN_DIR"/*; do
		[ -f "$binfile" ] || continue
		if ! file "$binfile" | grep -q 'Mach-O'; then
			continue
		fi
		copy_deps "$binfile"
	done

	for libfile in "$LIB_DIR"/*.dylib; do
		[ -f "$libfile" ] || continue
		if ! file "$libfile" | grep -q 'Mach-O'; then
			continue
		fi
		copy_deps "$libfile"
	done

	[ -f "$LIB_DIR/.new_dylib_copied" ] || break
done

for binfile in "$BIN_DIR"/*; do
	[ -f "$binfile" ] || continue
	if ! file "$binfile" | grep -q 'Mach-O'; then
		continue
	fi
	otool -L "$binfile" | tail -n +2 | awk '{print $1}' | while read -r dylib; do
		if [[ "$dylib" == /usr/lib/* || "$dylib" == /System/* ]]; then
			continue
		fi
		libname=$(basename "$dylib")
		install_name_tool -change "$dylib" "@executable_path/../lib/$libname" "$binfile"
	done
done

for libfile in "$LIB_DIR"/*.dylib; do
	[ -f "$libfile" ] || continue
	if ! file "$libfile" | grep -q 'Mach-O'; then
		continue
	fi
	echo "Patching $libfile"
	install_name_tool -id "@loader_path/$(basename "$libfile")" "$libfile"
	otool -L "$libfile" | tail -n +2 | awk '{print $1}' | while read -r dep; do
		if [[ "$dep" == /usr/lib/* || "$dep" == /System/* ]]; then
			continue
		fi
		depname=$(basename "$dep")
		if [[ -f "$LIB_DIR/$depname" ]]; then
			install_name_tool -change "$dep" "@loader_path/$depname" "$libfile"
		fi
	done
done

echo "Verification:"

for binfile in "$BIN_DIR"/*; do
	[ -f "$binfile" ] || continue
	if ! file "$binfile" | grep -q 'Mach-O'; then
		continue
	fi
	echo "== $binfile =="
	otool -L "$binfile"
done

for libfile in "$LIB_DIR"/*.dylib; do
	echo "== $libfile =="
	otool -L "$libfile"
done

if [ $# -gt 0 ]; then
	popd >/dev/null
fi

echo "Done!"
