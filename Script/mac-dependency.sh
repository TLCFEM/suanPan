#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# mac-dependency.sh
#
# Downloads and extracts prebuilt library dependencies for macOS runners,
# organizing them into the appropriate target directory based on the runner
# image name.
#
# Usage:
#   ./mac-dependency.sh <runner-image-name>
#
#   <runner-image-name>: The name of the GitHub Actions runner image
#                        (one of macos-14-large, macos-14-xlarge, macos-15-large, macos-15-xlarge).
#
# What it does:
#   1. Checks that exactly one argument is provided.
#   2. Sets the target directory:
#        - If the runner image name contains "xlarge", uses ../Libs/mac-arm.
#        - Otherwise, uses ../Libs/mac.
#   3. Removes and recreates the target directory.
#   4. Downloads and extracts:
#        - HDF5: Copies all .a files from the extracted lib directory.
#        - TBB:  Copies all .dylib files from the extracted tbb-install/lib directory.
#        - OpenBLAS: Copies all .dylib files.
#   5. Cleans up temporary directories.
#
# Requirements:
#   - wget
#   - tar
#
# -----------------------------------------------------------------------------

set -e

if [ "$#" -ne 1 ]; then
  echo "Usage: $0 <runner-image-name>"
  exit 1
fi

if ! command -v wget &>/dev/null; then
  echo "Error: wget is not installed. Please install wget and try again."
  exit 1
fi

RUNNER_IMAGE_NAME="$1"

if [[ "$RUNNER_IMAGE_NAME" == *"xlarge"* ]]; then
  TARGET_DIR="$(dirname "$0")/../Libs/mac-arm"
else
  TARGET_DIR="$(dirname "$0")/../Libs/mac"
fi

rm -rf "$TARGET_DIR"
mkdir -p "$TARGET_DIR"

TARBALL_URL="https://github.com/TLCFEM/prebuilds/releases/download/latest/HDF5-2.0.0-$RUNNER_IMAGE_NAME.tar.gz"
TMP_DIR="$(mktemp -d)"

wget -q -O "$TMP_DIR/archive.tar.gz" "$TARBALL_URL"
tar -xzf "$TMP_DIR/archive.tar.gz" -C "$TMP_DIR"

find "$TMP_DIR/lib" -name "*.a" -exec cp {} "$TARGET_DIR" \;

rm -rf "$TMP_DIR"

TARBALL_URL="https://github.com/TLCFEM/prebuilds/releases/download/latest/tbb-$RUNNER_IMAGE_NAME.tar.gz"
TMP_DIR="$(mktemp -d)"

wget -q -O "$TMP_DIR/archive.tar.gz" "$TARBALL_URL"
tar -xzf "$TMP_DIR/archive.tar.gz" -C "$TMP_DIR"

find "$TMP_DIR/tbb-install/lib" -name "*.dylib" -exec cp -P {} "$TARGET_DIR" \;

rm -rf "$TMP_DIR"

TARBALL_URL="https://github.com/TLCFEM/prebuilds/releases/download/latest/OpenBLAS-0.3.31-$RUNNER_IMAGE_NAME.tar.gz"
TMP_DIR="$(mktemp -d)"

wget -q -O "$TMP_DIR/archive.tar.gz" "$TARBALL_URL"
tar -xzf "$TMP_DIR/archive.tar.gz" -C "$TMP_DIR"

find "$TMP_DIR" -name "*.dylib" -exec cp -P {} "$TARGET_DIR" \;
# cp "$TMP_DIR/libopenblas.a" "$TARGET_DIR"

rm -rf "$TMP_DIR"

TARBALL_URL="https://github.com/TLCFEM/prebuilds/releases/download/latest/VTK-9.5.2-$RUNNER_IMAGE_NAME.tar.gz"
TMP_DIR="$(mktemp -d)"

wget -q -O "$TMP_DIR/archive.tar.gz" "$TARBALL_URL"
mkdir -p "$(dirname "$0")/../VTK"
tar -xzf "$TMP_DIR/archive.tar.gz" -C "$(dirname "$0")/../VTK"

rm -rf "$TMP_DIR"
