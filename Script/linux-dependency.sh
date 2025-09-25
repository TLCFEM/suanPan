#!/usr/bin/env bash

# -----------------------------------------------------------------------------
# linux-dependency.sh
#
# Downloads and extracts prebuilt library dependencies for Linux runners,
# organizing them into the appropriate target directory based on the runner
# image name.
#
# Usage:
#   ./linux-dependency.sh <runner-image-name>
#
#   <runner-image-name>: The name of the GitHub Actions runner image
#                        (one of ubuntu-22.04, ubuntu-22.04-arm).
#
# What it does:
#   1. Checks that exactly one argument is provided.
#   2. Sets the target directory:
#        - If the runner image name contains "arm", uses ../Libs/linux-arm.
#        - Otherwise, uses ../Libs/linux.
#   3. Removes and recreates the target directory.
#   4. Downloads and extracts:
#        - HDF5: Copies all .a files from the extracted lib directory.
#        - TBB:  Copies all .so files from the extracted tbb-install/lib directory.
#        - OpenBLAS: Copies libopenblas.a file.
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

if ! command -v wget &> /dev/null; then
    echo "Error: wget is not installed. Please install wget and try again."
    exit 1
fi

RUNNER_IMAGE_NAME="$1"

if [[ "$RUNNER_IMAGE_NAME" == *"arm"* ]]; then
    TARGET_DIR="$(dirname "$0")/../Libs/linux-arm"
else
    TARGET_DIR="$(dirname "$0")/../Libs/linux"
fi

if [ -d "$TARGET_DIR" ]; then
    find "$TARGET_DIR" -mindepth 1 ! -name '*aocl*' -delete
else
    mkdir -p "$TARGET_DIR"
fi

TARBALL_URL="https://github.com/TLCFEM/prebuilds/releases/download/latest/HDF5-1.14.6-$RUNNER_IMAGE_NAME.tar.gz"
TMP_DIR="$(mktemp -d)"

wget -q -O "$TMP_DIR/archive.tar.gz" "$TARBALL_URL"
tar -xzf "$TMP_DIR/archive.tar.gz" -C "$TMP_DIR"

find "$TMP_DIR/lib" -name "*.a" -exec cp {} "$TARGET_DIR" \;

rm -rf "$TMP_DIR"

TARBALL_URL="https://github.com/TLCFEM/prebuilds/releases/download/latest/tbb-$RUNNER_IMAGE_NAME.tar.gz"
TMP_DIR="$(mktemp -d)"

wget -q -O "$TMP_DIR/archive.tar.gz" "$TARBALL_URL"
tar -xzf "$TMP_DIR/archive.tar.gz" -C "$TMP_DIR"

find "$TMP_DIR/tbb-install/lib" -name "lib*so*" -exec cp -P {} "$TARGET_DIR" \;

rm -rf "$TMP_DIR"

TARBALL_URL="https://github.com/TLCFEM/prebuilds/releases/download/latest/OpenBLAS-0.3.30-$RUNNER_IMAGE_NAME-32.tar.gz"
TMP_DIR="$(mktemp -d)"

wget -q -O "$TMP_DIR/archive.tar.gz" "$TARBALL_URL"
tar -xzf "$TMP_DIR/archive.tar.gz" -C "$TMP_DIR"

cp "$TMP_DIR/libopenblas.a" "$TARGET_DIR"

rm -rf "$TMP_DIR"

TARBALL_URL="https://github.com/TLCFEM/prebuilds/releases/download/latest/VTK-9.5.2-$RUNNER_IMAGE_NAME.tar.gz"
TMP_DIR="$(mktemp -d)"

wget -q -O "$TMP_DIR/archive.tar.gz" "$TARBALL_URL"
tar -xzf "$TMP_DIR/archive.tar.gz" -C .

rm -rf "$TMP_DIR"
