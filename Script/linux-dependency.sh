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
#                        (one of ubuntu-24.04, ubuntu-24.04-arm).
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
#        - mimalloc: Copies all .so files from the extracted lib directory.
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

if ! command -v tar &>/dev/null; then
  echo "Error: tar is not installed. Please install tar and try again."
  exit 1
fi

if command -v wget &>/dev/null; then
  DOWNLOAD_TOOL="wget"
elif command -v curl &>/dev/null; then
  DOWNLOAD_TOOL="curl"
else
  echo "Error: neither wget nor curl is installed. Please install one and try again."
  exit 1
fi

fetch_archive() {
  local DOWNLOAD_TO="$1"
  local DOWNLOAD_FROM="$2"

  if [ "$DOWNLOAD_TOOL" = "wget" ]; then
    wget -q -O "$DOWNLOAD_TO" "$DOWNLOAD_FROM"
  else
    curl -fsSL -o "$DOWNLOAD_TO" "$DOWNLOAD_FROM"
  fi
}

RUNNER_IMAGE_NAME="$1"

if [[ "$RUNNER_IMAGE_NAME" == *"arm"* ]]; then
  TARGET_DIR="$(dirname "$0")/../Libs/linux-arm"
else
  TARGET_DIR="$(dirname "$0")/../Libs/linux"
fi

TARGET_DIR="$(realpath "$TARGET_DIR")"

if [ -d "$TARGET_DIR" ]; then
  find "$TARGET_DIR" -mindepth 1 ! -name "libflame.a" ! -name "libblis-mt.a" ! -name "libaoclutils.a" -delete
else
  mkdir -p "$TARGET_DIR"
fi

TARBALL_URL="https://github.com/TLCFEM/prebuilds/releases/download/latest/HDF5-2.1.1-$RUNNER_IMAGE_NAME.tar.gz"
TMP_DIR="$(mktemp -d)"

fetch_archive "$TMP_DIR/archive.tar.gz" "$TARBALL_URL"
tar -xzf "$TMP_DIR/archive.tar.gz" -C "$TMP_DIR"

find "$TMP_DIR/lib" -name "*.a" -exec cp {} "$TARGET_DIR" \;

rm -rf "$TMP_DIR"

TARBALL_URL="https://github.com/TLCFEM/prebuilds/releases/download/latest/tbb-$RUNNER_IMAGE_NAME.tar.gz"
TMP_DIR="$(mktemp -d)"

fetch_archive "$TMP_DIR/archive.tar.gz" "$TARBALL_URL"
tar -xzf "$TMP_DIR/archive.tar.gz" -C "$TMP_DIR"

find "$TMP_DIR/tbb-install/lib" -name "lib*so*" -exec cp -P {} "$TARGET_DIR" \;

rm -rf "$TMP_DIR"

TARBALL_URL="https://github.com/TLCFEM/prebuilds/releases/download/latest/mimalloc-$RUNNER_IMAGE_NAME.tar.gz"
TMP_DIR="$(mktemp -d)"

fetch_archive "$TMP_DIR/archive.tar.gz" "$TARBALL_URL"
tar -xzf "$TMP_DIR/archive.tar.gz" -C "$TMP_DIR"

find "$TMP_DIR/lib" -name "lib*so*" -exec cp -P {} "$TARGET_DIR" \;

rm -rf "$TMP_DIR"

TARBALL_URL="https://github.com/TLCFEM/prebuilds/releases/download/latest/OpenBLAS-0.3.33-$RUNNER_IMAGE_NAME-32.tar.gz"
TMP_DIR="$(mktemp -d)"

fetch_archive "$TMP_DIR/archive.tar.gz" "$TARBALL_URL"
tar -xzf "$TMP_DIR/archive.tar.gz" -C "$TMP_DIR"

cp "$TMP_DIR/libopenblas.a" "$TARGET_DIR"

rm -rf "$TMP_DIR"

TARBALL_URL="https://github.com/TLCFEM/prebuilds/releases/download/latest/VTK-9.6.1-$RUNNER_IMAGE_NAME.tar.gz"
TMP_DIR="$(mktemp -d)"

fetch_archive "$TMP_DIR/archive.tar.gz" "$TARBALL_URL"
mkdir -p "$(dirname "$0")/../VTK"
tar -xzf "$TMP_DIR/archive.tar.gz" -C "$(dirname "$0")/../VTK"

rm -rf "$TMP_DIR"

echo "Dependencies for $RUNNER_IMAGE_NAME have been downloaded and extracted to $TARGET_DIR."
