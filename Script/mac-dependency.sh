#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# mac-dependency.sh - Download macOS dependencies for CI/CD
#
# Downloads prebuilt library dependencies for macOS GitHub Actions runners
# and organizes them by runner image type (standard vs xlarge).
#
# Usage: ./mac-dependency.sh <runner-image-name>
#   runner-image-name: GitHub Actions runner image (macos-15-large, macos-15-xlarge, etc.)
#
# Downloads: HDF5, TBB, OpenBLAS, VTK
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
  TOKEN="mac-arm"
else
  TOKEN="mac"
fi

TARBALL_URL="https://github.com/TLCFEM/prebuilds/releases/download/latest/HDF5-2.2.0-$RUNNER_IMAGE_NAME.tar.gz"
TMP_DIR="$(mktemp -d)"

wget -q -O "$TMP_DIR/archive.tar.gz" "$TARBALL_URL"
tar -xzf "$TMP_DIR/archive.tar.gz" -C "$TMP_DIR"

TARGET_DIR="$(dirname "$0")/../Libs/$TOKEN"

rm -rf "$TARGET_DIR"
mkdir -p "$TARGET_DIR"

find "$TMP_DIR/lib" -name "*.a" -exec cp {} "$TARGET_DIR" \;

INCLUDE_DIR="$(dirname "$0")/../Include/hdf"

rm -rf "${INCLUDE_DIR:?}/*" "$INCLUDE_DIR-$TOKEN/*"

find "$TMP_DIR/include" -name "*.h" -exec cp {} "$INCLUDE_DIR" \;

rm -rf "$TMP_DIR"

TARBALL_URL="https://github.com/TLCFEM/prebuilds/releases/download/latest/tbb-v2023.1.0-$RUNNER_IMAGE_NAME.tar.gz"
TMP_DIR="$(mktemp -d)"

wget -q -O "$TMP_DIR/archive.tar.gz" "$TARBALL_URL"
tar -xzf "$TMP_DIR/archive.tar.gz" -C "$TMP_DIR"

rm -rf "$TARGET_DIR/libtbb*"
rm -rf "$(dirname "$0")/../Include/tbb/*"
find "$TMP_DIR/tbb-install/lib" -name "*.dylib" -exec cp -P {} "$TARGET_DIR" \;
cp -a "$TMP_DIR/tbb-install/include/." "$(dirname "$0")/../Include/tbb"

rm -rf "$TMP_DIR"

TARBALL_URL="https://github.com/TLCFEM/prebuilds/releases/download/latest/OpenBLAS-0.3.34-$RUNNER_IMAGE_NAME.tar.gz"
TMP_DIR="$(mktemp -d)"

wget -q -O "$TMP_DIR/archive.tar.gz" "$TARBALL_URL"
tar -xzf "$TMP_DIR/archive.tar.gz" -C "$TMP_DIR"

rm -rf "$TARGET_DIR/libopenblas*"
find "$TMP_DIR" -name "*.dylib" -exec cp -P {} "$TARGET_DIR" \;
# cp "$TMP_DIR/libopenblas.a" "$TARGET_DIR"

rm -rf "$TMP_DIR"

TARBALL_URL="https://github.com/TLCFEM/prebuilds/releases/download/latest/VTK-9.7.1-$RUNNER_IMAGE_NAME.tar.gz"
TMP_DIR="$(mktemp -d)"

wget -q -O "$TMP_DIR/archive.tar.gz" "$TARBALL_URL"
mkdir -p "$(dirname "$0")/../VTK"
tar -xzf "$TMP_DIR/archive.tar.gz" -C "$(dirname "$0")/../VTK"

rm -rf "$TMP_DIR"
