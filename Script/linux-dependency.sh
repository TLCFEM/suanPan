#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# linux-dependency.sh - Download Linux dependencies for CI/CD
#
# Downloads prebuilt library dependencies for Linux GitHub Actions runners
# and organizes them by architecture (x86_64 vs ARM).
#
# Usage: ./linux-dependency.sh <runner-image-name>
#   runner-image-name: GitHub Actions runner image (ubuntu-24.04, ubuntu-24.04-arm, etc.)
#
# Downloads: HDF5, TBB, mimalloc, OpenBLAS, VTK
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
  TOKEN="linux-arm"
else
  TOKEN="linux"
fi

TARGET_DIR="$(dirname "$0")/../Libs/$TOKEN"
TARGET_DIR="$(realpath "$TARGET_DIR")"

if [ -d "$TARGET_DIR" ]; then
  find "$TARGET_DIR" -mindepth 1 ! -name "libflame.a" ! -name "libblis-mt.a" ! -name "libaoclutils.a" -delete
else
  mkdir -p "$TARGET_DIR"
fi

TARBALL_URL="https://github.com/TLCFEM/prebuilds/releases/download/latest/HDF5-2.2.0-$RUNNER_IMAGE_NAME.tar.gz"
TMP_DIR="$(mktemp -d)"

fetch_archive "$TMP_DIR/archive.tar.gz" "$TARBALL_URL"
tar -xzf "$TMP_DIR/archive.tar.gz" -C "$TMP_DIR"

find "$TMP_DIR/lib" -name "*.a" -exec cp {} "$TARGET_DIR" \;

INCLUDE_DIR="$(dirname "$0")/../Include/hdf"

rm -rf "${INCLUDE_DIR:?}/*" "$INCLUDE_DIR-$TOKEN/*"

find "$TMP_DIR/include" -name "*.h" -exec cp {} "$INCLUDE_DIR" \;

rm -rf "$TMP_DIR"

TARBALL_URL="https://github.com/TLCFEM/prebuilds/releases/download/latest/tbb-v2023.1.0-$RUNNER_IMAGE_NAME.tar.gz"
TMP_DIR="$(mktemp -d)"

fetch_archive "$TMP_DIR/archive.tar.gz" "$TARBALL_URL"
tar -xzf "$TMP_DIR/archive.tar.gz" -C "$TMP_DIR"

rm -rf "$TARGET_DIR/libtbb*"
rm -rf "$(dirname "$0")/../Include/tbb/*"
find "$TMP_DIR/tbb-install/lib" -name "lib*so*" -exec cp -P {} "$TARGET_DIR" \;
cp -a "$TMP_DIR/tbb-install/include/." "$(dirname "$0")/../Include/tbb"

rm -rf "$TMP_DIR"

TARBALL_URL="https://github.com/TLCFEM/prebuilds/releases/download/latest/mimalloc-$RUNNER_IMAGE_NAME.tar.gz"
TMP_DIR="$(mktemp -d)"

fetch_archive "$TMP_DIR/archive.tar.gz" "$TARBALL_URL"
tar -xzf "$TMP_DIR/archive.tar.gz" -C "$TMP_DIR"

find "$TMP_DIR/lib" -name "lib*so*" -exec cp -P {} "$TARGET_DIR" \;

rm -rf "$TMP_DIR"

TARBALL_URL="https://github.com/TLCFEM/prebuilds/releases/download/latest/OpenBLAS-0.3.34-$RUNNER_IMAGE_NAME-32.tar.gz"
TMP_DIR="$(mktemp -d)"

fetch_archive "$TMP_DIR/archive.tar.gz" "$TARBALL_URL"
tar -xzf "$TMP_DIR/archive.tar.gz" -C "$TMP_DIR"

rm -rf "$TARGET_DIR/libopenblas*"
cp "$TMP_DIR/libopenblas.a" "$TARGET_DIR"

rm -rf "$TMP_DIR"

TARBALL_URL="https://github.com/TLCFEM/prebuilds/releases/download/latest/VTK-9.7.1-$RUNNER_IMAGE_NAME.tar.gz"
TMP_DIR="$(mktemp -d)"

fetch_archive "$TMP_DIR/archive.tar.gz" "$TARBALL_URL"
mkdir -p "$(dirname "$0")/../VTK"
tar -xzf "$TMP_DIR/archive.tar.gz" -C "$(dirname "$0")/../VTK"

rm -rf "$TMP_DIR"

echo "Dependencies for $RUNNER_IMAGE_NAME have been downloaded and extracted to $TARGET_DIR."
