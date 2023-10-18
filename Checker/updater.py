import re
import sys
import urllib.request
from os.path import abspath


def check_version(_major: int, _minor: int, _patch: int):
    html = urllib.request.urlopen(
        "https://github.com/TLCFEM/suanPan/releases").read()

    version = re.search(r'suanPan-v(\d)\.(\d)\.?(\d)?', html.decode('utf-8'))

    if not version:
        return

    new_major = int(version.group(1))
    new_minor = int(version.group(2))
    new_patch = int(version.group(3)) if version.group(3) else 0

    if 100 * new_major + 10 * new_minor + new_patch <= 100 * _major + 10 * _minor + _patch:
        return

    url = f'https://github.com/TLCFEM/suanPan/releases/download/{version.group(0)}/'

    result = input(f'New version {version.group(0)} available, download now? [y/N] ')

    if result == '' or result[0] != 'y' and result[0] != 'Y':
        return

    print('\nDownload the new version:')

    version_list = []
    if sys.platform.startswith('win32'):
        version_list = [
            "suanPan-win-mkl-vtk.exe",
            "suanPan-win-mkl-vtk.zip",
            "suanPan-win-openblas.7z",
            "suanPan-win-openblas-vtk.7z",
            "suanPan-win-openblas-no-avx.7z",
            "suanPan-win-openblas-vtk-no-avx.7z"
        ]
    elif sys.platform.startswith('linux'):
        version_list = [
            "suanPan-linux-mkl.tar.gz",
            "suanPan-linux-mkl-vtk.tar.gz",
            "suanPan-linux-openblas.tar.gz",
            "suanPan-linux-openblas-vtk.tar.gz",
            "suanPan-linux-mkl-no-avx.tar.gz",
            "suanPan-linux-mkl-vtk-no-avx.tar.gz",
            "suanPan-linux-openblas-no-avx.tar.gz",
            "suanPan-linux-openblas-vtk-no-avx.tar.gz"
        ]
    elif sys.platform.startswith('darwin'):
        version_list = [
            "suanPan-macos-openblas.tar.gz",
            "suanPan-macos-openblas-vtk.tar.gz"
        ]

    for index, item in enumerate(version_list):
        print(f'  [{index}] {item}')

    result = input("\nPlease select the version you want to download (leave empty to exit): ")

    if not result.isdigit():
        return

    result = int(result)

    if result < 0 or result >= len(version_list):
        return

    file_name = version_list[result]

    urllib.request.urlretrieve(url + file_name, file_name)

    print(f"\nDownloaded {abspath(file_name)}")
    print("You can manually extract the archive to overwrite the existing folder.")


if __name__ == '__main__':
    if len(sys.argv) == 3:
        major, minor, patch = int(sys.argv[1]), int(sys.argv[2]), 0
    elif len(sys.argv) == 4:
        major, minor, patch = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    else:
        major, minor, patch = 0, 0, 0

    check_version(major, minor, patch)
