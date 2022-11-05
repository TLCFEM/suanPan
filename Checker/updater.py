import re
import sys
import urllib.request


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

    result = input(f'New version {version.group(0)} available, download now? [Y/n] ')

    if result == '' or result[0] != 'y' and result[0] != 'Y':
        return

    print('\nDownload the new version:')

    version_list = []

    if sys.platform.startswith('win32'):
        version_list.append(('suanPan-win-mkl-vtk.exe', '(Installer)'))
        version_list.append(('suanPan-win-mkl-vtk.zip', '(Portable Archive)'))
        version_list.append(('suanPan-win-openblas-vtk.7z', '(Portable Archive)'))
    elif sys.platform.startswith('linux'):
        version_list.append(('suanPan-linux-mkl-vtk.tar.gz', '(Portable Archive)'))
        version_list.append(('suanPan-linux-mkl.tar.gz', '(Portable Archive)'))
        version_list.append(('suanPan-linux-openblas-vtk.tar.gz', '(Portable Archive)'))
        version_list.append(('suanPan-linux-openblas.tar.gz', '(Portable Archive)'))
        version_list.append((f'suanPan-{new_major}.{new_minor}.{new_patch}-1.x86_64.deb', '(Debian Installer)'))
        version_list.append((f'suanPan-{new_major}.{new_minor}.{new_patch}-1.x86_64.rpm', '(Red Hat Installer)'))
    elif sys.platform.startswith('darwin'):
        version_list.append(('suanPan-macos-openblas-vtk.tar.gz', '(Portable Archive)'))
        version_list.append(('suanPan-macos-openblas.tar.gz', '(Portable Archive)'))

    for index, item in enumerate(version_list):
        print(f'  [{index}] {item[0]} {item[1]}')

    result = input("\nPlease select the version you want to download (leave empty to exit): ")

    if not result.isdigit():
        return

    result = int(result)

    if result < 0 or result >= len(version_list):
        return

    file_name = version_list[result][0]

    urllib.request.urlretrieve(url + file_name, file_name)

    print(f"\nDownloaded {file_name}")
    print("You can manually extract the archive to overwrite the existing folder.")


if __name__ == '__main__':
    if len(sys.argv) == 3:
        major, minor, patch = int(sys.argv[1]), int(sys.argv[2]), 0
    elif len(sys.argv) == 4:
        major, minor, patch = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    else:
        major, minor, patch = 0, 0, 0

    check_version(major, minor, patch)
