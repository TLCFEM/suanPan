#  Copyright (C) 2017-2025 Theodore Chang
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <https://www.gnu.org/licenses/>.

import json
import re
import sys
from os.path import abspath
from urllib.request import urlopen


def check_version(_major: int, _minor: int, _patch: int):
    try:
        with urlopen(
            "https://api.github.com/repos/TLCFEM/suanPan/releases/latest", timeout=10
        ) as response:
            data = json.load(response)
            tag_name = data["tag_name"]
            assets = [asset["name"] for asset in data["assets"]]
    except Exception:
        return

    if not (version := re.search(r"suanPan-v(\d)\.(\d)\.?(\d)?", tag_name)):
        return

    new_major = int(version.group(1))
    new_minor = int(version.group(2))
    new_patch = int(version.group(3)) if version.group(3) else 0

    if (
        100 * new_major + 10 * new_minor + new_patch
        <= 100 * _major + 10 * _minor + _patch
    ):
        return

    result = input(f"New version {tag_name} available, download now? [y/N] ")

    if result == "" or result[0] != "y" and result[0] != "Y":
        return

    print("\nDownload the new version:")

    def _get_platform_list(token: str):
        return [
            x for x in assets if token in x and x.endswith((".zip", ".tar.gz", ".7z"))
        ]

    version_list = []
    if sys.platform.startswith("win32"):
        version_list = _get_platform_list("win")
    elif sys.platform.startswith("linux"):
        version_list = _get_platform_list("linux")
    elif sys.platform.startswith("darwin"):
        version_list = _get_platform_list("mac")

    version_list.sort()

    for index, item in enumerate(version_list):
        print(f"  [{index}] {item}")

    result = input(
        "\nPlease select the version you want to download (leave empty to exit): "
    )

    if not result.isdigit():
        return

    result = int(result)

    if result < 0 or result >= len(version_list):
        return

    try:
        file_name = version_list[result]
        url: str = f"https://github.com/TLCFEM/suanPan/releases/download/{tag_name}/{file_name}"
        with urlopen(url, timeout=30) as source, open(file_name, "wb") as destination:
            destination.write(source.read())
        print(f"\nDownloaded {abspath(file_name)}.")
        print("You can manually extract the archive to overwrite the existing folder.")
    except Exception:
        print("Download failed.")


if __name__ == "__main__":
    if len(sys.argv) == 3:
        major, minor, patch = int(sys.argv[1]), int(sys.argv[2]), 0
    elif len(sys.argv) == 4:
        major, minor, patch = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    else:
        major, minor, patch = 0, 0, 0

    check_version(major, minor, patch)
