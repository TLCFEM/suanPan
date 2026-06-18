#!/usr/bin/env python3


from pathlib import Path
from difflib import unified_diff
import sys


skip_files: tuple = (
    "DKTS3.supan",
    "Example.supan",
    "FEAST.BAND.supan",
    "FEAST.FULL.supan",
    "FEAST.SPARSE.supan",
    "MPC.supan",
    "MultilinearElastic1D.supan",
)


def skip(line: str):
    if "Time Wasted" in line:
        return True
    if "by tlc @" in line:
        return True

    return False


def readlines(file_path: Path):
    with file_path.open("r") as file:
        return [x for x in file if not skip(x)]


def compare_commits(current: Path, parent: Path):
    if not current.is_dir() or not parent.is_dir():
        return True

    error_flag: bool = False

    for current_path in current.rglob("*"):
        if not current_path.is_file():
            continue

        if current_path.name in skip_files:
            continue

        relative_path = current_path.relative_to(current)
        parent_path = parent / relative_path

        if not parent_path.exists():
            print(f"{parent_path} does not exist in the previous commit.")
            continue

        try:
            lines = list(
                unified_diff(
                    readlines(current_path),
                    readlines(parent_path),
                    fromfile=current_path.as_posix(),
                    tofile=parent_path.as_posix(),
                    lineterm="",
                )
            )

            if not lines:
                continue

            error_flag = True

            print(f"\n{'=' * 80}\nComparing: {relative_path}\n{'=' * 80}")

            for line in lines:
                print(line)

        except Exception as e:
            print(f"Error reading files: {e}.")

        return error_flag


if __name__ == "__main__":
    if len(sys.argv) > 2:
        if compare_commits(Path(sys.argv[1]), Path(sys.argv[2])):
            sys.exit(1)
