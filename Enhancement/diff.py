from difflib import unified_diff
from pathlib import Path


def compare(a: Path, b: Path):
    for file_a in a.rglob("*"):
        if not file_a.is_file():
            continue

        rel_path = file_a.relative_to(a)
        file_b = b / rel_path
        if not file_b.exists() or not file_b.is_file():
            raise RuntimeError(f"Missing in b: {file_b}")

        diff_line = []
        for line in unified_diff(
            file_a.read_text().splitlines(),
            file_b.read_text().splitlines(),
            fromfile=str(file_a),
            tofile=str(file_b),
        ):
            if not line.startswith(("+", "-")):
                continue
            if "Time Wasted:" in line:
                continue
            if "by tlc @" in line:
                continue
            if str(a) in line or str(b) in line:
                continue
            diff_line.append(line)

        if diff_line:
            print(f"\033[91mDifferences in {rel_path}:\033[0m")
            print("\n".join(diff_line))


if __name__ == "__main__":
    compare(
        Path(
            "/home/theodore/Downloads/Results-51acc636976ba8e5f25056f2021aa4c4db8ea780a2a2347afbdbf408cfe4a219/Results"
        ),
        Path(
            "/home/theodore/Downloads/Results-41986ef2cf258b9eddfcfb67b49cf29dcba83c8d4a39b545c689093444860dfc/Results"
        ),
    )
