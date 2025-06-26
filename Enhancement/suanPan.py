import hashlib
import os
import shutil
import sys


def pack(root_folder: str):
    dest = "tools"

    if os.path.exists(dest):
        shutil.rmtree(dest)

    os.mkdir(dest)

    with open(dest + "/VERIFICATION.txt", "w") as checksum:
        checksum.write(
            "VERIFICATION\n"
            "\tTo verify, please visit https://github.com/TLCFEM/suanPan/releases"
            " where the same archive is uploaded.\n"
            "\nDEPENDENCY\n"
        )

        for f in os.listdir(root_folder):
            if f.endswith(".dll"):
                full_path = os.path.join(root_folder, f)
                with open(full_path, "rb") as target:
                    checksum.write(
                        f"\t{f} : {hashlib.sha256(target.read()).hexdigest()}\n"
                    )
                shutil.copy(full_path, dest)

        checksum.write("\nEXECUTABLE\n")

        for f in ["suanPan.exe"]:
            full_path = os.path.join(root_folder, f)
            if os.path.exists(full_path):
                with open(full_path, "rb") as target:
                    checksum.write(
                        f"\t{f} : {hashlib.sha256(target.read()).hexdigest()}\n"
                    )
                shutil.copy(full_path, dest)

        checksum.write("\nTEXT\n")

        for f in [
            "AddAssociation.bat",
            "suanPan.sublime-completions",
            "suanPan.sublime-syntax",
            "LICENSE",
            "CHANGELOG.md",
            "README.md",
        ]:
            full_path = os.path.join(root_folder, f)
            if os.path.exists(full_path):
                checksum.write(f"\t{f}\n")
                shutil.copy(full_path, dest)

    with open("suanPan-win-mkl-vtk.sha256", "w") as hashfile, open(
        "suanPan-win-mkl-vtk.zip", "rb"
    ) as target:
        hashfile.write(hashlib.sha256(target.read()).hexdigest())

    shutil.copy("Enhancement/suanPan.nuspec", ".")

    os.system("choco pack")


if __name__ == "__main__":
    pack(sys.argv[1] if len(sys.argv) > 1 else "build/dist/bin")
