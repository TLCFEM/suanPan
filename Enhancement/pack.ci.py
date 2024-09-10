import hashlib
import os
import shutil
import sys
import zipfile


def pack(root_folder: str):
    fn = "tools"

    if os.path.exists(fn):
        shutil.rmtree(fn)

    os.mkdir(fn)

    verification = open(fn + "/VERIFICATION.txt", "w")

    verification.write("VERIFICATION\n")
    verification.write(
        "\tTo verify, please visit https://github.com/TLCFEM/suanPan/releases"
        " where the same archive is uploaded.\n"
    )

    print("[+] adding dependencies ...")

    verification.write("\nDEPENDENCY\n")

    for f in os.listdir(root_folder):
        if f.endswith(".dll"):
            full_path = os.path.join(root_folder, f)
            with open(full_path, "rb") as target:
                verification.write(
                    f"\t{os.path.basename(full_path)} : {hashlib.sha256(target.read()).hexdigest()}\n"
                )
            shutil.copy(full_path, fn)

    print("[+] adding executables ...")

    verification.write("\nEXECUTABLE\n")

    for f in ["suanPan.exe"]:
        full_path = os.path.join(root_folder, f)
        if os.path.exists(full_path):
            with open(full_path, "rb") as target:
                verification.write(
                    f"\t{f} : {hashlib.sha256(target.read()).hexdigest()}\n"
                )
            shutil.copy(full_path, fn)

    print("[+] adding text files ...")

    verification.write("\nTEXT\n")

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
            verification.write(f"\t{f}\n")
            shutil.copy(full_path, fn)

    verification.close()

    print("[+] packing archive ...")

    an = "suanPan-win-mkl-vtk"

    pwd = os.getcwd()

    os.chdir(fn)
    archive = zipfile.ZipFile("../" + an + ".zip", "w", zipfile.ZIP_DEFLATED)
    for f in os.listdir("."):
        archive.write(f)

    archive.close()

    os.chdir(pwd)

    with open(an + ".sha256", "w") as hashfile, open(an + ".zip", "rb") as target:
        hashfile.write(hashlib.sha256(target.read()).hexdigest())

    print("[+] packing chocolatey package ...")

    shutil.copy("Enhancement/suanpan.nuspec", ".")

    os.system("choco pack")


if __name__ == "__main__":
    # "build/dist/bin"
    pack(sys.argv[1])
