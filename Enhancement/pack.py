import hashlib
import os
import shutil
import zipfile

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

# intel mkl + tbb + msvc redistributable
for folder in [".", "../suanPan/Libs/vs", "../suanPan/Libs/vc"]:
    for f in os.listdir(folder):
        if f.endswith(".dll"):
            full_path = os.path.join(folder, f)
            base_name = os.path.basename(full_path)
            with open(full_path, "rb") as target:
                verification.write(
                    f"\t{base_name} : {hashlib.sha256(target.read()).hexdigest()}\n"
                )
            shutil.copy(full_path, fn)

print("[+] adding executables ...")

verification.write("\nEXECUTABLE\n")

for f in ["../suanPan/build/Release/suanPan.exe", "../suanPan/Checker/updater.exe"]:
    if os.path.exists(f):
        with open(f, "rb") as target:
            verification.write(
                f"\t{os.path.basename(f)} : {hashlib.sha256(target.read()).hexdigest()}\n"
            )
        shutil.copy(f, fn)

print("[+] adding text files ...")

verification.write("\nTEXT\n")

for f in [
    "../suanPan/Enhancement/AddAssociation.bat",
    "../suanPan/Enhancement/suanPan.sublime-completions",
    "../suanPan/Enhancement/suanPan.sublime-syntax",
    "../suanPan/LICENSE",
    "../suanPan/CHANGELOG.md",
    "../suanPan/README.md",
]:
    verification.write(f"\t{os.path.basename(f)}\n")
    shutil.copy(f, fn)

verification.close()

print("[+] packing archive ...")

an = "suanPan-win-mkl-vtk"
os.chdir(fn)
archive = zipfile.ZipFile("../" + an + ".zip", "w", zipfile.ZIP_DEFLATED)
for f in os.listdir("."):
    archive.write(f)

archive.close()

os.chdir("..")

with open(an + ".sha256", "w") as hashfile, open(an + ".zip", "rb") as target:
    hashfile.write(hashlib.sha256(target.read()).hexdigest())


print("[+] packing chocolatey package ...")

os.system("choco pack")

print("[+] compiling inno installer ...")

os.system(
    '"C:\\Program Files (x86)\\Inno Setup 6\\ISCC.exe" ..\\suanPan\\Enhancement\\suanPan.iss'
)
