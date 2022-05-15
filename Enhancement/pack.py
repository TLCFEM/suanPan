import os
import shutil
import zipfile
import hashlib

fn = 'tools'

if os.path.exists(fn):
    shutil.rmtree(fn)

os.mkdir(fn)

verification = open(fn + '/VERIFICATION.txt', 'w')

verification.write('VERIFICATION\n')
verification.write(
    '\tTo verify, please visit https://github.com/TLCFEM/suanPan/releases'
    ' where the same archive is uploaded.\n')

print('[+] adding dependencies ...')

verification.write('\nDEPENDENCY\n')

for f in os.listdir('.'):
    if f.endswith('.dll'):
        sha256 = hashlib.sha256(open(f, 'rb').read()).hexdigest()
        verification.write('\t' + f + ' : ' + sha256 + '\n')
        shutil.copy(f, fn)

print('[+] adding executables ...')

verification.write('\nEXECUTABLE\n')

fl = [
    '../suanPan-vs/Release/suanPan.exe',
    'updater.exe'
    ]
cl = [
    '\tsuanPan.exe : ',
    '\tupdater.exe : '
    ]

for f, c in zip(fl, cl):
    sha256 = hashlib.sha256(open(f, 'rb').read()).hexdigest()
    verification.write(c + sha256 + '\n')
    shutil.copy(f, fn)

print('[+] adding text files ...')

verification.write('\nTEXT\n')

fl = [
    '../suanPan/Enhancement/AddAssociation.bat',
    '../suanPan/Enhancement/suanPan.sublime-completions',
    '../suanPan/Enhancement/suanPan.sublime-syntax',
    '../suanPan/LICENSE',
    '../suanPan/CHANGELOG.md',
    '../suanPan/README.md'
    ]
cl = [
    '\tAddAssociation.bat\n',
    '\tsuanPan.sublime-completions\n',
    '\tsuanPan.sublime-syntax\n',
    '\tLICENSE\n',
    '\tCHANGELOG.md\n',
    '\tREADME.md\n'
    ]

for f, c in zip(fl, cl):
    verification.write(c)
    shutil.copy(f, fn)

verification.close()

print('[+] packing archive ...')

an = 'suanPan-win-mkl-vtk'
os.chdir(fn)
archive = zipfile.ZipFile('../' + an + '.zip', 'w', zipfile.ZIP_DEFLATED)
for f in os.listdir('.'):
    archive.write(f)

archive.close()

os.chdir('..')

hashfile = open(an + '.sha256', 'w')
hashfile.write(hashlib.sha256(open(an + '.zip', 'rb').read()).hexdigest())
hashfile.close()

print('[+] packing chocolatey package ...')

os.system('choco pack')

print('[+] compiling inno installer ...')

os.system('"C:\\Program Files (x86)\\Inno Setup 6\\ISCC.exe" ..\\suanPan\\Enhancement\\suanPan.iss')
