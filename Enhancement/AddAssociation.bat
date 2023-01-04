@echo off

assoc .sp=suanpanmodel
assoc .supan=suanpanmodel

set "program=%~dp0suanPan.exe"

if not exist "%program%" (
    echo suanPan.exe does not exist in current folder
    goto byebye
)

ftype suanpanmodel="%program%" "-f" "%%1"

set "program=%program:\=/%"

set "target=%appdata%\Sublime Text 3\Packages\User\"

if exist "%target%\" goto copyfile

set "target=%appdata%\Sublime Text\Packages\User\"

if exist "%target%\" goto copyfile

set target=

set /p target="Input path contains sublime_text.exe (leave empty to quit): "

if "%target%" == "" goto byebye

if not exist "%target%\sublime_text.exe" (
    echo sublime_text.exe does not exist in the given folder
    goto byebye
)

set "target=%target%\Data\Packages\User\"
set "target=%target:\=/%"

:copyfile

echo {"cmd":["%program%","-nc","-f","$file"],"selector":"source.supan","file_patterns":["*.supan","*.sp"]} > "%~dp0suanPan.sublime-build"
xcopy "%~dp0suanPan.sublime*" "%target%"

:byebye