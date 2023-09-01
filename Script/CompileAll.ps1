if (Test-Path -Path "CompileAll.ps1") {
    Set-Location -Path ..
}

$folders = Get-ChildItem -Directory -Name -Path . | Where-Object { $_ -like "cmake-build*" }
foreach ($folder in $folders) {
    Write-Output "Compiling $folder"
    Set-Location -Path $folder
    cmake --build .
    Set-Location -Path ..
}
