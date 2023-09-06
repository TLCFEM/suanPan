# This script collects all cmake-build* folders and compile them in parallel.
# It is mainly used to test different compilers and configurations.
# The cmake-build* naming convention is used by CLion.
# CLion also enables auto-reloading of CMakeLists.txt files, which means the cmake-build* folders are always up-to-date.
# One can then invoke this script such that
# ```ps
# .\Script\CompileAll.ps1
# ```

if (Test-Path -Path "CompileAll.ps1") { Set-Location -Path .. }

$folders = Get-ChildItem -Directory -Name -Path . | Where-Object { $_ -like "cmake-build*" } | ForEach-Object { Convert-Path $_ }
$folder_count = $folders.length
$per_core = [int][math]::Ceiling($env:NUMBER_OF_PROCESSORS / ($folder_count - 1) - 1)

$script = {
    param ( 
        [string]$folder,
        [int]$jobs
    )
    Set-Location $folder
    Write-Host "Compiling $folder with $jobs jobs."
    cmake --build . --parallel $jobs
}

foreach ($folder in $folders) { Start-Job -ScriptBlock $script -ArgumentList $folder, $per_core }

Wait-Job -Job (Get-Job) | Receive-Job
