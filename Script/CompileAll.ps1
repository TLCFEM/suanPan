# -----------------------------------------------------------------------------
# CompileAll.ps1 - Compile all CMake build folders in parallel (Windows)
#
# Discovers all cmake-build* and build/* folders and compiles them in parallel
# using multiple cores. Primarily used to test different compilers/configurations.
#
# Usage: .\CompileAll.ps1
# -----------------------------------------------------------------------------

if (Test-Path -Path "CompileAll.ps1") { Set-Location -Path .. }

$cmake_folders = Get-ChildItem -Directory -Path . -Filter "cmake-build*" | ForEach-Object { Convert-Path $_.FullName }
$build_folders = if (Test-Path -Path ".\build") {
    Get-ChildItem -Directory -Path ".\build" | ForEach-Object { Convert-Path $_.FullName }
}
else {
    @()
}
$folders = @($cmake_folders + $build_folders | Sort-Object -Unique)
$folder_count = $folders.Count
if ($folder_count -eq 0) {
    Write-Warning "No cmake-build* folders or build/* folders found."
    return
}
$per_core = [int][math]::Ceiling($env:NUMBER_OF_PROCESSORS / ($folder_count - 1) - 1)

$script = {
    param ( 
        [string]$folder,
        [int]$jobs
    )
    Set-Location $folder
    Write-Output "Compiling $folder with $jobs jobs."
    cmake --build . --parallel $jobs
}

foreach ($folder in $folders) { Start-Job -ScriptBlock $script -ArgumentList $folder, $per_core }

Wait-Job -Job (Get-Job) | Receive-Job
