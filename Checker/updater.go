package main

import (
	"fmt"
	"io"
	"net/http"
	"os"
	"path/filepath"
	"regexp"
	"runtime"
	"strconv"
	"strings"
	"time"
)

const URL = "https://github.com/TLCFEM/suanPan/releases"

func main() {
	client := &http.Client{
		Timeout: 2 * time.Second,
	}

	response, err := client.Get(URL)
	if err != nil {
		return
	}

	defer func(Body io.ReadCloser) {
		err := Body.Close()
		if err != nil {
			return
		}
	}(response.Body)

	html, err := io.ReadAll(response.Body)
	if err != nil {
		return
	}

	regex, _ := regexp.Compile(`suanPan-v(\d)\.(\d)\.?(\d?)`)

	newVersion := string(regex.Find(html))

	if len(os.Args) > 1 {
		fmt.Printf("Checking new version, delete/rename file updater/updater.exe or execute with '-nu' flag if not wanted.\n")

		number := regex.FindStringSubmatch(newVersion)

		newMajor, _ := strconv.Atoi(number[1])
		newMinor, _ := strconv.Atoi(number[2])
		newPatch := 0

		if len(number) == 3 {
			newPatch, _ = strconv.Atoi(number[3])
		}

		currentVersion, _ := strconv.Atoi(os.Args[1])

		if 100*newMajor+10*newMinor+newPatch <= currentVersion {
			return
		}
	} else {
		fmt.Printf("Downloading the latest version.\n")
	}

	_ = downloadLatestVersion(newVersion)

	return
}

func downloadLatestVersion(versionString string) error {
	cos := runtime.GOOS

	if cos != "windows" && cos != "linux" && cos != "darwin" {
		return nil
	}

	fmt.Printf("Found new version %s, you can download it using this updater, or using package managers.\nDo you want to download now? [y/N] ", versionString)
	var downloadSwitch string
	_, err := fmt.Scanln(&downloadSwitch)
	if err != nil {
		return err
	}

	if len(downloadSwitch) == 0 || (downloadSwitch[0] != 'y' && downloadSwitch[0] != 'Y') {
		return nil
	}

	fmt.Printf("\nPlease note the following:\n")
	fmt.Printf("  `mkl` uses oneMKL that has the best performance on Intel platforms.\n")
	fmt.Printf("      Please use `mkl` version if possible.\n")
	fmt.Printf("  `openblas` uses OpenBLAS that may show better performance on AMD platforms.\n")
	fmt.Printf("      If performance degradation is spotted on your platform, consider switch to `openblas` version.\n")
	fmt.Printf("  `vtk` uses VTK for visualisation.\n")
	fmt.Printf("      Visualisation may be useful when it comes to post-processing, but it requires OpenGL support. Please make sure the corresponding packages are installed.\n")
	fmt.Printf("  `no-avx` disables AVX.\n")
	fmt.Printf("      For CPUs that do not support AVX, please use this version.\n")
	fmt.Printf("\nDownload the new version:\n")

	regex, _ := regexp.Compile(`suanPan-v(\d)\.(\d)\.?(\d?)`)

	number := regex.FindStringSubmatch(versionString)

	newMajor, _ := strconv.Atoi(number[1])
	newMinor, _ := strconv.Atoi(number[2])
	newPatch := 0

	if len(number) == 3 {
		newPatch, _ = strconv.Atoi(number[3])
	}

	if cos == "windows" {
		fmt.Printf("    [0] suanPan-win-mkl-vtk.exe (Installer)\n")
		fmt.Printf("    [1] suanPan-win-mkl-vtk.zip (Portable Archive)\n")
		fmt.Printf("    [2] suanPan-win-openblas-vtk.7z (Portable Archive)\n")
		fmt.Printf("    [3] suanPan-win-openblas-vtk-no-avx.7z (Portable Archive)\n")
	} else if cos == "linux" {
		fmt.Printf("    [0] suanPan-linux-mkl-vtk.tar.gz (Portable Archive)\n")
		fmt.Printf("    [1] suanPan-linux-mkl-vtk-no-avx.tar.gz (Portable Archive)\n")
		fmt.Printf("    [2] suanPan-linux-mkl.tar.gz (Portable Archive)\n")
		fmt.Printf("    [3] suanPan-linux-openblas-vtk.tar.gz (Portable Archive)\n")
		fmt.Printf("    [4] suanPan-linux-openblas.tar.gz (Portable Archive)\n")
		fmt.Printf("    [5] suanPan-%d.%d.%d-1.x86_64.deb (Debian Installer)\n", newMajor, newMinor, newPatch)
		fmt.Printf("    [6] suanPan-%d.%d.%d-1.x86_64.rpm (Red Hat Installer)\n", newMajor, newMinor, newPatch)
	} else if cos == "darwin" {
		fmt.Printf("    [0] suanPan-macos-openblas-vtk.tar.gz (Portable Archive)\n")
		fmt.Printf("    [1] suanPan-macos-openblas-vtk-no-avx.tar.gz (Portable Archive)\n")
		fmt.Printf("    [2] suanPan-macos-openblas.tar.gz (Portable Archive)\n")
	}

	fmt.Printf("\nPlease select the version you want to download (leave empty to exit): ")
	downloadOption := 0
	_, err = fmt.Scanf("%d", &downloadOption)
	if err != nil {
		return err
	}

	link := URL + "/download/" + versionString
	fileName := ""
	if cos == "windows" {
		if 0 == downloadOption {
			fileName = "suanPan-win-mkl-vtk.exe"
		} else if 1 == downloadOption {
			fileName = "suanPan-win-mkl-vtk.zip"
		} else if 2 == downloadOption {
			fileName = "suanPan-win-openblas-vtk.7z"
		} else if 3 == downloadOption {
			fileName = "suanPan-win-openblas-vtk-no-avx.7z"
		}
	} else if cos == "linux" {
		if 0 == downloadOption {
			fileName = "suanPan-linux-mkl-vtk.tar.gz"
		} else if 1 == downloadOption {
			fileName = "suanPan-linux-mkl-vtk-no-avx.tar.gz"
		} else if 2 == downloadOption {
			fileName = "suanPan-linux-mkl.tar.gz"
		} else if 3 == downloadOption {
			fileName = "suanPan-linux-openblas-vtk.tar.gz"
		} else if 4 == downloadOption {
			fileName = "suanPan-linux-openblas.tar.gz"
		} else if 5 == downloadOption {
			fileName = fmt.Sprintf("suanPan-%d.%d.%d-1.x86_64.deb", newMajor, newMinor, newPatch)
		} else if 6 == downloadOption {
			fileName = fmt.Sprintf("suanPan-%d.%d.%d-1.x86_64.rpm", newMajor, newMinor, newPatch)
		}
	} else if cos == "darwin" {
		if 0 == downloadOption {
			fileName = "suanPan-macos-openblas-vtk.tar.gz"
		} else if 1 == downloadOption {
			fileName = "suanPan-macos-openblas-vtk-no-avx.tar.gz"
		} else if 2 == downloadOption {
			fileName = "suanPan-macos-openblas.tar.gz"
		}
	}

	if fileName == "" {
		return nil
	}

	link += "/" + fileName

	fmt.Printf("Downloading files...\n")

	response, err := http.Get(link)
	if err != nil {
		return err
	}
	defer func(Body io.ReadCloser) {
		err := Body.Close()
		if err != nil {
			return
		}
	}(response.Body)

	storage, err := os.Create(fileName)
	if err != nil {
		return err
	}
	defer func(storage *os.File) {
		err := storage.Close()
		if err != nil {
			return
		}
	}(storage)

	_, _ = io.Copy(storage, response.Body)

	absPath, err := filepath.Abs(fileName)
	if err != nil {
		return err
	}

	fmt.Printf("Downloaded %s\n", absPath)

	isArchive := false
	isArchive = isArchive || strings.HasSuffix(absPath, "zip")
	isArchive = isArchive || strings.HasSuffix(absPath, "tar.gz")
	isArchive = isArchive || strings.HasSuffix(absPath, "7z")

	if isArchive {
		fmt.Printf("You can manually extract the archive to overwrite the existing folder.")
	}

	return nil
}
