/*******************************************************************************
 * Copyright (C) 2017-2025 Theodore Chang
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 ******************************************************************************/

package main

import (
	"archive/tar"
	"compress/gzip"
	"fmt"
	"io"
	"net/http"
	"os"
	"os/exec"
	"path/filepath"
	"regexp"
	"runtime"
	"strconv"
	"strings"
	"time"
)

const URL = "https://github.com/TLCFEM/suanPan/releases"

func main() {
	if len(os.Args) > 1 && os.Args[1] == "--rename" && len(os.Args) > 2 {
		originalPath := os.Args[2]
		for i := 1; i <= 10; i++ {
			time.Sleep(time.Duration(100*i) * time.Millisecond)
			err := os.Rename(originalPath+".new", originalPath)
			if err == nil {
				break
			}
		}
	} else {
		fetch()
	}
}

func fetch() {
	client := &http.Client{
		Timeout: 2 * time.Second,
	}

	response, err := client.Get(URL)
	if err != nil {
		return
	}
	defer response.Body.Close()

	html, err := io.ReadAll(response.Body)
	if err != nil {
		return
	}

	regex, _ := regexp.Compile(`suanPan-v(\d)\.(\d)\.?(\d?)`)

	newVersion := string(regex.Find(html))

	fromMain := len(os.Args) > 1

	if fromMain {
		fmt.Printf("Checking new version, delete/rename file updater/updater.exe or execute with '-nu' flag if not wanted.\n")

		number := regex.FindStringSubmatch(newVersion)

		newMajor, _ := strconv.Atoi(number[1])
		newMinor, _ := strconv.Atoi(number[2])
		newPatch := 0

		if "" != number[3] {
			newPatch, _ = strconv.Atoi(number[3])
		}

		currentVersion, _ := strconv.Atoi(os.Args[1])

		if 100*newMajor+10*newMinor+newPatch <= currentVersion {
			return
		}
	} else {
		fmt.Printf("Downloading the latest version.\n")
	}

	_ = downloadLatestVersion(newVersion, fromMain)
}

func downloadLatestVersion(versionString string, fromMain bool) error {
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
	fmt.Printf("      Please use `mkl` version on Intel platforms.\n")
	fmt.Printf("  `aocl` uses AMD Optimizing CPU Libraries (AOCL) that has the best performance on AMD platforms.\n")
	fmt.Printf("      Please use `aocl` version on AMD platforms.\n")
	fmt.Printf("  `openblas` uses OpenBLAS that is general but does not have the best performance.\n")
	fmt.Printf("      Always prefer `mkl` and `aocl` versions if they are available.\n")
	fmt.Printf("  `vtk` uses VTK for visualisation.\n")
	fmt.Printf("      Visualisation may be useful when it comes to post-processing, but it requires OpenGL support. Please make sure the corresponding packages are installed.\n")
	fmt.Printf("  `no-avx` disables AVX2.\n")
	fmt.Printf("      For CPUs that do not support AVX2, please use this version.\n")
	fmt.Printf("\nDownload the new version:\n")

	var package_array []string

	switch cos {
	case "windows":
		package_array = []string{
			"suanPan-win-mkl-no-avx.zip",
			"suanPan-win-mkl-vtk-no-avx.zip",
			"suanPan-win-mkl-vtk.zip",
			"suanPan-win-mkl.zip",
			"suanPan-win-openblas-no-avx.7z",
			"suanPan-win-openblas-vtk-no-avx.7z",
			"suanPan-win-openblas-vtk.7z",
			"suanPan-win-openblas.7z",
		}
	case "linux":
		package_array = []string{
			"suanPan-linux-mkl-no-avx.tar.gz",
			"suanPan-linux-mkl-vtk-no-avx.tar.gz",
			"suanPan-linux-mkl-vtk.tar.gz",
			"suanPan-linux-mkl.tar.gz",
			"suanPan-linux-aocl-no-avx.tar.gz",
			"suanPan-linux-aocl-vtk-no-avx.tar.gz",
			"suanPan-linux-aocl-vtk.tar.gz",
			"suanPan-linux-aocl.tar.gz",
			"suanPan-linux-openblas-no-avx.tar.gz",
			"suanPan-linux-openblas-vtk-no-avx.tar.gz",
			"suanPan-linux-openblas-vtk.tar.gz",
			"suanPan-linux-openblas.tar.gz",
		}
	case "darwin":
		package_array = []string{
			"suanPan-macos-13-large-openblas-vtk.tar.gz",
			"suanPan-macos-13-large-openblas.tar.gz",
			"suanPan-macos-14-large-openblas-vtk.tar.gz",
			"suanPan-macos-14-large-openblas.tar.gz",
			"suanPan-macos-15-large-openblas-vtk.tar.gz",
			"suanPan-macos-15-large-openblas.tar.gz",
		}
	}

	for i, v := range package_array {
		fmt.Printf("    [%d] %s\n", i, v)
	}

	fmt.Printf("\nPlease select the version you want to download (leave empty to exit): ")
	downloadOption := 0
	_, err = fmt.Scanf("%d", &downloadOption)
	if err != nil {
		return err
	}

	link := URL + "/download/" + versionString
	fileName := ""
	if downloadOption < len(package_array) && downloadOption >= 0 {
		fileName = package_array[downloadOption]
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
	defer response.Body.Close()

	parentPath := filepath.Clean(filepath.Join(".", ".."))
	absPath, err := filepath.Abs(filepath.Join(parentPath, fileName))
	if err != nil {
		return err
	}

	storage, err := os.Create(absPath)
	if err != nil {
		return err
	}
	defer storage.Close()

	_, err = io.Copy(storage, response.Body)
	if err != nil {
		return err
	}

	fmt.Printf("Downloaded %s\n", absPath)

	isArchive := false
	isArchive = isArchive || strings.HasSuffix(absPath, "zip")
	isArchive = isArchive || strings.HasSuffix(absPath, "tar.gz")
	isArchive = isArchive || strings.HasSuffix(absPath, "7z")

	if isArchive {
		if cos != "linux" || fromMain {
			fmt.Printf("You can manually extract the archive to overwrite the existing folder.\n")
		} else {
			fmt.Printf("Do you want me to unpack the archive? [y/N] ")
			var unpackSwitch string
			_, err := fmt.Scanln(&unpackSwitch)
			if err != nil {
				return err
			}

			if len(unpackSwitch) == 0 || (unpackSwitch[0] != 'y' && unpackSwitch[0] != 'Y') {
				return nil
			}

			fmt.Printf("Overwriting the parent folder.\n")
			err = unpack(absPath, parentPath)
			if err != nil {
				return fmt.Errorf("Error unpacking the archive: %v\n", err)
			}

			selfPath, tmpPath, err := copySelf()
			if err != nil {
				return err
			}

			cmd := exec.Command(tmpPath, "--rename", selfPath)
			cmd.Stdout = os.Stdout
			cmd.Stderr = os.Stderr
			if err := cmd.Start(); err != nil {
				return err
			}
		}
	}

	return nil
}

func unpack(src, dest string) error {
	file, err := os.Open(src)
	if err != nil {
		return err
	}
	defer file.Close()

	gzipReader, err := gzip.NewReader(file)
	if err != nil {
		return err
	}
	defer gzipReader.Close()

	tarReader := tar.NewReader(gzipReader)

	for {
		header, err := tarReader.Next()
		if err == io.EOF {
			break
		}
		if err != nil {
			return err
		}

		targetPath := filepath.Join(dest, header.Name)
		if strings.HasSuffix(header.Name, "updater") {
			targetPath = targetPath + ".new"
		}

		switch header.Typeflag {
		case tar.TypeDir:
			if err := os.MkdirAll(targetPath, os.FileMode(header.Mode)); err != nil {
				return err
			}
		case tar.TypeReg:
			if err := os.MkdirAll(filepath.Dir(targetPath), os.ModePerm); err != nil {
				return err
			}
			target, err := os.OpenFile(targetPath, os.O_CREATE|os.O_RDWR, os.FileMode(header.Mode))
			if err != nil {
				return err
			}
			defer target.Close()
			if _, err := io.Copy(target, tarReader); err != nil {
				return err
			}
		}
	}

	return nil
}

func copySelf() (string, string, error) {
	selfPath, err := os.Executable()
	if err != nil {
		return "", "", err
	}
	sourceFile, err := os.Open(selfPath)
	if err != nil {
		return "", "", err
	}
	defer sourceFile.Close()

	tmpPath := filepath.Join(os.TempDir(), "updater")
	tmpFile, err := os.OpenFile(tmpPath, os.O_WRONLY|os.O_CREATE|os.O_TRUNC, 0755)
	if err != nil {
		return "", "", err
	}
	defer tmpFile.Close()

	_, err = io.Copy(tmpFile, sourceFile)

	return selfPath, tmpPath, err
}
