/*******************************************************************************
 * Copyright (C) 2017-2026 Theodore Chang
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
	"bufio"
	"compress/gzip"
	"encoding/json"
	"fmt"
	"io"
	"log"
	"net/http"
	"os"
	"os/exec"
	"path/filepath"
	"regexp"
	"runtime"
	"sort"
	"strconv"
	"strings"
	"time"
)

const URL = "https://github.com/TLCFEM/suanPan/releases"

func main() {
	if len(os.Args) > 2 && os.Args[1] == "--rename" {
		originalPath := os.Args[2]
		newPath := originalPath + ".new"
		for i := 1; i <= 10; i++ {
			time.Sleep(time.Duration(100*i) * time.Millisecond)
			err := os.Rename(newPath, originalPath)
			if err == nil {
				break
			}
		}
	} else if err := fetch(); err != nil {
		log.Printf("Error: %v", err)
	}
}

type Asset struct {
	Name string `json:"name"`
}

type Release struct {
	TagName string  `json:"tag_name"`
	Assets  []Asset `json:"assets"`
}

func fetch() error {
	client := &http.Client{
		Timeout: 2 * time.Second,
	}

	response, err := client.Get("https://api.github.com/repos/TLCFEM/suanPan/releases/latest")
	if err != nil {
		return err
	}
	defer response.Body.Close()

	var release Release
	if err := json.NewDecoder(response.Body).Decode(&release); err != nil {
		return err
	}

	fromMain := len(os.Args) > 1

	if fromMain {
		regex, _ := regexp.Compile(`suanPan-v(\d)\.(\d)\.?(\d?)`)
		number := regex.FindStringSubmatch(release.TagName)

		newMajor, _ := strconv.Atoi(number[1])
		newMinor, _ := strconv.Atoi(number[2])
		newPatch := 0

		if number[3] != "" {
			newPatch, _ = strconv.Atoi(number[3])
		}

		currentVersion, _ := strconv.Atoi(os.Args[1])

		if 100*newMajor+10*newMinor+newPatch <= currentVersion {
			return nil
		}
	} else {
		fmt.Printf("Downloading the latest version.\n")
	}

	return downloadLatestVersion(release, fromMain)
}

func isArchive(name string) bool {
	ext := strings.ToLower(name)
	return strings.HasSuffix(ext, ".zip") || strings.HasSuffix(ext, ".tar.gz") || strings.HasSuffix(ext, ".7z")
}

func getPackageList(release Release, platforms []string) []string {
	var packageList []string
	for _, asset := range release.Assets {
		for _, platform := range platforms {
			if strings.Contains(asset.Name, platform) && isArchive(asset.Name) {
				packageList = append(packageList, asset.Name)
				break
			}
		}
	}
	sort.Strings(packageList)
	return packageList
}

func downloadLatestVersion(release Release, fromMain bool) error {
	cos := runtime.GOOS

	if cos != "windows" && cos != "linux" && cos != "darwin" {
		return nil
	}

	fmt.Printf("Found new version %s, you can download it using this updater, or using package managers.\nDo you want to download now? [y/N] ", release.TagName)

	reader := bufio.NewReader(os.Stdin)
	input, err := reader.ReadString('\n')
	if err != nil {
		return err
	}

	downloadSwitch := strings.TrimSpace(input)

	if len(downloadSwitch) == 0 || (downloadSwitch[0] != 'y' && downloadSwitch[0] != 'Y') {
		return nil
	}

	fmt.Printf("\nPlease note the following:\n")
	fmt.Printf("  `amd64` represents x86_64 Intel architecture.\n")
	fmt.Printf("  `arm64` represents ARM64/AAarch64 architecture, not applicable to Windows.\n")
	fmt.Printf("  `mkl` uses oneMKL that has the best performance on Intel platforms.\n")
	fmt.Printf("      Please use `mkl` version on Intel platforms.\n")
	fmt.Printf("  `aocl` uses AMD Optimizing CPU Libraries (AOCL) that has the best performance on AMD platforms.\n")
	fmt.Printf("      Please use `aocl` version on AMD platforms.\n")
	fmt.Printf("  `openblas` uses OpenBLAS that is general but does not have the best performance.\n")
	fmt.Printf("      Always prefer `mkl` and `aocl` versions if they are available.\n")
	fmt.Printf("  `vtk` uses VTK for visualisation.\n")
	fmt.Printf("      Visualisation may be useful when it comes to post-processing, but it requires OpenGL support. Please make sure the corresponding packages are installed.\n")
	fmt.Printf("  `avx`/`no-avx` enables/disables AVX2, not applicable to `arm64` builds.\n")
	fmt.Printf("      For CPUs that do not support AVX2, please use the `no-avx` version.\n")
	fmt.Printf("  `mpi` enables distributed parallelism.\n")
	fmt.Printf("  `ilp64` enables 64-bit integer for indexing (default is 32-bit), not well tested, extensive testing is welcome.\n")
	fmt.Printf("\nDownload the new version:\n")

	var packageArray []string

	switch cos {
	case "windows":
		packageArray = getPackageList(release, []string{"win"})
	case "linux":
		packageArray = getPackageList(release, []string{"linux"})
	case "darwin":
		packageArray = getPackageList(release, []string{"macos"})
	}

	for i, v := range packageArray {
		fmt.Printf("    [%d] %s\n", i, v)
	}

	fmt.Printf("\nPlease select the version you want to download (leave empty to exit): ")

	reader = bufio.NewReader(os.Stdin)
	input, err = reader.ReadString('\n')
	if err != nil {
		return err
	}

	input = strings.TrimSpace(input)

	if input == "" {
		return nil
	}

	downloadOption, err := strconv.Atoi(input)
	if err != nil {
		return err
	}

	if downloadOption >= len(packageArray) || downloadOption < 0 {
		return nil
	}

	fileName := packageArray[downloadOption]
	link := URL + "/download/" + release.TagName + "/" + fileName

	fmt.Printf("Downloading files...\n")

	response, err := http.Get(link)
	if err != nil {
		return err
	}
	defer response.Body.Close()

	var parentPath string
	if cos != "windows" {
		parentPath = filepath.Clean(filepath.Join(".", ".."))
	} else {
		parentPath = filepath.Clean(".")
	}
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

	if isArchive(absPath) {
		if cos == "windows" || fromMain {
			fmt.Printf("You can manually extract the archive to overwrite the existing folder.\n")
		} else {
			fmt.Printf("Do you want me to unpack the archive? [y/N] ")

			reader := bufio.NewReader(os.Stdin)
			input, err := reader.ReadString('\n')
			if err != nil {
				return err
			}

			unpackSwitch := strings.TrimSpace(input)

			if len(unpackSwitch) == 0 || (unpackSwitch[0] != 'y' && unpackSwitch[0] != 'Y') {
				return nil
			}

			fmt.Printf("Overwriting the parent folder.\n")
			if err := unpack(absPath, parentPath); err != nil {
				return err
			}

			selfPath, tmpPath, err := copySelf()
			if err != nil {
				return err
			}

			cmd := exec.Command(tmpPath, "--rename", selfPath)
			cmd.Stdout = os.Stdout
			cmd.Stderr = os.Stderr
			cmd.Start()
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
