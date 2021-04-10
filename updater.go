package main

import (
	"io/ioutil"
	"net/http"
	"os"
	"regexp"
	"time"
)

func main() {
	url := "https://github.com/TLCFEM/suanPan/releases"

	client := &http.Client{
		Timeout: 2 * time.Second,
	}

	response, error := client.Get(url)
	if error != nil {
		removeExistingFile()
		return
	}

	defer response.Body.Close()

	html, error := ioutil.ReadAll(response.Body)
	if error != nil {
		removeExistingFile()
		return
	}

	regex, _ := regexp.Compile(`suanPan-v[0-9]\.[0-9]\.?[0-9]?`)

	removeExistingFile()

	versionFile, error := os.Create(".latest.version.sp")

	if error != nil {
		removeExistingFile()
		return
	}

	defer versionFile.Close()

	_, _ = versionFile.WriteString(string(regex.Find(html)))

	return
}

func removeExistingFile() {
	name := ".latest.version.sp"
	if _, error := os.Stat(name); error == nil {
		_ = os.Remove(name)
	}
}
