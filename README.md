# suanPan

<img src="Resource/suanPan-qr-ua.svg" width="150" align="middle"/><img src="Resource/suanPan-ua.svg" width="150" align="middle"/>

[![doi](https://zenodo.org/badge/DOI/10.5281/zenodo.1285221.svg)](https://doi.org/10.5281/zenodo.1285221)
[![license](https://img.shields.io/github/license/TLCFEM/suanPan.svg?color=44cc11)](https://www.gnu.org/licenses/gpl-3.0)
[![documentation](https://readthedocs.org/projects/suanpan-manual/badge/?version=latest)](https://suanpan-manual.readthedocs.io/?badge=latest)
[![release](https://img.shields.io/github/release-pre/TLCFEM/suanPan.svg?color=44cc11)](https://github.com/TLCFEM/suanPan/releases)
[![suanpan](https://snapcraft.io//suanpan/badge.svg)](https://snapcraft.io/suanpan)
[![choco](https://img.shields.io/chocolatey/v/suanpan?color=44cc11)](https://chocolatey.org/packages/suanpan)
[![choco](https://img.shields.io/chocolatey/dt/suanpan?color=44cc11&label=choco%20install)](https://chocolatey.org/packages/suanpan)
[![flathub](https://img.shields.io/flathub/downloads/io.github.tlcfem.suanPan?label=flathub%20install)](https://flathub.org/apps/io.github.tlcfem.suanPan)
[![github](https://img.shields.io/github/downloads/TLCFEM/suanPan/total.svg?color=44cc11&label=github%20download)](https://github.com/TLCFEM/suanPan/releases)
[![docker](https://img.shields.io/docker/pulls/tlcfem/suanpan?color=44cc11)](https://hub.docker.com/r/tlcfem/suanpan/tags)
[![dev](https://github.com/TLCFEM/suanPan/actions/workflows/dev-all.yml/badge.svg?branch=dev)](https://github.com/TLCFEM/suanPan/actions/workflows/dev-all.yml)
[![codecov](https://codecov.io/gh/TLCFEM/suanPan/branch/dev/graph/badge.svg?token=65BF9DF697)](https://codecov.io/gh/TLCFEM/suanPan)
[![codacy](https://app.codacy.com/project/badge/Grade/1ea08c43edf342a8b00b21e585e63503)](https://app.codacy.com/gh/TLCFEM/suanPan/dashboard?utm_source=gh&utm_medium=referral&utm_content=&utm_campaign=Badge_grade)
[![codefactor](https://www.codefactor.io/repository/github/tlcfem/suanpan/badge)](https://www.codefactor.io/repository/github/tlcfem/suanpan)
[![language](https://img.shields.io/github/languages/count/TLCFEM/suanPan.svg?color=44cc11)](https://github.com/TLCFEM/suanPan)
[![language](https://img.shields.io/github/languages/top/TLCFEM/suanPan.svg?color=44cc11&logo=c%2B%2B)](https://github.com/TLCFEM/suanPan)
[![size](https://img.shields.io/github/languages/code-size/TLCFEM/suanPan.svg?color=44cc11)](https://github.com/TLCFEM/suanPan)
[![issues](https://img.shields.io/github/issues/TLCFEM/suanPan.svg?color=44cc11)](https://github.com/TLCFEM/suanPan/issues)
[![fossa](https://app.fossa.com/api/projects/git%2Bgithub.com%2FTLCFEM%2FsuanPan.svg?type=shield)](https://app.fossa.com/projects/git%2Bgithub.com%2FTLCFEM%2FsuanPan?ref=badge_shield)

[![gplv3-or-later](https://www.gnu.org/graphics/gplv3-or-later.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html)

[![marketplace](https://img.shields.io/visual-studio-marketplace/v/tlc.suanpan?label=VS%20Code&color=44cc11)](https://marketplace.visualstudio.com/items?itemName=tlc.suanpan)

> [!IMPORTANT]
> - **Feature requests can be made via creating [new issues](https://github.com/TLCFEM/suanPan/issues/new/choose).**
> - Check out the VS Code [extension](https://marketplace.visualstudio.com/items?itemName=tlc.suanpan) for syntax highlighting and autocompletion.
> - Check out the [documentation](https://tlcfem.github.io/suanPan-manual/latest/SUMMARY/) for a summary of all available functionalities.
> - Please star ⭐ the project!

## Introduction

[🧮 **suanPan**](https://tlcfem.github.io/suanPan/) is a finite element method (FEM) simulation platform for applications in fields such as solid mechanics and civil/structural/seismic engineering.
**suanPan** is written in modern high-quality C++ code and is targeted to provide an efficient, concise, flexible and reliable FEM simulation platform.

**suanPan** is partially influenced by popular (non-)commercial FEA packages, such
as [ABAQUS UNIFIED FEA](https://www.3ds.com/products-services/simulia/products/abaqus/), [ANSYS](http://www.ansys.com/)
and [OpenSees](http://opensees.berkeley.edu/).

<p align="center"><img src="https://repobeats.axiom.co/api/embed/dc2e0bf6251aec5423157991015e9f2582618226.svg"></p>

## Features

The highlights of **suanPan** are

- ✅ ***fast***, memory and thread safe
- ✅ [shared memory](https://en.wikipedia.org/wiki/Shared_memory) parallelism support
- ✅ [distributed memory](https://en.wikipedia.org/wiki/Distributed_memory) parallelism support
- ✅ multi-architecture (`amd64` and `arm64`) [support](https://hub.docker.com/r/tlcfem/suanpan)
- ✅ rich collection of elements, materials, time integration methods, etc.
- ✅ open source and easy to extend
- ✅ modern language features

## Quick Start

Execute the application out-of-the-box in terminal on Linux using one of the following commands depending on how the application is obtained.
See details below.

```bash
# in folder bin/ for linux portable tarball
./suanPan.sh
# for linux packages and snap
suanPan
# for flatpak
flatpak run io.github.tlcfem.suanPan
```

Or on Windows,

```powershell
# in the folder containing suanPan.exe
.\suanPan.exe
```

First time users can use `overview` command to go through a quick introduction.

```text
+--------------------------------------------------------+
|             ____             suanPan is an open source |
|   ___ _   _|  _ \ __ _ _ __     FEM framework (64-bit) |
|  / __| | | | |_) / _` | '_ \           Canopus (3.7.0) |
|  \__ \ |_| |  __/ (_| | | | |        by tlc @ c34df242 |
|  |___/\__,_|_|   \__,_|_| |_|      all rights reserved |
|                                 10.5281/zenodo.1285221 |
+--------------------------------------------------------+
|  https://github.com/TLCFEM/suanPan                     |
|  https://tlcfem.github.io/suanPan-manual/latest        |
+--------------------------------------------------------+
|  https://bit.ly/vsc-sp                                 |
+--------------------------------------------------------+

suanPan ~<> overview
```

Sample models are available for almost all models/commands.
Please check the `Example` folder for details.
Further details can be seen [here](https://tlcfem.gitbook.io/suanpan-manual/tutorial/obtain) regarding how to run model files.

## Installation

> [!WARNING]
> Only the 64-bit version is compiled.
> It is assumed that [**AVX2**](https://en.wikipedia.org/wiki/Advanced_Vector_Extensions) is available thus if the program
> fails, please check if your CPU supports AVX2.
> Alternatively, you can try the `no-avx` version.

### Windows

> [!NOTE]
> The dependencies are bundled with the archive.
> One may also install the VC++ redistributable [package](https://aka.ms/vs/17/release/vc_redist.x64.exe).
> If the application prompts that some file, such as `msvcp140.dll`, is missing (unlikely), please install the redistributable package.

#### Binary Package

The archives of binaries are released under [Release](https://github.com/TLCFEM/suanPan/releases) page.

1. `suanpan-win-mkl-vtk.zip` is the portable archive.
2. `suanpan-win-mkl-vtk.exe` is the installer.

#### Chocolatey

The binaries, which are compiled with Intel MKL and VTK, are available
on [Chocolatey](https://chocolatey.org/packages/suanpan), please use the following command to install the package.

1. Follow the [instructions](https://chocolatey.org/install) to install Chocolatey.

2. Use the following command to install `suanPan`.

   ```ps
   choco install suanpan
   ```

3. It is recommended to use a modern terminal such as [Windows Terminal](https://github.com/microsoft/terminal) for better output display.

The Chocolatey repo available to you may not be up-to-date.
If the latest version is not available, please try alternatives, such as portable binaries or scoop.

[<p align="center"><img src="Resource/choco.gif"></p>](https://asciinema.org/a/684063)

#### Scoop

It is also possible to use [Scoop](https://scoop.sh/) to install the package.

1. Install [Scoop](https://scoop.sh/).

   ```ps
   Set-ExecutionPolicy RemoteSigned -scope CurrentUser
   iwr -useb get.scoop.sh | iex
   ```

2. Install `suanPan`.

   ```ps
   scoop install suanpan
   ```

### Linux

Linux users are recommended to obtain the binaries via `snap` or `flatpak`.

#### Snap

The snap supports visualisation via VTK and uses Intel MKL for linear algebra.
The `edge` channel is in sync with the `dev` branch.
The `stable` channel is in sync with the `master` branch.

[![snap](https://snapcraft.io/static/images/badges/en/snap-store-black.svg)](https://snapcraft.io/suanpan)

[<p align="center"><img src="Resource/snap.gif"></p>](https://asciinema.org/a/684061)

#### Flatpak

Flatpak is also available if preferred.
The `beta` channel is in sync with the `dev` branch.
The `stable` channel is in sync with the `master` branch.

<a href='https://flathub.org/apps/details/io.github.tlcfem.suanPan'><img width='200' alt='Download on Flathub' src='https://flathub.org/assets/badges/flathub-badge-en.svg'/></a>

```bash
# add repo
flatpak remote-add --if-not-exists flathub https://flathub.org/repo/flathub.flatpakrepo
# or the beta channel
# flatpak remote-add --if-not-exists flathub-beta https://flathub.org/beta-repo/flathub-beta.flatpakrepo
# install
flatpak install suanPan
# define alias
echo "alias suanpan=\"flatpak run io.github.tlcfem.suanPan\"" >> ~/.bashrc
```

### Docker

It is also possible to compile the package via docker, check the dockerfiles under the `Script` folder, for any
questions please open an issue.

One can directly pull the image.
Using [Docker Hub](https://hub.docker.com/r/tlcfem/suanpan).

```bash
docker pull tlcfem/suanpan
```

Using [GitHub Container Registry](https://github.com/TLCFEM/suanPan/pkgs/container/suanpan).

```bash
docker pull ghcr.io/tlcfem/suanpan
```

### Other Platforms

Precompiled binaries are provided via CI/CD on macOS, Windows, and Ubuntu.
Please download the file from the [release](https://github.com/TLCFEM/suanPan/releases) page.

A few flavors are available:

1. `vtk` --- visualisation support is enabled, with this you can record VTK files for postprocessing, however, OpenGL
   may be missing on server systems
2. `mkl` --- linear algebra operations are offloaded to MKL, which gives the optimal performance on Intel chips
3. `openblas` --- linear algebra operations are offloaded to OpenBLAS, which may outperform MKL on AMD platforms
4. `aocl` --- linear algebra operations are offloaded to AOCL, which is optimized for AMD platforms
5. `no-avx` --- AVX2 support is disabled, useful for older CPUs that do not support AVX2 instructions
6. `win-gcc` --- GCC is used to compile the binary
7. `win` --- MSVC is used to compile the binary

Advanced users can compile the program from source by themselves to enable GPU based solvers which require
an available [CUDA](https://docs.nvidia.com/cuda/cuda-toolkit-release-notes/) and/or [MAGMA](https://icl.utk.edu/magma/) library.

### Automation Related

#### VS Code

The VS Code extension is available [here](https://marketplace.visualstudio.com/items?itemName=tlc.suanpan).

#### Sublime Text

On Windows, a batch file named `AddAssociation.bat` is provided in the archive.
It provides file associations and prepares a proper working environment (build system, autocompletion, highlighting) with [Sublime Text](https://www.sublimetext.com/).
It also adds file associations with `.sp` and `.supan` files, please run the `AddAssociation.bat` file with administrator privilege.
[Sublime Text](https://www.sublimetext.com/) autocompletion and syntax highlighting files are also provided.
Please install Sublime Text first and execute the batch file with the administrator privilege.

On Linux, a script file named as `suanPan.sh` is provided.

```bash
./suanPan.sh --create-link
```

The above command adds Sublime Text autocompletion and syntax highlighting files to proper location if Sublime Text configuration folder is found.
It also adds a command alias `suanpan` to `~/.local/bin` and a desktop file to `~/.local/share/applications`.

## Dependency

Additional libraries used in **suanPan** are listed as follows.

- [**AMD Optimizing CPU Libraries (AOCL)**](https://www.amd.com/en/developer/aocl.html) version 5.1
- [**ARPACK**](https://github.com/opencollab/arpack-ng)
- [**Armadillo**](http://arma.sourceforge.net/) version 14.6.2
- [**CUDA**](https://docs.nvidia.com/cuda/cuda-toolkit-release-notes/) version 12.9
- [**Catch2**](https://github.com/catchorg/Catch2) version 3.9.1
- [**FEAST**](http://www.feast-solver.org/) version 4.0
- [**HDF5**](https://www.hdfgroup.org/solutions/hdf5/) version 1.14.5
- [**MAGMA**](https://icl.utk.edu/magma/) version 2.9.0
- [**OpenBLAS**](https://github.com/xianyi/OpenBLAS) version 0.3.30
- [**SPIKE**](http://www.spike-solver.org/) version 1.0
- [**SuperLU MT**](https://portal.nersc.gov/project/sparse/superlu/) version 4.0.0
- [**SuperLU**](https://portal.nersc.gov/project/sparse/superlu/) version 7.0.1
- [**TBB** Threading Building Blocks](https://github.com/oneapi-src/oneTBB) version 2022.1.0
- [**VTK**](https://vtk.org/) version 9.5.0
- [**argparse**](https://github.com/p-ranav/argparse)
- [**exprtk**](https://github.com/ArashPartow/exprtk) version 0.0.3
- [**ezp**](https://github.com/TLCFEM/ezp)
- [**fmt**](https://github.com/fmtlib/fmt) version 10.2.1
- [**magic_enum**](https://github.com/Neargye/magic_enum) version 0.9.7
- [**oneMKL**](https://software.intel.com/content/www/us/en/develop/tools/oneapi/components/onemkl.html) version 2025.2.0
- [**whereami**](https://github.com/gpakosz/whereami)
- **thread_pool** abridged version of [`thread-pool`](https://github.com/bshoshany/thread-pool)

Those libraries may depend on other libraries such as [zlib](https://zlib.net/) and [Szip](https://support.hdfgroup.org/doc_resource/SZIP/).
Additional tools may be used by **suanPan**, they are

- [**UPX** the Ultimate Packer for eXecutables](https://upx.github.io/)

## How To Compile

Please refer to the corresponding [page](https://tlcfem.github.io/suanPan-manual/latest/Basic/Compile/) in the manual for details.

## Happy Modelling

![an example simulation of particle collision](Resource/particle-collision.gif)

## Licence

[![FOSSA Status](https://app.fossa.com/api/projects/git%2Bgithub.com%2FTLCFEM%2FsuanPan.svg?type=large)](https://app.fossa.com/projects/git%2Bgithub.com%2FTLCFEM%2FsuanPan?ref=badge_large)
