# suanPan

<img src="Resource/suanPan-qr-ua.svg" width="150" align="middle"/><img src="Resource/suanPan-ua.svg" width="150" align="middle"/>

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1285221.svg)](https://doi.org/10.5281/zenodo.1285221)
[![license](https://img.shields.io/github/license/TLCFEM/suanPan.svg?color=44cc11)](https://www.gnu.org/licenses/gpl-3.0)
[![documentation](https://readthedocs.org/projects/suanpan-manual/badge/?version=latest)](https://suanpan-manual.readthedocs.io/?badge=latest)
[![release](https://img.shields.io/github/release-pre/TLCFEM/suanPan.svg?color=44cc11)](https://github.com/TLCFEM/suanPan/releases)
[![suanpan](https://snapcraft.io//suanpan/badge.svg)](https://snapcraft.io/suanpan)
[![Chocolatey](https://img.shields.io/chocolatey/v/suanpan?color=44cc11)](https://chocolatey.org/packages/suanpan)
[![Chocolatey](https://img.shields.io/chocolatey/dt/suanpan?color=44cc11&label=choco%20install)](https://chocolatey.org/packages/suanpan)
[![download](https://img.shields.io/github/downloads/TLCFEM/suanPan/total.svg?color=44cc11)](https://img.shields.io/github/downloads/TLCFEM/suanPan/total.svg?color=44cc11)
[![stable build](https://github.com/TLCFEM/suanPan/workflows/Stable%20Release/badge.svg?branch=master)](https://github.com/TLCFEM/suanPan/actions)
[![AppVeyor](https://img.shields.io/appveyor/ci/TLCFEM/suanPan/master.svg?label=master&logo=appveyor)](https://ci.appveyor.com/project/TLCFEM/suanpan/branch/master)
[![codecov](https://codecov.io/gh/TLCFEM/suanPan/branch/dev/graph/badge.svg?token=65BF9DF697)](https://codecov.io/gh/TLCFEM/suanPan)
[![codacy](https://app.codacy.com/project/badge/Grade/1ea08c43edf342a8b00b21e585e63503)](https://www.codacy.com/gh/TLCFEM/suanPan/dashboard?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=TLCFEM/suanPan&amp;utm_campaign=Badge_Grade)
[![CodeFactor](https://www.codefactor.io/repository/github/tlcfem/suanpan/badge)](https://www.codefactor.io/repository/github/tlcfem/suanpan)
[![Codiga](https://api.codiga.io/project/22357/score/svg)](https://app.codiga.io/public/project/22357/suanPan/dashboard)
[![language](https://img.shields.io/github/languages/count/TLCFEM/suanPan.svg?color=44cc11)](https://github.com/TLCFEM/suanPan)
[![language](https://img.shields.io/github/languages/top/TLCFEM/suanPan.svg?color=44cc11&logo=c%2B%2B)](https://github.com/TLCFEM/suanPan)
[![size](https://img.shields.io/github/languages/code-size/TLCFEM/suanPan.svg?color=44cc11)](https://img.shields.io/github/languages/code-size/TLCFEM/suanPan.svg?color=44cc11)
[![issues](https://img.shields.io/github/issues/TLCFEM/suanPan.svg?color=44cc11)](https://github.com/TLCFEM/suanPan/issues)
[![chat](https://badges.gitter.im/suanPan-dev/community.svg)](https://gitter.im/suanPan-dev/community?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge)
[![FOSSA Status](https://app.fossa.com/api/projects/git%2Bgithub.com%2FTLCFEM%2FsuanPan.svg?type=shield)](https://app.fossa.com/projects/git%2Bgithub.com%2FTLCFEM%2FsuanPan?ref=badge_shield)

## Introduction

[🧮 **suanPan**](https://tlcfem.github.io/suanPan/) is a finite element method (FEM) simulation platform for
applications in fields such as solid mechanics and civil/structural/seismic engineering. The name **suanPan** (in some
places such as suffix it is also abbreviated as **suPan**) comes from the term *Suan Pan* (算盤), which
is [Chinese abacus](https://en.wikipedia.org/wiki/Suanpan). **suanPan** is written in modern high quality C++ code and
is targeted to provide an efficient, concise, flexible and reliable FEM simulation platform.

**suanPan** is partially influenced by popular (non-)commercial FEA packages, such
as [ABAQUS UNIFIED FEA](https://www.3ds.com/products-services/simulia/products/abaqus/), [ANSYS](http://www.ansys.com/)
and [OpenSees](http://opensees.berkeley.edu/).

Please check documentation [here](https://tlcfem.gitbook.io/suanpan-manual/) and [here](http://suanpan-manual.rtfd.io/)
for command references. Please consider starring ⭐ the project!

## Features

The highlights of **suanPan** are

- **suanPan** is *fast*, both memory and thread safe.
- **suanPan** is designed based on the [shared memory](https://en.wikipedia.org/wiki/Shared_memory) model and supports
  parallelism on heterogeneous architectures, for example multithreaded CPU + optional GPU. The parallelism is available
  for both element state updating and global matrix assembling.
- **suanPan** is open source and easy to be extended to incorporate user-defined elements, materials, etc.
- **suanPan** separates the FEA model part from the linear algebra operation part, which significantly reduces the
  complexity and cost of development of new models.
- **suanPan** utilizes the new language features shipped with the latest standards (C++11 to C++20), such as new STL
  containers, smart pointers and many others.
- **suanPan** supports simple visualization supported by [VTK](https://vtk.org/).

## Quick Start

Sample models are available for almost all models/commands. Please check the `Example` folder for details. Further
details can be seen [here](https://tlcfem.gitbook.io/suanpan-manual/tutorial/obtain) regarding how to run model files.

## Installation

Only 64-bit version is compiled. It is assumed that [**AVX**](https://en.wikipedia.org/wiki/Advanced_Vector_Extensions)
is available thus if the program fails, please check if your CPU supports AVX.

### Windows

#### Binary Package

The archives of binaries are released under [Release](https://github.com/TLCFEM/suanPan/releases) page.

1. `suanpan-win-mkl-vtk.zip` is the portable version.
2. `suanpan-win-mkl-vtk.exe` is the installer.

#### Chocolatey

The binaries, which are compiled with Intel MKL and VTK, are available
on [Chocolatey](https://chocolatey.org/packages/suanpan), please use the following command to install the package.

1. Follow the [instructions](https://chocolatey.org/install) to install Chocolatey.

2. Use the following command to install `suanPan`.

   ```ps
   choco install suanpan
   ```

3. It is recommended to use a modern terminal such as [Windows Terminal](https://github.com/microsoft/terminal)
   and [Fluent Terminal](https://github.com/felixse/FluentTerminal) for better output display.

[![asciicast](https://asciinema.org/a/491350.svg)](https://asciinema.org/a/491350)

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

Linux's users are recommended to obtain the binaries via snap. The snap supports visualization via VTK and uses Intel
MKL for linear algebra.

[![Get it from the Snap Store](https://snapcraft.io/static/images/badges/en/snap-store-black.svg)](https://snapcraft.io/suanpan)

[![asciicast](https://asciinema.org/a/491330.svg)](https://asciinema.org/a/491330)

Alternatively, download the RPM (Fedora 35) or DEB (Ubuntu 22.04) package from the release page. The packages may not be
compatible with older distributions (due to different version of `libstdc++`). It is also possible to compile the
package via docker, check the dockerfiles under the `Script` folder, for any questions please open an issue.

### Other Platforms

Precompiled binaries are provided via CI/CD on MacOS, Windows and Ubuntu. Please download the file from
the [release](https://github.com/TLCFEM/suanPan/releases) page.

A few flavors are available:

1. `vtk` --- visualization support is enabled, with this you can record VTK files for postprocessing, however, OpenGL
   may be missing on server systems
2. `mkl` --- linear algebra operations are offloaded to MKL, which gives optimal performance on Intel chips
3. `openblas` --- linear algebra operations are offloaded to OpenBLAS, which may outperform MKL on AMD platforms

Advanced users can compile the program from source by themselves in order to enable GPU based solvers which require
available [CUDA](https://docs.nvidia.com/cuda/cuda-toolkit-release-notes/) library.

Since CI/CD uses `GCC 10.3.0` (on Linux) and `Clang 13.0.1` (on MacOS), it may be required to update/install
proper `libstdc++` (or `libc++`) version. The easiest way is to install the same compiler. For example, on Ubuntu,

```bash
# Ubuntu
sudo apt install gcc-10 g++-10 gfortran-10 libomp5
```

For VTK enabled versions, it may be necessary to install OpenGL.

```bash
# Ubuntu
sudo apt install libglu1-mesa-dev freeglut3-dev mesa-common-dev libglvnd-dev
```

### Automation Related

On Windows, a batch file named as `AddAssociation.bat` is provided in the archive. It provides file associations and
prepares a proper working environment (build system, autocompletion, highlighting)
with [Sublime Text](https://www.sublimetext.com/). It also adds file associations with `.sp` and `.supan` files, please
run the `AddAssociation.bat` file with administrator privilege. [Sublime Text](https://www.sublimetext.com/)
autocompletion and syntax highlighting files are also provided. Please install Sublime Text first and execute the batch
file with the administrator privilege.

On Linux, a script file named as `suanPan.sh` is provided. By executing

```bash
./suanPan.sh --create-link
```

It adds Sublime Text autocompletion and syntax highlighting files to proper location if Sublime Text configuration
folder is found. It also adds a command alias `suanpan` to `~/.local/bin` and a desktop file
to `~/.local/share/applications`.

## Dependency

Additional libraries used in **suanPan** are listed as follows.

- [**ARPACK**](https://www.caam.rice.edu/software/ARPACK/) version 0.96
- [**SPIKE**](http://www.spike-solver.org/) version 1.0
- [**FEAST**](http://www.feast-solver.org/) version 4.0
- [**SuperLU**](https://portal.nersc.gov/project/sparse/superlu/) version 5.3.0 and [**SuperLU MT**](https://portal.nersc.gov/project/sparse/superlu/) version 3.1
- [**OpenBLAS**](https://github.com/xianyi/OpenBLAS) version 0.3.15
- [**TBB** Threading Building Blocks](https://github.com/oneapi-src/oneTBB) version 2021.5.0
- [**HDF5**](https://www.hdfgroup.org/solutions/hdf5/) version 1.10.6
- [**MUMPS**](http://mumps.enseeiht.fr/) version 5.2.1
- [**METIS**](http://glaros.dtc.umn.edu/gkhome/metis/metis/overview) version 5.1.0
- [**VTK**](https://vtk.org/) version 9.1
- [**CUDA**](https://docs.nvidia.com/cuda/cuda-toolkit-release-notes/) version 11.7
- [**Armadillo**](http://arma.sourceforge.net/) version 11.0
- [**ensmallen**](https://ensmallen.org/) version 2.19.0
- [**oneMKL**](https://software.intel.com/content/www/us/en/develop/tools/oneapi/components/onemkl.html) version 2022.1.0
- [**Catch2**](https://github.com/catchorg/Catch2) version 2.13.9
- **thread_pool** abridged version of [`thread-pool`](https://github.com/bshoshany/thread-pool)

Those libraries may depend on other libraries such as [zlib](https://zlib.net/)
and [Szip](https://support.hdfgroup.org/doc_resource/SZIP/). Additional tools may be used by **suanPan**, they are

- [**UPX** the Ultimate Packer for eXecutables](https://upx.github.io/)

## How To Compile

Please refer to the corresponding [page](https://github.com/TLCFEM/suanPan-manual/blob/dev/docs/Tutorial/Compile.md) in
manual for details.

## Happy Modelling

![an example of simulation of particle collision](Resource/particle-collision.gif)

## Licence

[![FOSSA Status](https://app.fossa.com/api/projects/git%2Bgithub.com%2FTLCFEM%2FsuanPan.svg?type=large)](https://app.fossa.com/projects/git%2Bgithub.com%2FTLCFEM%2FsuanPan?ref=badge_large)
