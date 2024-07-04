# This Dockerfile is creating a base image for compiling suanPan.
# It preinstalls the MKL library, which is required as a fast linear algebra driver that does the heavy lifting.
# It installs the VTK library, which is required to generate VTK files for visualization.

FROM rockylinux:9

RUN dnf upgrade --refresh -y && dnf install -y libglvnd-devel gcc g++ gfortran rpm-build rpm-devel rpmdevtools cmake wget git

RUN wget -q https://registrationcenter-download.intel.com/akdlm/IRC_NAS/cdff21a5-6ac7-4b41-a7ec-351b5f9ce8fd/l_onemkl_p_2024.2.0.664_offline.sh && \
    sh ./l_onemkl_p_2024.2.0.664_offline.sh -f intel_tmp -a --silent --eula accept && rm -r intel_tmp && rm ./l_onemkl_p_2024.2.0.664_offline.sh

RUN mkdir vtk-build && cd vtk-build && \
    wget -q https://www.vtk.org/files/release/9.2/VTK-9.2.6.tar.gz && tar xf VTK-9.2.6.tar.gz && \
    cmake -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=OFF ./VTK-9.2.6 && \
    make install -j"$(nproc)" && cd .. && rm -r vtk-build

RUN useradd -u 1000 nonroot
USER nonroot
