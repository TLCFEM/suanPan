FROM nvidia/cuda:12.5.0-devel-rockylinux9

RUN dnf install -y dnf-plugins-core
RUN dnf config-manager --enable crb
RUN dnf install -y libglvnd-devel gfortran rpm-build rpm-devel rpmdevtools cmake wget git ninja-build

RUN wget -q https://registrationcenter-download.intel.com/akdlm/IRC_NAS/cdff21a5-6ac7-4b41-a7ec-351b5f9ce8fd/l_onemkl_p_2024.2.0.664_offline.sh && \
    sh ./l_onemkl_p_2024.2.0.664_offline.sh -f intel_tmp -a --silent --eula accept && rm -r intel_tmp && rm ./l_onemkl_p_2024.2.0.664_offline.sh

RUN mkdir vtk-build && cd vtk-build && \
    wget -q https://www.vtk.org/files/release/9.2/VTK-9.2.6.tar.gz && tar xf VTK-9.2.6.tar.gz && \
    cmake -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=OFF ./VTK-9.2.6 && \
    make install -j"$(nproc)" && cd .. && rm -r vtk-build

RUN useradd -u 1000 nonroot
USER nonroot
