# This Dockerfile is creating a base image for compiling suanPan.
# It is based on the official NVIDIA CUDA image.
# Compared to the Rocky.Pre.Dockerfile, it adds the CUDA toolkit.
# If you need the GPU related solvers, you should use this Dockerfile.

FROM nvidia/cuda:12.5.1-devel-rockylinux9

RUN dnf install -y dnf-plugins-core
RUN dnf config-manager --enable crb
RUN dnf install -y libglvnd-devel gfortran rpm-build rpm-devel rpmdevtools cmake wget git ninja-build

RUN wget -q https://registrationcenter-download.intel.com/akdlm/IRC_NAS/6e00e368-b61d-4f87-a409-9b510c022a37/l_onemkl_p_2024.2.1.105_offline.sh && \
    sh ./l_onemkl_p_2024.2.1.105_offline.sh -f intel_tmp -a --silent --eula accept && rm -r intel_tmp && rm ./l_onemkl_p_2024.2.1.105_offline.sh

RUN mkdir vtk-build && cd vtk-build && \
    wget -q https://www.vtk.org/files/release/9.2/VTK-9.2.6.tar.gz && tar xf VTK-9.2.6.tar.gz && \
    cmake -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=OFF ./VTK-9.2.6 && \
    make install -j"$(nproc)" && cd .. && rm -r vtk-build

RUN dnf install -y procps
RUN mkdir magma-build && cd magma-build && \
    wget -q https://icl.utk.edu/projectsfiles/magma/downloads/magma-2.8.0.tar.gz && tar xf magma-2.8.0.tar.gz && \
    source /opt/intel/oneapi/setvars.sh && cmake -DCMAKE_BUILD_TYPE=Release -DGPU_TARGET="Turing Ampere Hopper" -DBUILD_SHARED_LIBS=OFF -DFORTRAN_CONVENTION="-DADD_" -DBLA_STATIC=ON -DBLA_VENDOR="Intel10_64lp" ./magma-2.8.0 && \
    make install -j"$(nproc)" && cd .. && rm -r magma-build

ARG USERNAME=nonroot
ARG USER_UID=1000
ARG USER_GID=$USER_UID

RUN groupadd --gid $USER_GID $USERNAME \
    && useradd --uid $USER_UID --gid $USER_GID -m $USERNAME \
    && dnf install -y sudo \
    && echo $USERNAME ALL=\(root\) NOPASSWD:ALL > /etc/sudoers.d/$USERNAME \
    && chmod 0440 /etc/sudoers.d/$USERNAME

USER $USERNAME