# This Dockerfile is creating a base image for compiling suanPan.
# It is based on the official NVIDIA CUDA image.
# Compared to the Rocky.Pre.Dockerfile, it adds the CUDA toolkit.
# If you need the GPU related solvers, you should use this Dockerfile.

FROM nvidia/cuda:12.9.1-devel-rockylinux9

RUN echo "[oneAPI]" > /etc/yum.repos.d/oneAPI.repo && \
    echo "name=Intel oneAPI repository" >> /etc/yum.repos.d/oneAPI.repo && \
    echo "baseurl=https://yum.repos.intel.com/oneapi" >> /etc/yum.repos.d/oneAPI.repo && \
    echo "enabled=1" >> /etc/yum.repos.d/oneAPI.repo && \
    echo "gpgcheck=1" >> /etc/yum.repos.d/oneAPI.repo && \
    echo "repo_gpgcheck=1" >> /etc/yum.repos.d/oneAPI.repo && \
    echo "gpgkey=https://yum.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB" >> /etc/yum.repos.d/oneAPI.repo

RUN dnf install -y dnf-plugins-core && \
    dnf config-manager --enable crb && \
    dnf install -y libglvnd-devel gfortran rpm-build rpm-devel rpmdevtools cmake wget git ninja-build intel-oneapi-mkl-devel procps sudo

RUN mkdir vtk-build && cd vtk-build && \
    wget -q https://www.vtk.org/files/release/9.5/VTK-9.5.1.tar.gz && tar xf VTK-9.5.1.tar.gz && \
    cmake -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=OFF ./VTK-9.5.1 && \
    make install -j"$(nproc)" && cd .. && rm -r vtk-build

RUN mkdir magma-build && cd magma-build && \
    wget -q https://icl.utk.edu/projectsfiles/magma/downloads/magma-2.9.0.tar.gz && tar xf magma-2.9.0.tar.gz && \
    source /opt/intel/oneapi/setvars.sh && cmake -DCMAKE_BUILD_TYPE=Release -DGPU_TARGET="Turing Ampere Hopper" -DBUILD_SHARED_LIBS=OFF -DFORTRAN_CONVENTION="-DADD_" -DBLA_STATIC=ON -DBLA_VENDOR="Intel10_64lp" ./magma-2.9.0 && \
    make install -j"$(nproc)" && cd .. && rm -r magma-build

ARG USERNAME=nonroot
ARG USER_UID=1000
ARG USER_GID=$USER_UID

RUN groupadd --gid $USER_GID $USERNAME \
    && useradd --uid $USER_UID --gid $USER_GID -m $USERNAME \
    && echo $USERNAME ALL=\(root\) NOPASSWD:ALL > /etc/sudoers.d/$USERNAME \
    && chmod 0440 /etc/sudoers.d/$USERNAME

USER $USERNAME
