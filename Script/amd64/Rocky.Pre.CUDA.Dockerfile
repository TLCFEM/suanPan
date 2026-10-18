# amd64/Rocky.Pre.CUDA.Dockerfile - Build dependencies for CUDA-enabled suanPan
# 
# Extends Rocky.Pre.Dockerfile with NVIDIA CUDA toolkit and cuDNN support
# for GPU-accelerated computations. Required for CUDA-based solvers.
# 
# Build: docker build -f amd64/Rocky.Pre.CUDA.Dockerfile -t tlcfem/suanpan-env:cuda .

FROM nvidia/cuda:13.3.1-devel-rockylinux9

RUN dnf install -y dnf-plugins-core && \
    dnf config-manager --enable crb && \
    dnf install -y libglvnd-devel gfortran rpm-build rpm-devel rpmdevtools cmake wget git sudo && \
    dnf clean all

# changes of mkl sparse utility names force to install an older version
RUN wget -q https://registrationcenter-download.intel.com/akdlm/IRC_NAS/6a17080f-f0de-41b9-b587-52f92512c59a/intel-onemkl-2025.3.1.11_offline.sh && \
    bash ./intel-onemkl-2025.3.1.11_offline.sh -a --silent --eula accept && \
    rm intel-onemkl-2025.3.1.11_offline.sh

RUN mkdir magma-build && cd magma-build && \
    wget -q https://icl.utk.edu/projectsfiles/magma/downloads/magma-2.10.0.tar.gz && tar xf magma-2.10.0.tar.gz && \
    source /opt/intel/oneapi/setvars.sh && cmake -DCMAKE_BUILD_TYPE=Release -DGPU_TARGET="Turing Ampere Hopper" -DBUILD_SHARED_LIBS=OFF -DFORTRAN_CONVENTION="-DADD_" -DBLA_STATIC=ON -DBLA_VENDOR="Intel10_64lp" ./magma-2.10.0 && \
    make install -j"$(nproc)" && cd .. && rm -r magma-build

RUN mkdir vtk-build && cd vtk-build && \
    wget -q https://www.vtk.org/files/release/9.7/VTK-9.7.1.tar.gz && tar xf VTK-9.7.1.tar.gz && \
    cmake -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=OFF ./VTK-9.7.1 && \
    make install -j"$(nproc)" && cd .. && rm -r vtk-build

ARG USERNAME=nonroot
ARG USER_UID=1000
ARG USER_GID=$USER_UID

RUN groupadd --gid $USER_GID $USERNAME \
    && useradd --uid $USER_UID --gid $USER_GID -m $USERNAME \
    && echo $USERNAME ALL=\(root\) NOPASSWD:ALL > /etc/sudoers.d/$USERNAME \
    && chmod 0440 /etc/sudoers.d/$USERNAME

USER $USERNAME
