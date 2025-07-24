# This Dockerfile is creating a base image for compiling suanPan.
# It preinstalls the MKL library, which is required as a fast linear algebra driver that does the heavy lifting.
# It installs the VTK library, which is required to generate VTK files for visualization.

FROM rockylinux:9

RUN echo "[oneAPI]" > /etc/yum.repos.d/oneAPI.repo && \
    echo "name=Intel oneAPI repository" >> /etc/yum.repos.d/oneAPI.repo && \
    echo "baseurl=https://yum.repos.intel.com/oneapi" >> /etc/yum.repos.d/oneAPI.repo && \
    echo "enabled=1" >> /etc/yum.repos.d/oneAPI.repo && \
    echo "gpgcheck=1" >> /etc/yum.repos.d/oneAPI.repo && \
    echo "repo_gpgcheck=1" >> /etc/yum.repos.d/oneAPI.repo && \
    echo "gpgkey=https://yum.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB" >> /etc/yum.repos.d/oneAPI.repo

RUN dnf install -y libglvnd-devel gcc g++ gfortran rpm-build rpm-devel rpmdevtools cmake wget git intel-oneapi-mkl-devel

RUN mkdir vtk-build && cd vtk-build && \
    wget -q https://www.vtk.org/files/release/9.5/VTK-9.5.0.tar.gz && tar xf VTK-9.5.0.tar.gz && \
    cmake -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=OFF ./VTK-9.5.0 && \
    make install -j"$(nproc)" && cd .. && rm -r vtk-build

# ARG USERNAME=nonroot
# ARG USER_UID=1000
# ARG USER_GID=$USER_UID

# RUN groupadd --gid $USER_GID $USERNAME \
#     && useradd --uid $USER_UID --gid $USER_GID -m $USERNAME \
#     && dnf install -y sudo \
#     && echo $USERNAME ALL=\(root\) NOPASSWD:ALL > /etc/sudoers.d/$USERNAME \
#     && chmod 0440 /etc/sudoers.d/$USERNAME

# USER $USERNAME
