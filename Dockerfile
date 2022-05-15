FROM ubuntu:focal

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y locales && rm -rf /var/lib/apt/lists/* \
    && localedef -i en_US -c -f UTF-8 -A /usr/share/locale/locale.alias en_US.UTF-8

ENV LANG en_US.utf8

RUN apt-get update -y && apt-get upgrade -y && apt-get install -y wget gnupg

RUN wget https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB
RUN apt-key add GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB
RUN echo "deb https://apt.repos.intel.com/oneapi all main" | tee /etc/apt/sources.list.d/oneAPI.list
RUN rm GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB

RUN apt-get update -y && apt-get upgrade -y && apt-get install -y gcc g++ gfortran cmake git intel-oneapi-mkl-devel-2021.2.0 libxt-dev freeglut3-dev libxcursor-dev

RUN mkdir vtk-build && cd vtk-build && \
    wget https://www.vtk.org/files/release/9.1/VTK-9.1.0.tar.gz && tar xf VTK-9.1.0.tar.gz && \
    cmake -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=OFF ./VTK-9.1.0 && \
    make install -j8 && \
    cd .. && rm -r vtk-build

RUN git clone -b master https://github.com/TLCFEM/suanPan.git && \
    cd suanPan && mkdir build && cd build && \
    cmake -DCMAKE_BUILD_TYPE=Release -DBUILD_MULTITHREAD=ON -DUSE_HDF5=ON -DUSE_EXTERNAL_VTK=ON -DUSE_MKL=ON -DMKLROOT=/opt/intel/oneapi/mkl/latest/ -DUSE_INTEL_OPENMP=OFF -DLINK_DYNAMIC_MKL=OFF -DCMAKE_INSTALL_PREFIX=suanPan-linux-mkl-vtk .. && \
    make install -j8 && \
    tar czf /suanPan-linux-mkl-vtk.tar.gz suanPan-linux-mkl-vtk && \
    cd suanPan-linux-mkl-vtk/bin && ./suanPan.sh -v && \
    cd / && ls -al && rm -r suanPan