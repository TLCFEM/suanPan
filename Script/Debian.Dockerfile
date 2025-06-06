FROM debian:12 AS build

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y locales && rm -rf /var/lib/apt/lists/* \
    && localedef -i en_US -c -f UTF-8 -A /usr/share/locale/locale.alias en_US.UTF-8

ENV LANG=en_US.utf8

RUN apt-get update -y && apt-get install -y wget gnupg

RUN wget -q https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB
RUN apt-key add GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB
RUN rm GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB
RUN echo "deb https://apt.repos.intel.com/oneapi all main" | tee /etc/apt/sources.list.d/oneAPI.list

RUN apt-get update -y && apt-get install -y gcc g++ gfortran cmake git intel-oneapi-mkl-devel libxt-dev freeglut3-dev libxcursor-dev file dpkg-dev

RUN mkdir vtk-build && cd vtk-build && \
    wget -q https://www.vtk.org/files/release/9.4/VTK-9.4.2.tar.gz && tar xf VTK-9.4.2.tar.gz && \
    cmake -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=OFF ./VTK-9.4.2 && \
    make install -j"$(nproc)" && cd .. && rm -r vtk-build

RUN git clone --recurse-submodules -b dev --depth 1 https://github.com/TLCFEM/suanPan.git && cd suanPan && mkdir build && cd build && \
    cmake -DCMAKE_BUILD_TYPE=Release -DSP_BUILD_PARALLEL=ON -DSP_ENABLE_HDF5=ON -DSP_ENABLE_VTK=ON -DSP_ENABLE_MKL=ON -DSP_ENABLE_IOMP=OFF -DSP_ENABLE_SHARED_MKL=OFF -DBUILD_PACKAGE=DEB .. && \
    make package -j"$(nproc)" && cp suanPan*.deb / && \
    cd / && rm -r suanPan

FROM debian:12

COPY --from=build /suanPan*.deb /

RUN apt-get install ./suanPan*.deb -y && rm ./suanPan*.deb

RUN ln -s /usr/bin/suanPan /usr/bin/suanpan
RUN ln -s /usr/bin/suanPan /usr/bin/sp

RUN sp -v
