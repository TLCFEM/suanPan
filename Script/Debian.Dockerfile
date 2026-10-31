# Debian.Dockerfile - Build suanPan on Debian 12
# 
# Compiles suanPan with Intel oneAPI MKL and creates a Debian DEB package.
# Includes VTK for visualization support.
# 
# Build: docker build -f Debian.Dockerfile -t tlcfem/suanpan-debian:latest .
FROM debian:12 AS build

RUN apt-get update -y && apt-get install -y --no-install-recommends wget gnupg ca-certificates gcc g++ gfortran cmake git libxt-dev freeglut3-dev libxcursor-dev file dpkg-dev && rm -rf /var/lib/apt/lists/*

RUN wget -q https://registrationcenter-download.intel.com/akdlm/IRC_NAS/97415287-271a-4d14-967f-72e2d9ceeb1e/intel-onemkl-2026.2.0.890_offline.sh && sh ./intel-onemkl-2026.2.0.890_offline.sh -a --silent --eula accept && rm ./intel-onemkl-2026.2.0.890_offline.sh

RUN mkdir vtk-build && cd vtk-build && \
    wget -q https://www.vtk.org/files/release/9.7/VTK-9.7.1.tar.gz && tar xf VTK-9.7.1.tar.gz && \
    cmake -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=OFF ./VTK-9.7.1 && \
    make install -j"$(nproc)" && cd .. && rm -r vtk-build

RUN git clone --recurse-submodules -b dev --depth 1 https://github.com/TLCFEM/suanPan.git && cd suanPan && mkdir build && cd build && \
    cmake \
    -DBUILD_PACKAGE=DEB \
    -DCMAKE_BUILD_TYPE=Release \
    -DSP_BUILD_PARALLEL=ON \
    -DSP_ENABLE_AVX2=OFF \
    -DSP_ENABLE_HDF5=ON \
    -DSP_ENABLE_IOMP=OFF \
    -DSP_ENABLE_MKL=ON \
    -DSP_ENABLE_SHARED_MKL=OFF \
    -DSP_ENABLE_VTK=ON \
    -DCMAKE_INSTALL_PREFIX=/suanPan-debian12 \
    .. && \
    make package -j"$(nproc)" && cp suanPan*.deb / && \
    tar -C /suanPan-debian12 -czf /suanPan-debian12.tar.gz . \
    cd / && rm -r suanPan

FROM debian:12

COPY --from=build /suanPan*.deb /
COPY --from=build /suanPan-debian12.tar.gz /

RUN apt-get update -y && apt-get install -y --no-install-recommends ./suanPan*.deb && rm ./suanPan*.deb && rm -rf /var/lib/apt/lists/*

RUN ln -s /usr/bin/suanPan /usr/bin/suanpan
RUN ln -s /usr/bin/suanPan /usr/bin/sp

RUN sp -v
