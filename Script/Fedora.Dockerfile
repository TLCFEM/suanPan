FROM fedora:41 AS build

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

RUN git clone --recurse-submodules -b dev --depth 1 https://github.com/TLCFEM/suanPan.git && \
    cd suanPan && mkdir build && cd build && cmake -DCMAKE_BUILD_TYPE=Release -DSP_BUILD_PARALLEL=ON -DSP_ENABLE_HDF5=ON -DSP_ENABLE_VTK=ON -DSP_ENABLE_MKL=ON -DSP_ENABLE_IOMP=OFF -DSP_ENABLE_SHARED_MKL=OFF -DBUILD_PACKAGE=RPM .. && \
    make package -j"$(nproc)" && cp suanPan*.rpm / && \
    cd / && rm -r suanPan

FROM fedora:41

COPY --from=build /suanPan*.rpm /

RUN dnf install ./suanPan*.rpm -y && rm ./suanPan*.rpm

RUN ln -s /usr/bin/suanPan /usr/bin/suanpan
RUN ln -s /usr/bin/suanPan /usr/bin/sp

RUN sp -v
