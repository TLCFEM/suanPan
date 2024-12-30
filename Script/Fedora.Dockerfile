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
    wget -q https://www.vtk.org/files/release/9.4/VTK-9.4.1.tar.gz && tar xf VTK-9.4.1.tar.gz && \
    cmake -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=OFF ./VTK-9.4.1 && \
    make install -j"$(nproc)" && cd .. && rm -r vtk-build

RUN git clone -b dev --depth 1 https://github.com/TLCFEM/suanPan.git && \
    cd suanPan && mkdir build && cd build && cmake -DCMAKE_BUILD_TYPE=Release -DBUILD_MULTITHREAD=ON -DUSE_HDF5=ON -DUSE_VTK=ON -DUSE_MKL=ON -DMKLROOT=/opt/intel/oneapi/mkl/latest/ -DUSE_INTEL_OPENMP=OFF -DLINK_DYNAMIC_MKL=OFF -DBUILD_PACKAGE=RPM .. && \
    make package -j"$(nproc)" && cp suanPan*.rpm / && \
    cd / && rm -r suanPan

FROM fedora:41

COPY --from=build /suanPan*.rpm /

RUN dnf install ./suanPan*.rpm -y && rm ./suanPan*.rpm

RUN ln -s /usr/bin/suanPan /usr/bin/suanpan
RUN ln -s /usr/bin/suanPan /usr/bin/sp

RUN sp -v
