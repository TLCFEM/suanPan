FROM fedora:40 AS build

RUN echo "[oneAPI]" > /etc/yum.repos.d/oneAPI.repo && \
    echo "name=Intel oneAPI repository" >> /etc/yum.repos.d/oneAPI.repo && \
    echo "baseurl=https://yum.repos.intel.com/oneapi" >> /etc/yum.repos.d/oneAPI.repo && \
    echo "enabled=1" >> /etc/yum.repos.d/oneAPI.repo && \
    echo "gpgcheck=1" >> /etc/yum.repos.d/oneAPI.repo && \
    echo "repo_gpgcheck=1" >> /etc/yum.repos.d/oneAPI.repo && \
    echo "gpgkey=https://yum.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB" >> /etc/yum.repos.d/oneAPI.repo

RUN dnf upgrade --refresh -y && dnf install -y libglvnd-devel gcc g++ gfortran rpm-build rpm-devel rpmdevtools cmake wget git intel-oneapi-mkl-devel

RUN mkdir vtk-build && cd vtk-build && \
    wget -q https://www.vtk.org/files/release/9.2/VTK-9.2.6.tar.gz && tar xf VTK-9.2.6.tar.gz && \
    cmake -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=OFF ./VTK-9.2.6 && \
    make install -j"$(nproc)" && cd .. && rm -r vtk-build

RUN git clone -b dev --depth 1 https://github.com/TLCFEM/suanPan.git
RUN cd suanPan && mkdir build && cd build && cmake -DCMAKE_BUILD_TYPE=Release -DBUILD_MULTITHREAD=ON -DUSE_HDF5=ON -DUSE_VTK=ON -DUSE_MKL=ON -DMKLROOT=/opt/intel/oneapi/mkl/latest/ -DUSE_INTEL_OPENMP=OFF -DLINK_DYNAMIC_MKL=OFF -DCMAKE_INSTALL_PREFIX=suanPan-linux-mkl-vtk -DBUILD_PACKAGE=RPM ..
RUN cd suanPan/build && make install -j"$(nproc)" && make package
RUN cd suanPan/build && cp suanPan*.rpm / && \
    tar czf /suanPan-linux-mkl-vtk.tar.gz suanPan-linux-mkl-vtk && \
    cd suanPan-linux-mkl-vtk/bin && ./suanPan.sh -v && \
    cd / && ls -al && rm -r suanPan

FROM fedora:40 AS runtime

COPY --from=build /suanPan*.rpm /

RUN dnf install ./suanPan*.rpm -y && rm ./suanPan*.rpm

RUN ln -s /usr/bin/suanPan /usr/bin/suanpan
RUN ln -s /usr/bin/suanPan /usr/bin/sp

RUN sp -v
