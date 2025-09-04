FROM tlcfem/suanpan-env:latest AS build

RUN dnf install -y intel-oneapi-mpi-devel procps-ng

RUN source /opt/intel/oneapi/setvars.sh && \
    git clone --recurse-submodules --branch dev --depth 1 https://github.com/TLCFEM/suanPan.git && cd suanPan && mkdir build && cd build && \
    cmake -DCMAKE_BUILD_TYPE=Release -DSP_BUILD_PARALLEL=ON -DSP_ENABLE_HDF5=ON -DSP_ENABLE_VTK=ON -DSP_ENABLE_AVX2=OFF -DSP_ENABLE_MKL=ON -DSP_ENABLE_IOMP=OFF -DSP_ENABLE_SHARED_MKL=OFF -DSP_ENABLE_64BIT_INDEXING=ON -DSP_ENABLE_MPI=ON -DBUILD_PACKAGE=RPM .. && make package -j"$(nproc)" && cp suanPan*.rpm / && \
    cd / && rm -r suanPan

FROM rockylinux:9

RUN echo "[oneAPI]" > /etc/yum.repos.d/oneAPI.repo && \
    echo "name=Intel oneAPI repository" >> /etc/yum.repos.d/oneAPI.repo && \
    echo "baseurl=https://yum.repos.intel.com/oneapi" >> /etc/yum.repos.d/oneAPI.repo && \
    echo "enabled=1" >> /etc/yum.repos.d/oneAPI.repo && \
    echo "gpgcheck=1" >> /etc/yum.repos.d/oneAPI.repo && \
    echo "repo_gpgcheck=1" >> /etc/yum.repos.d/oneAPI.repo && \
    echo "gpgkey=https://yum.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB" >> /etc/yum.repos.d/oneAPI.repo

COPY --from=build /suanPan*.rpm /

RUN dnf install ./suanPan*.rpm intel-oneapi-mpi procps-ng -y && rm ./suanPan*.rpm

RUN ln -s /usr/bin/suanPan /usr/bin/suanpan
RUN ln -s /usr/bin/suanPan /usr/bin/sp
