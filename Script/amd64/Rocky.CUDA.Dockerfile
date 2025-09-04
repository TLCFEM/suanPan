FROM tlcfem/suanpan-env:cuda AS build

RUN git clone --recurse-submodules --branch dev --depth 1 https://github.com/TLCFEM/suanPan.git && cd suanPan && mkdir build && cd build && \
    cmake -DCMAKE_BUILD_TYPE=Release -DSP_BUILD_PARALLEL=ON -DSP_ENABLE_HDF5=ON -DSP_ENABLE_AVX2=OFF -DSP_ENABLE_VTK=ON -DSP_ENABLE_MKL=ON -DSP_ENABLE_IOMP=OFF -DSP_ENABLE_SHARED_MKL=OFF -DSP_ENABLE_CUDA=ON -DBUILD_PACKAGE=RPM .. && make package -j"$(nproc)" && cp suanPan*.rpm / && \
    cd / && rm -r suanPan

FROM nvidia/cuda:12.9.1-runtime-rockylinux9

COPY --from=build /suanPan*.rpm /

RUN dnf install ./suanPan*.rpm -y && rm ./suanPan*.rpm

RUN ln -s /usr/bin/suanPan /usr/bin/suanpan
RUN ln -s /usr/bin/suanPan /usr/bin/sp

RUN sp -v
