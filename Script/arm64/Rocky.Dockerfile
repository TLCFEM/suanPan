FROM tlcfem/suanpan-env:latest AS build

# part 3: suanpan
RUN git clone --recurse-submodules --depth 1 --branch dev https://github.com/TLCFEM/suanPan.git && cd suanPan && mkdir build && cd build && \
    cmake -DCMAKE_BUILD_TYPE=Release -DSP_ENABLE_SYSLIB=ON -DSP_OPENBLAS_PATH=/OpenBLAS -DSP_BUILD_PARALLEL=ON -DSP_ENABLE_HDF5=ON -DSP_ENABLE_VTK=ON -DSP_ENABLE_AVX2=OFF .. && make -j"$(nproc)"

FROM rockylinux:9

RUN dnf install -y epel-release && crb enable && dnf install -y libgomp hdf5 libX11

COPY --from=build /tbb-build /tbb-build

RUN find /tbb-build -name "libtbb*.so*" -exec cp {} /usr/local/lib64 \; && rm -r tbb-build

COPY --from=build /suanPan/build/suanPan /usr/local/bin/suanPan

RUN ln -s /usr/local/bin/suanPan /usr/local/bin/suanpan
RUN ln -s /usr/local/bin/suanPan /usr/local/bin/sp

RUN sp -v
