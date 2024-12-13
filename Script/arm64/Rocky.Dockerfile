FROM tlcfem/suanpan-env:latest AS build

# part 3: suanpan
RUN git clone --depth 1 --branch dev https://github.com/TLCFEM/suanPan.git

RUN cd suanPan && mkdir build && cd build && cmake -DCMAKE_BUILD_TYPE=Release -DUSE_SYS_LIB=ON -DCUSTOM_OPENBLAS=/OpenBLAS -DBUILD_MULTITHREAD=ON -DUSE_HDF5=ON -DUSE_VTK=ON -DUSE_AVX2=OFF .. && make -j"$(nproc)"

FROM rockylinux:9 AS runtime

RUN dnf install -y epel-release && crb enable
RUN dnf install -y libgomp hdf5 libX11

COPY --from=build /tbb-build /tbb-build

RUN find /tbb-build -name "libtbb*.so*" -exec cp {} /usr/local/lib64 \;

RUN rm -rf tbb-build

COPY --from=build /suanPan/build/suanPan /usr/local/bin/suanPan

RUN ln -s /usr/local/bin/suanPan /usr/local/bin/suanpan
RUN ln -s /usr/local/bin/suanPan /usr/local/bin/sp

RUN sp -v
