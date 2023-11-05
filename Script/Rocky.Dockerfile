FROM tlcfem/suanpan-env:latest as build

RUN git clone -b dev --depth 1 https://github.com/TLCFEM/suanPan.git
RUN cd suanPan && mkdir build && cd build && cmake -DCMAKE_BUILD_TYPE=Release -DBUILD_MULTITHREAD=ON -DUSE_HDF5=ON -DUSE_VTK=ON -DUSE_MKL=ON -DMKLROOT=/opt/intel/oneapi/mkl/latest/ -DUSE_INTEL_OPENMP=OFF -DLINK_DYNAMIC_MKL=OFF -DCMAKE_INSTALL_PREFIX=suanPan-linux-mkl-vtk -DBUILD_PACKAGE=RPM ..
RUN cd suanPan/build && make install -j"$(nproc)" && make package
RUN cd suanPan/build && cp suanPan*.rpm / && \
    tar czf /suanPan-linux-mkl-vtk.tar.gz suanPan-linux-mkl-vtk && \
    cd suanPan-linux-mkl-vtk/bin && ./suanPan.sh -v && \
    cd / && ls -al && rm -r suanPan

FROM rockylinux:9 as runtime

COPY --from=build /suanPan*.rpm /

RUN dnf install ./suanPan*.rpm -y

RUN suanPan -v