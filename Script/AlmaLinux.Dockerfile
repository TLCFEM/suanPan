FROM almalinux:9 as build

RUN dnf upgrade --refresh -y && dnf install -y libglvnd-devel gcc g++ gfortran rpm-build rpm-devel rpmdevtools cmake wget git

RUN wget -q https://registrationcenter-download.intel.com/akdlm/irc_nas/19138/l_onemkl_p_2023.0.0.25398_offline.sh
RUN sh ./l_onemkl_p_2023.0.0.25398_offline.sh -a --silent --eula accept && rm ./l_onemkl_p_2023.0.0.25398_offline.sh

RUN mkdir vtk-build && cd vtk-build && \
    wget -q https://www.vtk.org/files/release/9.2/VTK-9.2.6.tar.gz && tar xf VTK-9.2.6.tar.gz && \
    cmake -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=OFF ./VTK-9.2.6 && \
    make install -j"$(nproc)" && cd .. && rm -r vtk-build

RUN git clone -b dev --depth 1 https://github.com/TLCFEM/suanPan.git
RUN cd suanPan && mkdir build && cd build && cmake -DCMAKE_BUILD_TYPE=Release -DBUILD_MULTITHREAD=ON -DUSE_HDF5=ON -DUSE_EXTERNAL_VTK=ON -DUSE_MKL=ON -DMKLROOT=/opt/intel/oneapi/mkl/latest/ -DUSE_INTEL_OPENMP=OFF -DLINK_DYNAMIC_MKL=OFF -DCMAKE_INSTALL_PREFIX=suanPan-linux-mkl-vtk -DBUILD_PACKAGE=RPM ..
RUN cd suanPan/build && make install -j"$(nproc)" && make package
RUN cd suanPan/build && cp suanPan*.rpm / && \
    tar czf /suanPan-linux-mkl-vtk.tar.gz suanPan-linux-mkl-vtk && \
    cd suanPan-linux-mkl-vtk/bin && ./suanPan.sh -v && \
    cd / && ls -al && rm -r suanPan

FROM almalinux:9 as runtime

COPY --from=build /suanPan*.rpm /

RUN dnf upgrade --refresh -y && dnf install ./suanPan*.rpm -y

RUN suanPan -v