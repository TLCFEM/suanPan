FROM rockylinux:9

RUN dnf install -y epel-release && crb enable
RUN dnf install -y gcc g++ gfortran cmake wget git hdf5-devel libglvnd-devel

# part 1: openblas
# change the following configurations to match your needs
RUN git clone --depth 1 --branch v0.3.30 https://github.com/OpenMathLib/OpenBLAS.git openblas-build && cd openblas-build && \
    make TARGET=ARMV8 DYNAMIC_ARCH=1 BINARY=64 USE_THREAD=1 USE_OPENMP=1 NUM_THREADS=20 NO_SHARED=1 NO_CBLAS=1 NO_LAPACKE=1 GEMM_MULTITHREAD_THRESHOLD=64 && \
    cd .. && mkdir OpenBLAS && cp openblas-build/*.a OpenBLAS && rm -r openblas-build

# part 2: tbb
RUN mkdir tbb-build && cd tbb-build && \
    git clone --depth 1 --branch v2022.1.0 https://github.com/oneapi-src/oneTBB.git && \
    cmake -DCMAKE_BUILD_TYPE=Release -DTBB_TEST=OFF ./oneTBB && cmake --build . --target install --config Release --parallel "$(nproc)"

# part 3: vtk
RUN mkdir vtk-build && cd vtk-build && \
    wget -q https://www.vtk.org/files/release/9.5/VTK-9.5.1.tar.gz && tar xf VTK-9.5.1.tar.gz && \
    cmake -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=OFF ./VTK-9.5.1 && \
    cmake --build . --target install --config Release --parallel "$(nproc)" && \
    cd .. && rm -r vtk-build
