FROM rockylinux:9

RUN dnf install -y epel-release && crb enable
RUN dnf install -y gcc g++ gfortran cmake wget git hdf5-devel

# part 1: openblas
RUN git clone --depth 1 --branch v0.3.28 https://github.com/OpenMathLib/OpenBLAS.git

# change the following configurations to match your needs
RUN cd OpenBLAS && make TARGET=ARMV8 DYNAMIC_ARCH=1 BINARY=64 USE_THREAD=1 USE_OPENMP=1 NUM_THREADS=20 NO_SHARED=1 NO_CBLAS=1 NO_LAPACKE=1 GEMM_MULTITHREAD_THRESHOLD=64

# part 2: tbb
RUN git clone --depth 1 --branch v2021.12.0 https://github.com/oneapi-src/oneTBB.git
RUN mkdir tbb-build && cd tbb-build && cmake -DCMAKE_BUILD_TYPE=Release -DTBB_TEST=OFF ../oneTBB && cmake --build . --target install --config Release --parallel "$(nproc)"
