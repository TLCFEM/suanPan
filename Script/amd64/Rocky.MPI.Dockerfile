FROM tlcfem/suanpan-env:latest AS build

RUN source /opt/intel/oneapi/setvars.sh && \
    git clone --recurse-submodules --branch dev --depth 1 https://github.com/TLCFEM/suanPan.git && cd suanPan && mkdir build && cd build && \
    cmake -DCMAKE_BUILD_TYPE=Release -DSP_BUILD_PARALLEL=ON -DSP_ENABLE_HDF5=ON -DSP_ENABLE_VTK=ON -DSP_ENABLE_AVX2=OFF -DSP_ENABLE_MKL=ON -DSP_ENABLE_IOMP=OFF -DSP_ENABLE_SHARED_MKL=OFF -DSP_ENABLE_64BIT_INDEXING=ON -DSP_ENABLE_MPI=ON -DBUILD_PACKAGE=RPM ..

ENV QA_RPATHS=0x0002

RUN source /opt/intel/oneapi/setvars.sh && \
    cd suanPan/build && \
    make package -j"$(nproc)" && cp suanPan*.rpm / && \
    cd / && rm -r suanPan

WORKDIR /tester

RUN echo -e '#include <mpi.h>\n#include <stdio.h>\n\nint main(int argc, char** argv) {\n    MPI_Init(&argc, &argv);\n    int world_size, world_rank;\n    MPI_Comm_size(MPI_COMM_WORLD, &world_size);\n    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);\n    printf("Hello world from rank %d out of %d processors\\n", world_rank, world_size);\n    MPI_Finalize();\n    return 0;\n}' > mpi_tester.c

RUN source /opt/intel/oneapi/setvars.sh && mpicc mpi_tester.c -o mpi_tester

FROM almalinux:10

RUN echo "[oneAPI]" > /etc/yum.repos.d/oneAPI.repo && \
    echo "name=Intel oneAPI repository" >> /etc/yum.repos.d/oneAPI.repo && \
    echo "baseurl=https://yum.repos.intel.com/oneapi" >> /etc/yum.repos.d/oneAPI.repo && \
    echo "enabled=1" >> /etc/yum.repos.d/oneAPI.repo && \
    echo "gpgcheck=1" >> /etc/yum.repos.d/oneAPI.repo && \
    echo "repo_gpgcheck=1" >> /etc/yum.repos.d/oneAPI.repo && \
    echo "gpgkey=https://yum.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB" >> /etc/yum.repos.d/oneAPI.repo

RUN dnf install -y intel-oneapi-mpi procps-ng openssh-server openssh-clients && dnf clean all

COPY --from=build /suanPan*.rpm /

RUN dnf install ./suanPan*.rpm -y && rm ./suanPan*.rpm

RUN ln -s /usr/bin/suanPan /usr/bin/suanpan
RUN ln -s /usr/bin/suanPan /usr/bin/sp

COPY --from=build /tester/mpi_tester /usr/bin/mpi_tester

RUN ssh-keygen -A

# RUN useradd -m runner
# USER runner

RUN echo "source /opt/intel/oneapi/setvars.sh > /dev/null 2>&1" >> ~/.bashrc

RUN ssh-keygen -t rsa -f ~/.ssh/id_rsa -N "" && \
    cp ~/.ssh/id_rsa.pub ~/.ssh/authorized_keys && \
    echo -e "Host *\n    StrictHostKeyChecking no\n    UserKnownHostsFile /dev/null" > ~/.ssh/config

VOLUME ["/model"]

WORKDIR /model

CMD ["/usr/sbin/sshd", "-D"]
