FROM nobodyxu/intel-mkl:2020-ubuntu-bionic

WORKDIR /root
RUN apt-get update && DEBIAN_FRONTEND="noninteractive" apt-get install -y \
    build-essential libssl-dev \
    wget \
    python3.8 python3-pip python3.8-distutils \
    git \
    htop \
    screen \
    tar

# install cmake
WORKDIR /root
RUN wget https://github.com/Kitware/CMake/releases/download/v3.20.2/cmake-3.20.2.tar.gz
RUN tar -zxvf cmake-3.20.2.tar.gz
WORKDIR /root/cmake-3.20.2
RUN ./bootstrap && make && make install

RUN wget https://bootstrap.pypa.io/get-pip.py && python3.8 get-pip.py

RUN ln -s /usr/bin/python3.8 /usr/bin/python && \
    ln -s /usr/bin/pip3 /usr/bin/pip && \
    pip install --upgrade pip

WORKDIR /root
# install libcint
RUN git clone https://github.com/sunqm/libcint.git
WORKDIR /root/libcint
RUN git checkout v4.0.7 && mkdir build
WORKDIR /root/libcint/build
RUN cmake -DWITH_F12=1 -DWITH_RANGE_COULOMB=1 -DWITH_COULOMB_ERF=1 \
    -DCMAKE_INSTALL_PREFIX:PATH=/opt -DCMAKE_INSTALL_LIBDIR:PATH=lib ..
RUN make && make install

# install libxc
WORKDIR /root
RUN git clone https://gitlab.com/libxc/libxc.git
WORKDIR /root/libxc
RUN git checkout 4.3.4 && mkdir -p build
WORKDIR /root/libxc/build
RUN cmake -DCMAKE_BUILD_TYPE=RELEASE -DBUILD_SHARED_LIBS=1 \
    -DENABLE_FORTRAN=0 -DDISABLE_KXC=0 -DDISABLE_LXC=1 \
    -DCMAKE_INSTALL_PREFIX:PATH=/opt -DCMAKE_INSTALL_LIBDIR:PATH=lib ..
RUN make && make install

# install xcfun
WORKDIR /root
RUN git clone https://github.com/sunqm/xcfun.git
WORKDIR /root/xcfun
RUN git checkout cmake-3.5 && mkdir build
WORKDIR /root/xcfun/build
RUN cmake -DCMAKE_BUILD_TYPE=RELEASE -DBUILD_SHARED_LIBS=1 -DXC_MAX_ORDER=3 -DXCFUN_ENABLE_TESTS=0 \
    -DCMAKE_INSTALL_PREFIX:PATH=/opt -DCMAKE_INSTALL_LIBDIR:PATH=lib ..
RUN make && make install


WORKDIR /root
ENV LD_LIBRARY_PATH /opt/intel/compilers_and_libraries/linux/mkl/lib/intel64:/opt/lib
RUN pip install numpy scipy pandas h5py==2.10.0