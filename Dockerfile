FROM ubuntu:bionic

ENV DEBIAN_FRONTEND=noninteractive

# use older version of Boost to avoid compilation problems with smina
RUN apt-get update && \
    apt-get -y install \
        cmake \
        cython \
        doxygen \
        freeglut3-dev \
        git \
        gfortran \
        libboost-date-time1.62-dev \
        libboost-filesystem1.62-dev \
        libboost-iostreams1.62-dev \
        libboost-program-options1.62-dev \
        libboost-serialization1.62-dev \
        libboost-system1.62-dev \
        libboost-thread1.62-dev \
        libboost-timer1.62-dev \
        libeigen3-dev \
        libfftw3-dev \
        libfreetype6-dev \
        libglew-dev \
        libglm-dev \
        libpng-dev \
        libpython-dev \
        libxml2-dev \
        python-pip \
        python-pmw \
        swig \
        wget

# Intel MKL; see https://github.com/eddelbuettel/mkl4deb
RUN wget -O- https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS-2019.PUB | apt-key add - && \
    sh -c 'echo deb https://apt.repos.intel.com/mkl all main > /etc/apt/sources.list.d/intel-mkl.list' && \
    apt-get update && \
    apt-get -y install intel-mkl-64bit-2018.2-046 && \
    update-alternatives --install /usr/lib/x86_64-linux-gnu/libblas.so     libblas.so-x86_64-linux-gnu      /opt/intel/mkl/lib/intel64/libmkl_rt.so 50 && \
    update-alternatives --install /usr/lib/x86_64-linux-gnu/libblas.so.3   libblas.so.3-x86_64-linux-gnu    /opt/intel/mkl/lib/intel64/libmkl_rt.so 50 && \
    update-alternatives --install /usr/lib/x86_64-linux-gnu/liblapack.so   liblapack.so-x86_64-linux-gnu    /opt/intel/mkl/lib/intel64/libmkl_rt.so 50 && \
    update-alternatives --install /usr/lib/x86_64-linux-gnu/liblapack.so.3 liblapack.so.3-x86_64-linux-gnu  /opt/intel/mkl/lib/intel64/libmkl_rt.so 50 && \
    echo "/opt/intel/lib/intel64"     >  /etc/ld.so.conf.d/mkl.conf && \
    echo "/opt/intel/mkl/lib/intel64" >> /etc/ld.so.conf.d/mkl.conf && \
    ldconfig && \
    echo "MKL_THREADING_LAYER=GNU" >> /etc/environment

# Pip install some dependencies
RUN pip install intel-scipy intel-numpy pandas

# OpenMM (no CUDA or OpenCL support enabled for now)
RUN git clone --depth 1 https://github.com/pandegroup/openmm.git && \
    mkdir openmm/build && \
    cd openmm/build && \
    cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/usr .. && \
    make -j`nproc` && \
    make install && \
    make PythonInstall && \
    cd / && \
    rm -rf openmm

# PDBFixer
RUN git clone --depth 1 https://github.com/pandegroup/pdbfixer.git && \
    cd pdbfixer && \
    #grep -v openmm setup.py > setup-fixed.py && \
    python setup.py install && \
    cd / && \
    rm -rf pdbfixer

# MDtraj
RUN pip install mdtraj

# OpenBabel; need unreleased version for SMINA!
RUN git clone --depth 1 https://github.com/openbabel/openbabel.git && \
    mkdir openbabel/build && \
    cd openbabel/build && \
    cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/usr -DPYTHON_BINDINGS=ON -DRUN_SWIG=ON .. && \
    make -j`nproc` && \
    make install && \
    cd / && \
    rm -rf openbabel

# SMINA
RUN git clone --depth 1 https://github.com/mwojcikowski/smina.git && \
    cd smina/build/linux/release && \
    make OPENBABEL_INCLUDE=/usr/include/openbabel-2.0 -j`nproc` && \
    cp smina server tosmina /usr/local/bin && \
    cd / && \
    rm -rf smina

# RCD
RUN wget -qO- http://chaconlab.org`wget -qO- http://chaconlab.org/modeling/rcd/rcd-download|grep txz|grep uk-button| cut -f4 -d\"` | tar Jxf - && \
    cp RCD_v1.40_Linux_20190228/bin/rcd /usr/local/bin && \
    rm -rf RCD_v1.40_Linux_20190228

# Pymol
RUN git clone --depth 1 https://github.com/schrodinger/pymol-open-source.git && \
    cd pymol-open-source && \
    python setup.py install --use-msgpackc=no --prefix=/root/pymol && \
    cd / && \
    rm -rf pymol-open-source

# Autodock Vina
RUN wget -qO- http://vina.scripps.edu/download/autodock_vina_1_1_2_linux_x86.tgz | tar zxf - && \
    mv autodock_vina_1_1_2_linux_x86/bin/vina* /usr/local/bin/ && \
    rm -rf /autodock_vina_1_1_2_linux_x86

# APE-Gen
WORKDIR /home/apegen
COPY . /home/apegen
ENTRYPOINT ["bash"]
