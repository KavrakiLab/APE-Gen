FROM ubuntu:bionic
ENV DEBIAN_FRONTEND=noninteractive

# The following line is important for singularity
#RUN mkdir /scratch /work /home1 /gpfs /corral-repl /corral-tacc /data

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
        wget \
        vim


RUN ["/bin/bash", "-c", "wget http://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh"]
RUN chmod 0755 Miniconda3-latest-Linux-x86_64.sh
RUN ./Miniconda3-latest-Linux-x86_64.sh -b -p /conda
ENV PATH "/conda/bin:$PATH"

# environment.yml has pdbfixer, numpy, mdtraj and openmm
ADD environment.yml /environment.yml
RUN conda env create -f /environment.yml
RUN echo "source activate apegen" >> ~/.bashrc
ENV PATH "/conda/envs/apegen/bin:$PATH"

# These are okay outside environment.yml because they install to PATH
# Perhaps better if inside yml?
RUN conda install -c bioconda smina
RUN conda install -c schrodinger pymol
RUN conda install -c bioconda autodock-vina
RUN conda install -c salilab modeller
RUN conda install -c conda-forge biopython

# RCD
RUN conda install -c intel mkl
ENV LD_LIBRARY_PATH "$LD_LIBRARY_PATH:/conda/lib/"
RUN wget -qO- http://chaconlab.org`wget -qO- http://chaconlab.org/modeling/rcd/rcd-download|grep txz|grep uk-button| cut -f4 -d\"` | tar Jxf - && \
    cp RCD_v1.40_Linux_20190228/bin/rcd /usr/local/bin

RUN git clone https://github.com/KavrakiLab/APE-Gen.git

ENTRYPOINT ["bash"]
