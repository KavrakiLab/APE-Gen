#FROM continuumio/miniconda3
#RUN conda create -n env python=3
#RUN echo "source activate env" >> ~/.bashrc
#ENV PATH /opt/conda/envs/env/bin:$PATH

FROM centos:6

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
        wget \
        vim


RUN ["/bin/bash", "-c", "wget http://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh -O /home/jra8/miniconda.sh"]
RUN chmod 0755 /home/jra8/miniconda.sh
RUN ["/bin/bash", "-c", "/home/jra8/miniconda.sh -b -p /home/jra8/conda"]
ENV PATH="/home/jra8/conda/bin:$PATH"
RUN mkdir -p /home/jra8/singularity
WORKDIR /home/jra8/singularity

# environment.yml has pdbfixer, numpy, mdtraj and openmm
ADD environment.yml /home/jra8/singularity/environment.yml
RUN conda env create -f /home/jra8/singularity/environment.yml
RUN echo "source activate $(head -1 /home/jra8/singularity/environment.yml | cut -d' ' -f2)" > ~/.bashrc
ENV PATH /opt/conda/envs/$(head -1 /home/jra8/singularity/environment.yml | cut -d' ' -f2)/bin:$PATH

# These are okay outside environment.yml because they install to PATH
# Perhaps better if inside yml?
RUN conda install -c bioconda smina
RUN conda install -c schrodinger pymol

# To call APE-Gen as module
RUN echo "export PYTHONPATH=$PYTHONPATH:/home/jra8/singularity" >> ~/.bashrc

# To get nglview to work inside jupyter notebook
#RUN echo "jupyter-nbextension enable --py --sys-prefix widgetsnbextension" >> ~/.bashrc
#RUN echo "jupyter-nbextension enable nglview --py --sys-prefix" >> ~/.bashrc

# RCD
RUN conda install -c intel mkl
RUN echo "export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/conda/lib/" >> ~/.bashrc
RUN wget -qO- http://chaconlab.org`wget -qO- http://chaconlab.org/modeling/rcd/rcd-download|grep txz|grep uk-button| cut -f4 -d\"` | tar Jxf - && \
    cp RCD_v1.40_Linux_20190228/bin/rcd /usr/local/bin
RUN echo "export PATH=/home/jra8/singularity/RCD_v1.40_Linux_20190228/bin:$PATH"

# Autodock Vina
RUN wget -qO- http://vina.scripps.edu/download/autodock_vina_1_1_2_linux_x86.tgz | tar zxf - 
RUN echo "export PATH=/home/jra8/singularity/autodock_vina_1_1_2_linux_x86/bin:$PATH"

RUN git clone https://github.com/KavrakiLab/APE-Gen.git

# APE-Gen
#WORKDIR /home/apegen
COPY . /home/jra8/singularity/
ENTRYPOINT ["bash"]
