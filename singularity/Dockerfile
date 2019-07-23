#FROM continuumio/miniconda3
#RUN conda create -n env python=3
#RUN echo "source activate env" >> ~/.bashrc
#ENV PATH /opt/conda/envs/env/bin:$PATH

FROM ubuntu:bionic
ENV DEBIAN_FRONTEND=noninteractive

RUN mkdir /scratch /work /home1 /gpfs /corral-repl /corral-tacc /data

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


RUN ["/bin/bash", "-c", "wget http://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh"]
RUN chmod 0755 Miniconda3-latest-Linux-x86_64.sh
#RUN ["/bin/bash", "-c", "Miniconda-latest-Linux-x86_64.sh -b -p /home/jra8/conda"]
RUN ./Miniconda3-latest-Linux-x86_64.sh -b -p /conda
ENV PATH="/conda/bin:$PATH"
#RUN mkdir -p /home/jra8/singularity
#WORKDIR /home/jra8/singularity

# environment.yml has pdbfixer, numpy, mdtraj and openmm
ADD environment.yml /environment.yml
RUN conda env create -f /environment.yml
RUN echo "source activate apegen" >> ~/.bashrc
ENV PATH /conda/envs/apegen/bin:$PATH

# These are okay outside environment.yml because they install to PATH
# Perhaps better if inside yml?
RUN conda install -c bioconda smina
RUN conda install -c schrodinger pymol

# To call APE-Gen as module
#RUN echo "export PYTHONPATH=$PYTHONPATH:/home/jra8/singularity" >> ~/.bashrc

# To get nglview to work inside jupyter notebook
#RUN echo "jupyter-nbextension enable --py --sys-prefix widgetsnbextension" >> ~/.bashrc
#RUN echo "jupyter-nbextension enable nglview --py --sys-prefix" >> ~/.bashrc

# RCD
RUN conda install -c intel mkl
#RUN echo "export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/conda/lib/" >> ~/.bashrc
ENV LD_LIBRARY_PATH $LD_LIBRARY_PATH:/conda/lib/
RUN wget -qO- http://chaconlab.org`wget -qO- http://chaconlab.org/modeling/rcd/rcd-download|grep txz|grep uk-button| cut -f4 -d\"` | tar Jxf - && \
    cp RCD_v1.40_Linux_20190228/bin/rcd /usr/local/bin
#RUN echo "export PATH=/RCD_v1.40_Linux_20190228/bin:$PATH" >> ~/.bashrc

# Autodock Vina
RUN wget -qO- http://vina.scripps.edu/download/autodock_vina_1_1_2_linux_x86.tgz | tar zxf - && \
    cp /autodock_vina_1_1_2_linux_x86/bin/vina_split /usr/local/bin
#RUN echo "export PATH=/autodock_vina_1_1_2_linux_x86/bin:$PATH" >> ~/.bashrc


RUN git clone https://github.com/KavrakiLab/APE-Gen.git

# APE-Gen
#WORKDIR /home/apegen
#COPY . /home/jra8/singularity/
ENTRYPOINT ["bash"]
