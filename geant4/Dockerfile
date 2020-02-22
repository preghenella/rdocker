FROM cern/cc7-base

USER root

RUN yum -y update && \
    yum -y install \
    git cmake gcc-c++ gcc binutils libX11-devel libXpm-devel libXft-devel libXext-devel \
    giflib libtiff \
    emacs-nox cmake3 make boost-devel xerces-c-devel python-devel \
    centos-release-scl devtoolset-3-gcc devtoolset-3-gcc-c++ \
    mesa-libGLU mesa-libGLU-devel libXmu-devel qt5-qtbase-devel libxkbcommon-x11

WORKDIR /root

RUN mkdir /app

### ROOT

ENV ROOT_DIR /app/ROOT-6.18.02

RUN curl -kOL https://root.cern/download/root_v6.18.02.Linux-centos7-x86_64-gcc4.8.tar.gz && \
    tar zxvf root_v6.18.02.Linux-centos7-x86_64-gcc4.8.tar.gz && \
    rm root_v6.18.02.Linux-centos7-x86_64-gcc4.8.tar.gz && \
    mv root $ROOT_DIR

### GEANT4

ENV GEANT4_DIR /app/Geant4-10.5.1

RUN curl -kOL https://github.com/Geant4/geant4/archive/v10.5.1.tar.gz && \
    tar zxvf v10.5.1.tar.gz && \
    rm -rf v10.5.1.tar.gz && \
    mkdir -p geant4-10.5.1/build && \
    cd geant4-10.5.1/build && \
    cmake3 ../. -DGEANT4_INSTALL_DATA=ON -DGEANT4_USE_OPENGL_X11=ON -DOpenGL_GL_PREFERENCE=LEGACY -DGEANT4_USE_QT=ON -DCMAKE_INSTALL_PREFIX=$GEANT4_DIR && \
    make -j2 install && \
    cd - && \
    rm -rf geant4-10.5.1

### HEPMC3

ENV HEPMC3_DIR /app/HepMC3-3.1.2

RUN curl -kOL http://hepmc.web.cern.ch/hepmc/releases/HepMC3-3.1.2.tar.gz && \
    tar zxvf HepMC3-3.1.2.tar.gz && \
    rm HepMC3-3.1.2.tar.gz && \
    cd HepMC3-3.1.2 && \
    cmake3 . -DROOT_DIR=$ROOT_DIR -DHEPMC3_INSTALL_INTERFACES=ON -DCMAKE_INSTALL_PREFIX=$HEPMC3_DIR && \
    make -j2 install && \
    cd - && \
    rm -rf HepMC3-3.1.2


ENV PATH $PATH:$ROOT_DIR/bin:$GEANT4_DIR/bin:$HEPMC3_DIR/bin
ENV LD_LIBRARY_PATH $LD_LIBRARY_PATH:$ROOT_DIR/lib:$GEANT4_DIR/lib64:$HEPMC3_DIR/lib64

### THE ENVIRONMENT

RUN echo -e "\nsource geant4.sh" >> .bash_profile

### READY

CMD ["bash"]
