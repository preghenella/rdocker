FROM rootproject/root-ubuntu16

USER root

RUN apt-get update && \
    apt-get -y install \
    emacs-nox \
    libboost-all-dev \
    swig

ENV FAST_ROOT /app
RUN mkdir -p /app

COPY bin     /app/bin
COPY patches /app/patches

RUN /bin/bash -c "source /app/bin/fastenv && /app/bin/fastbuild"

RUN /bin/bash -c "source /app/bin/fastenv && lhapdf install CT10"

COPY utils    /app/utils
COPY examples /app/examples

ARG username=fastuser
RUN userdel -r builder && useradd --create-home --home-dir /home/${username} ${username}
ENV HOME /home/${username}

RUN echo "root:fastroot" | chpasswd

RUN chmod +x /app/bin/fastenv

USER ${username}
WORKDIR /home/${username}

CMD ["/app/bin/fastenv"]
