FROM ubuntu:noble

WORKDIR /workspace

RUN apt update && apt -y dist-upgrade && DEBIAN_FRONTEND=noninteractive apt install -y locales \
octave octave-doc octave-dev octave-control octave-image octave-io octave-optim octave-signal octave-statistics octave-communications

ENV LC_ALL en_US.UTF-8
ENV LANG en_US.UTF-8
ENV LANGUAGE en_US.UTF-8