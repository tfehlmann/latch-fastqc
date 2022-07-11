FROM 812206152185.dkr.ecr.us-west-2.amazonaws.com/latch-base:9a7d-main

SHELL ["/bin/bash", "--login", "-c"]
ENV BASH_ENV ~/.bashrc


RUN apt-get update && apt-get install -y wget
RUN wget -qO- https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj bin/micromamba

COPY envs /envs

RUN bin/micromamba shell init -s bash -p ~/micromamba &&\
    source ~/.bashrc && \
    micromamba activate &&\
    micromamba install -f /envs/base.yml

# STOP HERE:
# The following lines are needed to ensure your build environement works
# correctly with latch.
ARG tag
ENV FLYTE_INTERNAL_IMAGE $tag
RUN python3 -m pip install --upgrade latch
COPY wf /root/wf
WORKDIR /root
