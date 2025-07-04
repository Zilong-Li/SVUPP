FROM ubuntu:24.04

# Prevent interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive

# Set environment variables for pyenv
ENV PYENV_ROOT="/opt/pyenv"
ENV PATH="$PYENV_ROOT/bin:$PYENV_ROOT/shims:$PATH"

RUN apt-get update
RUN apt-get install --assume-yes software-properties-common
RUN add-apt-repository universe
RUN apt-get update
RUN apt-get install --assume-yes build-essential libbz2-dev libcurl4-openssl-dev libssl-dev zlib1g-dev liblzma-dev libdeflate-dev libncurses5-dev curl autoconf git
RUN apt-get install --assume-yes make libreadline-dev libsqlite3-dev libncursesw5-dev xz-utils tk-dev libxml2-dev libxmlsec1-dev libffi-dev liblzma-dev

RUN git clone https://github.com/pyenv/pyenv.git /opt/pyenv
RUN export PYENV_ROOT="/opt/pyenv"
RUN export PATH="$PYENV_ROOT/bin:$PYENV_ROOT/shims:$PATH"
RUN eval "$(pyenv init -)" && \
    pyenv install 3.10.12 && \
    pyenv global 3.10.12 && \
    python -m pip install --no-cache-dir numpy==2.1.3 scikit-learn==1.6.1

RUN mkdir /metadata
RUN mkdir /repos
WORKDIR /repos

RUN git clone --branch v0.0.2 https://github.com/Zilong-Li/cuteSV
RUN cd cuteSV && \
    python setup.py install && \
    ln -s /opt/pyenv/versions/3.10.12/bin/cuteSV /usr/local/bin/ && \
    echo `git rev-parse --short HEAD` > /metadata/cutesv.git.version && \
    cd ..

RUN apt-get install --assume-yes tabix samtools bcftools
RUN apt-get install -y r-base r-cran-devtools

RUN git clone https://github.com/rwdavies/QUILT
RUN cd QUILT && \
    Rscript ./scripts/install-r-dependencies.R && \
    Rscript ./scripts/build-and-install.R && \
    cp *.R /usr/local/bin && \
    echo `git rev-parse --short HEAD` > /metadata/quilt.git.version && \
    cd ..

RUN apt-get remove --assume-yes git software-properties-common build-essential
RUN apt-get autoremove --assume-yes
RUN apt-get clean --assume-yes
