# Base image with COBRAme and StressME tools
FROM queensysbio/stressme:v1.1

USER root
WORKDIR /opt

# Install build tools and Python 3 deps
RUN apk --no-cache --force add \
    build-base \
    python3 \
    python3-dev \
    py3-pip \
    musl-dev \
    cython \
    wget \
    tar \
    libtool \
    autoconf \
    automake \
    openssh \
    bash

# Set python and pip defaults to version 3
RUN ln -sf python3 /usr/bin/python && ln -sf pip3 /usr/bin/pip

# Install OpenMPI from source
ENV OPENMPI_VERSION=4.1.5

RUN wget https://download.open-mpi.org/release/open-mpi/v4.1/openmpi-${OPENMPI_VERSION}.tar.gz && \
    tar -xzf openmpi-${OPENMPI_VERSION}.tar.gz && \
    cd openmpi-${OPENMPI_VERSION} && \
    ./configure --prefix=/usr/local && \
    make -j$(nproc) && \
    make install && \
    cd .. && rm -rf openmpi-${OPENMPI_VERSION}*

ENV PATH=/usr/local/bin:$PATH
ENV LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH

# Install DynamicME
COPY dynamicme /opt/dynamicme
WORKDIR /opt/dynamicme
RUN python setup.py install

# Set working directory for your project
WORKDIR /app

# Copy your project files
COPY . /app

# Install Python dependencies (Python 3)
RUN pip install --upgrade pip
RUN pip install -r requirements.txt

# Expose Jupyter Notebook port
EXPOSE 8888

CMD ["jupyter", "notebook", "--ip=0.0.0.0", "--port=8888", "--no-browser", "--allow-root", "--NotebookApp.token=''", "--NotebookApp.password=''"]
