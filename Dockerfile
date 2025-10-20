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
    bash \
    libxml2-dev \
    json-c-dev \
    libstdc++ \
    libstdc++-dev

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

# # Download and extract FreeSASA and Gemmi
# WORKDIR /opt
# RUN wget https://github.com/mittinatten/freesasa/archive/refs/tags/2.1.2.tar.gz -O freesasa-2.1.2.tar.gz && \
#     wget https://github.com/project-gemmi/gemmi/archive/refs/tags/v0.5.7.tar.gz -O gemmi-0.5.7.tar.gz && \
#     tar -xzf freesasa-2.1.2.tar.gz && \
#     tar -xzf gemmi-0.5.7.tar.gz && \
#     mkdir -p freesasa-2.1.2/third-party/gemmi && \
#     cp -r gemmi-0.5.7/* freesasa-2.1.2/third-party/gemmi/

# Build and install FreeSASA with mmCIF support
# WORKDIR /opt/freesasa-2.1.2
# RUN autoreconf -i && \
#     ./configure && \
#     sed -i 's/-lc++/-lstdc++/g' src/Makefile && \
#     make && \
#     make install

# Clean up sources
# WORKDIR /opt
# RUN rm -rf freesasa-2.1.2* gemmi-0.5.7* 

# Install MSMS (molecular surface calculation tool)
# WORKDIR /opt
# RUN wget https://ccsb.scripps.edu/msms/download/933/ -O msms.tar.gz && \
#     tar -xzf msms.tar.gz && \
#     find . -name "msms.x86_64Linux2.2.6.1" -exec cp {} /usr/local/bin/msms \; && \
#     chmod +x /usr/local/bin/msms && \
#     rm -rf msms* msms.tar.gz



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
