FROM debian:latest

# 1. Install essential dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    git \
    wget \
    curl \
    libgl1-mesa-dev \
    libxt-dev \
    libxrender-dev \
    libx11-dev \
    libxext-dev \
    libfreetype6-dev \
    libexpat1-dev \
    libjpeg-dev \
    libpng-dev \
    zlib1g-dev \
    libtiff-dev \
    libhdf5-dev \
    ninja-build \
    make \
 && apt-get clean

# 2. Download and extract VTK v9.4.2 from GitHub
WORKDIR /opt
RUN wget https://github.com/Kitware/VTK/archive/refs/tags/v9.4.2.tar.gz && \
    tar -xzf v9.4.2.tar.gz && \
    rm v9.4.2.tar.gz

# 3. Configure and build VTK with HDF support
WORKDIR /opt/VTK-9.4.2
RUN mkdir build && cd build && cmake -GNinja .. \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX=/usr/local \
    -DVTK_BUILD_ALL_MODULES=OFF \
    -DVTK_GROUP_ENABLE_Qt=NO \
    -DVTK_MODULE_ENABLE_VTK_IOHDF=YES \
    -DVTK_MODULE_ENABLE_VTK_IOImage=YES \
    -DVTK_MODULE_ENABLE_VTK_CommonCore=YES \
    -DVTK_MODULE_ENABLE_VTK_CommonDataModel=YES \
    -DVTK_MODULE_ENABLE_VTK_FiltersCore=YES \
    -DVTK_MODULE_ENABLE_VTK_FiltersGeneral=YES \
    -DVTK_MODULE_ENABLE_VTK_FiltersSources=YES \
    -DVTK_MODULE_ENABLE_VTK_IOGeometry=YES \
    -DVTK_MODULE_ENABLE_VTK_FiltersPoints=YES \
    -DVTK_MODULE_ENABLE_VTK_ImagingCore=YES \
    && ninja && ninja install

# 4. Set environment variables
ENV LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH
ENV PATH=/usr/local/bin:$PATH
ENV CMAKE_PREFIX_PATH=/usr/local

# 5. Working directory for your C++ project
WORKDIR /workspace

# 6. Default shell
CMD ["/bin/bash"]
