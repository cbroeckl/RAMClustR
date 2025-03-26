FROM rocker/r-ver:devel

RUN apt-get update
# Install system dependencies for R packages (you can add more if needed)
RUN apt-get install -y \
            libcurl4-openssl-dev \
            libssl-dev \
            libxml2-dev \
            libfontconfig1-dev \
            libharfbuzz-dev \
            libcairo2-dev \
            libsqlite3-dev \
            libgit2-dev \
            libv8-dev \
            libprotobuf-dev \
            protobuf-compiler \
            libfribidi-dev \
            libfreetype6-dev \
            libpango1.0-dev \
            pkg-config \
            libtiff-dev \
            libpng-dev \
            libjpeg-dev \
            libnetcdf-dev \
            libhdf5-dev

# Install RAMClustR from source
COPY . /RAMClustR
WORKDIR /RAMClustR

# Install R packages
RUN R -e "install.packages('ps')"
RUN R -e "install.packages('ragg')"
RUN R -e "install.packages(c('miniUI', 'pkgdown'))"
RUN R -e "install.packages('devtools')"
RUN R -e 'devtools::install_deps(dependencies = TRUE)'
RUN R -e "BiocManager::install(c('ncdf4', 'mzR', 'MSnbase', 'xcms'))"


# Default command to run R
CMD [ "R"]