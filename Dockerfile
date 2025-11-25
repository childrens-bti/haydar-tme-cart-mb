FROM rocker/tidyverse:4.4.0
LABEL maintainer="Sam Chen (schen8@childrensnational.org)"
WORKDIR /rocker-build/

RUN apt-get update && apt-get install -y --no-install-recommends apt-utils dialog

# Add curl, bzip2 and some dev libs
RUN apt-get update -qq && apt-get -y --no-install-recommends install \
    curl \
    bzip2 \
    zlib1g \
    libbz2-dev \
    liblzma-dev \
    libreadline-dev \
    libglpk40 \
    libglpk-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libmagick++-dev

# libmagick++-dev is needed for coloblindr to install
RUN apt-get -y --no-install-recommends install \
    libgdal-dev \
    libudunits2-dev \
    libmagick++-dev \
    libgsl27 \
    libgsl-dev

# Required for installing pdftools, which is a dependency of gridGraphics
RUN apt-get -y --no-install-recommends install \
    libpoppler-cpp-dev

# Install java
RUN apt-get update && apt-get -y --no-install-recommends install \
   default-jdk \
   libxt6

# Enable github copilot
RUN mkdir -p /etc/rstudio && echo "copilot-enabled=1" >> /etc/rstudio/rsession.conf

# ---- R packages (scRNA-seq only) ----
# Use CRAN mirror explicitly; Seurat v5 depends on SeuratObject (v5), sctransform, etc.
RUN R -e 'options(Ncpus = max(1, parallel::detectCores()-1), repos = c(CRAN="https://cloud.r-project.org")); \
          install.packages(c( \
            "Seurat",         \
            "sctransform",    \
            "SeuratObject",   \
            "rprojroot",      \
            "tidyverse",      \
            "Matrix",         \
            "devtools",       \
            "remotes",        \
            "uwot",           \
            "RcppAnnoy",       \
            "scCustomize",     \
            "SoupX", \
            "instantiate", \
            "msigdbr" \
          ))'

# ---- Bioconductor packages ----
RUN R -e 'if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org"); \
          BiocManager::install(c("biomaRt", "SingleR", "ComplexHeatmap", "dittoSeq", "DropletUtils", "Nebulosa", "celldex", \
          "fgsea", "AUCell", "MAST"), ask = FALSE, update = TRUE)'

## install GitHub packages
RUN R -e "remotes::install_github('clauswilke/colorblindr', ref = '1ac3d4d62dad047b68bb66c06cee927a4517d678', dependencies = TRUE)"
RUN R -e "remotes::install_github('thomasp85/patchwork')"
RUN R -e 'remotes::install_github("chris-mcginnis-ucsf/DoubletFinder")'
RUN R -e 'remotes::install_github("immunogenomics/presto")'
RUN R -e 'remotes::install_github("MangiolaLaboratory/sccomp")'

## finish sccomp installation
RUN R -e 'install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev/", getOption("repos")))'
RUN R -e 'cmdstanr::check_cmdstan_toolchain(fix = TRUE)'
RUN R -e 'cmdstanr::install_cmdstan()'

WORKDIR /rocker-build/

ADD Dockerfile .