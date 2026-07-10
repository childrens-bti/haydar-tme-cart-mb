FROM rocker/tidyverse:4.4.0
LABEL maintainer="Sam Chen (schen8@childrensnational.org), Bicna Song (bsong@childrensnational.org)"
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
RUN R -e 'options(Ncpus = max(1, parallel::detectCores()-1), repos = c(CRAN="https://cloud.r-project.org")); \
          remotes_urls <- c("https://cloud.r-project.org/src/contrib/remotes_2.5.0.tar.gz", "https://cloud.r-project.org/src/contrib/Archive/remotes/remotes_2.5.0.tar.gz"); \
          remotes_installed <- FALSE; \
          for (url in remotes_urls) { \
            remotes_installed <- tryCatch({ install.packages(url, repos = NULL, type = "source"); TRUE }, error = function(e) FALSE); \
            if (remotes_installed) break; \
          }; \
          stopifnot(remotes_installed, as.character(packageVersion("remotes")) == "2.5.0"); \
          remotes::install_version("BiocManager", version = "1.30.23", upgrade = "never")'

RUN R -e 'options(Ncpus = max(1, parallel::detectCores()-1), repos = c(CRAN="https://cloud.r-project.org")); \
          cran_packages <- c( \
            Seurat = "5.4.0",         \
            sctransform = "0.4.3",    \
            SeuratObject = "5.3.0",   \
            rprojroot = "2.1.1",      \
            tidyverse = "2.0.0",      \
            Matrix = "1.7-4",         \
            devtools = "2.4.5",       \
            uwot = "0.2.4",           \
            RcppAnnoy = "0.0.23",     \
            scCustomize = "3.2.4",    \
            SoupX = "1.6.2",          \
            instantiate = "0.2.3",    \
            msigdbr = "25.1.1",       \
            ggthemes = "5.2.0",       \
            future = "1.69.0",        \
            future.apply = "1.20.2",  \
            igraph = "2.2.2"          \
          ); \
          mapply(function(pkg, version) remotes::install_version(pkg, version = version, upgrade = "never"), names(cran_packages), cran_packages); \
          failed <- names(cran_packages)[!vapply(names(cran_packages), function(pkg) requireNamespace(pkg, quietly = TRUE) && identical(packageDescription(pkg)$Version, cran_packages[[pkg]]), logical(1))]; \
          if (length(failed)) stop("CRAN package version check failed: ", paste(failed, collapse = ", "))'

# ---- Bioconductor packages ----
RUN R -e 'options(Ncpus = max(1, parallel::detectCores()-1)); \
          BiocManager::install(version = "3.20", ask = FALSE); \
          bioc_packages <- c( \
            Biobase = "2.66.0",             \
            biomaRt = "2.62.1",             \
            SingleR = "2.8.0",              \
            ComplexHeatmap = "2.22.0",      \
            dittoSeq = "1.18.0",            \
            DropletUtils = "1.26.0",        \
            BiocParallel = "1.40.2",        \
            Nebulosa = "1.16.0",            \
            celldex = "1.16.0",             \
            fgsea = "1.32.4",               \
            AUCell = "1.28.0",              \
            MAST = "1.32.0",                \
            UCell = "2.10.1",               \
            DESeq2 = "1.46.0",              \
            EnhancedVolcano = "1.24.0",     \
            slingshot = "2.14.0"            \
          ); \
          BiocManager::install(names(bioc_packages), ask = FALSE, update = FALSE); \
          failed <- names(bioc_packages)[!vapply(names(bioc_packages), function(pkg) requireNamespace(pkg, quietly = TRUE) && identical(packageDescription(pkg)$Version, bioc_packages[[pkg]]), logical(1))]; \
          if (!identical(as.character(BiocManager::version()), "3.20")) stop("Bioconductor version check failed"); \
          if (length(failed)) stop("Bioconductor package version check failed: ", paste(failed, collapse = ", "))'

RUN R -e 'options(Ncpus = max(1, parallel::detectCores()-1), repos = c(CRAN="https://cloud.r-project.org")); \
          remotes::install_version("scGate", version = "1.7.2", upgrade = "never"); \
          stopifnot(requireNamespace("scGate", quietly = TRUE), identical(packageDescription("scGate")$Version, "1.7.2"))'

RUN R -e 'options(Ncpus = max(1, parallel::detectCores()-1), repos = c(CRAN="https://cloud.r-project.org")); \
          remotes::install_version("NMF", version = "0.27", upgrade = "never"); \
          stopifnot(requireNamespace("NMF", quietly = TRUE), identical(packageDescription("NMF")$Version, "0.27"))'

## install GitHub packages
# colorblindr 0.1.0
RUN R -e "remotes::install_github('clauswilke/colorblindr', ref = '1ac3d4d62dad047b68bb66c06cee927a4517d678', dependencies = TRUE, upgrade = 'never')"
# patchwork 1.3.2.9000
RUN R -e "remotes::install_github('thomasp85/patchwork', ref = '6b1d88ce1da1c5cae3818d984edf80dc0bb3de8c', upgrade = 'never')"
# DoubletFinder 2.0.6
RUN R -e 'remotes::install_github("chris-mcginnis-ucsf/DoubletFinder", ref = "1b244d8f0d54b4b1cb4365639931bbb16f01e1cd", upgrade = "never")'
# presto 1.0.0
RUN mkdir -p /root/.R && printf 'CXX11 = %s\nCXX11STD = -std=gnu++14\n' "$(R CMD config CXX14)" > /root/.R/Makevars
RUN R -e 'remotes::install_github("immunogenomics/presto", ref = "7636b3d0465c468c35853f82f1717d3a64b3c8f6", upgrade = "never")'
# sccomp 2.1.30
RUN R -e 'remotes::install_github("MangiolaLaboratory/sccomp", ref = "18a4be86f6497e96389339620ff0c3e7faa71492", upgrade = "never")'
# miloR 2.9.1
RUN R -e 'remotes::install_github("MarioniLab/miloR", ref = "ff744bbb5d793163b59f28483c1ad05192fddc15", upgrade = "never")'

## finish sccomp installation
# cmdstanr 0.9.0
RUN R -e 'remotes::install_github("stan-dev/cmdstanr", ref = "da99e2ba954658bdad63bffb738c4444c33a4e0e", upgrade = "never")'
RUN R -e 'cmdstanr::check_cmdstan_toolchain(fix = TRUE)'
RUN R -e 'cmdstanr::install_cmdstan()'

# Required for hdf5r (dependency of SeuratExtend)
RUN apt-get update && apt-get -y --no-install-recommends install \
    libhdf5-dev patch

RUN R -e 'remotes::install_version("hdf5r", version = "1.3.10", type = "source", upgrade = "never")'
# SeuratExtend 1.2.10
RUN R -e 'remotes::install_github("huayc09/SeuratExtend", ref = "f567d9c22a43ac538aedca0e4630421a23bd568f", dependencies = TRUE, upgrade = "never")'

RUN R -e 'expected_versions <- c( \
            remotes = "2.5.0", BiocManager = "1.30.23", \
            Seurat = "5.4.0", sctransform = "0.4.3", SeuratObject = "5.3.0", \
            rprojroot = "2.1.1", tidyverse = "2.0.0", Matrix = "1.7-4", \
            devtools = "2.4.5", \
            uwot = "0.2.4", RcppAnnoy = "0.0.23", scCustomize = "3.2.4", \
            SoupX = "1.6.2", instantiate = "0.2.3", msigdbr = "25.1.1", \
            ggthemes = "5.2.0", scGate = "1.7.2", NMF = "0.27", future = "1.69.0", \
            future.apply = "1.20.2", igraph = "2.2.2", \
            Biobase = "2.66.0", biomaRt = "2.62.1", SingleR = "2.8.0", \
            ComplexHeatmap = "2.22.0", dittoSeq = "1.18.0", DropletUtils = "1.26.0", \
            BiocParallel = "1.40.2", Nebulosa = "1.16.0", celldex = "1.16.0", fgsea = "1.32.4", \
            AUCell = "1.28.0", MAST = "1.32.0", UCell = "2.10.1", \
            DESeq2 = "1.46.0", EnhancedVolcano = "1.24.0", slingshot = "2.14.0", \
            colorblindr = "0.1.0", patchwork = "1.3.2.9000", DoubletFinder = "2.0.6", \
            presto = "1.0.0", sccomp = "2.1.30", miloR = "2.9.1", \
            cmdstanr = "0.9.0", hdf5r = "1.3.10", SeuratExtend = "1.2.10" \
          ); \
          expected_sha <- c( \
            colorblindr = "1ac3d4d62dad047b68bb66c06cee927a4517d678", \
            patchwork = "6b1d88ce1da1c5cae3818d984edf80dc0bb3de8c", \
            DoubletFinder = "1b244d8f0d54b4b1cb4365639931bbb16f01e1cd", \
            presto = "7636b3d0465c468c35853f82f1717d3a64b3c8f6", \
            sccomp = "18a4be86f6497e96389339620ff0c3e7faa71492", \
            miloR = "ff744bbb5d793163b59f28483c1ad05192fddc15", \
            cmdstanr = "da99e2ba954658bdad63bffb738c4444c33a4e0e", \
            SeuratExtend = "f567d9c22a43ac538aedca0e4630421a23bd568f" \
          ); \
          version_failed <- names(expected_versions)[!vapply(names(expected_versions), function(pkg) requireNamespace(pkg, quietly = TRUE) && identical(packageDescription(pkg)$Version, expected_versions[[pkg]]), logical(1))]; \
          sha_failed <- names(expected_sha)[!vapply(names(expected_sha), function(pkg) identical(packageDescription(pkg)$RemoteSha, expected_sha[[pkg]]), logical(1))]; \
          if (!identical(as.character(BiocManager::version()), "3.20")) stop("Bioconductor version check failed"); \
          if (length(version_failed)) stop("Package version check failed: ", paste(version_failed, collapse = ", ")); \
          if (length(sha_failed)) stop("GitHub SHA check failed: ", paste(sha_failed, collapse = ", "))'

WORKDIR /rocker-build/

ADD Dockerfile .
