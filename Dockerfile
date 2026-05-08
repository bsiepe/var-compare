# R version tag sets the CRAN snapshot date automatically
# (packages are pinned to the last day of R 4.3.2's lifecycle)
FROM rocker/verse:4.3.2

WORKDIR /home/rstudio

# CRAN packages
# rstan/BGGM require C++ compilation and take a long time to build
RUN install2.r --error --skipinstalled --ncpus -1 \
    doParallel doRNG foreach \
    BGGM graphicalVAR mlVAR mgm \
    BFpack bain bayestestR \
    qgraph igraph \
    matrixcalc mvtnorm corpcor glasso glmnet \
    cowplot patchwork ggridges ggh4x ggokabeito GGally ggdist \
    sysfonts showtext \
    rrapply reshape2 Hmisc \
    lme4 lavaan emmeans effectsize parameters insight psych \
    imputeTS gt gtExtras magick \
    here rmarkdown

# tsnet was not on CRAN at the time of the study
RUN installGithub.r bsiepe/tsnet

COPY . /home/rstudio/
