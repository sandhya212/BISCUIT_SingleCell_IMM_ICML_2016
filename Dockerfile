##### 
#  BISCUIT (Bayesian Inference for Single-cell ClUstering and ImpuTing) 
#  http://proceedings.mlr.press/v48/prabhakaran16.pdf 
#  https://github.com/sandhya212/BISCUIT_SingleCell_IMM_ICML_2016 
##### 

FROM ubuntu:latest 
RUN rm /bin/sh && ln -s /bin/bash /bin/sh
MAINTAINER Steve Tsang <mylagimail2004@yahoo.com> 
RUN apt-get clean && apt-get update

RUN apt-get install --yes \
 build-essential \
 git \
 wget \
 nano \
 r-base

# Get latest source from releases 
ENV SRC /opt 
ENV BIN /usr/local/bin 

RUN echo 'install.packages(c("MCMCpack","mvtnorm","ellipse","coda","Matrix","Rtsne","gtools","foreach","doParallel","doSNOW","snow","lattice","MASS","bayesm","robustbase","chron","mnormt","schoolmath","devtools","RColorBrewer"), repos="http://cran.us.r-project.org", dependencies=TRUE)' > /tmp/packages.R \
    && Rscript /tmp/packages.R

WORKDIR $SRC
RUN git clone https://github.com/sandhya212/BISCUIT_SingleCell_IMM_ICML_2016
COPY Dockerfile /opt/Dockerfile

## Add Example 
WORKDIR $SRC/BISCUIT_SingleCell_IMM_ICML_2016/
RUN wget https://storage.googleapis.com/linnarsson-lab-www-blobs/blobs/cortex/expression_mRNA_17-Aug-2014.txt
# set the command 
#CMD ["R"]
