## 21st Dec 2016
## BISCUIT R implementation
## Start_file with user inputs
## 
## Code author SP


###
###

############## packages required ##############

library(MCMCpack)
library(mvtnorm)
library(ellipse)
library(coda)
library(Matrix)
library(Rtsne)
library(gtools)
library(foreach)
library(doParallel)
library(doSNOW)
library(snow)
library(lattice)
library(MASS)
library(bayesm)
library(robustbase)
library(chron)
library(mnormt)
library(schoolmath)
library(RColorBrewer)

#############################################


input_file_name <- "expression_mRNA_17-Aug-2014.txt";

input_data_tab_delimited <- TRUE; #set to TRUE if the input data is tab-delimited

is_format_genes_cells <-  TRUE; #set to TRUE if input data has rows as genes and columns as cells

choose_cells <- 3000; #comment if you want all the cells to be considered

choose_genes <- 150; #comment if you want all the genes to be considered

gene_batch <- 50; #number of genes per batch, therefore num_batches = choose_genes (or numgenes)/gene_batch. Max value is 150

num_iter <- 20; #number of iterations, choose based on data size.

num_cores <- detectCores() - 4; #number of cores for parallel processing. Ensure that detectCores() > 1 for parallel processing to work, else set num_cores to 1.

z_true_labels_avl <- TRUE; #set this to TRUE if the true labels of cells are available, else set it to FALSE. If TRUE, ensure to populate 'z_true' with the true labels in 'BISCUIT_process_data.R'

num_cells_batch <- 1000; #set this to 1000 if input number of cells is in the 1000s, else set it to 100.

alpha <- 1; #DPMM dispersion parameter. A higher value spins more clusters whereas a lower value spins lesser clusters.

output_folder_name <- "output"; #give a name for your output folder.

## call BISCUIT
source("BISCUIT_main.R")

