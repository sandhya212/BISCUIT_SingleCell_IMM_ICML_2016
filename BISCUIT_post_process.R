## 24th May 2017
## BISCUIT postprocessing
##
## Code author SP
##

###
###

if(num_gene_batches ==1){
    
    source("BISCUIT_parallel_impute_onegenebatch.R")
    source("BISCUIT_extras_onegenebatch.R")
 
}else{
    
    source("BISCUIT_parallel_impute.R")
    source("BISCUIT_extras.R")

}

