## 21st Dec 2016
## BISCUIT main and helper functions 
## main()
## Code author SP
##

###
###

############## helper functions ##############

centralize.mat <- function(M){
    n <- nrow(M)
    Q <- matrix(-1/n, nrow=n, ncol = n)
    diag(Q) <-  diag(Q)+1
    M <- Q %*% M %*% Q
    M
}

######
center_colmeans <- function(x) {
    xcenter = colMeans(x)
    x - rep(xcenter, rep.int(nrow(x), ncol(x)))
}

######
norm_vec <- function(x) sqrt(sum(x^2))

######

getmode <- function(v) {
    uniqv <- unique(v)
    uniqv[which.max(tabulate(match(v, uniqv)))]
}
######

######### Projecting X

project.data <- function(data, dim_data){
    
    S <- data %*% t(data)
    Sc <- centralize.mat(S)
    Sc <- 0.5*(Sc + t(Sc))
    eig <- eigen(Sc)
    w <- which(eig$values>0.01)
    if(length(w) < dim_data){
        w <- c(1:dim_data)
    }
    
    sq_diag <- diag(sqrt((eig$values)[w]))
    sq_diag[is.na(sq_diag)] <- 0.001;
    data_pca <- (eig$vectors)[,w] %*% sq_diag

}

############


######### Projecting X

fiedler.vector <- function(data){
    
    s.eigen <- eigen(data) # the eigenvalues are in decreasing order so just extract the 2nd last one.
    return(s.eigen$vectors[,(ncol(data) -1)])
    
}

############


####Fix the number of parallel subprocesses

sub_batch <- function(num_gene_batches){
    flag1 <- TRUE
    if(num_gene_batches==1){
        num_gene_sub_batches <- 1
    }else{
        if(num_gene_batches %% 10 == 0 | (is.prim(num_gene_batches) & (num_gene_batches > 10)) ){
            num_gene_sub_batches <- 10
        }else{
            while(flag1==TRUE){
                for(count1 in 9:1){
                    if(num_gene_batches %% count1 ==0){
                        num_gene_sub_batches <- count1;
                        flag1 <- FALSE
                        break
                    }
                }
            }
        }
    }
    return(sub_batch <- num_gene_sub_batches)
}

##########MDS scaling

mds.tau <- function(H)
{
    n <- nrow(H)
    P <- diag(n) - 1/n
    return(-0.5 * P %*% H %*% P)
}


###### getting distinguishable colours for clusters #####
##ref: http://stackoverflow.com/questions/15282580/how-to-generate-a-number-of-most-distinctive-colors-in-r
##

qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))) #len = 74
num_col <- 40
#pie(rep(1,num_col), col=(col_vector[1:num_col]))
col_palette <- col_vector[1:num_col]; # or sample if you wish

###output directory creation

if( dir.exists(paste0(getwd(),"/", output_folder_name))){
    file.rename(paste0(getwd(),"/", output_folder_name),paste0(getwd(),"/","BISCUIT_previous_run","/"))
}

if(! dir.exists(paste0(getwd(),"/",output_folder_name))){
     dir.create(paste0(getwd(),"/",output_folder_name,"/"))
     dir.create(paste0(getwd(),"/",output_folder_name,"/plots/"))
     dir.create(paste0(getwd(),"/",output_folder_name,"/plots/Inferred_labels/"))
     dir.create(paste0(getwd(),"/",output_folder_name,"/plots/Inferred_labels_per_step_per_batch/"))
     dir.create(paste0(getwd(),"/",output_folder_name,"/plots/Inferred_alphas_betas/"))
     dir.create(paste0(getwd(),"/",output_folder_name,"/plots/Inferred_Sigmas/"))
     dir.create(paste0(getwd(),"/",output_folder_name,"/plots/Inferred_means/"))
     dir.create(paste0(getwd(),"/",output_folder_name,"/plots/extras/"))
}

############## Run BISCUIT ##############

start_time_overall <- Sys.time()

#1) Prepare the input data. Explain what is input and what has to be the output.
source("BISCUIT_process_data.R");


#2) Main MCMC engine. Do not change anything. This runs in parallel where each parallel run takes in a matrix X of all cells and a gene batch i.e. dim(X) is numcells x gene_batch.
source("BISCUIT_IMM_Gibbs_MCMC_parallel.R")

#3) Postprocess MCMC chains from multiple parallel runs
source("BISCUIT_post_MCMC_genesplit_merge.R")

#4) Compute imputed data based on inferred variables and generate plots
source("BISCUIT_post_process.R")
########################################

#print(Sys.time() - start_time_overall)
curr_time <- Sys.time()
print(curr_time - start_time_overall)
write(paste('Overall run time: ',curr_time - start_time_overall),file=f1, append=TRUE)



