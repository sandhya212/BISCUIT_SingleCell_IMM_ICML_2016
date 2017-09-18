## 4th April 2017
## Findng the best group of cells across gene-splits.
## by building a global confusion matrix.
##
## 17th April 2017
## parallelised the global confusion matrix creation.
##
## 18-19th April 2017
##
## Construct the MAP esitmates for momets and then resample
##
## 1st May 2017
## propagating weights for weighted CM
##
## Code author SP


###
###
###


################

#Construct the global Confusion matrix.
print("Computing the global confusion matrix")
print("Monitor log_CM.txt in outputs folder and debug_CM.txt")

strt_CM <- Sys.time()


divsr <- numcells %/% num_cells_batch
dividend <- numcells %% num_cells_batch


CM_local <- function(df){
    
    df_indicator <- paste0("Batch in process is: ",df);
    write.table(df_indicator,file=paste0(getwd(),"/",output_folder_name,"/log_CM.txt"),append=TRUE,sep="");
    
    cell_j_vals <- matrix(0,length(((1+(num_cells_batch*(df-1))):(num_cells_batch*df))),numcells)
    write.table("Cells selected are ",file=paste0(getwd(),"/",output_folder_name,"/log_CM.txt"),append=TRUE,sep="");
    write.table(((1+(num_cells_batch*(df-1))):(num_cells_batch*df)),file=paste0(getwd(),"/",output_folder_name,"/log_CM.txt"),append=TRUE,sep="");
    
    rownames(cell_j_vals) <- as.character(((1+(num_cells_batch*(df-1))):(num_cells_batch*df)))
    colnames(cell_j_vals) <- as.character(1:numcells)
    for ( cell_i in((1+(num_cells_batch*(df-1))):(num_cells_batch*df))){
        print(paste("comparing cell",cell_i));
        
        for(cell_j in (cell_i):numcells){
            
            for (i in 1:length(results.all.MCMC)){
                for(j in 1:num_gene_sub_batches){
                    if(results.all.MCMC[[i]][[j]][[1]][cell_i] == results.all.MCMC[[i]][[j]][[1]][cell_j]){
                        #global_Conf_Mat1[cell_i,cell_j] = global_Conf_Mat1[cell_i,cell_j] + 1;
                        #print(i)
                        cell_j_vals[as.character(cell_i),as.character(cell_j)] <- cell_j_vals[as.character(cell_i),as.character(cell_j)] + 1*results.all.MCMC[[i]][[j]][[6]];
                    }
                }
            }
        }
    }
    return(cell_j_vals)
}


# Start the cluster and register with doSNOW
cl <- makeCluster(num_cores, type = "SOCK",outfile="debug_CM.txt") #opens multiple socket connections
clusterExport(cl, "CM_local")
registerDoSNOW(cl)

# Call parallel processes to build CM

global_Conf_Mat <- foreach (df = 1:(divsr), .combine="rbind") %dopar%{
    local_CM_batch <- CM_local(df)
    return(local_CM_batch)
}

print("Time for creating the upper triangular global confusion matrix is ")
print(Sys.time()-strt_CM);

stopCluster(cl)

##Taking care of the remaining cells that fell off the parallel bins

if (dividend !=0){
    cell_j_vals <- matrix(0,dividend,numcells)
    
    df_indicator <- paste0("Batch in process is: ",divsr+1);
    write.table(df_indicator,file=paste0(getwd(),"/",output_folder_name,"/log_CM.txt"),append=TRUE,sep="");
    
    write.table("Cells selected are ",file=paste0(getwd(),"/",output_folder_name,"/log_CM.txt"),append=TRUE,sep="");
    write.table(((num_cells_batch*divsr+1):(num_cells_batch*divsr+dividend)),file=paste0(getwd(),"/",output_folder_name,"/log_CM.txt"),append=TRUE,sep="");
    
    rownames(cell_j_vals) <- as.character((num_cells_batch*divsr+1):(num_cells_batch*divsr+dividend))
    colnames(cell_j_vals) <- as.character(1:numcells)
    for ( cell_i in (num_cells_batch*divsr+1):(num_cells_batch*divsr+dividend)){
        #print(paste("comparing cell",cell_i));
        f <- paste("comparing cell ",cell_i)
        write.table(f,file=paste0(getwd(),"/debug_CM.txt"),append=TRUE,sep="");
        
        for(cell_j in (cell_i):numcells){
            
            for (i in 1:length(results.all.MCMC)){
                for(j in 1:num_gene_sub_batches){
                    if(results.all.MCMC[[i]][[j]][[1]][cell_i] == results.all.MCMC[[i]][[j]][[1]][cell_j]){
                        #global_Conf_Mat1[cell_i,cell_j] = global_Conf_Mat1[cell_i,cell_j] + 1;
                        #print(i)
                        cell_j_vals[as.character(cell_i),as.character(cell_j)] <- cell_j_vals[as.character(cell_i),as.character(cell_j)] + 1*results.all.MCMC[[i]][[j]][[6]];
                    }
                }
            }
        }
    }
    global_Conf_Mat  <- rbind(global_Conf_Mat , cell_j_vals)
}


print("Time for creating the overall upper triangular global confusion matrix is ")
print(Sys.time()-strt_CM);

CM <- global_Conf_Mat + t(global_Conf_Mat);
diag(CM) <- 0;

#average CM
ave_CM <- CM / num_gene_batches;


f <- paste0(getwd(),"/",output_folder_name,"/plots/extras/ave_global_Conf_mat.txt");
write.matrix(ave_CM, file=f, sep="\t")



####
####
####Using the CM to do a k-means

num_clusters_per_batch <- rep(0,num_gene_batches);
batch <- 1
for (i in 1:length(results.all.MCMC)){
        for(j in 1:num_gene_sub_batches){
            inferred_labels <- results.all.MCMC[[i]][[j]][[1]];
            uqe_inferred_labels <- unique(inferred_labels);
            num_clusters_per_batch[batch] <- length(uqe_inferred_labels);
            batch <- batch + 1;
        }
}

print(paste("Number of clusters over gene splits: ",num_clusters_per_batch));

mean_num_clusters <- round(max(num_clusters_per_batch));
print(paste("Mean number of clusters over gene splits: ",mean_num_clusters));

#
print("Cluster the global CM by kmeans")
kmeans_CM <- kmeans(ave_CM,mean_num_clusters,nstart=1)
z_inferred_final <- kmeans_CM$cluster
f <- paste0(getwd(),"/",output_folder_name,"/plots/Inferred_labels/Conf_mat_final_inferred_kmeans.pdf");
pdf(file=f);
plot(X_tsne_all$Y[,1],X_tsne_all$Y[,2],col = col_palette[1*(z_inferred_final)],  main="t-SNE of X (inferred labels) based on kmeans of CM");
dev.off();

rm(CM)
rm(ave_CM)

###########################
## Construct the MAP estimates cluster moments
##
## Since the splits are in equal sizes, we can create the MAP estimates for the moments by concatenating moments via law of total expectation and law of total moments from the MCMC chains.
## mu_MAP and Sigma_MAP are formed by calculating the means and covariances separately for each batch and averaging those vectors/matrices (with arithmetic mean) to form one total mean/covariance per cluster.
##
print("Construct the cluster moments based on k-means of the CM")
MAP.non.parallel <- FALSE
if (MAP.non.parallel){ # cycles per cluster
mu_MAP <- matrix(0, gene_batch*num_gene_batches, mean_num_clusters); #or numgenes?
Sigma_MAP <- list();


for(i in 1:mean_num_clusters){
    print(paste0("cluster is ", i));
    cell_ind <- which(z_inferred_final == i);
    #mu_temp <- matrix(0,gene_batch*num_gene_batches, length(cell_ind));
    mu_temp <- rep(0,gene_batch*num_gene_batches);
    
    #Sigma_temp <- list();
    Sigma_temp <- matrix(0,gene_batch*num_gene_batches,gene_batch*num_gene_batches);
    
    for (j in 1:length(cell_ind)){
        mu_temp_cellj <- rep(0,gene_batch*num_gene_batches);
        
        
        p <- paste0("cell is ", cell_ind[j]);
        print(p);
        write.table(p,file=paste0(getwd(),"/",output_folder_name,"/log_CM.txt"),append=TRUE,sep="");
        count <- 1;

        Sigma_temp_per_cell <- matrix(0,gene_batch*num_gene_batches,gene_batch*num_gene_batches)
        for (r1 in 1:length(results.all.MCMC)){
            for(r2 in 1:num_gene_sub_batches){
                inferred_label_per_batch <- results.all.MCMC[[r1]][[r2]][[1]][cell_ind[j]];
                
                pos.mu.temp <- (1+(gene_batch*(count-1))):(gene_batch*count)
                mu_temp_cellj[pos.mu.temp] <- results.all.MCMC[[r1]][[r2]][[4]][,as.character(inferred_label_per_batch)];
                
                
                Sigma_temp_per_cell[pos.mu.temp,pos.mu.temp] <- results.all.MCMC[[r1]][[r2]][[5]][,,as.character(inferred_label_per_batch)];
                count <- count + 1;
            }
        }
        #Sigma_temp[[j]] <- Sigma_temp_per_cell
        Sigma_temp <- Sigma_temp + Sigma_temp_per_cell
        rm(Sigma_temp_per_cell)
        mu_temp <- sum(mu_temp,mu_temp_cellj);
        rm(mu_temp_cellj);
        
    }
    
    mu_MAP[,i] <- mu_temp/length(cell_ind);
    Sigma_MAP[[i]] <- Sigma_temp/length(cell_ind);
    diag(Sigma_MAP[[i]]) <- diag(Sigma_MAP[[i]]) + 0.001;
    #diag(Sigma_MAP[[i]]) <- diag(Sigma_MAP[[i]]) + 0.001
}

}else{ #cycles per cluster in parallel. If there are space /vector allocation issues, switch MAP.non.parallel to TRUE
    
    
    MAP_local <- function(i){
        
        write.table(paste0("Cluster in process is: ",i),file=paste0(getwd(),"/",output_folder_name,"/log_CM.txt"),append=TRUE,sep="");
        
        cell_ind <- which(z_inferred_final == i);
        mu_temp <- matrix(0,gene_batch*num_gene_batches, length(cell_ind));
        #Sigma_temp <- matrix(0,gene_batch*num_gene_batches,gene_batch*num_gene_batches); # very bad idea as you need this one for each cell!
        Sigma_temp <- list();
        for (j in 1:length(cell_ind)){
            p <- paste0("cell is ", cell_ind[j]);
            print(p);
            write.table(p,file=paste0(getwd(),"/",output_folder_name,"/log_CM.txt"),append=TRUE,sep="");
            count <- 1
            Sigma_temp_per_cell <- matrix(0,gene_batch*num_gene_batches,gene_batch*num_gene_batches);
            for (r1 in 1:length(results.all.MCMC)){
                print(paste("r1", r1))
                for(r2 in 1:num_gene_sub_batches){
                    print(paste("r2", r2))
                    
                    inferred_label_per_batch <- results.all.MCMC[[r1]][[r2]][[1]][cell_ind[j]];
                    
                    pos.mu.temp <- (1+(gene_batch*(count-1))):(gene_batch*count)
                    mu_temp[pos.mu.temp,j] <- results.all.MCMC[[r1]][[r2]][[4]][,as.character(inferred_label_per_batch)];
                    
                    
                    Sigma_temp_per_cell[pos.mu.temp,pos.mu.temp] <- results.all.MCMC[[r1]][[r2]][[5]][,,as.character(inferred_label_per_batch)];
                    count <- count + 1;
                    print(paste("count is ",count))
                }
            }
            Sigma_temp[[j]] <- Sigma_temp_per_cell;
            
        }
        mu_MAP_per_k <- rowMeans(mu_temp);
        
        Sigma_MAP_per_k <- Reduce('+',Sigma_temp)/length(cell_ind);
        diag(Sigma_MAP_per_k) <- diag(Sigma_MAP_per_k) + 0.001;
        
        
        return(list(mu_MAP_per_k, Sigma_MAP_per_k))
        
    }
    
    print("Compute MAP estimates of moments per cluster")
    print("Monitor debug_CM.txt for parallel process diagnostics and log_CM.txt for MAP progression")
    
    # Start the cluster and register with doSNOW
    cl <- makeCluster(num_cores, type = "SOCK",outfile="debug_CM.txt") #opens multiple socket connections
    clusterExport(cl, "MAP_local")
    registerDoSNOW(cl)
    
    # Call parallel processes per cluster to construct MAP estimates for mean and covariances
    strt_MAP <- Sys.time();
    
    MAP_estimates <- foreach (i = 1:mean_num_clusters) %dopar%{
        local_MAP_batch <- MAP_local(i)
        return(local_MAP_batch)
    }
    
    stopCluster(cl)
    
    print("Time for creating the MAP estimates")
    print(Sys.time()-strt_MAP)
    
    mu_MAP <- matrix(0, gene_batch*num_gene_batches, mean_num_clusters); #or numgenes?
    Sigma_MAP <- list();
    for (k in 1:mean_num_clusters){
        mu_MAP[,k] <- MAP_estimates[[k]][[1]]
        Sigma_MAP[[k]] <- MAP_estimates[[k]][[2]]
    }

} # end of parallel MAP building

## Use mu_MAP and Sigma_MAP to sample for moment estimates
###########################
## Construct the cluster moments using riwish and rmnorm
## Increase draws for averaged moments
print("Resample moments per cluster")

mu_final <- matrix(0, gene_batch*num_gene_batches, mean_num_clusters); #or numgenes?
Sigma_final <- list();
draws <- 1;

for( i in 1:mean_num_clusters){
    
        t <- paste0("Sigma_",i)
        t1 <- paste0("Heatmap for Sigma_",i)
        
        print(paste0("cluster is ", i));
        
        Sigma_MAP_temp <- list();
        
        for ( dr in 1:draws){
            Sigma_MAP_temp[[dr]] <- riwish(numgenes+2,forceSymmetric(Sigma_MAP[[i]]));
            diag(Sigma_MAP_temp[[dr]]) <- diag(Sigma_MAP_temp[[dr]]) + 0.001
        }
        
        Sigma_final[[i]] <- Reduce('+', Sigma_MAP_temp)/draws;
        
        
        f <- paste0(getwd(),"/",output_folder_name,"/plots/Inferred_Sigmas/Sigma_final_",i,".txt");
        write.matrix(Sigma_final[[i]],file=f,sep="\t")


           
        if(draws == 1){
             mu_final[,i] <- rmnorm(draws,mu_MAP[,i],Sigma_final[[i]])
        }else{
             mu_final[,i] <- colMeans(rmnorm(draws,mu_MAP[,i],Sigma_final[[i]]))
        }
    
}

####collect mu_final and Sigma_final

f <- paste0(getwd(),"/",output_folder_name,"/plots/Inferred_means/mu_final.txt");
write.matrix(mu_final,file=f,sep="\t")





















