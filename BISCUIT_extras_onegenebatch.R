## 24th May 2017
## Post processing BISCUIT clusters in the case of a single gene split
##
## Code author SP
##
## 27th Feb 2018
## Adding code for a cluster probability matrix (cell x cluster)
##
##################




print("Recomputing cluster probabilities of each cell to its cluster")
cluster_prob <- rep(0,numcells);
mean_shift_per_K <- matrix(0,gene_batch*num_gene_batches,final_num_K);

# Use this for subsetting the columns (i.e. clusters)
cnames <- dimnames(mu_final)[[2]]
cnames <- cnames[!cnames == "dummy"]
colnames(mean_shift_per_K) <- cnames; #as.character(uq_z);
mean_rot_per_K <- Sigma_final

for (i in 1:final_num_K){
        # mean_shift_per_K[,uq_z[i]] <- mean_alpha_inferred_per_K[uq_z[i],] * mu_final[,uq_z[i]];
        # mean_shift_per_K[,uq_z[i]] <-  mu_final[,uq_z[i]];
        mean_shift_per_K[,cnames[i]] <-  mu_final[,cnames[i]];

        mean_rot_per_K[,,cnames[i]] <- matrix(forceSymmetric(Sigma_final[,,cnames[i]]),numgenes,numgenes);
        #mean_rot_per_K[,,uq_z[i]] <- forceSymmetric(Sigma_final[,,uq_z[i]]);
}

###above we removed the dummy column from the moments

for (i in 1:numcells){
    print(i);
    g <- as.character(z_inferred_final[i]);
    #cluster_prob[i] <- dmnorm(X_all[i,],mean_shift_per_K[,g],mean_rot_per_K[[which(uq_z==g)]])[1];
    cluster_prob[i] <- dmnorm(X_all[i,],mean_shift_per_K[,g],mean_rot_per_K[,,g])[1]
}

## max for that cell (below) is not tallying with z_inferred_final for that cell (above block)


cluster_prob <- cbind(z_inferred_final_plot[1:numcells],cluster_prob)
colnames(cluster_prob) <- c("z_inferred", "cluster_probability")

## but we need to use the below block for downstream sans max(z_inferred_final) per cell (above)
# cluster probability matrix (cell x cluster)
cluster_prob_mat <- matrix(0,numcells, final_num_K);
colnames(cluster_prob_mat) <- cnames;

for (i in 1:numcells){
    print(i);
    for (k in 1:final_num_K){
        cluster_prob_mat[i,cnames[k]] <- dmnorm(X_all[i,],mean_shift_per_K[,cnames[k]],mean_rot_per_K[,,cnames[k]])[1]
    }
}



## write to text file.

f <- paste0(getwd(),"/",output_folder_name,"/plots/extras/cluster_probabilities.csv")
write.csv(cluster_prob, file=f);

f <- paste0(getwd(),"/",output_folder_name,"/plots/extras/cluster_posterior_probabilities.csv")
write.csv(cluster_prob_mat, file=f);

write.csv(X_tsne_all_global_norm$Y,file=paste0(getwd(),"/",output_folder_name,"/plots/extras/X_tsne_all_global_norm.csv"))
write.csv(X_tsne_all$Y,file=paste0(getwd(),"/",output_folder_name,"/plots/extras/X_tsne_all.csv"))

if(z_true_labels_avl){
    
    f <- paste0(getwd(),"/",output_folder_name,"/plots/Inferred_labels/Final_true_inferred_labels_preimputed_X.pdf");
    pdf(file=f);
    par(mfrow=c(2,1))
    plot(X_tsne_all$Y[,1],X_tsne_all$Y[,2],col = col_palette[1*(z_true)],  main="t-SNE of pre-imputed X (true labels)");
    plot(X_tsne_all$Y[,1],X_tsne_all$Y[,2],col = col_palette[1*(z_inferred_final_plot)],  main="t-SNE of pre-imputed X (inferred labels)");
    dev.off()
    
    f <- paste0(getwd(),"/",output_folder_name,"/plots/Inferred_labels/Final_true_inferred_labels_imputed_X.pdf");
    pdf(file=f);
    par(mfrow=c(2,1))
    plot(Y_tsne$Y[,1],Y_tsne$Y[,2],col = col_palette[1*(z_true)],  main="t-SNE of imputed X (true labels)");
    plot(Y_tsne$Y[,1],Y_tsne$Y[,2],col = col_palette[1*(z_inferred_final_plot)],  main="t-SNE of imputed X (inferred labels)");
    dev.off()
    
    f <- paste0(getwd(),"/",output_folder_name,"/plots/Inferred_labels/Final_true_inferred_labels_globalnorm_X.pdf");
    pdf(file=f);
    par(mfrow=c(2,1))
    plot(X_tsne_all_global_norm$Y[,1],X_tsne_all_global_norm$Y[,2],col = col_palette[1*(z_true)],  main="t-SNE of global normalised X (true labels)");
    plot(X_tsne_all_global_norm$Y[,1],X_tsne_all_global_norm$Y[,2],col = col_palette[1*(z_inferred_final_plot)],  main="t-SNE of global normalised X (inferred labels)");
    dev.off()
    
}

## write all input parameters to text file

f1 <- paste0(getwd(),"/",output_folder_name,"/plots/extras/Input_parms_used.txt")
write('********************************************** ',file=f1, append=TRUE)
write('***** Input parameters and settings used ***** ',file=f1, append=TRUE)
write('********************************************** ',file=f1, append=TRUE)
write('',file=f1, append=TRUE)

write(paste('Number of cells chosen: ',numcells),file=f1, append=TRUE)
write(paste('Number of genes chosen: ',numgenes),file=f1, append=TRUE)
write(paste('Number of MCMC iterations: ',num_iter),file=f1, append=TRUE)
write(paste('Number of cells per batch: ',num_cells_batch),file=f1, append=TRUE)
write(paste('Number of genes per batch: ',gene_batch),file=f1, append=TRUE)
write(paste('Dispersion parameter, alpha: ',alpha),file=f1, append=TRUE)
write(paste('Number of parallel gene batches: ',num_gene_batches),file=f1, append=TRUE)
write(paste('Number of parallel gene subbatches: ',num_gene_sub_batches),file=f1, append=TRUE)
write(paste('Number of clusters: ',final_num_K),file=f1, append=TRUE)
write(paste('Time for MCMC iterations: ',MCMC_time-strt),file=f1, append=TRUE)


## write all output values to text file

f1 <- paste0(getwd(),"/",output_folder_name,"/plots/extras/Output_values.txt")
write('********************************************** ',file=f1, append=TRUE)
write('***** Output values ***** ',file=f1, append=TRUE)
write('********************************************** ',file=f1, append=TRUE)
write('',file=f1, append=TRUE)

write(paste('Number of clusters via DPMM ',final_num_K),file=f1, append=TRUE)
write(paste('cluster distribution ',table(z_inferred_final_plot)),file=f1, append=TRUE)

