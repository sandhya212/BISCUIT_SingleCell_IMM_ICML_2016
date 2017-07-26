## 16th May 2017
## Post processing BISCUIT clusters
##
## Code author SP
##################



print("Recomputing cluster probabilities of each cell to its cluster")
cluster_prob <- rep(0,numcells);
mean_shift_per_K <- matrix(0,gene_batch*num_gene_batches,final_num_K);
mean_rot_per_K <- list();

    for (i in 1:final_num_K){
        mean_shift_per_K[,i] <- mean_alpha_inferred_per_K[i] * mu_final[,i];
        mean_rot_per_K[[i]] <- matrix(forceSymmetric(mean_beta_inferred_per_K[i] * Sigma_final[[i]]))
    }

###

for (i in 1:numcells){
        g <- z_inferred_final[i];
        cluster_prob[i] <- dmnorm(X_all[i,],mean_shift_per_K[,g],mean_rot_per_K[[g]])[1]
}


cluster_prob <- cbind(z_inferred_final[1:numcells],cluster_prob)
colnames(cluster_prob) <- c("z_inferred", "cluster_probability")

## write to text file.

f <- paste0(getwd(),"/output/plots/extras/cluster_probabilities.csv")
write.csv(cluster_prob, file=f);

if(z_true_labels_avl){
    
    f <- paste0(getwd(),"/output/plots/Inferred_labels/Final_true_inferred_labels_preimputed_X.pdf");
    pdf(file=f);
    par(mfrow=c(2,1))
    plot(X_tsne_all$Y[,1],X_tsne_all$Y[,2],col = col_palette[1*(z_true)],  main="t-SNE of pre-imputed X (true labels)");
    plot(X_tsne_all$Y[,1],X_tsne_all$Y[,2],col = col_palette[1*(z_inferred_final)],  main="t-SNE of pre-imputed X (inferred labels)");
    dev.off()
    
    f <- paste0(getwd(),"/output/plots/Inferred_labels/Final_true_inferred_labels_imputed_X.pdf");
    pdf(file=f);
    par(mfrow=c(2,1))
    plot(Y_tsne$Y[,1],Y_tsne$Y[,2],col = col_palette[1*(z_true)],  main="t-SNE of imputed X (true labels)");
    plot(Y_tsne$Y[,1],Y_tsne$Y[,2],col = col_palette[1*(z_inferred_final)],  main="t-SNE of imputed X (inferred labels)");
    dev.off()
    
    f <- paste0(getwd(),"/output/plots/Inferred_labels/Final_true_inferred_labels_globalnorm_X.pdf");
    pdf(file=f);
    par(mfrow=c(2,1))
    plot(X_tsne_all_global_norm$Y[,1],X_tsne_all_global_norm$Y[,2],col = col_palette[1*(z_true)],  main="t-SNE of global normalised X (true labels)");
    plot(X_tsne_all_global_norm$Y[,1],X_tsne_all_global_norm$Y[,2],col = col_palette[1*(z_inferred_final)],  main="t-SNE of global normalised X (inferred labels)");
    dev.off()
    
}

## write all input parameters to text file

f1 <- paste0(getwd(),"/output/plots/extras/Input_parms_used.txt")
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
write(paste('Number of clusters per batch: ',num_clusters_per_batch),file=f1, append=TRUE)
write(paste('Mean number of clusters: ',mean_num_clusters),file=f1, append=TRUE)
write(paste('Time for MCMC iterations: ',MCMC_time-strt),file=f1, append=TRUE)
