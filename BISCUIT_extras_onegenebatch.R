## 24th May 2017
## Post processing BISCUIT clusters in the case of a single gene split
##
## Code author SP
##################




print("Recomputing cluster probabilities of each cell to its cluster")
cluster_prob <- rep(0,numcells);
mean_shift_per_K <- matrix(0,gene_batch*num_gene_batches,final_num_K);
colnames(mean_shift_per_K) <- as.character(uq_z);


mean_rot_per_K <- list();

    for (i in 1:final_num_K){
        mean_shift_per_K[,uq_z[i]] <- mean_alpha_inferred_per_K[uq_z[i],] * mu_final[,uq_z[i]];
        mean_rot_per_K[[which(uq_z==uq_z[i])]] <- matrix(forceSymmetric(mean_beta_inferred_per_K[uq_z[i],] * Sigma_final[,,uq_z[i]]))
    }

###

for (i in 1:numcells){
        g <- z_inferred_final[i];
        cluster_prob[i] <- dmnorm(X_all[i,],mean_shift_per_K[,g],mean_rot_per_K[[which(uq_z==g)]])[1]
}


cluster_prob <- cbind(z_inferred_final_plot[1:numcells],cluster_prob)
colnames(cluster_prob) <- c("z_inferred", "cluster_probability")

## write to text file.

f <- paste0(getwd(),"/output/plots/extras/cluster_probabilities.csv")
write.csv(cluster_prob, file=f);

if(z_true_labels_avl){
    
    f <- paste0(getwd(),"/output/plots/Inferred_labels/Final_true_inferred_labels_preimputed_X.pdf");
    pdf(file=f);
    par(mfrow=c(2,1))
    plot(X_tsne_all$Y[,1],X_tsne_all$Y[,2],col = col_palette[1*(z_true)],  main="t-SNE of pre-imputed X (true labels)");
    plot(X_tsne_all$Y[,1],X_tsne_all$Y[,2],col = col_palette[1*(z_inferred_final_plot)],  main="t-SNE of pre-imputed X (inferred labels)");
    dev.off()
    
    f <- paste0(getwd(),"/output/plots/Inferred_labels/Final_true_inferred_labels_imputed_X.pdf");
    pdf(file=f);
    par(mfrow=c(2,1))
    plot(Y_tsne$Y[,1],Y_tsne$Y[,2],col = col_palette[1*(z_true)],  main="t-SNE of imputed X (true labels)");
    plot(Y_tsne$Y[,1],Y_tsne$Y[,2],col = col_palette[1*(z_inferred_final_plot)],  main="t-SNE of imputed X (inferred labels)");
    dev.off()
    
    f <- paste0(getwd(),"/output/plots/Inferred_labels/Final_true_inferred_labels_globalnorm_X.pdf");
    pdf(file=f);
    par(mfrow=c(2,1))
    plot(X_tsne_all_global_norm$Y[,1],X_tsne_all_global_norm$Y[,2],col = col_palette[1*(z_true)],  main="t-SNE of global normalised X (true labels)");
    plot(X_tsne_all_global_norm$Y[,1],X_tsne_all_global_norm$Y[,2],col = col_palette[1*(z_inferred_final_plot)],  main="t-SNE of global normalised X (inferred labels)");
    dev.off()
    
}
