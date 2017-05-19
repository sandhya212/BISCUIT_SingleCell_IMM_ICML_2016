## 16th May 2017
## Post processing BISCUIT clusters
##
## Code author SP
##################


#print("Computing soft cluster probabilities of each cell to total inferred clusters")
#soft_cluster_prob <- matrix(0,numcells,final_num_K+1);

#for (i in 1:numcells){
#    g <- z_inferred_final[i];
#    print(paste("cell  is: ", i));
#    for (k in 1:final_num_K){
        
#        text6 <- (paste('soft cluster assignments for cluster: ', k));
#        write(text6,file=paste0(getwd(),"/output/log.txt"),append=TRUE,sep="");
#       print(paste("k is: ", k));
#       soft_cluster_prob[i,k] <- dmnorm(Y_rt_final[i,],alpha_inferred_final[i] * mu_final[,k],beta_inferred_final[i] * matrix(forceSymmetric(Sigma_final[[k]])))[1]

#    }
#    soft_cluster_prob[i,final_num_K+1] <- which(soft_cluster_prob[i,1:final_num_K]==max(soft_cluster_prob[i,1:final_num_K]))
# }

#soft_cluster_prob <- cbind(z_inferred_final[1:numcells],soft_cluster_prob)

## write to text file.
#f <- paste0(getwd(),"/output/plots/extras/soft_cluster_probabilities.txt")
#write.matrix(soft_cluster_prob, file=f,sep="/t");


print("Recomputing cluster probabilities of each cell to its cluster")
cluster_prob <- matrix(0,numcells,2);

for (i in 1:numcells){
    g <- z_inferred_final[i];
    print(paste("cell  is: ", i));
    cluster_prob[i,2] <- dmnorm(Y_rt_final[i,],alpha_inferred_final[i] * mu_final[,g],beta_inferred_final[i] * matrix(forceSymmetric(Sigma_final[[g]])))[1]
    cluster_prob[i,1] <- g;
}
## write to text file.

f <- paste0(getwd(),"/output/plots/extras/cluster_probabilities.txt")
write.matrix(cluster_prob, file=f,sep="/t");


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
