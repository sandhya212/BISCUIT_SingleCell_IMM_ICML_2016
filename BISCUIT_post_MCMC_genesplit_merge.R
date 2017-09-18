## 7th April 2017
## Merging the gene splits
## via Set cover approximation (NP complete)
##
##
## Code author SP



######
# Combine the runs/gene_batches


if(num_gene_batches > 1){
    alpha_inferred_final <- matrix(0,numcells, 1);
    beta_inferred_final <- matrix(0,numcells,1);
    
    for(a in 1:length(results.all.MCMC)){
        for(b in 1:num_gene_sub_batches){
                for (j in 1:N){
                    alpha_inferred_final[j] <- alpha_inferred_final[j] + as.numeric(results.all.MCMC[[a]][[b]][[2]])[j];
                    
                    beta_inferred_final[j] <- beta_inferred_final[j] + as.numeric(results.all.MCMC[[a]][[b]][[3]])[j];
                }
        }
    }
    
    print("Merging gene splits")
    #source("Confusion_matrix_stitch_old_nonparallel_MAP.R")
    source("Confusion_matrix_stitch_parallel_MAP.R")

}else{
    
   
    mu_final <- results.all.MCMC[[num_gene_batches]][[num_gene_sub_batches]][[4]]; #numgenes x K (c1,c2..)
    Sigma_final <- results.all.MCMC[[num_gene_batches]][[num_gene_sub_batches]][[5]]; #numgenes x numgenes X K (c1,c2..)
    
    
    alpha_inferred_final <- as.numeric(results.all.MCMC[[num_gene_batches]][[num_gene_sub_batches]][[2]]);
    beta_inferred_final <- as.numeric(results.all.MCMC[[num_gene_batches]][[num_gene_sub_batches]][[3]]);
   
    
    z_inferred_final_plot <- as.numeric(results.all.MCMC[[num_gene_batches]][[num_gene_sub_batches]][[1]]); # for classes that are missing z_inferred_final will be appropriately corrected. So the values in results.all.MCMC will not be in sync after as.numeric is applied. This is beneficial for plotting.
    
    z_inferred_final <- results.all.MCMC[[num_gene_batches]][[num_gene_sub_batches]][[1]]; # true labels as per logical subscripting.
    
    

    f <- paste0(getwd(),"/",output_folder_name,"/plots/Inferred_labels/Final_inferred_one_split.pdf");
    pdf(file=f);
    plot(X_tsne_all$Y[,1],X_tsne_all$Y[,2],col = col_palette[1*(z_inferred_final_plot)],  main="t-SNE of X (inferred labels)");
    dev.off();
}
########

alpha_inferred_final <- alpha_inferred_final/num_gene_batches;
beta_inferred_final <- beta_inferred_final/num_gene_batches;

## write alphas and betas to .csv

f <- paste0(getwd(),"/",output_folder_name,"/plots/Inferred_alphas_betas/Final_alphas.csv")
write.csv(alpha_inferred_final, file=f);
f <- paste0(getwd(),"/",output_folder_name,"/plots/Inferred_alphas_betas/Final_betas.csv")
write.csv(beta_inferred_final, file=f);


filename=paste0(getwd(),"/",output_folder_name,"/plots/Inferred_alphas_betas/Final_alpha_beta.pdf")
pdf(file=filename)
par(mfrow=c(1,2))
plot(alpha_inferred_final, type="l", main="alpha final spread");
plot(beta_inferred_final, type="l", main="beta final spread");
dev.off()


total_clusters <- length(which(tabulate(z_inferred_final)!=0));
total_singleton_clusters <- length(which(tabulate(z_inferred_final)==1));



##
min_a <- min(alpha_inferred_final)
max_a <- max(alpha_inferred_final)
mean_a <- mean(min_a,max_a)

f=paste0(getwd(),"/",output_folder_name,"/plots/Inferred_alphas_betas/Inferred_alpha_log_libsize.pdf")
pdf(file=f)
plot(sort(log_lib_size),alpha_inferred_final[order(log_lib_size)], main="alpha final spread",ylim=c(min_a-mean_a,max_a+mean_a),xlab ="log(library_size)",ylab="Inferred alphas");
dev.off()

##
min_b <- min(beta_inferred_final)
max_b <- max(beta_inferred_final)
mean_b <- mean(min_b,max_b)
diff_b <- max_b - min_b

f=paste0(getwd(),"/",output_folder_name,"/plots/Inferred_alphas_betas/Inferred_beta_log_libsize.pdf")
pdf(file=f)
plot(sort(log_lib_size), beta_inferred_final[order(log_lib_size)], main="beta final spread",ylim=c(min_b-diff_b,max_b+diff_b),xlab ="log(library_size)",ylab="Inferred betas");
dev.off()






