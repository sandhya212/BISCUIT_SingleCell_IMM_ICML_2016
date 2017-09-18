## 20th November 2015
## Reverse transform X to Y based on inferred parameters
##
## 3rd April 2017
## coding imputing based on gene-split code
##
## 24th April 2017
## Computing A based on inferred betas alone
##
## 25th April 2017
## Parallelising the code
##
## 27th April 2017
## Removed the Sigma-based imputing
##
## Code author SP
##################


X_std_all <- X_all;

final_num_K <- length(unique(z_inferred_final));

mean_alpha_inferred <- mean(alpha_inferred_final)
mean_beta_inferred <- mean(beta_inferred_final)

mean_alpha_inferred_per_K <- rep(0,final_num_K)
mean_beta_inferred_per_K <- rep(0,final_num_K)

A_rt <- matrix(0,final_num_K,numgenes)

##
#### compute cluster-based means of alphas and betas
##
print("Computing predictor matrix A and cluster-based means of alphas and betas")
for (i in 1:final_num_K){
    
    #cell_ids <- which(z_inferred_final == final_K[i]);
    cell_ids <- which(z_inferred_final == i);
    mean_alpha_inferred_per_K[i] <- median(alpha_inferred_final[cell_ids])
    mean_beta_inferred_per_K[i] <- median(beta_inferred_final[cell_ids])
    
    A_rt[i,] <- rep(1/sqrt(mean_beta_inferred_per_K[i]),numgenes);
}




print("Computing imputed data based on inferred alphas, betas, mu_k, Sigma_k and z")

Impute_batches <- function(df){
    
    
    df_indicator <- paste0("Imputed batch in process is: ",df);
    write.table(df_indicator,file=paste0(getwd(),"/",output_folder_name,"/log_CM.txt"),append=TRUE,sep="");



    Y_rt <- matrix(0,length(((1+(num_cells_batch*(df-1))):(num_cells_batch*df))), numgenes)
    



    rownames(Y_rt) <- as.character(((1+(num_cells_batch*(df-1))):(num_cells_batch*df)))
    colnames(Y_rt) <- as.character(1:numgenes)
   
    
    
    for (cell_ind in((1+(num_cells_batch*(df-1))):(num_cells_batch*df))){
        #print(paste("imputing cell",cell_ind));
    
        g <- z_inferred_final[cell_ind];
        
        b <- (diag(numgenes) - mean_alpha_inferred_per_K[g]*diag(A_rt[g,])) %*%(mu_final[,g])
        Y_rt[as.character(cell_ind),] <- t(diag(A_rt[g,]) %*% X_std_all[cell_ind,] + b)

    }
    return(list(Y_rt))

}




strt_imp <- Sys.time()


divsr <- numcells %/% num_cells_batch
dividend <- numcells %% num_cells_batch



# Start the cluster and register with doSNOW
cl <- makeCluster(num_cores, type = "SOCK",outfile="debug_CM.txt") #opens multiple socket connections
clusterExport(cl, "Impute_batches")
registerDoSNOW(cl)

# Call parallel processes to build imputed matrix


global_imputation <- foreach (df = 1:(divsr)) %dopar%{
    local_impute_batch <- Impute_batches(df)
    return(list(local_impute_batch[[1]]))
}

print("Time for imputing is ")
print(Sys.time()-strt_imp)

stopCluster(cl)

# Create the Y_rt matrix

Y_rt_final <- rep(0, numgenes)


for (df in 1:(divsr)){

    Y_rt_final <- rbind(Y_rt_final,global_imputation[[df]][[1]]);

}
    


##Taking care of the remaining cells that fell off the parallel bins

if (dividend !=0){
   
   
    df_indicator <- paste0("Imputed Batch in process is: ",divsr+1);
    write.table(df_indicator,file=paste0(getwd(),"/",output_folder_name,"/log_CM.txt"),append=TRUE,sep="");
    
    write.table("Cells selected are ",file=paste0(getwd(),"/",output_folder_name,"/log_CM.txt"),append=TRUE,sep="");
    write.table(((num_cells_batch*divsr+1):(num_cells_batch*divsr+dividend)),file=paste0(getwd(),"/",output_folder_name,"/log_CM.txt"),append=TRUE,sep="");
    Y_rt <- matrix(0, dividend, numgenes)

    rownames(Y_rt) <- as.character((num_cells_batch*divsr+1):(num_cells_batch*divsr+dividend))
    colnames(Y_rt) <- as.character(1:numgenes)


    for ( cell_ind in (num_cells_batch*divsr+1):(num_cells_batch*divsr+dividend)){
        #print(paste("imputing cell",cell_ind));
        
        g <- z_inferred_final[cell_ind];
        
        b <- (diag(numgenes) - mean_alpha_inferred_per_K[g]*diag(A_rt[g,])) %*%(mu_final[,g])
        Y_rt[as.character(cell_ind),] <- t(diag(A_rt[g,]) %*% X_std_all[cell_ind,] + b)
        

    }
    
    

    Y_rt_final <- rbind(Y_rt_final,Y_rt);

    
}



##strip off first row of Y_rt_final
Y_rt_final <- Y_rt_final[-1,]

##write the imputed matrix, still in logspace
f <- paste0(getwd(),"/",output_folder_name,"/plots/extras/Imputed_Y_logspace.txt");
write.matrix(Y_rt_final,file=f,sep="\t")

##write the imputed matrix, in countspace
Y_rt_count_space <- exp(Y_rt_final) - 0.1;
f <- paste0(getwd(),"/",output_folder_name,"/plots/extras/Imputed_Y_countspace.txt");
write.matrix(Y_rt_count_space,file=f,sep="\t")


## plotting tSNE of Y
print('Computing t-sne projection of the imputed data')


######


Y_tsne <- Rtsne(Y_rt_final,check_duplicates = FALSE);

#rm(Y_rt_final)
f <- paste0(getwd(),"/",output_folder_name,"/plots/Inferred_labels/Final_inferred_labels_imputed_X.pdf");
pdf(file=f);
plot(Y_tsne$Y[,1],Y_tsne$Y[,2],col = col_palette[1*(z_inferred_final)],  main="t-SNE of imputed X (inferred labels)");
dev.off()

f <- paste0(getwd(),"/",output_folder_name,"/plots/Inferred_labels/Final_inferred_labels_prepost_imputed_X.pdf");
pdf(file=f);
par(mfrow=c(2,1))
plot(X_tsne_all$Y[,1],X_tsne_all$Y[,2],col = col_palette[1*(z_inferred_final)],  main="t-SNE of pre-imputed X (inferred labels)");
plot(Y_tsne$Y[,1],Y_tsne$Y[,2],col = col_palette[1*(z_inferred_final)],  main="t-SNE of imputed X (inferred labels)");
dev.off()


f <- paste0(getwd(),"/",output_folder_name,"/plots/Inferred_labels/Final_inferred_labels_globalnorm_post_imputed_X.pdf");
pdf(file=f);
par(mfrow=c(2,1))
plot(X_tsne_all_global_norm$Y[,1],X_tsne_all_global_norm$Y[,2],col = col_palette[1*(z_inferred_final)],  main="t-SNE of global normalised X (inferred labels)");
plot(Y_tsne$Y[,1],Y_tsne$Y[,2],col = col_palette[1*(z_inferred_final)],  main="t-SNE of imputed X (inferred labels)");
dev.off()

f <- paste0(getwd(),"/",output_folder_name,"/plots/Inferred_labels/Final_inferred_labels_globalnorm_X.pdf");
pdf(file=f);
plot(X_tsne_all_global_norm$Y[,1],X_tsne_all_global_norm$Y[,2],col = col_palette[1*(z_inferred_final)],  main="t-SNE of global normalised X (inferred labels)");
dev.off()



##
##collecting the pre and post imputed tSNE coordinates
write.matrix(X_tsne_all$Y, file=paste0(working_path,"/",output_folder_name,"/plots/extras/pre_imputed_tSNE_coord.txt"),sep = "\t")

write.matrix(Y_tsne$Y, file=paste0(working_path,"/",output_folder_name,"/plots/extras/post_imputed_tSNE_coord.txt"),sep = "\t")




