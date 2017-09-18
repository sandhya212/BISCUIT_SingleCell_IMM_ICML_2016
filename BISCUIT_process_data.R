## 21st Dec 2016
## BISCUIT preparing input data
## (example code for Zeisel et al mouse cortex data
## Downloaded from : https://storage.googleapis.com/linnarsson-lab-www-blobs/blobs/cortex/expression_mRNA_17-Aug-2014.txt
##
## Code author SP

## 28th Dec 2016
## choose genes if needed

## 14th March 2017
## pre-processing data
##

## 18th April 2017
## added logic to split genes based on co-variances rather than variance alone
##
##
##############################################################################################################
#### Modify this portion to accommodate your input data.                                                  ####
##############################################################################################################


print("Loading Data")

########
if (input_file_name=="expression_mRNA_17-Aug-2014.txt"){
    full.data = data.frame(read.csv(input_file_name, header=TRUE, sep="\t", quote="",stringsAsFactors = TRUE)); # genes x cells
    dim(full.data)
    colnames(full.data); #names(full.data),
    rownames(full.data);
    
    gene_names <- full.data[11:nrow(full.data),1];
    
    #creating the cellsXgenes data
    num_rows <- length(c(11:nrow(full.data)))
    num_cols <- length(c(3:ncol(full.data)))
    
    #full.data.1 <- as.matrix(full.data[11:nrow(full.data),3:ncol(full.data)]);
    full.data.1 <- matrix(as.numeric(as.matrix(full.data[11:nrow(full.data),3:ncol(full.data)])),num_rows,num_cols);
    dim(full.data.1);
    
    #getting the true labels
    if (z_true_labels_avl){
            z_true <- as.numeric(factor(unlist(full.data[8,][3:(dim(full.data)[2])])));
    }
    
    
    full.data.1 <- t(full.data.1); #cells x genes
    
    
}else{ #this assumes the input data has both column and row names.
    #full.data <- data.frame(read.csv(input_file_name, header=TRUE, row.names=1, sep=",",stringsAsFactors = TRUE));
    
    if(input_data_tab_delimited == TRUE){
        full.data <- data.frame(read.csv(input_file_name, header=TRUE, row.names=1, sep="\t",stringsAsFactors = TRUE));
    }else{ #comma-separated input data
        full.data <- data.frame(read.csv(input_file_name, header=TRUE, row.names=1, sep=",",stringsAsFactors = TRUE));
    }
    
    if(is_format_genes_cells == TRUE){
        full.data <- t(full.data) #cellsxgenes
    }
    
       
    dim(full.data)
    gene_names <- colnames(full.data); #gene_names
    rownames(full.data); # #cell_names
    
    #creating the cellsXgenes data
    full.data.1 <- as.matrix(full.data);
    
}




dim(full.data.1);
full.data.1[is.na(full.data.1)] <- 0;




#save(full.data.1,file="full.data.1.RData")

#Idea 1 to get meaningful genes: choose genes based on highest global co-expression
stddev.genes <- apply(full.data.1,2,sd);## find std dev of genes
f <- paste0(getwd(),"/",output_folder_name,"/plots/Stddev_genes.pdf")
pdf(file=f);
plot(sort(stddev.genes,decreasing=T),typ='o')
dev.off();
#full.data.2 <- full.data.1[,order(stddev.genes,decreasing=TRUE)];
#gene_names <- gene_names[order(stddev.genes,decreasing=T)]; 

#Idea 2 to get meaningful genes: perform disparity check between rows and columns
emp.cov <- cov(full.data.1);
diag(emp.cov) <- 0;
rowsums.emp.cov <- rowSums(emp.cov);
colsums.emp.cov <- colSums(emp.cov);
gene.disparity <- rowsums.emp.cov + colsums.emp.cov
f <- paste0(getwd(),"/",output_folder_name,"/plots/Disparity_genes.pdf")
pdf(file=f);
plot(sort(gene.disparity,decreasing=T),typ='o')
dev.off();
#full.data.2 <- full.data.1[,order(gene.disparity,decreasing=TRUE)];
#gene_names <- gene_names[order(gene.disparity,decreasing=T)];

#Idea 3 to get meaningful genes: use Fiedler vector to split genes (graph theoretic partition)
print("Calculating the Fiedler vector of the data")


if (input_file_name=="expression_mRNA_17-Aug-2014.txt"){
    load(file="fvec.RData")
}else{
    L.mat <- diag(colsums.emp.cov) - emp.cov; #ensure L.mat is singular p.s.d
    f.vec <- fiedler.vector(L.mat);
    save(f.vec,file="fvec.RData");
}


#load(file="fvec.RData")
full.data.2 <- full.data.1[,order(f.vec,decreasing=TRUE)]
gene_names <- gene_names[order(f.vec,decreasing=TRUE)];

#choose cells

lib_size <- rowSums(full.data.2);

f <- paste0(getwd(),"/",output_folder_name,"/plots/lib_size_hist.pdf")
pdf(file=f);
par(mfrow=c(2,1))
plot(sort(lib_size,decreasing=T),typ='l')
hist(lib_size,breaks=300)
dev.off();

if (input_file_name=="MERGEDtumors_subsetgenes_counts.csv"){
    full.data.2 <- full.data.2[which(lib_size>400),]
}

##log transform data

print("Ensuring entire data is numeric and then log transforming it")

X_1 <- log(full.data.2+1);

##
if(exists("choose_cells")){
    numcells <- choose_cells;
}else{
    numcells <- dim(X_1)[1];
}
print(paste("numcells is", numcells))

if(exists("choose_genes")){
    numgenes <- choose_genes;
}else{
    numgenes <- dim(X_1)[2];
}
print(paste("numgenes is", numgenes))
###

##write the genes used in this run into a file
##$$$$
f <- paste0(getwd(),"/",output_folder_name,"/plots/Inferred_means/Genes_selected.csv")
write.csv(gene_names[1:numgenes], file=f);
f <- paste0(getwd(),"/",output_folder_name,"/plots/Inferred_Sigmas/Genes_selected.csv")
write.csv(gene_names[1:numgenes], file=f);


#############################################################################################
#############################################################################################
#############################################################################################

tot_numgenes <- dim(full.data.1)[2];
rm(full.data.1)
num_gene_batches <- floor(numgenes/gene_batch);
print(paste0('Number of gene batches is ', num_gene_batches));

numgenes <- num_gene_batches * gene_batch;

num_gene_sub_batches <- sub_batch(num_gene_batches);
print(paste0('Number of gene subbatches is ', num_gene_sub_batches));


#############################################################################################
#############################################################################################
#############################################################################################
# preparing the dataset as per user-defined number of cells and genes


    
full.data.3 <- full.data.2[1:numcells,1:numgenes];

## Ensure data is numeric
print("Ensuring user-specified data is numeric")

X_all <- matrix(as.numeric(full.data.3),nrow=numcells,ncol=numgenes);
X_all <- log(X_all + 0.1); # log normalisation. + 0.1 to account for zero entries in X that cannot be log transformed.

###
# Visualisation

## centering X and plotting
X_c_all<- center_colmeans(X_all);
N <- numcells;
D <- numgenes;
n <- rep(0,N)
for( i in 1:N){
    n[i] <- norm_vec(X_c_all[i,])
}
#X_c_norm_all <- X_c_all/max(n)

print('Computing t-sne projection of the data')
X_tsne_all <- Rtsne(X_all,check_duplicates = FALSE);


## plotting standardised X
#X_std_all <- project.data(X_all,D);


##Global normalised data

log_lib_size <- rowSums(X_all);
X_all_global_norm <- X_all/(log_lib_size + 0.0001);
X_tsne_all_global_norm <- Rtsne(X_all_global_norm,check_duplicates = FALSE);

#rm(X_all)

