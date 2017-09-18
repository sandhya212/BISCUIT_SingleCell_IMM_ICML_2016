## 1st March 2017
## MCMC engine using NIW likelihood
##
## 7th April 2017
## Parallelising the gene splits
##
## 1st May 2017
## Gathering weights per gene batch
##
## Code author SP


####### Main MCMC Engine #######################
####### DO NOT change anything in IMM.MCMC() ###########

IMM.MCMC <- function(r){
    
    ##import process_data_rerun.R##
    #############################################################################################
    #############################################################################################
    #############################################################################################
    full.data.subsample <- full.data.2[1:numcells,(1+(gene_batch*(r-1))):(gene_batch*r)];
    
    r_indicator <- paste0("Batch in process is: ",r);
    write.table(r_indicator,file=paste0(getwd(),"/",output_folder_name,"/log.txt"),append=TRUE,sep="");
    write.table("Genes selected are ",file=paste0(getwd(),"/",output_folder_name,"/log.txt"),append=TRUE,sep="");
    write.table((1+(gene_batch*(r-1))):(gene_batch*r),file=paste0(getwd(),"/",output_folder_name,"/log.txt"),append=TRUE,sep="");
    ## Ensure data is numeric
    write.table("Ensuring data is numeric",file=paste0(getwd(),"/",output_folder_name,"/log.txt"),append=TRUE,sep="");
    
    #X <- matrix(as.numeric(retina.data.subsample),nrow=numcells,ncol=gene_batch);
    X <- full.data.subsample;
    X <- log(X+1); # log normalisation. +1 to account for zero entries in X that cannot be log transformed.
    
    ###std of X
    stddev.genes.batch <- sd(X);
    # Visualisation
    
    ## centering X and plotting
    X_c<- center_colmeans(X);
    N <- numcells;
    D <- gene_batch;
    n <- rep(0,N)
    for( i in 1:N){
        n[i] <- norm_vec(X_c[i,])
    }
    X_c_norm <- X_c/max(n)
    
    
    write.table("Computing projection of data",file=paste0(getwd(),"/",output_folder_name,"/log.txt"),append=TRUE,sep="");
    ## plotting standardised X
    #X_std <- project.data(X,D);


    ## plotting tSNE of X
    write.table('Computing t-sne projection of the data',file=paste0(getwd(),"/",output_folder_name,"/log.txt"),append=TRUE,sep="");
    X_tsne <- Rtsne(X,check_duplicates = FALSE,perplexity=10);



    ##import BISCUIT_ADVI_init_a_b_4.R##
    #############################################################################################
    #############################################################################################
    #############################################################################################
    
    write.table("Computing empirical mean and covariance of data",file=paste0(getwd(),"/",output_folder_name,"/log.txt"),append=TRUE,sep="");
    ##Hyperprior layer
    #mu_dprime <- colMeans(X_c_norm);
    #Sigma_dprime <- cov(X_c_norm);
    
    #mu_dprime <- colMeans(X_std);
    #Sigma_dprime <- cov(X_std);
    
    #alpha0 <- 10 #1/rgamma(1,1,scale=1); ## CHECK
    
    #sigma_prime  <- rep(1,D);
    
    write.table("Estimating initial values of alphas and betas", file=paste0(getwd(),"/",output_folder_name,"/log.txt"),append=TRUE,sep="");

    
    ###antilog of X
    X_counts <- exp(X_c_norm);
    lib_size <- rowSums(X_counts);
    med_libsize <- median(lib_size);
    mode_libsize <- getmode(lib_size);
    val1 <- min(med_libsize, mode_libsize);
    alpha_j_init <- lib_size/val1;
    
    
    # safe alpha
    min_a <- min(min(lib_size),min(lib_size)/val1);
    alpha_j_init[which(alpha_j_init==min(alpha_j_init))] <- min_a;
    
    min_alpha <- min(alpha_j_init);
    max_alpha <- max(alpha_j_init);
    sum_alpha_j <- sum(alpha_j_init);
   
    
    #alpha priors
    v = 1*mean(alpha_j_init); #0
    delta = 1000*var(alpha_j_init); # 1
    
    
    ###betas
    
    beta_j_init <- rep(1,N); #apply(X_counts,1,var); #beta_j; #
    min_beta <- min(beta_j_init);
    max_beta <- max(beta_j_init);
    sum_beta_j <- sum(beta_j_init);
    
    #beta priors (skewed inv-Gamma distribution)
    omega <- 1; #mean(beta_j_init);
    theta <- 1; #var(beta_j_init);
    
    ####Printing
    
    l <- paste0(getwd(),"/",output_folder_name,"/plots/Inferred_alphas_betas/Initial_alpha_beta_spread(batch)",r,".pdf");
    pdf(file=l)
    par(mfrow=c(1,2))
    plot(alpha_j_init, type="l", main="alpha init spread");
    plot(beta_j_init, type="l", main="beta init spread");
    dev.off()


    ####"Gibbs_DPGMM_BISCUIT_alpha_beta.R" for r = 1 and "Gibbs_DPGMM_BISCUIT_alpha_beta_propagate.R" for r > 1####
    #####################################################################################################
    #####################################################################################################
    #####################################################################################################
   
   
        
        
        path <- working_path;
        Rprof(filename = paste(path,"/Rprof.out",sep=""), append = FALSE, interval = 0.02, memory.profiling=TRUE)
        
        #choose the dimension of the model
        d <- dim(X)[2];
        #x <- X_std;
        #x <- X_c_norm;
        #x <- X_c
        x <- X

        
        N <- dim(x)[1] #total number of observations in all classes
        
        
        startTime <- Sys.time()
        
        
        #set initial values for the hyperparameters:
        # nu, Delta, mu_0, kappa, alpha
        nu <- d+1
        
        Delta <- matrix(0,d,d)
        #Delta <- Sigma_dprime;
        for(i in 1:d){Delta[i,i]<- 5}
        
        mu0 <- rep(1.5,d)
        #mu0 <- mu_dprime;
        
        kappa <- 1
        #alpha <- 0.5; # now moved to start_file.R
        
        #Choose initial values for :K, mu, Sigma, pi, C , N_k
        
        K <- 2
        Kit <- K
        
        mu <- rep(mu0,K+1);
        dim(mu) <- c(d,K+1)
        dimnames(mu)<-list( 1:d,c("c1","c2","dummy")) #dummy is here to ensure that mu and sigma have at least two classes otherwise dim(mu) becomes null
        
        cnames <- dimnames(mu)[[2]]
        cnames <- cnames[!cnames == "dummy"]
        
        Sigma <- rep( matrix(1,d,d),K+1)
        #Sigma <- rep(Sigma_dprime,K+1)
        dim(Sigma)<- c(d,d,K+1)
        dimnames(Sigma)<-list( 1:d,1:d,c("c1","c2","dummy"))
        
        diag22<- c(rep( c(1, rep(0,d)),d-1),1)
        dim(diag22) <- c(d,d)
        for (k in 1:K){
            Sigma[,,k] <- diag22
        }
        Pi<- rep( 1/K, K)
        Class <- rep(c("c1","c2"),floor(N/K)+1) #initialisation
        length(Class) <- N #cuts Class to the right size
        Class <- as.factor(Class)
        Nk<-table(Class)
        
        #production of the data frame
        mydata <- data.frame(observ=x,class=Class)
        
        
        
        # Update the parameters of the NIW distribution after observing data
        
        updateNIW <- function(fctx,fctNc,fctnu,fctkappa,fctmu0,fctDelta)
        {
            fctnup <- fctnu + fctNc
            fctkappap <- fctkappa + fctNc
            sumxk <- rep(0,d)
            sumxkmatrix <- matrix(0,d,d)
            if (is.null(dim(fctx))){
                sumxk <- fctx
                sumxkmatrix <- fctx %o% fctx
            }else{
                for(j in 1:length(fctx[,1])){
                    sumxk <- sumxk + fctx[j,]
                    sumxkmatrix <- sumxkmatrix + fctx[j,] %o% fctx[j,]
                }
            }
            fctmu0p <- (fctkappa*fctmu0 + sumxk) / fctkappap
            fctDeltap <-(fctDelta + sumxkmatrix + fctkappa * fctmu0 %o% fctmu0 - fctkappap * fctmu0p %o% fctmu0p )
            c(fctnup,fctkappap,fctmu0p,fctDeltap)
        }
        Nc <- Nk[1]
        
        
        inferred_parm2 <- rep(0,N*(1+1+d+d*d))
        dim(inferred_parm2) <- c(N,1+1+d+d*d)
        
        ############## Infer Scale parameters ##############
        
        omega_p <- omega + d/2;
        
        beta_j_inferred <- beta_j_init;
        alpha_j_inferred <- alpha_j_init;
        delta_psq_j <- matrix(0,1,N);
        v_p_j <- matrix(0,1,N);
        delta_sq <- delta;
        
        
        
        write.table("Computing initial NIW moments",file=paste0(getwd(),"/",output_folder_name,"/log.txt"),append=TRUE,sep="");
        
        for(i in 1:N){
            #print(i)
            inferred_parm2[i,] <- updateNIW(x[i,],1,nu,kappa,mu0,Delta)
        }
        
        
        timedmnorm <- 0
        timedmnorm2 <- 0
        mydata$class <- Class
        
        #####Gibbs sampling########
        
        steps <- num_iter
        alpha_beta_kickin <- 5
        step <- 0
        
        while(step <steps){ #begin of sampling loop
            step <- step +1
            
            text1 <- (paste("step is ", step));
            print(text1);
            #write(text1, file=paste0(getwd(),"/output/log.txt"),append=TRUE,sep="");
            
            #evaluate the dmnorm grouped by classes
            timedmnormstart2 <- Sys.time()
            
            # if (r == 1 | (r > 1 & step > 1)){
            Q <- rep(0,N*K)
            dim(Q) <- c(N,K)
            dimnames(Q)<- list(1:N,cnames)
            #}
            
            Q <- data.frame(Q,Class)
            for(j in 1:K){
                for(k in 1:K){
                    text2 <- (paste("cluster is", k));
                    print(text2);
                    #write(text2, file=paste0(getwd(),"/output/log.txt"),append=TRUE,sep="");
                    
                    diag(Sigma[,,cnames[k]]) <- diag(Sigma[,,cnames[k]]) + 0.001;
                    temp_Sigma <- matrix(forceSymmetric(Sigma[,,cnames[k]]),d,d);
                    
                    Q[Q$Class==cnames[j],cnames[k]] <- dmnorm(as.matrix(mydata[mydata$class==cnames[j],1:d]),mu[,cnames[k]],temp_Sigma) ##revert to this
                }
            }
            timedmnormend2 <- Sys.time()
            timedmnorm2 <- timedmnorm2 + difftime(timedmnormstart2,timedmnormend2,tz="",units="mins")
            cnamespreviousstep <- cnames
            
            
            #go through the observations and reassign classes
            for (i in 1:N){
                text3 <- paste("cell is ", i);
                print(text3);
                #write(text3, file=paste0(getwd(),"/output/log.txt"),append=TRUE,sep="");
                
                if(step > alpha_beta_kickin){
                    ##infer alpha and beta
                    z <- as.character(Class[i]);
                    theta_p <- abs(theta + 0.5*(x[i,]-alpha_j_inferred[i]*mu[,z])%*%Sigma[,,z]%*%t(t(x[i,]-alpha_j_inferred[i]*mu[,z]))/(d*d))
                    beta_j_inferred[i] <- rinvgamma(1,omega_p,theta_p)  ##
                    
                    
                    ## setting the sampled betas within constraints
                    if(beta_j_inferred[i] < min_beta){
                        beta_j_inferred[i] <- min_beta + rinvgamma(1,omega_p,1/omega_p)
                    }
                    
                    if(beta_j_inferred[i] > max_beta){
                        beta_j_inferred[i] <- max_beta - rinvgamma(1,omega_p,1/omega_p)
                    }
                    
                    
                    
                    
                    #beta_j_inferred[cell_ind] <- rinvchisq(1,omega_p,theta_p)
                    
                    A <- 1/sqrt(beta_j_inferred[i]*abs(diag(Sigma[,,z]))) # changed to diag
                    A_mu_k <- A%*%t(mu[,z])/(d*d)
                    delta_xsq <- 1/sum(A_mu_k)
                    delta_psq <-  abs(1/ ( (1/delta_xsq) + 1/(delta_sq)) )
                    
                    
                    A_x_j <- A%*%x[i,]/(d*d)
                    v_x <- delta_xsq * A_x_j;
                    v_p <- abs(delta_psq*(v_x/delta_xsq + v/(delta_sq))) ## bad idea
                    
                    alpha_j_inferred[i] <- rnorm(1,v_p,delta_psq) ##
                   
                    delta_psq_j[i] <- delta_psq;
                    v_p_j[i] <- v_p
                    
                    ## setting the sampled alphas within constraints
                    if(alpha_j_inferred[i] < min_alpha){
                       alpha_j_inferred[i] <- min_alpha + rnorm(1,0,0.1)
                    }
                    
                    if(alpha_j_inferred[i] > max_alpha){
                        alpha_j_inferred[i] <- max_alpha - rnorm(1,0,0.1)
                    }
                    
                    
                }
                
                # end alpha and beta
                
                qi <- rep(0,K)
                q0 <- 0
                names(qi)<- cnames
                inferred_parm <- inferred_parm2[i,]
                inferred_parm5 <- inferred_parm[-(1:(d+2))]
                dim(inferred_parm5) <- c(d,d)
                tdeg <- inferred_parm[[1]]-d +1
                tmu <- inferred_parm[3:(3+d-1)]
                tsigma <- solve(inferred_parm5) * (inferred_parm[[2]]+1) / (inferred_parm[[2]]*tdeg)
                diag(tsigma) <- diag(tsigma)+ 0.001;
                
                #make symmetric psd for d >=10
                tsigma[lower.tri(tsigma)] <- 0;
                diag_tsigma <- diag(tsigma);
                te <- tsigma + t(tsigma);
                diag(te) <- diag(te) - diag_tsigma;
              
              
                if(step <= alpha_beta_kickin){
                    q0<- alpha * mnormt::dmt(x[i,],tmu,te,tdeg) #revert to this
                    
                }else{
                    q0<- alpha * mnormt::dmt(x[i,],alpha_j_inferred[i]*tmu,beta_j_inferred[i]*te,tdeg) #revert to this
                    
                }
                
                
                
                timedmnormstart <- Sys.time()
                for (k in 1:K){
                    
                    text6 <- (paste('cluster is ', k));
                    write(text6,file=paste0(getwd(),"/",output_folder_name,"/log.txt"),append=TRUE,sep="");
                    print(paste("k is: ", k));
                    if( cnames[k] %in% cnamespreviousstep){
                        qi[cnames[k]] <- Nk[cnames[k]]*Q[i,cnames[k]] ; #revert to this
                        
                        
                    }else{
                        ##added
                        diag(Sigma[,,cnames[k]]) <- diag(Sigma[,,cnames[k]]) + 0.001;
                        temp_Sigma <- matrix(forceSymmetric(Sigma[,,cnames[k]]),d,d);
                        
                        if(step <= alpha_beta_kickin){
                            qi[cnames[k]] <- Nk[cnames[k]]* (dmnorm(x[i,],mu[,cnames[k]],temp_Sigma)[1])
                            
                            
                           
                        }else{
                            qi[cnames[k]] <- Nk[cnames[k]]* (dmnorm(x[i,],alpha_j_inferred[i] * mu[,cnames[k]],beta_j_inferred[i] * temp_Sigma)[1])
                            
                            
                        }
                        
                    }
                }
                timedmnormend <- Sys.time()
                timedmnorm <- timedmnorm + difftime(timedmnormstart,timedmnormend,tz="",units="mins")
                
                ##to take care of log values
                
                #q0 <- exp(q0);
                #qi <- exp(qi);
                ################
                c <- sum(qi) + q0
                q0 <- q0/c
                qi <- qi /c
                
                #c <- exp(q0) + exp(sum(qi))
                #q0 <- exp(q0) /c
                #qi <- exp(qi) /c
                
                ClassOld <- as.character(Class)
                NkOld <- Nk
                
                Classtemp <- sample(c("new",cnames),1, replace = TRUE,c(q0,qi) )#resample the class indicator parameters and sum(qi)+q0 = 1one
                
                
                if( Classtemp=="new") #0 means i starts a new class
                {
                    if (NkOld[ClassOld[i]]==1)# to take care of a singleton cluster that changes to a new class; in effect it is not a new class.
                    {
                        Kit <- Kit +1
                        Class <- factor(Class, levels= c( levels(Class), paste("c",as.character(Kit),sep="")))
                        Class[i] <- paste("c",as.character(Kit),sep="")
                        Sigma[,,ClassOld[i]] <- riwish(inferred_parm[[1]],inferred_parm5)
                        mu[,ClassOld[i]] <- rmnorm(1,inferred_parm[3:(3+d-1)], Sigma[,,ClassOld[i]]/inferred_parm[[2]])
                        for(s in 1:length(dimnames(mu)[[2]]))
                        {
                            if(dimnames(mu)[[2]][s]==ClassOld[i]) {dimnames(mu)[[2]][s] <- paste("c",as.character(Kit),sep="")}
                        }
                        dimnames(Sigma)[[3]] <- dimnames(mu)[[2]]
                        Nk <- table(Class)
                        cnames <- dimnames(mu)[[2]]
                        cnames <- cnames[!cnames == "dummy"]
                        
                    } else
                    {
                        K<-K+1
                        Kit <- Kit +1
                        Class <- factor(Class, levels= c( levels(Class), paste("c",as.character(Kit),sep="")))
                        Class[i] <- paste("c",as.character(Kit),sep="")
                        Sigmatemp <- riwish(inferred_parm[[1]],inferred_parm5)
                        mutemp <-  rmnorm(1,inferred_parm[3:(3+d-1)], Sigma[,,ClassOld[i]]/inferred_parm[[2]])
                        dimnamesmu <- dimnames(mu)
                        
                        mu <- cbind(mu,mutemp);
                        dimnames(mu) <- list( dimnamesmu[[1]],c( dimnamesmu[[2]],as.character(Class[i])))
                        dimnamessig <- dimnames(Sigma)
                        Sigma <- c(Sigma, Sigmatemp)
                        dim(Sigma)<- c(d,d,K+1)
                        dimnames(Sigma) <- list( dimnamessig[[1]], dimnamessig[[2]], c(dimnamessig[[3]],as.character(Class[i])))
                        Nk <- table(Class)
                        cnames <- dimnames(mu)[[2]]
                        cnames <- cnames[!cnames == "dummy"]
                    }
                } else
                {
                    Class[i] <- Classtemp
                    Nk <- table(Class)
                    if ((NkOld[ClassOld[i]]==1) && (as.character(Class[i])!=ClassOld[i]) )
                    {
                        mu <- mu[ , !(colnames(mu) %in% ClassOld[i])]
                        Sigma <- Sigma[,,!(dimnames(Sigma)[[3]] %in% ClassOld[i])]
                        cnames <- dimnames(mu)[[2]]
                        cnames <- cnames[!cnames == "dummy"]
                        K<- K-1      #do not decrease Kit !
                        Class <- factor(Class) #Removes the corresponding level in the class factor
                        Nk <- table(Class)
                    }
                    
                }
            }#end of the loop on the observations
            
            ClassChar <- as.character(Class)
            mydata$class <- Class
            Nk <- table(Class)
            
            
            #for each class resample the parameters
            median_alpha <-rep(1,K)
            median_beta <-rep(1,K)
            
            for(k in 1:K){
                
                
                if(step >= alpha_beta_kickin){
                    median_alpha[k] <- median(alpha_j_inferred[which(mydata$class==cnames[k])])
                    median_beta[k] <- median(beta_j_inferred[which(mydata$class==cnames[k])])
                }
                
                inferred_parm <-  updateNIW(as.matrix(mydata[mydata$class==cnames[k],1:d]),Nk[cnames[k]],nu,kappa,mu0,Delta)
                
                inferred_parm5 <- inferred_parm[-(1:(d+2))]
                dim(inferred_parm5) <- c(d,d)
                Sigma[,,cnames[k]] <-  median_beta[k] * riwish(inferred_parm[[1]],inferred_parm5)
                mu[,cnames[k]] <- median_alpha[k] * rmnorm(1,inferred_parm[3:(3+d-1)], Sigma[,,cnames[k]]/inferred_parm[[2]])
                
                
            }
            
            
            
            f <- paste0(getwd(),"/",output_folder_name,"/plots/Inferred_labels_per_step_per_batch/inferred_labels_(step)_",step,"_(batch)_",r,".pdf")
            pdf(file=f)
            plot(X_tsne$Y[,1],X_tsne$Y[,2],col = col_palette[1*(as.numeric(factor(Class)))],  main=paste("t-SNE of X (true labels) step", step," batch ",r));
            dev.off();
            
        } #end of sampling loop######
        
        
        endTime <- Sys.time()
        diffTime <-  difftime(startTime,endTime,tz="",units="mins")                 
        
        
        
        
        
        
        mydata <- data.frame(cbind(x,Class))
        
        
        
        #x11();
        f <- paste0(getwd(),"/",output_folder_name,"/plots/Inferred_labels/inferred_labels_(batch)_",r,".pdf")
        pdf(file=f)
        
        plot(X_tsne$Y[,1],X_tsne$Y[,2],col = col_palette[1*(as.numeric(factor(Class)))],   main=paste("t-SNE of X (inferred labels) batch ",r));
        dev.off()
        
        ##plot inferred alphas
        #x11();
        g <- paste0(getwd(),"/",output_folder_name,"/plots/Inferred_alphas_betas/alpha_beta_(batch)_",r,".pdf")
        pdf(file=g)
        par(mfrow=c(4,1));
        
        plot(alpha_j_init[order(as.numeric(factor(Class)))], type="l", main="alpha init spread based on z inferred");
        plot(alpha_j_inferred[order(as.numeric(factor(Class)))], type="l", main="alpha inferred spread based on z inferred");
        plot(beta_j_init[order(as.numeric(factor(Class)))], type="l", main="beta init spread based on z inferred");
        plot(beta_j_inferred[order(as.numeric(factor(Class)))], type="l", main="beta inferred spread based on z inferred");
        dev.off();
        
        
        #collect R summary
        #summaryRprof(filename = paste(c(path,"Rprof.out"),sep=""))
        
    
    return(list(z_inferred <- Class, alpha_j_inferred <- alpha_j_inferred, beta_j_inferred <- beta_j_inferred, mu <- mu, Sigma <- Sigma, sd_per_batch <- stddev.genes.batch))
    
    
}


# The main function called once each loop
main.fun <- function(r)
{
    mcmc <- IMM.MCMC(r);
    return(mcmc)
}




##
##
##


# Start the cluster and register with doSNOW
cl <- makeCluster(num_cores, type = "SOCK",outfile="debug.txt") #opens multiple socket connections
clusterExport(cl, c("main.fun", "IMM.MCMC"))
registerDoSNOW(cl)



###Call MCMC per gene split, in parallel
print("Monitor log.txt and outputs/plots/ folder for outputs")
strt <- Sys.time()

results.all.MCMC <- list()
results.per.MCMC <- list()
#num_gene_sub_batches <- 2

runs.all <- floor(num_gene_batches/num_gene_sub_batches)
print(paste0("floor(num_gene_batches/num_gene_sub_batches): ", runs.all))

print("MCMC begins")
print("Begin parallel processing of gene splits");

for(r_all in 1:runs.all){
    
    print(paste("Beginning of batch ",r_all));

    results.per.MCMC <- foreach(r=((1+(num_gene_sub_batches*(r_all-1))):(num_gene_sub_batches*r_all)),.packages=c("Rtsne","lattice","MASS","bayesm","robustbase","chron","mnormt","MCMCpack", "coda","Matrix", "mvtnorm")) %dopar% {
        values <- list()
        values <- main.fun(r)
        return(list(unlist(values[[1]]),unlist(values[[2]]),unlist(values[[3]]),values[[4]],values[[5]],values[[6]])) ##z, alpha, beta, mu and Sigma
    }
    results.all.MCMC[[r_all]] <- results.per.MCMC
    print(paste("End of batch ",r_all));
}
    
print("End of parallel runs");


stopCluster(cl)
#print(Sys.time()-strt)
MCMC_time <- Sys.time()
print(MCMC_time-strt)

















