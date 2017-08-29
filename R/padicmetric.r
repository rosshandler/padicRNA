####################################
####################################
# padicRNA R package
# Ivan Imaz
# August-2017
####################################
####################################

# Internal functions for ranking expression values
#
# Parameters:
# i: sample iterator
# data: data matrix
# IDs: feature names (eg. genes)
# k: unassigned vector ranks.
# ranked_data: unsorted ranked data.
## Output: matrix of ranks (K) ordered as the input matrix

ranking <- function(i,data,IDs,k){
    data_tmp    <- data.table::data.table(cbind(sample=data[,i],IDs))
    ranked      <- data.table::setorder(data_tmp,-sample)
    ranked_data <- data.frame(cbind(ranked,k))
    rm(data_tmp)
    return(ranked_data)
}

rankingSort <- function(i,ranked_data,IDs){
    data_tmp <- ranked_data[match(IDs,ranked_data[,i-1]), c(i-2,i-1,i)]
    return(data_tmp)
}

# Internal functions for computing v values
#
# Parameters:
# i: sample iterator
# data: data matrix
# K: ranks matrix.
# nfeatures: number of features.
## Output: matrix of v values

vcomp <- function(i,data,K,nfeatures){
    # solve linear-log regression
    x <- log(K[,i], base = exp(1))
    y <- as.numeric(data[,i])
    Sxx <- sum( (x - mean(x))^2 ) / nfeatures
    Syy <- sum( (y - mean(y))^2 ) / nfeatures
    Sxy <- sum( (x - mean(x)) * (y - mean(y)) ) / nfeatures
    b <- Sxy / Sxx
    a <- mean(y) - b * mean(x)
    
    l <- log((a/y), base <- exp(1)) # Inf whether gene expression=0
    r <- l / x                      # Inf whether ranking=1 (log(K))
    D <- exp(r/b)
    inf_index   <- which(x==0)
    mean_signal <- mean(y[-inf_index])
    pvalue      <- D^mean_signal
    
    v <- log(y, base = exp(1)) / log(pvalue, base = exp(1))
    return(v)
}

# Internal functions for evaluating logarithmic approximation
#
# Parameters:
# i: sample iterator
# data: data matrix
# K: ranks matrix.
# nfeatures: number of features.
## Output: R^2 values of fitting for each item

regeval <- function(i,data,K,nfeatures){
    # solve linear-log regression
    x <- log(K[,i], base = exp(1))
    y <- as.numeric(data[,i])
    Sxx <- sum( (x - mean(x))^2 ) / nfeatures
    Syy <- sum( (y - mean(y))^2 ) / nfeatures
    Sxy <- sum( (x - mean(x)) * (y - mean(y)) ) / nfeatures
    b <- Sxy / Sxx
    a <- mean(y) - b * mean(x)
    
    predicted <- a + b*x
    R2 <- 1 - (sum((y-predicted )^2)/sum((y-mean(y))^2))
    return(R2)
}


#' Computes v metric to create a matrix of p-adic transformed data.
#'
#' @param data (numeric) data matrix with features as rows and items as columns.
#'
#' @param reg.eval (logical) Compute R^2 to evaluate fitting of logarithmic approximation. Default is TRUE.
#'
#' @param ncores (Integer) Number of cores. Default is 1.
#'
#' @return matrix of v values (p-adic metric):
#'
#' @examples
#'
#' test ...
#'
#' @author Ivan Imaz \email{ii236@@cam.ac.uk}
#'
#' @export

vmetric <- function(data,reg_eval=TRUE,ncores=1){
    
    k   <- 1:nrow(data)
    IDs <- rownames(data)
    nitems    <- ncol(data)
    nfeatures <- nrow(data)
    
    ### ranking
    if(ncores>1){
        doMC::registerDoMC(cores=ncores)
        #doParallel::registerDoParallel(cores=ncores)
        message("Number of Cores: ",foreach::getDoParWorkers())
        message("Ranking data...")
        `%dopar%` <- foreach::`%dopar%`
        ranked_data <- foreach::foreach(i = 1:nitems,.combine=cbind) %dopar% {
            ranking(i,data,IDs,k)
        }
    }else{
        message("Ranking data...")
        `%do%` <- foreach::`%do%`
        ranked_data <- foreach::foreach(i = 1:nitems,.combine=cbind) %do% {
            ranking(i,data,IDs,k)
        }
    }
    ### re-sort ranking
    ranking_indx <- 3*(1:(ncol(ranked_data)/3))
    sorted_ranked_data <- foreach::foreach(i = ranking_indx,.combine=cbind) %dopar% {
        rankingSort(i,ranked_data,IDs)
    }
    
    # extract ranking matrix
    K  <- matrix(as.numeric(unlist(sorted_ranked_data[, ranking_indx])),
    nrow=nfeatures,ncol=nitems)
    colnames(K) <- paste(rep("K",nitems),1:nitems,sep="_sample")
    rownames(K) <- IDs
    
    message("Done")
    
    #### p-adic metric computation (V)
    message("Computing v metric...")
    
    if(reg_eval){
        R2 <- foreach::foreach(i = 1:nitems,.inorder=FALSE) %dopar% {
            regeval(i,data,K,nfeatures)
        }
        R2 <- round(as.numeric(unlist(R2)),digits=3)
        message(paste("Rsquared range of linear log regression: [",min(R2),",",max(R2),"]",sep=""))
    }
    V <- foreach::foreach(i = 1:nitems,.combine=cbind) %dopar% {
        vcomp(i,data,K,nfeatures)
    }
    colnames(V) <- colnames(data)
    rownames(V) <- IDs
    
    message("Done")
    
    return(-V)
}

#' Applies dimensionality reduction methods to p-adic transformed data from scRNAseq.
#'
#' @param data (numeric) data matrix with features as rows and items as columns.
#'
#' @param log2Transform (logica) Apply log2 transformation to scRNAseq data.
#'
#' @param reg.eval (logical) Compute R^2 to evaluate fitting of logarithmic approximation. Default is TRUE.
#'
#' @param method (character) Method or list of methods for application of dimensionality reduction.
#'
#' @param ndim (numeric) Number of components to extract from dimensionaity reduction methods.
#'
#' @param tsne.iter (numeric) Number of iterations for tsne or Rtsne.
#'
#' @param ncores (Integer) Number of cores. Default is 1.
#'
#' @return ....
#'
#' @examples
#'
#' test ...
#'
#' @author Ivan Imaz \email{ii236@@cam.ac.uk}
#'
#' @export

padicDimred <- function(data,log2Transform=FALSE,reg_eval=TRUE,method="none",ndims=2,tsne.iter=1000,ncores=1){
    
    if(logTransform){
        data <- log2(data+1)
    }
    
    rdim_padic <- list()
    
    V <- vmetric(data=data,reg_eval=reg_eval,ncores=ncores)
    rdim_padic$metric <- V
    logV <- log(abs(V)+1)
    rm(V)
    
    if(length(grep("$tsne",method)) > 0 & length(intersect(method,"none"))> 0){
        rdim_padic$tsne <- tsne::tsne(t(logV),k=ndims,max_iter=tsne_iter)
    }
    
    if(length(grep("$Rtsne",method)) > 0 & length(intersect(method,"none"))> 0){
        rdim_padic$tsne <- Rtsne::Rtsne(t(logV),k=ndims,max_iter=tsne_iter)
    }
    
    if(length(grep("$pca",method)) > 0 & length(intersect(method,"none"))> 0){
        rdim_padic$pca <- prcomp(t(logV))
    }
    
    if(length(grep("$dmap",method)) > 0 & length(intersect(method,"none"))> 0){
        rdim_padic$dmap <- roots::diffuseMat(logV,ndims = ndims)
    }
    return(rdim_padic)
}

