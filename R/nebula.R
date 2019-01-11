#' Nonparametric empirical Bayes classifier using latent annotations: chi-square test statistics; training
#'
#' Assumes that chi-square test statistics for each SNP are available from another study. Treats the true control and case minor allele frequencies and the chi-square non-centrality parameters as random triples from a bivariate prior distribution G, and estimates the optimal Bayesian classifier given G. Nonparametric maximum likelihood is used as a plug-in estimator for G.
#' 
#' @param pi0,pi1 p x 1 vectors of control and case minor allele frequencies, respectively; IMPORTANT: must be relative to the same allele in both cases and controls
#' @param n0,n1 number of controls and number of cases, respectively
#' @param T p x 1 vector of chi-square test statistics
#' @param d if a single number, G is estimated on a d x d x d grid; if a three-component vector (d0,d1,dt), G is estimated on a d0 x d1 x dt grid
#' @param maxit maximum number of EM iterations
#' @param tol error tolerance
#' @param verbose TRUE to print the error attained by each EM iteration
#'
#' @return
#' \item{Pi0}{grid points for estimating the distribution of the control minor allele frequencies}
#' \item{Pi1}{grid points for estimating the distribution of the case minor allele frequencies}
#' \item{Lam}{grid points for estimating the distribution of the non-centrality parameter}
#' \item{D0}{conditional density matrix for controls}
#' \item{D1}{conditional density matrix for cases}
#' \item{DT}{conditional density matrix for test statistics}
#' \item{g}{estimated mixing probability mass function}
#' \item{P}{proportion of cases}
#'
#' @examples
#' p <- 1000; ## number of snps
#' I <- rep(0,p); I[1:10] <- 1; ## which snps are causal
#' set.seed(1); pi0 <- runif(p,0.1,0.5); ## control minor allele frequencies
#' set.seed(1); ors <- runif(sum(I),-1,1); ## odds ratios
#' pi1 <- pi0;
#' pi1[I==1] <- expit(ors+logit(pi0[I==1]));
#' set.seed(1); lam <- rep(0,p); lam[I==1] <- rchisq(sum(I==1),1,25); ## ncps
#' ## training data
#' n0 <- 100; ## number of controls
#' X0 <- t(replicate(n0,rbinom(p,2,pi0))); ## controls
#' n1 <- 50; ## number of cases
#' X1 <- t(replicate(n1,rbinom(p,2,pi1))); ## cases
#' T <- rchisq(p,1,lam); ## chi-square statistics
#' nebula <- nebula.chisq.train(colMeans(X0)/2,colMeans(X1)/2,n0,n1,T,d=c(20,25,30));
#' par(mfrow=c(1,3));
#' contour(nebula$Pi0,nebula$Pi1,apply(nebula$g,c(1,2),sum));
#' points(pi0,pi1);
#' contour(nebula$Pi0,nebula$Lam,apply(nebula$g,c(1,3),sum));
#' points(pi0,lam);
#' contour(nebula$Pi1,nebula$Lam,apply(nebula$g,c(2,3),sum));
#' points(pi1,lam);
#'
#' @import stats
#' @export

nebula.chisq.train <- function(pi0,pi1,n0,n1,T,d=25,maxit=200,tol=1e-4,verbose=FALSE){
    if(sum(is.na(pi0))>0||sum(is.na(pi1))>0){
        stop("Missing values in training data.");
    }
    ## check for MAF
    if(min(pi0)==0||min(pi1)==0||max(pi0)==1||max(pi1)==1){
        stop("At least one SNP has MAF=0.");
    }
    if(sum(is.na(T))>0){
        stop("Missing values in T.");
    }
    if(sum(T<=0)>0){
        warning("There is least one non-positive T. Replacing with smallest positive T divided by 10.");
        T[T<=0] <- min(T[T>0])/10; ## otherwise dchisq(0,1,ncp)=Inf
    }
    if(length(d)==0||length(d)>3||length(d)==2){
        stop("d can only have length 1 or 3.");
    }
    ## grids
    if(length(d)==1){
        d0 <- d; d1 <- d; dt <- d;
    } else {
        d0 <- d[1]; d1 <- d[2]; dt <- d[3];
    }
    Pi0 <- seq(min(pi0),max(pi0),length=d0);
    Pi1 <- seq(min(pi1),max(pi1),length=d1);
    Lam <- seq(min(T),max(T),length=dt);
    ## calculate density matrices
    D0 <- outer(pi0*n0*2,Pi0,function(x,y){ dbinom(x,2*n0,y); });
    D1 <- outer(pi1*n1*2,Pi1,function(x,y){ dbinom(x,2*n1,y); });
    DT <- outer(T,Lam,function(x,y){ dchisq(x,df=1,ncp=y); });
    g <- tri.npmle(D0,D1,DT,maxit,tol,verbose);
    return(list(D0=D0,D1=D1,DT=DT,Pi0=Pi0,Pi1=Pi1,Lam=Lam,g=g,P=n1/(n1+n0)));
}

#' Nonparametric empirical Bayes classifier using latent annotations: chi-square test statistics; prediction
#'
#' @param newX n x p matrix of additively coded genotypes to be predicted; IMPORTANT: must be coded relative to the same allele as in the cases and controls
#' @param nebula output of nebula.chisq.train()
#' @param P prevalence of cases in the testing set; if NULL, P is taken from the train object
#' @param cores number of cores to use
#' @return
#' \item{ll}{2 x p matrix of log-likelihoods, first row is from controls}
#' \item{score}{risk score}
#' \item{class}{predicted class, 0=control, 1=case}
#'
#' @examples
#' p <- 1000; ## number of snps
#' I <- rep(0,p); I[1:10] <- 1; ## which snps are causal
#' set.seed(1); pi0 <- runif(p,0.1,0.5); ## control minor allele frequencies
#' set.seed(1); ors <- runif(sum(I),-1,1); ## odds ratios
#' pi1 <- pi0;
#' pi1[I==1] <- expit(ors+logit(pi0[I==1]));
#' set.seed(1); lam <- rep(0,p); lam[I==1] <- rchisq(sum(I==1),1,50); ## ncps
#' ## training data
#' n0 <- 100; ## number of controls
#' X0 <- t(replicate(n0,rbinom(p,2,pi0))); ## controls
#' n1 <- 50; ## number of cases
#' X1 <- t(replicate(n1,rbinom(p,2,pi1))); ## cases
#' T <- rchisq(p,1,lam); ## chi-square statistics
#' nebula <- nebula.chisq.train(colMeans(X0)/2,colMeans(X1)/2,n0,n1,T,d=c(10,12,14));
#' ## testing data
#' newX <- rbind(t(replicate(n0,rbinom(p,2,pi0))),
#'               t(replicate(n1,rbinom(p,2,pi1))));
#' newY <- c(rep(0,n0),rep(1,n1));
#' Yhat <- nebula.chisq.predict(newX,nebula);
#' mean(abs(newY-Yhat$class));
#' 
#' @import stats
#' @import iterators
#' @import parallel
#' @export

nebula.chisq.predict <- function(newX,nebula,P=NULL,cores=1){
    if(sum(complete.cases(newX))!=nrow(newX)){
        warning("New data contains missing values.");
    }
    ## pre-calculate some quantities
    p <- nrow(nebula$D0);
    d0 <- ncol(nebula$D0);
    d1 <- ncol(nebula$D1);
    dt <- ncol(nebula$DT);
    tmp.ctrl <- matrix(0,nrow=p,ncol=d0);
    tmp.case <- matrix(0,nrow=p,ncol=d1);
    ## integrate across u3
    for(u3 in 1:dt){
        ## integrate across u2, fix u1
        tmp.ctrl <- tmp.ctrl+tcrossprod(nebula$D1,nebula$g[,,u3])*nebula$DT[,u3];
        ## integrate across u1, fix u2
        tmp.case <- tmp.case+(nebula$D0%*%nebula$g[,,u3])*nebula$DT[,u3];
      }
    tmp0 <- nebula$D0*tmp.ctrl;
    tmp1 <- nebula$D1*tmp.case;
    Pi0 <- nebula$Pi0; Pi1 <- nebula$Pi1;
    ## calculate scores
    if(cores>1){
        ## parallel apply
        maxcores <- detectCores();
        if(maxcores<cores){
            cores <- maxcores;
            warning(paste("Number of cores set to maximum: ",maxcores,".",sep=""));
        }
        cl <- makeCluster(cores);
        clusterExport(cl,list("tmp0","Pi0","tmp1","Pi1",
                              "p","d0","d1"),envir=environment());
        ## parallel apply
        ll <- parApply(cl,newX,MARGIN=1,function(r){
            D0 <- outer(r,Pi0,function(x,y){ dbinom(x,2,y); });
            ll0 <- sum(log(.rowSums(D0*tmp0,p,d0)));
            D1 <- outer(r,Pi1,function(x,y){ dbinom(x,2,y); });
            ll1 <- sum(log(.rowSums(D1*tmp1,p,d1)));
            return(c(ll0,ll1));
        }); ## 2 x nrow(newX)
        stopCluster(cl);
    } else {
        ll <- apply(newX,1,function(r){
            D0 <- outer(r,Pi0,function(x,y){ dbinom(x,2,y); });
            ll0 <- sum(log(.rowSums(D0*tmp0,p,d0)));
            D1 <- outer(r,Pi1,function(x,y){ dbinom(x,2,y); });
            ll1 <- sum(log(.rowSums(D1*tmp1,p,d1)));
            return(c(ll0,ll1));
        }); ## 2 x nrow(newX)
    }
    if(is.null(P)){
        P <- nebula$P;
    }
    score <- colSums(ll*c(-1,1));
    class <- as.numeric(score>=-logit(P));
    return(list(ll=ll,score=score,class=class));
}

#' Nonparametric empirical Bayes classifier using latent annotations: binary indicators; training
#'
#' Assumes that binary indicators for each SNP are available; e.g. indicate whether the SNP is an eQTL. Treats the true control and case minor allele frequencies for SNPs with indicators equal to 0 and 1 as random triples from bivariate prior distributions G0 and G1, and estimates the optimal Bayesian classifier given G0 and G1. Nonparametric maximum likelihood is used as a plug-in estimator for G0 and G1.
#' 
#' @param pi0,pi1 p x 1 vectors of control and case minor allele frequencies, respectively; IMPORTANT: must be relative to the same allele in both cases and controls
#' @param n0,n1 number of controls and number of cases, respectively
#' @param I p x 1 vector of binary indicators
#' @param d if a single number, G0 and G1 are estimated on d x d grids; if a two-component vector (d0,d1), G0 and G1 are estimated on d0 x d1 grids
#' @param maxit maximum number of EM iterations
#' @param tol error tolerance
#' @param verbose TRUE to print the error attained by each EM iteration
#'
#' @return
#' \item{neb0}{output of neb.train using only SNPs with I==0}
#' \item{neb1}{output of neb.train using only SNPs with I==1}
#' \item{I}{binary indicator}
#'
#' @examples
#' p <- 1000; ## number of snps
#' I <- rep(0,p); I[1:10] <- 1; ## which snps are causal
#' set.seed(1); pi0 <- runif(p,0.1,0.5); ## control minor allele frequencies
#' set.seed(1); ors <- runif(sum(I),-1,1); ## odds ratios
#' pi1 <- pi0;
#' pi1[I==1] <- expit(ors+logit(pi0[I==1]));
#' set.seed(1); lam <- rep(0,p); lam[I==1] <- rchisq(sum(I==1),1,25); ## ncps
#' ## training data
#' n0 <- 100; ## number of controls
#' X0 <- t(replicate(n0,rbinom(p,2,pi0))); ## controls
#' n1 <- 50; ## number of cases
#' X1 <- t(replicate(n1,rbinom(p,2,pi1))); ## cases
#' T <- rchisq(p,1,lam); ## chi-square statistics
#' nebula <- nebula.bin.train(colMeans(X0)/2,colMeans(X1)/2,n0,n1,I,d=c(20,25));
#' par(mfrow=c(1,2));
#' contour(nebula$neb0$Pi0,nebula$neb0$Pi1,apply(nebula$neb0$g,c(1,2),sum));
#' points(pi0[I==0],pi1[I==0]);
#' contour(nebula$neb1$Pi0,nebula$neb1$Pi1,apply(nebula$neb1$g,c(1,2),sum));
#' points(pi0[I==1],pi1[I==1]);
#' 
#' @import stats
#' @export

nebula.bin.train <- function(pi0,pi1,n0,n1,I,d=25,maxit=200,tol=1e-4,verbose=FALSE){
    if(verbose){
        cat("I==0\n");
    }
    neb0 <- neb.train(pi0[I==0],pi1[I==0],n0,n1,
                      d=d,maxit=maxit,tol=tol,verbose=verbose);
    if(verbose){
        cat("I==1\n");
    }
    neb1 <- neb.train(pi0[I==1],pi1[I==1],n0,n1,
                      d=d,maxit=maxit,tol=tol,verbose=verbose);
    return(list(neb0=neb0,neb1=neb1,I=I));
}

#' Nonparametric empirical Bayes classifier using latent annotations: binary indicators; prediction
#'
#' @param newX n x p matrix of additively coded genotypes to be predicted; IMPORTANT: must be coded relative to the same allele as in the cases and controls
#' @param nebula output of nebula.chisq.train()
#' @param P prevalence of cases in the testing set; if NULL, P is taken from the train object
#' @param cores number of cores to use
#' 
#' @return
#' \item{ll}{2 x p matrix of log-likelihoods, first row is from controls}
#' \item{score}{risk score}
#' \item{class}{predicted class, 0=control, 1=case}
#'
#' @examples
#' p <- 1000; ## number of snps
#' I <- rep(0,p); I[1:10] <- 1; ## which snps are causal
#' set.seed(1); pi0 <- runif(p,0.1,0.5); ## control minor allele frequencies
#' set.seed(1); ors <- runif(sum(I),-1,1); ## odds ratios
#' pi1 <- pi0;
#' pi1[I==1] <- expit(ors+logit(pi0[I==1]));
#' set.seed(1); lam <- rep(0,p); lam[I==1] <- rchisq(sum(I==1),1,25); ## ncps
#' ## training data
#' n0 <- 100; ## number of controls
#' X0 <- t(replicate(n0,rbinom(p,2,pi0))); ## controls
#' n1 <- 50; ## number of cases
#' X1 <- t(replicate(n1,rbinom(p,2,pi1))); ## cases
#' T <- rchisq(p,1,lam); ## chi-square statistics
#' nebula <- nebula.bin.train(colMeans(X0)/2,colMeans(X1)/2,n0,n1,I,d=c(20,25));
#' ## testing data
#' newX <- rbind(t(replicate(n0,rbinom(p,2,pi0))),
#'               t(replicate(n1,rbinom(p,2,pi1))));
#' newY <- c(rep(0,n0),rep(1,n1));
#' Yhat <- nebula.bin.predict(newX,nebula);
#' mean(abs(newY-Yhat$class));
#' 
#' @import stats
#' @import iterators
#' @import parallel
#' @export

nebula.bin.predict <- function(newX,nebula,P=NULL,cores=1){
    if(sum(complete.cases(newX))!=nrow(newX)){
        stop("New data contains missing values.");
    }
    if(is.null(P)){
        P <- nebula$neb0$P;
    }
    pred0 <- neb.predict(newX[,nebula$I==0],nebula$neb0,P,cores);
    pred1 <- neb.predict(newX[,nebula$I==1],nebula$neb1,P,cores);
    ## classify
    ll <- pred0$ll+pred1$ll;
    score <- colSums(ll*c(-1,1));
    class <- as.numeric(score>=-logit(P));
    return(list(ll=ll,score=score,class=class));
}

#' Nonparametric empirical Bayes classifier using latent annotations: chi-square test statistics and binary indicators; training
#'
#' Assumes that chi-square test statistics for each SNP are available from another study, and binary indicators for each SNP are available as well. Treats the true control and case minor allele frequencies and the chi-square non-centrality parameters as random triples from a bivariate prior distribution G0 for SNPs with indicators 0, and from G1 for SNPs with indicators equal to 1. Estimates the optimal Bayesian classifier given G. Nonparametric maximum likelihood is used as a plug-in estimator for G.
#' 
#' @param pi0,pi1 p x 1 vectors of control and case minor allele frequencies, respectively; IMPORTANT: must be relative to the same allele in both cases and controls
#' @param n0,n1 number of controls and number of cases, respectively
#' @param T p x 1 vector of chi-square test statistics
#' @param I p x 1 vector of binary indicators
#' @param d if a single number, G0 and G1 are estimated on d x d x d grids; if a three-component vector (d0,d1,dt), G0 and G1 are estimated on d0 x d1 x dt grids
#' @param maxit maximum number of EM iterations
#' @param tol error tolerance
#' @param verbose TRUE to print the error attained by each EM iteration
#' 
#' @return
#' \item{nebula0}{output of nebula.chisq.train using only SNPs with I==0}
#' \item{nebula1}{output of nebula.chisq.train using only SNPs with I==1}
#' \item{I}{binary indicator}
#' 
#' @examples
#' p <- 1000; ## number of snps
#' I <- rep(0,p); I[1:10] <- 1; ## which snps are causal
#' set.seed(1); pi0 <- runif(p,0.1,0.5); ## control minor allele frequencies
#' set.seed(1); ors <- runif(sum(I),-1,1); ## odds ratios
#' pi1 <- pi0;
#' pi1[I==1] <- expit(ors+logit(pi0[I==1]));
#' set.seed(1); lam <- rep(0,p); lam[I==1] <- rchisq(sum(I==1),1,25); ## ncps
#' ## training data
#' n0 <- 100; ## number of controls
#' X0 <- t(replicate(n0,rbinom(p,2,pi0))); ## controls
#' n1 <- 50; ## number of cases
#' X1 <- t(replicate(n1,rbinom(p,2,pi1))); ## cases
#' T <- rchisq(p,1,lam); ## chi-square statistics
#' nebula <- nebula.chisq.bin.train(colMeans(X0)/2,colMeans(X1)/2,n0,n1,T,I,d=c(10,15,20));
#' par(mfrow=c(2,3));
#' contour(nebula$nebula0$Pi0,nebula$nebula0$Pi1,apply(nebula$nebula0$g,c(1,2),sum));
#' points(pi0[I==0],pi1[I==0]);
#' contour(nebula$nebula0$Pi0,nebula$nebula0$Lam,apply(nebula$nebula0$g,c(1,3),sum));
#' points(pi0[I==0],lam[I==0]);
#' contour(nebula$nebula0$Pi1,nebula$nebula0$Lam,apply(nebula$nebula0$g,c(2,3),sum));
#' points(pi1[I==0],lam[I==0]);
#' contour(nebula$nebula1$Pi0,nebula$nebula1$Pi1,apply(nebula$nebula1$g,c(1,2),sum));
#' points(pi0[I==1],pi1[I==1]);
#' contour(nebula$nebula1$Pi0,nebula$nebula1$Lam,apply(nebula$nebula1$g,c(1,3),sum));
#' points(pi0[I==1],lam[I==1]);
#' contour(nebula$nebula1$Pi1,nebula$nebula1$Lam,apply(nebula$nebula1$g,c(2,3),sum));
#' points(pi1[I==1],lam[I==1]);
#'
#' @import stats
#' @export

nebula.chisq.bin.train <- function(pi0,pi1,n0,n1,T,I,d=25,maxit=200,tol=1e-4,verbose=FALSE){
    if(verbose){
        cat("I==0\n");
    }
    nebula0 <- nebula.chisq.train(pi0[I==0],pi1[I==0],n0,n1,T[I==0],
                                  d=d,maxit=maxit,tol=tol,verbose=verbose);
    if(verbose){
        cat("I==1\n");
    }
    nebula1 <- nebula.chisq.train(pi0[I==1],pi1[I==1],n0,n1,T[I==1],
                                  d=d,maxit=maxit,tol=tol,verbose=verbose);
    return(list(nebula0=nebula0,nebula1=nebula1,I=I));
}

#' Nonparametric empirical Bayes classifier using latent annotations: chi-square test statistics and binary indicators; prediction
#'
#' @param newX n x p matrix of additively coded genotypes to be predicted; IMPORTANT: must be coded relative to the same allele as in the cases and controls
#' @param nebula output of nebula.chisq.train()
#' @param P prevalence of cases in the testing set; if NULL, P is taken from the train object
#' @param cores number of cores to use
#' @return
#' \item{ll}{2 x p matrix of log-likelihoods, first row is from controls}
#' \item{score}{risk score}
#' \item{class}{predicted class, 0=control, 1=case}
#' 
#' @examples
#' p <- 1000; ## number of snps
#' I <- rep(0,p); I[1:10] <- 1; ## which snps are causal
#' set.seed(1); pi0 <- runif(p,0.1,0.5); ## control minor allele frequencies
#' set.seed(1); ors <- runif(sum(I),-1,1); ## odds ratios
#' pi1 <- pi0;
#' pi1[I==1] <- expit(ors+logit(pi0[I==1]));
#' set.seed(1); lam <- rep(0,p); lam[I==1] <- rchisq(sum(I==1),1,25); ## ncps
#' ## training data
#' n0 <- 100; ## number of controls
#' X0 <- t(replicate(n0,rbinom(p,2,pi0))); ## controls
#' n1 <- 50; ## number of cases
#' X1 <- t(replicate(n1,rbinom(p,2,pi1))); ## cases
#' T <- rchisq(p,1,lam); ## chi-square statistics
#' nebula <- nebula.chisq.bin.train(colMeans(X0)/2,colMeans(X1)/2,n0,n1,T,I,d=c(10,12,14));
#' ## testing data
#' newX <- rbind(t(replicate(n0,rbinom(p,2,pi0))),
#'               t(replicate(n1,rbinom(p,2,pi1))));
#' newY <- c(rep(0,n0),rep(1,n1));
#' Yhat <- nebula.chisq.bin.predict(newX,nebula);
#' mean(abs(newY-Yhat$class));
#' 
#' @return
#' \item{ll}{2 x p matrix of log-likelihoods, first row is from controls}
#' \item{score}{risk score}
#' \item{class}{predicted class, 0=control, 1=case}
#'
#' @import stats
#' @import iterators
#' @import parallel
#' @export

nebula.chisq.bin.predict <- function(newX,nebula,P=NULL,cores=1){
    if(sum(complete.cases(newX))!=nrow(newX)){
        stop("New data contains missing values.");
    }
    if(is.null(P)){
        P <- nebula$nebula0$P;
    }
    pred0 <- nebula.chisq.predict(newX[,nebula$I==0],nebula$nebula0,P,cores);
    pred1 <- nebula.chisq.predict(newX[,nebula$I==1],nebula$nebula1,P,cores);
    ## classify
    ll <- pred0$ll+pred1$ll;
    score <- colSums(ll*c(-1,1));
    class <- as.numeric(score>=-logit(P));
    return(list(ll=ll,score=score,class=class));
}

#' Nonparametric empirical Bayes classifier using latent annotations: wrapper function; training
#'
#' @param pi0,pi1 p x 1 vectors of control and case minor allele frequencies, respectively; IMPORTANT: must be relative to the same allele in both cases and controls
#' @param n0,n1 number of controls and number of cases, respectively
#' @param T p x 1 vector of chi-square test statistics
#' @param I p x 1 vector of binary indicators
#' @param d if a single number, G0 and G1 are estimated on d x d x d grids; if a three-component vector (d0,d1,dt), G0 and G1 are estimated on d0 x d1 x dt grids
#' @param maxit maximum number of EM iterations
#' @param tol error tolerance
#' @param verbose TRUE to print the error attained by each EM iteration
#'
#' @return
#' \item{type}{1=given neither T nor I; 2=given T but not I; 3=not given T but given I; 4=given both T and I}
#' \item{nebula}{trained classifier}
#'
#' @examples
#' p <- 1000; ## number of snps
#' I <- rep(0,p); I[1:10] <- 1; ## which snps are causal
#' set.seed(1); pi0 <- runif(p,0.1,0.5); ## control minor allele frequencies
#' set.seed(1); ors <- runif(sum(I),-1,1); ## odds ratios
#' pi1 <- pi0;
#' pi1[I==1] <- expit(ors+logit(pi0[I==1]));
#' set.seed(1); lam <- rep(0,p); lam[I==1] <- rchisq(sum(I==1),1,50); ## ncps
#' ## training data
#' n0 <- 100; ## number of controls
#' X0 <- t(replicate(n0,rbinom(p,2,pi0))); ## controls
#' n1 <- 50; ## number of cases
#' X1 <- t(replicate(n1,rbinom(p,2,pi1))); ## cases
#' T <- rchisq(p,1,lam); ## chi-square statistics
#' nebula1 <- nebula.train(colMeans(X0)/2,colMeans(X1)/2,n0,n1,d=c(10,15));
#' nebula2 <- nebula.train(colMeans(X0)/2,colMeans(X1)/2,n0,n1,T=T,d=c(10,15,20));
#' nebula3 <- nebula.train(colMeans(X0)/2,colMeans(X1)/2,n0,n1,I=I,d=c(10,15));
#' nebula4 <- nebula.train(colMeans(X0)/2,colMeans(X1)/2,n0,n1,T=T,I=I,d=c(10,15,20));
#'
#' @import stats
#' @export

nebula.train <- function(pi0,pi1,n0,n1,T=NULL,I=NULL,d=25,maxit=200,tol=1e-4,verbose=FALSE){
    if(is.null(T[1])&&is.null(I[1])){
        return(list(type=1,
                    nebula=neb.train(pi0,pi1,n0,n1,d,maxit,tol,verbose)));
    } else if(!is.null(T[1])&&is.null(I[1])){
        return(list(type=2,
                    nebula=nebula.chisq.train(pi0,pi1,n0,n1,T,d,maxit,tol,verbose)));
    } else if(is.null(T[1])&&!is.null(I[1])){
        return(list(type=3,
                    nebula=nebula.bin.train(pi0,pi1,n0,n1,I,d,maxit,tol,verbose)));
    } else if(!is.null(T[1])&&!is.null(I[1])){
        return(list(type=4,
                    nebula=nebula.chisq.bin.train(pi0,pi1,n0,n1,T,I,d,maxit,tol,verbose)));
    }
}

#' Nonparametric empirical Bayes classifier using latent annotations: wrapper function; predict
#'
#' @param newX n x p matrix of additively coded genotypes to be predicted; IMPORTANT: must be coded relative to the same allele as in the cases and controls
#' @param nebula output of nebula.chisq.train()
#' @param P prevalence of cases in the testing set; if NULL, P is taken from the train object
#' @param cores number of cores to use
#'
#' @return
#' \item{ll}{2 x p matrix of log-likelihoods, first row is from controls}
#' \item{score}{risk score}
#' \item{class}{predicted class, 0=control, 1=case}
#' 
#' @examples
#' p <- 1000; ## number of snps
#' I <- rep(0,p); I[1:10] <- 1; ## which snps are causal
#' set.seed(1); pi0 <- runif(p,0.1,0.5); ## control minor allele frequencies
#' set.seed(1); ors <- runif(sum(I),-1,1); ## odds ratios
#' pi1 <- pi0;
#' pi1[I==1] <- expit(ors+logit(pi0[I==1]));
#' set.seed(1); lam <- rep(0,p); lam[I==1] <- rchisq(sum(I==1),1,25); ## ncps
#' ## training data
#' n0 <- 100; ## number of controls
#' X0 <- t(replicate(n0,rbinom(p,2,pi0))); ## controls
#' n1 <- 50; ## number of cases
#' X1 <- t(replicate(n1,rbinom(p,2,pi1))); ## cases
#' T <- rchisq(p,1,lam); ## chi-square statistics
#' nebula1 <- nebula.train(colMeans(X0)/2,colMeans(X1)/2,n0,n1,d=c(5,7));
#' nebula2 <- nebula.train(colMeans(X0)/2,colMeans(X1)/2,n0,n1,T=T,d=c(5,7,9));
#' nebula3 <- nebula.train(colMeans(X0)/2,colMeans(X1)/2,n0,n1,I=I,d=c(5,7));
#' nebula4 <- nebula.train(colMeans(X0)/2,colMeans(X1)/2,n0,n1,T=T,I=I,d=c(5,7,9));
#' ## testing data
#' newX <- rbind(t(replicate(n0,rbinom(p,2,pi0))),
#'               t(replicate(n1,rbinom(p,2,pi1))));
#' newY <- c(rep(0,n0),rep(1,n1));
#' Yhat1 <- nebula.predict(newX,nebula1);
#' Yhat2 <- nebula.predict(newX,nebula2);
#' Yhat3 <- nebula.predict(newX,nebula3);
#' Yhat4 <- nebula.predict(newX,nebula4);
#' mean(abs(newY-Yhat1$class));
#' mean(abs(newY-Yhat2$class));
#' mean(abs(newY-Yhat3$class));
#' mean(abs(newY-Yhat4$class));
#' 
#' @import stats
#' @import iterators
#' @import parallel
#' @export

nebula.predict <- function(newX,nebula,P=NULL,cores=1){
    switch(nebula$type,
           "1"=return(neb.predict(newX,nebula$nebula,P,cores)),
           "2"=return(nebula.chisq.predict(newX,nebula$nebula,P,cores)),
           "3"=return(nebula.bin.predict(newX,nebula$nebula,P,cores)),
           "4"=return(nebula.chisq.bin.predict(newX,nebula$nebula,P,cores)));
}
