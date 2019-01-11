#' Nonparametric empirical Bayes classifier without annotations; training
#'
#' Treats the control and case minor allele frequencies as random tuples from a bivariate prior distribution G and then estimates the optimal Bayesian classifier given G. Nonparametric maximum likelihood is used as a plug-in estimator for G.
#' 
#' @param pi0,pi1 p x 1 vectors of control and case minor allele frequencies, respectively; IMPORTANT: must be relative to the same allele in both cases and controls
#' @param n0,n1 number of controls and number of cases, respectively
#' @param d if a single number, G is estimated on a d x d grid; if a two-component vector (d0,d1), G is estimated on a d0 x d1 grid
#' @param maxit maximum number of EM iterations
#' @param tol error tolerance
#' @param verbose TRUE to print the error attained by each EM iteration
#'
#' @return
#' \item{Pi0}{grid points for estimating the distribution of the control minor allele frequencies}
#' \item{Pi1}{grid points for estimating the distribution of the case minor allele frequencies}
#' \item{D0}{conditional density matrix for controls}
#' \item{D1}{conditional density matrix for cases}
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
#' ## training data
#' n0 <- 100; ## number of controls
#' X0 <- t(replicate(n0,rbinom(p,2,pi0))); ## controls
#' n1 <- 50; ## number of cases
#' X1 <- t(replicate(n1,rbinom(p,2,pi1))); ## cases
#' neb <- neb.train(colMeans(X0)/2,colMeans(X1)/2,n0,n1,d=c(20,25));
#' contour(neb$Pi0,neb$Pi1,neb$g);
#' points(pi0,pi1);
#'
#' @import stats
#' @export

neb.train  <- function(pi0,pi1,n0,n1,d=25,maxit=200,tol=1e-4,verbose=FALSE){
    if(sum(is.na(pi0))>0||sum(is.na(pi1))>0){
        stop("Missing values in training data.");
    }
    if(length(d)==0||length(d)>2){
        stop("d can only length 1 or 2.");
    }
    ## check for MAF
    if(min(pi0)==0||min(pi1)==0||max(pi0)==1||max(pi1)==1){
        stop("At least one SNP has MAF=0.");
    }
    ## grids
    if(length(d)==1){
        d0 <- d; d1 <- d;
    } else {
        d0 <- d[1]; d1 <- d[2];
    }
    Pi0 <- seq(min(pi0),max(pi0),length=d0);
    Pi1 <- seq(min(pi1),max(pi1),length=d1);
    ## density matrices
    D0 <- outer(pi0*n0*2,Pi0,function(x,y){ dbinom(x,2*n0,y); });
    D1 <- outer(pi1*n1*2,Pi1,function(x,y){ dbinom(x,2*n1,y); });
    g <- bi.npmle(D0,D1,maxit,tol,verbose);
    return(list(Pi0=Pi0,Pi1=Pi1,D0=D0,D1=D1,g=g,P=n1/(n1+n0)));
}

#' Nonparametric empirical Bayes classifier without annotations; prediction
#' 
#' @param newX n x p matrix of additively coded genotypes to be predicted; IMPORTANT: must be coded relative to the same allele as in the cases and controls
#' @param neb output of neb.train()
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
#' ## training data
#' n0 <- 100; ## number of controls
#' X0 <- t(replicate(n0,rbinom(p,2,pi0))); ## controls
#' n1 <- 50; ## number of cases
#' X1 <- t(replicate(n1,rbinom(p,2,pi1))); ## cases
#' neb <- neb.train(colMeans(X0)/2,colMeans(X1)/2,n0,n1,d=c(20,25));
#' ## testing data
#' newX <- rbind(t(replicate(n0,rbinom(p,2,pi0))),
#'               t(replicate(n1,rbinom(p,2,pi1))));
#' newY <- c(rep(0,n0),rep(1,n1));
#' Yhat <- neb.predict(newX,neb);
#' mean(abs(newY-Yhat$class));
#'
#' @import stats
#' @import iterators
#' @import parallel
#' @export

neb.predict <- function(newX,neb,P=NULL,cores=1){
    if(sum(complete.cases(newX))!=nrow(newX)){
        stop("New data contains missing values.");
    }
    ## pre-calculate some quantities
    p <- ncol(newX);
    d0 <- nrow(neb$g); d1 <- ncol(neb$g);
    tmp0 <- neb$D0*tcrossprod(neb$D1,neb$g); ## p x d0
    tmp1 <- neb$D1*(neb$D0%*%neb$g); ## p x d1
    Pi0 <- neb$Pi0; Pi1 <- neb$Pi1;
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
        ll <- parApply(cl,newX,MARGIN=1,function(x){
            ## case Y=0
            D0 <- outer(x,neb$Pi0,function(x,y){ dbinom(x,2,y); });
            ll0 <- sum(log(.rowSums(D0*tmp0,p,d0)));
            ## case Y=1
            D1 <- outer(x,neb$Pi1,function(x,y){ dbinom(x,2,y); });
            ll1 <- sum(log(.rowSums(D1*tmp1,p,d1)));
            return(c(ll0,ll1));
        }); ## 2 x nrow(newX)
        stopCluster(cl);
    } else {
        ll <- apply(newX,1,function(x){
            ## case Y=0
            D0 <- outer(x,neb$Pi0,function(x,y){ dbinom(x,2,y); });
            ll0 <- sum(log(.rowSums(D0*tmp0,p,d0)));
            ## case Y=1
            D1 <- outer(x,neb$Pi1,function(x,y){ dbinom(x,2,y); });
            ll1 <- sum(log(.rowSums(D1*tmp1,p,d1)));
            return(c(ll0,ll1));
        }); ## 2 x nrow(newX)
    }
    ## classify
    if(is.null(P)){
        P <- neb$P;
    }
    score <- colSums(ll*c(-1,1));
    class <- as.numeric(score>=-logit(P));
    return(list(ll=ll,score=score,class=class));
}
