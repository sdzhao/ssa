#' Polygenic risk score (given only allele frequencies); training
#'
#' Equivalent to maximum likelihood naive Bayes classifier. The discriminant function is
#' \deqn{\sum_j\hat{\beta}_jX_j,}
#' where \eqn{X_j} is the additively coded genotype of SNP j.
#' 
#' @param pi0,pi1 p x 1 vectors of control and case minor allele frequencies, respectively; IMPORTANT: must be relative to the same allele in both cases and controls
#' @param n0,n1 number of controls and number of cases, respectively
#'
#' @return
#' \item{pi0}{minor allele frequencies in controls}
#' \item{pi1}{minor allele frequencies in cases}
#' \item{P}{proportion of cases}
#'
#' @examples
#' p <- 1000; ## number of snps
#' I <- rep(0,p); I[1:10] <- 1; ## which snps are causal
#' set.seed(1); pi0 <- runif(p,0.1,0.5); ## control minor allele frequencies
#' set.seed(1); ors <- runif(sum(I),-1,1); ## odds ratios
#' pi1 <- pi0;
#' pi1[I==1] <- expit(ors+logit(pi0[I==1]));
#' n0 <- 100; ## number of controls
#' X0 <- t(replicate(n0,rbinom(p,2,pi0))); ## controls
#' n1 <- 50; ## number of cases
#' X1 <- t(replicate(n1,rbinom(p,2,pi1))); ## cases
#' prs.train(colMeans(X0)/2,colMeans(X1)/2,n0,n1);
#'
#' @import stats
#' @export

prs.train <- function(pi0,pi1,n0,n1){
    if(sum(is.na(pi0))>0||sum(is.na(pi1))>0){
        stop("Missing values in training data.");
    }
    ## check for MAF
    if(min(pi0)==0||min(pi1)==0||max(pi0)==1||max(pi1)==1){
        stop("At least one SNP has MAF=0.");
    }
    return(list(pi0=pi0,pi1=pi1,P=n1/(n0+n1)));
}

#' Polygenic risk score (given only allele frequencies); training with CV
#'
#' Uses CV to select how many SNPs to include. SNPs are ordered by the magnitude of their estimated (ang possibly weighted) allelic log-odds ratio. The discriminant function is
#' \deqn{\sum_j\hat{\beta}_jI\left(\left\vert\frac{\hat{\beta}_j}{w_j}\right\vert>\lambda\right)X_j,}
#' where \eqn{X_j} is the additively coded genotype of SNP j.
#' 
#' @param X0,X1 n x p vectors of control and case genotypes, additively coded; IMPORTANT: coding must be relative to the same allele in both cases and controls
#' @param K number of folds for CV
#' @param w p x 1 weight vector
#' @param nlambda number of thresholds to tune over
#' @param verbose if TRUE, report current fold of CV
#'
#' @return
#' \item{pi0}{minor allele frequencies in controls of kept SNPs}
#' \item{pi1}{minor allele frequencies in cases of kept SNPs}
#' \item{w}{weight vector for kept SNPs}
#' \item{lambda}{lambda cutoff}
#' \item{P}{proportion of cases}
#'
#' @examples
#' p <- 1000; ## number of snps
#' I <- rep(0,p); I[1:10] <- 1; ## which snps are causal
#' set.seed(1); pi0 <- runif(p,0.1,0.5); ## control minor allele frequencies
#' set.seed(1); ors <- runif(sum(I),-1,1); ## odds ratios
#' pi1 <- pi0;
#' pi1[I==1] <- expit(ors+logit(pi0[I==1]));
#' n0 <- 100; ## number of controls
#' X0 <- t(replicate(n0,rbinom(p,2,pi0))); ## controls
#' n1 <- 50; ## number of cases
#' X1 <- t(replicate(n1,rbinom(p,2,pi1))); ## cases
#' prs.train.cv(X0,X1,K=3,w=1,nlambda=100,verbose=TRUE);
#'
#' @import stats
#' @export

prs.train.cv <- function(X0,X1,K=3,w=1,nlambda=100,verbose=FALSE){
    n0 <- nrow(X0); n1 <- nrow(X1); p <- ncol(X0);
    P <- n1/(n0+n1);
    pi1 <- colMeans(X1,na.rm=TRUE)/2;
    pi0 <- colMeans(X0,na.rm=TRUE)/2;
    b <- logit(pi1)-logit(pi0);
    lambdas <- seq(0,max(abs(b/w)),length=nlambda);
    ## cv
    ind0 <- sample(1:n0,n0,replace=FALSE);
    ind1 <- sample(1:n1,n1,replace=FALSE);
    inc0 <- floor(n0/K);
    inc1 <- floor(n1/K);
    errs.cv <- rep(0,nlambda);
    for(k in 1:K){
        if(verbose){
            cat("Fold ",k," out of ",K,"...\n",sep="");
        }
        pi1.cv <- colMeans(X1[-ind1[((k-1)*inc1+1):(k*inc1)],],
                           na.rm=TRUE)/2;
        pi0.cv <- colMeans(X0[-ind0[((k-1)*inc0+1):(k*inc0)],],
                           na.rm=TRUE)/2;
        b.cv <- logit(pi1.cv)-logit(pi0.cv);
        errs <- sapply(lambdas,function(l){
            ## thresholds at lambda
            ## CV might result in some infinite log-odds ratios
            thresh <- which(abs(b.cv/w)>l&abs(b.cv/w)!=Inf);
            if(length(thresh)==0){
                scores <- rep(0,inc0+inc1);
            } else {
                b.cv <- b.cv[thresh];
                pi1.cv <- pi1.cv[thresh];
                pi0.cv <- pi0.cv[thresh];
                newX.cv <- rbind(matrix(X0[ind0[((k-1)*inc0+1):(k*inc0)],thresh],
                                        ncol=length(thresh)),
                                 matrix(X1[ind1[((k-1)*inc1+1):(k*inc1)],thresh],
                                        ncol=length(thresh)));
                scores <- 2*sum(log((1-pi1.cv)/(1-pi0.cv)))+newX.cv%*%b.cv;
            }
            newY.cv <- c(rep(0,inc0),rep(1,inc1));
            classes <- as.numeric(scores>-logit(inc1/(inc0+inc1)));
            return(mean(abs(newY.cv-classes)));
        });
        errs.cv <- errs.cv+errs;
    }
    opt <- which.min(errs.cv);
    return(list(pi0=pi0,pi1=pi1,w=w,lambda=lambdas[opt],P=P));
}

#' Polygenic risk score; prediction
#'
#' @param newX n x p matrix of additively coded genotypes to be predicted; IMPORTANT: must be coded relative to the same allele as in the cases and controls
#' @param prs output of prs.train()
#' @param P prevalence of cases in the testing set; if NULL, P is taken from the train object
#'
#' @return
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
#' prs <- prs.train(colMeans(X0)/2,colMeans(X1)/2,n0,n1);
#' ## testing data
#' newX <- rbind(t(replicate(n0,rbinom(p,2,pi0))),
#'               t(replicate(n1,rbinom(p,2,pi1))));
#' newY <- c(rep(0,n0),rep(1,n1));
#' Yhat <- prs.predict(newX,prs);
#' mean(abs(newY-Yhat$class));
#' 
#' @import stats
#' @export

prs.predict <- function(newX,prs,P=NULL){
    if(sum(complete.cases(newX))!=nrow(newX)){
        stop("New data contains missing values.");
    }
    if(is.null(P)){
        P <- prs$P;
    }
    b <- logit(prs$pi1)-logit(prs$pi0);
    score <- 2*sum(log((1-prs$pi1)/(1-prs$pi0)))+
        newX%*%b;
    return(list(score=score,class=as.numeric(score>=-logit(P))));
}

#' Polygenic risk score; prediction for classifier trained with CV
#'
#' @param newX n x p matrix of additively coded genotypes to be predicted; IMPORTANT: must be coded relative to the same allele as in the cases and controls
#' @param prs output of prs.train.cv()
#' @param P prevalence of cases in the testing set; if NULL, P is taken from the train object
#'
#' @return
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
#' prs <- prs.train.cv(X0,X1,K=3,w=1,nlambda=100,verbose=TRUE);
#' ## testing data
#' newX <- rbind(t(replicate(n0,rbinom(p,2,pi0))),
#'               t(replicate(n1,rbinom(p,2,pi1))));
#' newY <- c(rep(0,n0),rep(1,n1));
#' Yhat <- prs.predict.cv(newX,prs);
#' mean(abs(newY-Yhat$class));
#' 
#' @import stats
#' @export

prs.predict.cv <- function(newX,prs,P=NULL){
    if(sum(complete.cases(newX))!=nrow(newX)){
        stop("New data contains missing values.");
    }
    if(is.null(P)){
        P <- prs$P;
    }
    b <- logit(prs$pi1)-logit(prs$pi0);
    keep <- which(abs(b/prs$w)>prs$lambda);
    if(length(keep)==0){
        score <- rep(0,nrow(newX));
    } else {
        newX <- matrix(newX[,keep],ncol=length(keep));
        score <- 2*sum(log((1-prs$pi1[keep])/(1-prs$pi0[keep])))+
            newX%*%b[keep];
    }
    return(list(score=score,class=as.numeric(score>=-logit(P))));
}
