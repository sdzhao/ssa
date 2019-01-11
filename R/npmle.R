#' Univariate NPMLE
#'
#' General nonparametric maximum likelihood estimation for a univariate mixing distribution, implemented using EM. Assumes that the observed data are \eqn{X_i} with marginal likelihood
#' \deqn{\int f(X_i;u)dG(u),}
#' where \eqn{G} is the mixing distribution to be estimated. Suppose there are p observations and \eqn{G} is to be estimated on a grid of d points.
#'
#' @param D p x d matrix of conditional density values, where the \eqn{ij}th entry is \eqn{f(X_i;u_j)}
#' @param maxit maximum number of EM iterations
#' @param tol error tolerance
#' @param verbose TRUE to print the error attained by each EM iteration
#'
#' @return
#' \item{g}{d x 1 vector of probability masses at each grid point}
#'
#' @examples
#' ## generate parameters from mixing distribution
#' p <- 1000;
#' set.seed(1); theta <- rnorm(p);
#' ## generate observed variables
#' X <- rnorm(p,theta,1);
#' ## set grid points
#' d <- 25;
#' Theta <- seq(min(X),max(X),length=d);
#' ## calculate D matrix
#' D <- outer(X,Theta,function(x,y){
#'   dnorm(x,y,1);
#' });
#' ## fit npmle
#' g <- uni.npmle(D);
#' plot(Theta,g/sum(g*(Theta[2]-Theta[1])),type="l");
#' lines(Theta,dnorm(Theta),col="red");
#' 
#' @import stats
#' @export

uni.npmle <- function(D,maxit=200,tol=1e-4,verbose=FALSE){
    p <- nrow(D); d <- ncol(D);
    g <- rep(1/d,d);
    oldll <- 0;
    for(it in 1:maxit){
        tmp <- D%*%g;
        newll <- sum(log(tmp));
        err <- abs(newll-oldll)/abs(oldll);
        if(err<=tol){
            break;
        } else {
            oldll <- newll;
            g <- crossprod(D,1/tmp)*g/p;
            if(verbose) cat("Iter",it,", error",err,"\n");
        }
    }
    if(it==maxit){
        warning("Maximum number of iterations reached.");
    }
    return(g);
}

#' Bivariate NPMLE
#'
#' General nonparametric maximum likelihood estimation for a bivariate mixing distribution, implemented using EM. Assumes that the observed data are tuples \eqn{(X_{1i},X_{2i})} with marginal likelihood
#' \deqn{\int f_1(X_{1i};u_1)f_2(X_{2i};u_2)dG(u_1,u_2),}
#' where \eqn{G} is the mixing distribution to be estimated. Suppose there are p observed tuples and \eqn{G} is to be estimated on a grid of d1 x d2 points.
#'
#' @param D1 p x d1 matrix of conditional density values, where the \eqn{ij}th entry is \eqn{f_1(X_{1i};u_{1j})}.
#' @param D2 p x d2 matrix of conditional density values, where the \eqn{ij}th entry is \eqn{f_2(X_{2i};u_{2j})}.
#' @param maxit maximum number of EM iterations
#' @param tol error tolerance
#' @param verbose TRUE to print the error attained by each EM iteration
#' 
#' @return
#' \item{g}{d1 x d2 matrix of probability masses at each grid point}
#'
#' @examples
#' ## generate parameters from mixing distribution
#' p <- 1000;
#' set.seed(1); theta1 <- rnorm(p); theta2 <- -theta1+rnorm(p);
#' ## generate observed variables
#' X1 <- rnorm(p,theta1,1); X2 <- rnorm(p,theta2,1);
#' ## set grid points
#' d1 <- 25; d2 <- 30;
#' Theta1 <- seq(min(X1),max(X1),length=d1);
#' Theta2 <- seq(min(X2),max(X2),length=d2);
#' ## calculate D matrices
#' D1 <- outer(X1,Theta1,function(x,y){
#'   dnorm(x,y,1);
#' });
#' D2 <- outer(X2,Theta2,function(x,y){
#'   dnorm(x,y,1);
#' });
#' ## fit npmle
#' g <- bi.npmle(D1,D2);
#' contour(Theta1,Theta2,g);
#' points(theta1,theta2);
#'
#' @import stats
#' @export

bi.npmle <- function(D1,D2,maxit=200,tol=1e-4,verbose=FALSE){
    p <- nrow(D1); d1 <- ncol(D1); d2 <- ncol(D2);
    g <- matrix(1/d1/d2,nrow=d1,ncol=d2);
    oldll <- 0;
    for(it in 1:maxit){
        tmp <- .rowSums(D1*tcrossprod(D2,g),p,d1);
        newll <- sum(log(tmp));
        err <- abs(newll-oldll)/abs(oldll);
        if(err<=tol){
            break;
        } else {
            oldll <- newll;
            g <- crossprod(D1,D2/tmp)*g/p;
            if(verbose) cat("Iter",it,", error",err,"\n");
        }
    }
    if(it==maxit){
        warning("Maximum number of iterations reached.");
    }
    return(g);
}

#' Trivariate NPMLE
#'
#' General nonparametric maximum likelihood estimation for a trivariate mixing distribution, implemented using EM. Assumes that the observed data are triples \eqn{(X_{1i},X_{2i},X_{3i})} with marginal likelihood
#' \deqn{\int f_1(X_{1i};u_1)f_2(X_{2i};u_2)f_3(X_{3i};u_3)dG(u_1,u_2,u_3),}
#' where \eqn{G} is the mixing distribution to be estimated. Suppose there are p observed tuples and \eqn{G} is to be estimated on a grid of d1 x d2 x d3 points.
#'
#' @param D1 p x d1 matrix of conditional density values, where the \eqn{ij}th entry is \eqn{f_1(X_{1i};u_{1j})}.
#' @param D2 p x d2 matrix of conditional density values, where the \eqn{ij}th entry is \eqn{f_2(X_{2i};u_{2j})}.
#' @param D3 p x d3 matrix of conditional density values, where the \eqn{ij}th entry is \eqn{f_3(X_{3i};u_{3j})}.
#' @param maxit maximum number of EM iterations
#' @param tol error tolerance
#' @param verbose TRUE to print the error attained by each EM iteration
#' 
#' @return
#' \item{g}{d1 x d2 x d3 array of probability masses at each grid point}
#'
#' @examples
#' ## generate parameters from mixing distribution
#' p <- 1000;
#' set.seed(1);
#' theta1 <- rnorm(p);
#' theta2 <- -theta1+rnorm(p);
#' theta3 <- 0.5*theta1+theta2+rnorm(p);
#' ## generate observed variables
#' X1 <- rnorm(p,theta1,1);
#' X2 <- rnorm(p,theta2,1);
#' X3 <- rnorm(p,theta3,1);
#' ## set grid points
#' d1 <- 15; d2 <- 20; d3 <- 25;
#' Theta1 <- seq(min(X1),max(X1),length=d1);
#' Theta2 <- seq(min(X2),max(X2),length=d2);
#' Theta3 <- seq(min(X3),max(X3),length=d3);
#' ## calculate D matrices
#' D1 <- outer(X1,Theta1,function(x,y){
#'   dnorm(x,y,1);
#' });
#' D2 <- outer(X2,Theta2,function(x,y){
#'   dnorm(x,y,1);
#' });
#' D3 <- outer(X3,Theta3,function(x,y){
#'   dnorm(x,y,1);
#' });
#' ## fit npmle
#' g <- tri.npmle(D1,D2,D3);
#' par(mfrow=c(1,3));
#' contour(Theta1,Theta2,apply(g,c(1,2),sum));
#' points(theta1,theta2);
#' contour(Theta1,Theta3,apply(g,c(1,3),sum));
#' points(theta1,theta3);
#' contour(Theta2,Theta3,apply(g,c(2,3),sum));
#' points(theta2,theta3);
#'
#' @import stats
#' @export

tri.npmle <- function(D1,D2,D3,maxit=200,tol=1e-4,verbose=FALSE){
    ## calculate
    p <- nrow(D1);
    d1 <- ncol(D1); d2 <- ncol(D2); d3 <- ncol(D3);
    g <- array(1/d1/d2/d3,dim=c(d1,d2,d3));
    oldll <- 0;
    for(it in 1:maxit){
        tmp <- matrix(0,nrow=p,ncol=d1);
        for(u3 in 1:d3){
            tmp <- tmp+tcrossprod(D2,g[,,u3])*D3[,u3];
        }
        tmp <- .rowSums(D1*tmp,p,d1);
        newll <- sum(log(tmp));
        err <- abs(newll-oldll)/abs(oldll);
        if(err<=tol){
            break;
        } else {
            oldll <- newll;
            g.tmp <- apply(D3,2,function(col3){
                apply(D2,2,function(col2){
                    .colSums(D1*(col2*col3/tmp),p,d1);
                });
            });
            g <- array(as.vector(g.tmp)/p,dim=c(d1,d2,d3))*g;
            if(verbose) cat("Iter",it,", error",err,"\n");
        }
    }
    if(it==maxit){
        warning("Maximum number of iterations reached.");
    }
    return(g);
}
