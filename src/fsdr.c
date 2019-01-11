#include <R.h>
#include <Rinternals.h>

// for any test statistic with stochastically larger alternative
// T1 = test stat in study 1
// T2 = test stat in study 2
// entries of T1 are sorted from largest to smallest
// T2ord = order(T2,decreasing=T)
// maxprod = maximum possible product, used to initialize minprod
// m1 = restrict threshold for study 1 to be larger than the m1^th largest T1
// m2 = same as m1, but for the T2

SEXP fsdr(SEXP R_T1,SEXP R_T2,SEXP R_T2ord,SEXP R_maxprod,SEXP R_m1,SEXP R_m2,SEXP R_alpha)
{
  // import arguments
  double *T1,*T2,*t;
  int *T2ord;
  T1 = REAL(R_T1); T2 = REAL(R_T2);
  T2ord = INTEGER(R_T2ord);
  int m = length(R_T1);
  double maxprod = *REAL(R_maxprod);
  int m1 = *INTEGER(R_m1);
  int m2 = *INTEGER(R_m2);
  double alpha = *REAL(R_alpha);
  
  // returned threshold pair
  SEXP R_t = PROTECT(allocVector(REALSXP,2));
  t = REAL(R_t);
  t[0] = -99; t[1] = -99; // initialize
  
  // R = number of rejections for curreint pair
  // maxR = max number of rejections observed so far
  // fsdr = fsdr incurred by current pair
  // prod = product of current pair
  // minprod = min product of thresholds observed so far
  int i,j,R,maxR;
  double fsdr,prod,minprod;
  minprod = maxprod; maxR = 0;
  for(i=0;i<m1;i++) // loop through each observed T1 (largest to smallest)
  {
    R = 0;
    for(j=0;j<m2;j++) // loop through each observed T2 (largest to smallest)
    {
      if(T1[T2ord[j]]>=T1[i]){ R++; }
      fsdr = (double)(i+1)*(j+1)/m;
      if(R>0){ fsdr /= R; }
      prod = T1[i]*T2[T2ord[j]];
      
      if(fsdr<=alpha){ // the current pair is feasible
	if(R>maxR){ // update max R, reset min prod, update threshold
	  //printf("T1:%f,T2:%f,R:%d\n",T1[i],T2[T2ord[j]],R);
	  maxR = R;
	  minprod = prod;
	  t[0] = T1[i]; t[1] = T2[T2ord[j]];
	}
	if(R==maxR){ // only update threshold if prod < minprod
	  if(prod<minprod){
	    //printf("T1:%f,T2:%f,R:%d\n",T1[i],T2[T2ord[j]],R);
	    minprod = prod;
	    t[0] = T1[i]; t[1] = T2[T2ord[j]];
	  }
	}
	// do nothing if the pair is feasible but R<maxR
      }
      
    }
  }
  
  UNPROTECT(1);
  return(R_t);
}

// for p-values
// T1 and T2 are sorted from smallest to largest
SEXP fsdr_p(SEXP R_T1,SEXP R_T2,SEXP R_T2ord,SEXP R_m1,SEXP R_m2,SEXP R_alpha)
{
  // import arguments
  double *T1,*T2,*t;
  int *T2ord;
  T1 = REAL(R_T1); T2 = REAL(R_T2);
  T2ord = INTEGER(R_T2ord);
  int m = length(R_T1);
  int m1 = *INTEGER(R_m1);
  int m2 = *INTEGER(R_m2);
  double alpha = *REAL(R_alpha);
  //double attained; // attained FSDR
  
  // returned threshold pair
  SEXP R_t = PROTECT(allocVector(REALSXP,2));
  t = REAL(R_t);
  t[0] = 0; t[1] = 0; // initialize
  
  // R = number of rejections for curreint pair
  // maxR = max number of rejections observed so far
  // fsdr = fsdr incurred by current pair
  // prod = product of current pair
  // maxprod = max product of thresholds observed so far
  int i,j,R,maxR;
  double fsdr,prod,maxprod;
  maxprod = 0; maxR = 0;
  for(i=0;i<m1;i++) // loop through each observed T1
  {
    R = 0;
    for(j=0;j<m2;j++) // loop through each observed T2
    {
      if(T1[T2ord[j]]<=T1[i]){ R++; }
      fsdr = (double)(i+1)*(j+1)/m;
      if(R>0){ fsdr /= R; }
      prod = T1[i]*T2[T2ord[j]];
      
      if(fsdr<=alpha){ // the current pair is feasible
	if(R>maxR){ // update max R, reset min prod, update threshold
	  //printf("T1:%f,T2:%f,R:%d\n",T1[i],T2[T2ord[j]],R);
	  maxR = R;
	  maxprod = prod;
	  t[0] = T1[i]; t[1] = T2[T2ord[j]];
	  //attained = fsdr;
	}
	if(R==maxR){ // only update threshold if prod < minprod
	  if(prod>maxprod){
	    //printf("T1:%f,T2:%f,R:%d\n",T1[i],T2[T2ord[j]],R);
	    maxprod = prod;
	    t[0] = T1[i]; t[1] = T2[T2ord[j]];
	    //attained = fsdr;
	  }
	}
	// do nothing if the pair is feasible but R<maxR
      }
      
    }
  }
  
  //printf("fsdr:%f\n",attained);
  UNPROTECT(1);
  return(R_t);
}

// fsdr_all reports all thresholds in \mathcal{T}\cup\mathcal{P}
SEXP fsdr_all(SEXP R_T1,SEXP R_T2ord,SEXP R_m1,SEXP R_m2,SEXP R_alpha)
{
  // import arguments
  double *T1,*t;
  int *T2ord;
  T1 = REAL(R_T1);
  T2ord = INTEGER(R_T2ord);
  int m = length(R_T1);
  int m1 = *INTEGER(R_m1);
  int m2 = *INTEGER(R_m2);
  double alpha = *REAL(R_alpha);
  
  // calculate max number of rejections that have fdr<=alpha
  // R = number of rejections
  int i,j,R,maxR;
  double fdr;
  maxR = 0;
  for(i=0;i<m1;i++) // loop through each observed T1 (largest to smallest)
  {
    R = 0;
    for(j=0;j<m2;j++) // loop through each observed T2 (largest to smallest)
    {
      if(T1[T2ord[j]]>=T1[i]){ R++; }
      fdr = (double)(i+1)*(j+1)/m;
      if(R>0){ fdr /= R; }
      if(fdr<=alpha&&R>maxR){
	maxR = R;
      }
    }
  }
  
  // inefficient code, fix later using stack, then allocate SEXP and copy
  
  // how many thresholds attain maxfdr
  // recalculate thresholds because there might be too many to store
  int nt = 0;
  for(i=0;i<m1;i++){
    R = 0;
    for(j=0;j<m2;j++){
      if(T1[T2ord[j]]>=T1[i]){ R++; }
      fdr = (double)(i+1)*(j+1)/m;
      if(R>0){ fdr /= R; }
      if(R>=maxR&&fdr<=alpha){
	nt++;
	//printf("maxR:%d,fdr:%f,nt:%d\n",maxR,fdr,nt);
      }
    }
  }
  
  // calculate optimal thresholds
  if(nt==0||maxR==0){
    SEXP R_t = PROTECT(allocVector(REALSXP,1));
    t = REAL(R_t);
    t[0] = -99;
    UNPROTECT(1);
    return(R_t);
  } else {
    SEXP R_t = PROTECT(allocVector(REALSXP,2*nt));
    t = REAL(R_t);
    int wt = 0; // which threshold
    for(i=0;i<m1;i++){
      R = 0;
      for(j=0;j<m2;j++){
	if(T1[T2ord[j]]>=T1[i]){ R++; }
	fdr = (double)(i+1)*(j+1)/m;
	if(R>0){ fdr /= R; }
	if(R>=maxR&&fdr<=alpha){
	  t[wt*2] = T1[i]; t[wt*2+1] = T2ord[j];
	  wt++;
	}
      }
    }
    UNPROTECT(1);
    return(R_t);
  }
}
