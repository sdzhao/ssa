#include <R.h>
#include <Rinternals.h>
// this code won't work with tied data

SEXP ldd(SEXP R_U,SEXP R_Vord,SEXP R_nU,SEXP R_nV)
// entries of U,V are sorted from smaller U to largest U
// entries of Fv are the eCDF value of the corresponding entries of V
// Vord = order(V)
{
  int i,j,rankv;
  // rankv is the rank of V among the U<=U[i]
  double temp,Fuv,Fu,Fv;
  SEXP R_D = PROTECT(allocVector(REALSXP, 3)); // last two components are optimal thresholds
  
  // pointers for easy reference
  double *U,*D;
  int *Vord;
  U = REAL(R_U);
  D = REAL(R_D);
  Vord = INTEGER(R_Vord);

  // allow user to only take max in the lower tails of U and V
  int n = length(R_U);
  int nU = INTEGER(R_nU)[0];
  int nV = INTEGER(R_nV)[0];
  
  // if i==n-1, Fu=1, rankv equals j+1, and Fuv-Fu*Fv = 0
  // if j==n-1, Fv=1, rankv equals rank of U, and Fuv-Fu*Fv = 0
  // just loop up to n-1 so denominator doesn't blow up
  D[0] = 0; D[1] = 0; D[2] = 0;
  for(i=0;i<(nU-1);i++) // loop through each observed U (smallest to largest)
  {
    rankv = 0;
    for(j=0;j<(nV-1);j++) // loop through each observed V (smallest to largest)
    {
      if(U[Vord[j]]<=U[i]){ rankv++; }
      // printf("%f,%d,%f,%d,%d;\t",U[Vord[j]],Vord[j],U[i],i,rankv);
      Fuv = (double)rankv/n;
      Fu = (double)(i+1)/n;
      Fv = (double)(j+1)/n;
      // will never have Fu*Fv = 1
      temp = fabs(Fuv-Fu*Fv)/sqrt(Fu*Fv-Fu*Fv*Fu*Fv);
      if(temp>D[0]){ D[0] = temp; D[1] = i+1; D[2] = j+1;}
      //printf("%f,%d,%f,%f,%f,%f,%f,%f\n",U[i],j,Fuv,Fu,Fv,sqrt(Fu*Fv-Fu*Fv*Fu*Fv),sqrt(Fu*(1-Fu)*Fv*(1-Fv)),temp);
    }
  }
  
  UNPROTECT(1);
  return(R_D);
}
