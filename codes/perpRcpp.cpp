#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double perp(IntegerVector d,IntegerVector w,IntegerVector x,NumericMatrix lognu,NumericMatrix lnZ){
  int kk,ii,di,wi,xc;
  double probsum;
  double prp=0;
  int K=lnZ.nrow();
  int Ns=d.size();

  for(ii=0;ii<Ns;ii++){
    di=d[ii]; wi=w[ii]; xc=x[ii];
    probsum=0;
    for(kk=0;kk<K;kk++){
      probsum+=lognu(kk,wi)*lnZ(kk,di);
    }
    prp+=xc*log(probsum);
  }
  return prp;
}
