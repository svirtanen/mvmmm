#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List hybridVBGibbs(IntegerVector zinit,IntegerVector w,IntegerVector d,NumericMatrix b,NumericVector alpha,IntegerVector Nl,int sweeps,int burnin){
  int ii,wi,di,t,kk,iter;
  double probsum,r,maxprobs;
  int Ns=d.size();
  int K=b.nrow();
  int D=b.ncol();
  int M=Nl.size();
  NumericVector alp(K);
  NumericVector probs(K);
  IntegerMatrix N(K,M);
  IntegerMatrix avgN(K,M);
  IntegerMatrix Gam(K,D);
  IntegerVector z(Ns);

  double sumalpha,nom,div;
  int rep;
  int reps=10;

  // collect counts in N
  for(ii=0;ii<Ns;ii++){
    z[ii]=zinit[ii];
    t=z[ii]; wi=w[ii]; di=d[ii];
    N(t,di)+=1;
  }

  for(iter=0;iter<sweeps;iter++){
    for(ii=0;ii<Ns;ii++){
      t=z[ii];/* current topic */
      wi=w[ii];/* current word */
      di=d[ii]; /* current document */
      N(t,di)-=1; // remove current count
      for(kk=0;kk<K;kk++){
        probs[kk]=b(kk,wi)+log(N(kk,di)+alpha[kk]);
        if(kk==0){
          maxprobs=probs[kk];
        }
        if(kk>0){
          if(probs[kk]>maxprobs){
            maxprobs=probs[kk];
          }
        }
      }
      probsum=0;
      for(kk=0;kk<K;kk++){
        probs[kk]=exp(probs[kk]-maxprobs);
        probsum+=probs[kk];
        probs[kk]=probsum;
      }
      r = R::runif(0,1);
      for(kk=0;kk<K;kk++){
        if(r < (probs[kk]/probsum) ){
          t = kk;
          break;
        }
      }
      z[ii] = t;
      N(t,di)+=1;
      // discard burnin and average topic proportions
      if(iter>=burnin){
        avgN(t,di)+=1; Gam(t,wi)+=1;
      }
    }
  }
  // update using Minka's iterations
  sumalpha=0;
  for(kk=0;kk<K;kk++){
    alp[kk] = alpha[kk];
    sumalpha+=alp[kk];
  }
  for(rep=0;rep<reps;rep++){
    for(kk=0;kk<K;kk++){
      nom=-M*R::digamma(alp[kk]);
      div=-M*R::digamma(sumalpha);
      for(di=0;di<M;di++){
        nom+=R::digamma(N(kk,di)+alp[kk]);
        /*div+=digamma(Nl[di]+sumalpha); approx below is accurate for arg>1 */
        div+=log(Nl[di]+sumalpha-0.5);
      }
      sumalpha-=alp[kk];
      alp[kk]*=nom/div;
      sumalpha+=alp[kk];
    }
  }
  return List::create(
    _["z"]=z,
    _["N"]=N,
    _["avgN"]=avgN,
    _["Gam"]=Gam,
    _["alpha"]=alp
    );
}
