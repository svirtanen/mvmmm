#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List hybridVBGibbsTopics(IntegerVector zinit,IntegerVector w,IntegerVector d,NumericMatrix lnZ,double gam, int D,int sweeps,int burnin){
 int ii,wi,di,t,kk,iter;
 double probsum,r,maxprobs;
 int Ns=d.size();
 int K=lnZ.nrow();
 int M=lnZ.ncol();
 NumericVector probs(K);
 IntegerMatrix N(K,M);
 IntegerMatrix Gam(K,D);
 IntegerMatrix avgGam(K,D);
 IntegerVector z(Ns);
 IntegerVector Ntot(K);
 double gamD=gam*D;

 double alp,sumalpha,nom,div;
 int rep;
 int reps=10;

// collect counts in N
for(ii=0;ii<Ns;ii++){
  z[ii]=zinit[ii];
  t=z[ii]; wi=w[ii]; di=d[ii];
  Gam(t,wi)+=1; Ntot[t]+=1;
}

for(iter=0;iter<sweeps;iter++){
  for(ii=0;ii<Ns;ii++){
    t=z[ii];/* current topic */
    wi=w[ii];/* current word */
    di=d[ii]; /* current document */
    Gam(t,wi)-=1; Ntot[t]-=1;
    for(kk=0;kk<K;kk++){
      probs[kk]=log(Gam(kk,wi)+gam)-log(Ntot[kk]+gamD)+lnZ(kk,di);
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
    Gam(t,wi)+=1;
    Ntot[t]+=1;
    // Then discard burnin and average topic proportions
    if(iter>=burnin){
      N(t,di)+=1;
      avgGam(t,wi)+=1;
    }
  }
}

// update using Minka's iterations
alp=gam;
sumalpha=D*alp;
for(rep=0;rep<reps;rep++){
  nom=-D*K*R::digamma(alp);
  div=-D*K*R::digamma(sumalpha);
  for(kk=0;kk<K;kk++){
    for(di=0;di<D;di++){
      nom+=R::digamma(Gam(kk,di)+alp);
    }
    /*div+=digamma(Nl[di]+sumalpha); approx below is accurate for arg>1 */
    div+=D*log(Ntot[kk]+sumalpha-0.5);
  }
  alp*=nom/div;
  sumalpha=D*alp;
}

return List::create(
  _["z"]=z,
  _["N"]=N,
  _["Gam"]=Gam,
  _["avgGam"]=avgGam,
  _["gam"]=alp
  );
}
