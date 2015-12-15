#
#
#

getDefaults <- function(){
  opts <- vector("list")
  opts$batches <- 200 # number of iterations for training
  opts$batch <- 10 # Note, this is the number of latent variables to jointly optimise in batches for improved performance (MUST be <= M)
  opts$burnin <- 50 # burnin for prediction
  opts$sweeps <- opts$burnin + 100 # set the number of posterior (predictive) samples for prediction

  opts$beta.tau <- 1 # assume Gamma(1,beta.tau) prior for beta
  opts$alpha.tau <- 1 # -the same- for alpha

  return(opts)
}

# Set of auxiliary functions

E.l2mat <- function(x,const,AdivB,U,R,Reg){
  L <- matrix(x,ncol=R)
  E <- sum(const*L) - sum(AdivB*exp( tcrossprod(L,U) )) - Reg*sum(L[,-1]^2)/2
  return(-E)
}

gradE.l2mat <- function(x,const,AdivB,U,R,Reg){
  L <- matrix(x,ncol=R)
  grad <- const - (AdivB*exp(tcrossprod(L,U)))%*%U
  L[,1] <- 0
  grad <- grad - Reg*L
  return(-as.vector(grad))
}


E.V <- function(x,const,AdivB,U,R,Reg){
  L <- matrix(x,ncol=R)
  E <- sum(const*L) - sum(AdivB*exp( tcrossprod(L,U) )) - Reg*sum(L^2)/2
  return(-E)
}

gradE.V <- function(x,const,AdivB,U,R,Reg){
  L <- matrix(x,ncol=R)
  grad <- const - (AdivB*exp(tcrossprod(L,U)))%*%U - Reg*L
  return(-as.vector(grad))
}


E.lu2mat <- function(x,AdivB,AdivB2,L,V,cnstGrad,R){
  U <- matrix(x,ncol=R)
  E <- sum(U*cnstGrad) - sum(AdivB2*exp(tcrossprod(V,U))) - crossprod(x)/2
  for(m in 1:length(L)){
    E <- E - sum(AdivB[[m]]*exp(tcrossprod(L[[m]],U)))
  }
  return(-E)
}

gradE.lu2mat <- function(x,AdivB,AdivB2,L,V,cnstGrad,R){
  U <- matrix(x,ncol=R)
  grad <- cnstGrad - crossprod(AdivB2*exp(tcrossprod(V,U)),V) - U
  for(m in 1:length(L)){
    grad <- grad - crossprod(AdivB[[m]]*exp(tcrossprod(L[[m]],U)),L[[m]])
  }
  return(-as.vector(grad))
}

E.betas2 <- function(beta,par){
  beta <- exp(beta)
  E <- beta*par$const - par$M*lgamma(beta)
  return(-E)
}
gradE.betas2 <- function(beta,par){
  beta <- exp(beta)
  grad <- (par$const - par$M*digamma(beta))*beta
  return(-grad)
}

mvmmm <- function(ses,dat,K,R,D,opts){
  # Inputs
  # ses: app launches
  #  - list of size S. Each element contains a matrix, where the first row contains items of the app vocabulary of size M and the second indicates counts, correspondingly.
  #    The items are pointers (1 to M) for the app vocabulary.
  # dat: co-occurring app text groups
  #  - list of size V (number of views). Each element is an M-length list that follows a similar format to ses (defined above).
  # K is a vector of component numbers for dat and ses, respectively. We assume the V text groups use the same number of components, K[1].
  # R is the dimensionality of the latent variables
  # D is a V-dimensional vector for vocabulary sizes for TEXT GROUPS (dat).

  # Outputs
  # list
  #  - nu, (unnormalised) topics for dat (list of size V). lognu is log nu useful for prediction purposes.
  #  - nu2 (unnormalised) patterns for ses. lognu2 is log nu2
  #  - beta, gamma shape parameter for dat
  #  - beta2, gamma shape parameter for ses
  #  - th, normalised expected text topic proportions (dat)
  #  - V, projection from latent space to ses
  #  - L, projection from latent space to dat (list)
  #  - Lmu, mean of the projection (List)
  #  - U, the latent variables
  #  - lbounds, dat likelihoods (wrt. iterations)
  #  - lbounds, ses likelihoods (wrt. iterations)
  #  - alpha, ses Dirichlet concentration parameter
  #  - gam, dat Dirichlet concentration parameter
  #  - sumEnu, normalisation constant for ses topics (useful for prediction)

  Mods <- length(dat)

  batches <- opts$batches
  beta <- rep(1,Mods)
  beta2 <- 1
  gam <- rep(1,Mods)

  M <- length(dat[[1]]) # number of apps, also dimensionality of topics for ses
  S <- length(ses) # number of users

  if(opts$batch>M){
    stop("The batch size, opts$batch, must be smaller than M.")
  }

  # Initialize parameters for dat
  N <- A <- B <- th <- lnZ <- sumEZ <- Gam <- lognu <- nu <- L <- Lmu <- AdivB <- vector("list",length=Mods)
  for(m in 1:Mods){
  N[[m]] <- matrix(0,K[1],M)
  A[[m]] <- B[[m]] <- matrix(0,K[1],M) + 10
  AdivB[[m]] <- matrix(0,K[1],M)
  th[[m]] <- A[[m]]/B[[m]]
  lnZ[[m]] <- log(th[[m]])
  sumEZ[[m]] <- colSums(th[[m]])

  Gam[[m]] <- matrix(abs(rnorm(D[m]*K[1])),K[1],D[m]) + 10
  lognu[[m]] <- Gam[[m]]
  nu[[m]] <- lognu[[m]]

  L[[m]] <- matrix(rnorm(R*K[1])/1e2,K[1],R)
  Lmu[[m]] <- rnorm(K[1])/1e2
  }

  # Initialize parameters for sessions data (ses)
  alpha <- rep(1,K[2])
  A2 <- matrix( abs(rnorm(K[2]*M)),K[2],M) + 10
  B2 <- matrix( abs(rnorm(K[2]*M)),K[2],M) + 10
  nu2 <- A2/B2
  lognu2 <- log(nu2)
  sumEnu <- rowSums(nu2)
  Gam2 <- matrix(0,K[2],M) + 10

  lbounds <- lbounds.ses <- vector()

  U <- matrix(rnorm(R*M)/1e2,M,R)
  V <- matrix(rnorm(R*K[2])/1e2,K[2],R)
  Xmu <- colSums(U)
  idU <- rep(1,M)

  control.opt <- list(maxit=1,factr=1e7)
  id <- rep(1,K[1]); id2 <- rep(1,K[2])

  Nm <- obs <- obs.ind <- z <- dall <- vector("list",length=Mods)
  for(m in 1:Mods){
    Nm[[m]] <- unlist(lapply(dat[[m]],function(x){sum(x[2,])}))
    obs[[m]] <- as.double(Nm[[m]]!=0); obs.ind[[m]] <- which(obs[[m]]==1)
    dat[[m]] <- as.integer(unlist(lapply(dat[[m]],function(x){sample(rep(x[1,],x[2,]))}))-1)
    z[[m]] <- as.integer(sample.int(K[1],size=sum(Nm[[m]]),replace=TRUE)-1)
    dall[[m]] <- as.integer(rep(1:M,Nm[[m]])-1)
  }

  Nsums <- as.integer(unlist(lapply(ses,function(x){sum(x[2,])})))
  ses <- as.integer(unlist(lapply(ses,function(x){sample(rep(x[1,],x[2,]))}))-1)
  z2 <- as.integer(sample.int(K[2],size=sum(Nsums),replace=TRUE)-1)
  dall2 <- as.integer(rep(1:S,Nsums)-1)

  par.beta <- list(); par.beta$M <- 0; par.beta$const <- 0
  lbound <- 0

  newgam <- gam

  minibatches <- floor(M/opts$batch); bend <- c(c(1:minibatches)*opts$batch,M); bstart <- c(1,bend[-length(bend)]+1)

  for(batch.count in 1:batches){
    # SESSION DATA
    lbound <- 0
    lbound.ses <- 0

    sweeps <- as.integer(2); burnin <- as.integer(1)
    res <- hybridVBGibbs(z2,ses,dall2,lognu2-log(sumEnu),alpha,Nsums,sweeps,burnin)
    z2 <- res$z; newalpha <- res$alpha; Gam2 <- res$Gam/burnin

    # update session topics
    A2 <- Gam2 + beta2
    B2 <- rowSums(Gam2)/sumEnu + exp(tcrossprod(V,U))
    nu2 <- A2/B2
    lognu2 <- log(nu2)
    lbound.ses <- lbound.ses + sum(Gam2*lognu2) - crossprod( log(rowSums(nu2)),rowSums(Gam2))
    sumEnu <- rowSums(nu2)

    # REVIEW DATA
    for(m in 1:Mods){
    sweeps <- as.integer(2); burnin <- as.integer(1)
    res <- hybridVBGibbsTopics(z[[m]],dat[[m]],dall[[m]],lnZ[[m]],gam[m],D[m],sweeps,burnin)

    z[[m]] <- res$z; N[[m]] <- res$N/burnin; newgam[m] <- res$gam
    nu[[m]] <- res$avgGam/burnin+gam[m]; nu[[m]] <- nu[[m]]/rowSums(nu[[m]]); lognu[[m]] <- log(nu[[m]])
    lbound <- lbound + sum(res$Gam*lognu[[m]])

    A[[m]] <- N[[m]] + beta[m]; # note some are empty
    B[[m]] <- outer(id,Nm[[m]]/sumEZ[[m]]) + exp(tcrossprod(L[[m]],U)+Lmu[[m]]);
    th[[m]] <- A[[m]]/B[[m]]
    lnZ[[m]] <- log(th[[m]])
    sumEZ[[m]] <- colSums(th[[m]])
    }

    if(batch.count>1){
      alpha <- newalpha
      gam <- newgam

      for(bs in 1:length(bstart)){
        tmp <- bstart[bs]:bend[bs]
        cnstGrad <- outer(rep(1,length(tmp)),beta2*colSums(V))

        for(m in 1:Mods){
          cnstGrad <- cnstGrad + outer(obs[[m]][tmp],beta[m]*colSums(L[[m]]))
          AdivB[[m]] <- th[[m]][,tmp]*exp(Lmu[[m]])*outer(rep(1,K[1]),obs[[m]][tmp])
        }

        x.tmp <- try(optim(as.vector(U[tmp,]),fn=E.lu2mat,gr=gradE.lu2mat,AdivB=AdivB,
                           AdivB2=nu2[,tmp],L=L,V=V,
                           cnstGrad=cnstGrad,R=R, method="L-BFGS-B", control=control.opt)$par,silent=TRUE)
        if(inherits(x.tmp,"try-error")){ print("Problems in U") }else{ U[tmp,] <- matrix(x.tmp,ncol=R) }
      }

      for(m in 1:Mods){
      Xmu <- colSums(U[obs.ind[[m]],]);
      x.tmp <- try(optim(as.vector(cbind(Lmu[[m]],L[[m]])),fn=E.l2mat,gr=gradE.l2mat,const=outer(rep(beta[m],K[1]),c(sum(obs[[m]]),Xmu)),
                         AdivB=th[[m]][,obs.ind[[m]]],U=cbind(idU,U)[obs.ind[[m]],],R=R+1,Reg=R, method="L-BFGS-B", control=control.opt)$par,silent=TRUE)
      if(inherits(x.tmp,"try-error")){ print("Problems in L") }else{ L[[m]] <- matrix(x.tmp,ncol=R+1); Lmu[[m]] <- L[[m]][,1]; L[[m]] <- L[[m]][,-1] }
      par.beta$const <- sum(L[[m]]%*%Xmu + sum(obs[[m]])*Lmu[[m]]) + sum(lnZ[[m]][,obs.ind[[m]]])
      par.beta$M <- sum(obs[[m]])*K[1]
      beta.tmp <- try(optim(log(beta[m]),fn=E.betas2,gr=gradE.betas2,par.beta,
                           method="L-BFGS-B",control=control.opt)$par,silent=TRUE)
      if(!inherits(beta.tmp,"try-error")){ beta[m] <- exp(beta.tmp) }else{ print("Beta problems") }
      }

      # update V
      Xmu <- colSums(U)
      x.tmp <- try(optim(as.vector(V),fn=E.V,gr=gradE.V,const=outer(rep(beta2,K[2]),Xmu),
                         AdivB=nu2,U=U,R=R,Reg=R,method="L-BFGS-B", control=control.opt)$par,silent=TRUE)
      if(inherits(x.tmp,"try-error")){ print("Problems in V") }else{
        V <- matrix(x.tmp,ncol=R); # need to average...? yes, i guess
      }

      # update beta2
      par.beta$const <- sum(V%*%Xmu) + sum(lognu2);
      par.beta$M <- M*K[2]
      beta.tmp <- try(optim(log(beta2),fn=E.betas2,gr=gradE.betas2,par.beta,
                           method="L-BFGS-B",control=control.opt)$par,silent=TRUE)
      if(!inherits(beta.tmp,"try-error")){ beta2 <- exp(beta.tmp) }else{ print("Beta2 problems") }
    }

    lbounds[batch.count] <- lbound
    lbounds.ses[batch.count] <- lbound.ses

    #print(paste("Iter",batch.count,"LB",lbound,"LB-Ses",lbound.ses))
  }

  for(m in 1:Mods){
    th[[m]] <- apply(th[[m]],2,function(x){x/sum(x)}) # normalised expected text topic proportions
  }
  lbounds[1] <- NA # the first likelihood is not meaningful

  return(list(nu=nu,lognu=lognu,nu2=nu2,lognu2=lognu2,beta=beta,beta2=beta2,th=th,V=V,Lmu=Lmu,L=L,U=U,lbounds=lbounds,
              alpha=alpha,sumEnu=sumEnu,lbounds.ses=lbounds.ses,gam=gam))
}

mvmmmPredict <- function(model,dat,opts){
  # Predictive modeling, now dat here corresponds to ses!
  #   - infer and output topic proportions, $th, for new groups (in dat)
  M <- length(dat)

  nu2 <- model$nu2 # unnormalised topics for ses
  lognu2 <- model$lognu2 # log unnormalised topics
  sumEnu <- model$sumEnu # normalisation for topics
  alpha <- model$alpha # dirichlet topic concentration parameter

  K <- nrow(nu2) # number of topics
  D <- ncol(nu2) # the size of app vocabulary (topics)

  id <- rep(1,K)

  Nsums <- as.integer(unlist(lapply(dat,function(x){sum(x[2,])})))
  wall <- as.integer(unlist(lapply(dat,function(x){sample(rep(x[1,],x[2,]))}))-1)
  z <- as.integer(sample.int(K,size=sum(Nsums),replace=TRUE)-1)
  dall <- as.integer(rep(1:M,Nsums)-1)

  sweeps <- as.integer(opts$sweeps); burnin <- as.integer(opts$burnin)
  res <- hybridVBGibbs(z,wall,dall,lognu2-log(sumEnu),alpha,Nsums,sweeps,burnin)

  N <- res$avgN/(sweeps-burnin) + alpha
  th <- apply(N,2,function(x){x/sum(x)})

  return(list(th=th))
}

predLogLike <- function(dat,lognu,lnZ){
  # compute likelihood for dat under the model with;
  #   - lnZ expected normalised topic proportions
  #   - lognu expected normalised topics
  M <- length(dat)
  K <- as.integer(nrow(lognu))
  dall <- as.integer(rep(1:M,unlist(lapply(dat,ncol)))-1)
  wall <- as.integer(unlist(lapply(dat,function(x){x[1,]}))-1)
  xall <- as.integer(unlist(lapply(dat,function(x){x[2,]})))

  prp <- perp(dall,wall,xall,lognu,lnZ)

  return(prp)
}

#
