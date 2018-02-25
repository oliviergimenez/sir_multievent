---
title: Benchmarking comparing the implementation of SIR multi-event capture-recapture
  models in TMB vs. native R
output:
  html_document:
    keep_md: yes
---

# Data simulation

We first build a function to simulate data from a SIR multi-event capture-recapture model:

```r
simul <- function(
n.occasions = 5, # number of occasions
n.states = 4, # number of states: S, I and R plus dead
n.obs = 4, # number of events: 1) detected and diagnosed as S, 2) as I, 3) as R, 4) detected and undiagnosed
phiS = 0.96, # survival susceptible
phiI = 0.53, # survival infected
phiR = 0.84, # survival recovered
betaSI = 0.60, # infection prob 
pS = 0.95, # detection susceptible
pI = 0.95, # detection infected
pR = 0.95, # detection recovered
deltaS = 0.1, # ass susceptible
deltaI = 0.1, # ass infected
deltaR = 0.1, # ass recovered
unobservable = NA){
 
# number of individuals marked in each state
marked <- matrix(NA, ncol = n.states, nrow = n.occasions) 
marked[,1] <- rep(12, n.occasions) 
marked[,2] <- rep(24, n.occasions) 
marked[,3] <- rep(12, n.occasions) 
marked[,4] <- rep(0, n.occasions) 
tot <- marked[1,1] + marked[1,2] + marked[1,3] + marked[1,4]
parsim <- c(marked[1,1]/tot, marked[1,2]/tot, phiS, phiI, phiR, betaSI, pS, pI, pR, deltaS, deltaI, deltaR) 

# Create the structure of the SIR model  
# Define matrices with survival, transition and recapture probabilities 
# These are 4-dimensional matrices, with 
# Dimension 1: state of departure 
# Dimension 2: state of arrival 
# Dimension 3: individual 
# Dimension 4: time 
    
# 1. State process matrix SIR transition matrix given survival of individuals
totrel <- sum(marked)*(n.occasions-1) 
PSI.STATE <- array(NA, dim=c(n.states, n.states, totrel, n.occasions-1)) 
    for (i in 1:totrel){ 
      for (t in 1:(n.occasions-1)){ 
        PSI.STATE[,,i,t] <- matrix(c( 
          phiS*(1-betaSI), phiS*(betaSI), 0, 1-phiS, 
          0, 0, phiI, 1-phiI, 
          0, 0, phiR, 1-phiR, 
          0, 0, 0, 1 ), nrow = n.states, byrow = TRUE) 
      } #t 
    } #i 
    
# 2.Observation process matrix 
PSI.OBS <- array(NA, dim=c(n.states, n.obs, totrel, n.occasions-1)) 
    for (i in 1:totrel){ 
      for (t in 1:(n.occasions-1)){ 
        PSI.OBS[,,i,t] <- matrix(c( 
          pS, 0, 0, 1-pS, 
          0, pI, 0, 1-pI, 
          0, 0, pR, 1-pR, 
          0, 0, 0, 1), nrow = n.states, byrow = TRUE) 
      } #t 
    } #i 
        
# Unobservable: number of state that is unobservable 
n.occasions <- dim(PSI.STATE)[4] + 1 
CH <- CH.TRUE <- matrix(NA, ncol = n.occasions, nrow = sum(marked)) 

# Define a vector with the occasion of marking 
mark.occ <- matrix(0, ncol = dim(PSI.STATE)[1], nrow = sum(marked)) 
g <- colSums(marked) 
  
for (s in 1:dim(PSI.STATE)[1]){ 
    if (g[s]==0) next # avoid error message if nothing to replace 
    mark.occ[(cumsum(g[1:s])-g[s]+1)[s]:cumsum(g[1:s])[s],s] <- 
      rep(1:n.occasions, marked[1:n.occasions,s]) # repeat the occasion t a certain time i.e the number of observation in state s 
} #s 
  
for (i in 1:sum(marked)){ 
    for (s in 1:dim(PSI.STATE)[1]){ 
      if (mark.occ[i,s]==0) next 
      first <- mark.occ[i,s] 
      CH[i,first] <- s 
      CH.TRUE[i,first] <- s 
    } #s 
    for (t in (first+1):n.occasions){ 
      # Multinomial trials for state transitions 
      if (first==n.occasions) next 
      state <- which(rmultinom(1, 1, PSI.STATE[CH.TRUE[i,t-1],,i,t-1])==1) 
      CH.TRUE[i,t] <- state 
      # Multinomial trials for observation process 
      event <- which(rmultinom(1, 1, PSI.OBS[CH.TRUE[i,t],,i,t-1])==1) 
      CH[i,t] <- event 
    } #t 
} #i 

# Replace the NA and the highest state number (dead) in the file by 0 
CH[is.na(CH)] <- 0 
CH[CH==dim(PSI.STATE)[1]] <- 0 
CH[CH==unobservable] <- 0 
id <- numeric(0) 
for (i in 1:dim(CH)[1]){ 
    z <- min(which(CH[i,]!=0)) 
    ifelse(z==dim(CH)[2], id <- c(id,i), id <- c(id)) 
} 

# capture histories to be used 
CH <- CH[-id,]
# capture histories with perfect observation 
CH.TRUE <- CH.TRUE[-id,] 

# To artificially generate uncertainty on states 
# we alter the raw capture-recapture data from file CH
    
# nb of capture occasions
ny <- ncol(CH)

# nb of individuals
nind <- nrow(CH)
    
titi2 <- CH
for (i in 1:nind)
{
      # deltaS<- max(0.01, deltaS*(1-deltaScor))
      # deltaR<- max(0.01, deltaR*(1-deltaRcor))
      # 
      # deltaScor<-0
      # deltaRcor<-0
      
      for (j in 1:ny){
        # 1 seen and ascertained A (with probability .2)
        # 2 seen and ascertained B (with probability .7)
        # 3 seen but not ascertained (A with probability .8 + B with probability .3)
        # 0 not seen

        if (CH[i,j] == 1)
        {
          temp <- rbinom(1,size=1,prob=deltaS) 
          if (temp == 1) titi2[i,j] <- 1 # if ascertained NB, event = 1
          if (temp == 0) titi2[i,j] <- 4 # if not ascertained, event = 4
        }
        
        if (CH[i,j] == 2) 
        {
          temp <- rbinom(1,size=1,prob=deltaI)
          if (temp == 1) titi2[i,j] <- 2 # if ascertained B, event = 2
          if (temp == 0) titi2[i,j] <- 4 # if not ascertained, event = 3
          #uncertain<-which(titi2[i,]==4)
        }
        
        if (CH[i,j] == 3) 
        {
          temp <- rbinom(1,size=1,prob=deltaR)
          if (temp == 1) titi2[i,j] <- 3 # if ascertained B, event = 3
          if (temp == 0) titi2[i,j] <- 4 # if not ascertained, event = 4
          #uncertain<-which(titi2[i,]==4)
        }
      }
     }
     
return(titi2)
}
```

# Likelihood function in native `R`

First a short function to protect the log from "exploding"

```r
logprot <- function(v){
  eps <- 2.2204e-016
  u <- log(eps) * (1+vector(length=length(v)))
  index <- (v>eps)
  u[index] <- log(v[index])
  u
}
```

Then the likelihood:

```r
devMULTIEVENT <- function(b,data,eff,e,garb,nh,km1){
  
  # b parameters
  # data encounter histories
  # eff counts
  # e vector of dates of first captures
  # garb vector of initial states
  # nh nb ind
  # km1 nb of recapture occasions (nb of capture occ - 1)
  
  # logit link for all parameters

  par <- plogis(b) # we could use the multinomial logit link instead
  piS <- par[1]
  piI <- par[2]
  
  phiS <- par[3]
  phiI <- par[4]
  phiR <- par[5]
  betaSI <- par[6]
  
  pS <- par[7]
  pI <- par[8]
  pR <- par[9]
  
  
  deltaS <- par[10]
  deltaI <- par[11]
  deltaR <- par[12]
  
  
  # prob of obs (rows) cond on states (col)
  B1 <- matrix(c(
  1-pS,pS,0,0,
  1-pI,0,pI,0,
  1-pR,0,0,pR,
  1,0,0,0),nrow=4,ncol=4,byrow=T)
  B2 = matrix(c(
    1,0,0,0,0,
    0,deltaS,0,0,1-deltaS,
    0,0,deltaI,0,1-deltaI,
    0, 0, 0,deltaR, 1-deltaR),nrow=4,ncol=5,byrow=T)
  B = t(B1 %*% B2)

  # first encounter
  BE1 = matrix(c(
  0,1,0,0,
  0,0,1,0,
  0,0,0,1,
  1,0,0,0),nrow=4,ncol=4,byrow=T)
  BE2 = matrix(c(
    1,0,0,0,0,
    0,deltaS,0,0,1-deltaS,
    0,0,deltaI,0,1-deltaI,
    0, 0, 0,deltaR, 1-deltaR),nrow=4,ncol=5,byrow=T)
  BE = t(BE1 %*% BE2)
  
  # prob of states at t+1 given states at t
  A1 <- matrix(c(
  phiS,0,0,1-phiS,
  0,phiI,0,1-phiI,
  0,0,phiR,1-phiR,
  0,0,0,1),nrow=4,ncol=4,byrow=T)
  A2 <- matrix(c(
  1-betaSI,betaSI,0,0,
  0,0,1,0,
  0,0,1,0,
  0,0,0,1),nrow=4,ncol=4,byrow=T)
  A <- A1 %*% A2

  # init states
  PI <- c(piS,piI,1-piS-piI,0)
  
  # likelihood
  l <- 0
  for (i in 1:nh) # loop on ind
  {
    ei <- e[i] # date of first det
    oe <- garb[i] + 1 # init obs
    evennt <- data[,i] + 1 # add 1 to obs to avoid 0s in indexing
    ALPHA <- PI*BE[oe,]
    for (j in (ei+1):(km1+1)) # cond on first capture
    {
      if ((ei+1)>(km1+1)) {break} 
      ALPHA <- (ALPHA %*% A)*B[evennt[j],]
    }
    l <- l + logprot(sum(ALPHA)) #*eff[i]
  }
  l <- -l
  l
}
```

# Likelihood function in `TMB`

First, create the model template:

```r
tmb_model <- "
// multi-event SIR capture-recapture model
#include <TMB.hpp>

//template<class Type>
//matrix<Type> multmat(array<Type> A, matrix<Type> B) {
//  int nrowa = A.rows();
//  int ncola = A.cols(); 
//  int ncolb = B.cols(); 
//  matrix<Type> C(nrowa,ncolb);
//  for (int i = 0; i < nrowa; i++)
//  {
//    for (int j = 0; j < ncolb; j++)
//    {
//      C(i,j) = Type(0);
//      for (int k = 0; k < ncola; k++)
//        C(i,j) += A(i,k)*B(k,j);
//    }
//  }
//  return C;
//}
//
/* implement the vector - matrix product */
template<class Type>
vector<Type> multvecmat(array<Type>  A, matrix<Type>  B) {
  int nrowb = B.rows();
  int ncolb = B.cols(); 
  vector<Type> C(ncolb);
  for (int i = 0; i < ncolb; i++)
  {
	    C(i) = Type(0);
      for (int k = 0; k < nrowb; k++){
        C(i) += A(k)*B(k,i);
    }
  }
  return C;
}

template<class Type>
Type objective_function<Type>::operator() () {
  
  // b = parameters
  PARAMETER_VECTOR(b);
  
  // ch = capture-recapture histories (individual format)
  // fc = date of first capture
  // fs = state at first capture
  DATA_IMATRIX(ch);
  DATA_IVECTOR(fc);
  DATA_IVECTOR(fs);
  
  // OBSERVATIONS
  // 0 = non-detected
  // 1 = seen and ascertained as non-breeder
  // 2 = seen and ascertained as breeder
  // 3 = not ascertained
  //   
  // STATES
  // 1 = alive non-breeder
  // 2 = alive breeder
  // 3 = dead
  //   
  // PARAMETERS
  // phiNB  survival prob. of non-breeders
  // phiB  survival prob. of breeders
  // pNB  detection prob. of non-breeders
  // pB  detection prob. of breeders
  // psiNBB transition prob. from non-breeder to breeder
  // psiBNB transition prob. from breeder to non-breeder
  // piNB prob. of being in initial state non-breeder
  // deltaNB prob to ascertain the breeding status of an individual encountered as non-breeder
  // deltaB prob to ascertain the breeding status of an individual encountered as breeder
  //   
  // logit link for all parameters
  // note: below, we decompose the state and obs process in two steps composed of binomial events, 
  // which makes the use of the logit link appealing; 
  // if not, a multinomial (aka generalised) logit link should be used

  int km = ch.rows();
  int nh = ch.cols();  
  int npar = b.size();
  vector<Type> par(npar);
  for (int i = 0; i < npar; i++) {
    par(i) = Type(1.0) / (Type(1.0) + exp(-b(i)));
  }
  Type piS = par(0); // careful, indexing starts at 0 in Rcpp!
  Type piI = par(1); // careful, indexing starts at 0 in Rcpp!
  Type phiS = par(2);
  Type phiI = par(3);
  Type phiR = par(4);
  Type betaSI = par(5);
  Type pS = par(6);
  Type pI = par(7);
  Type pR = par(8);
  Type deltaS = par(9);
  Type deltaI = par(10);
  Type deltaR = par(11);
  
  // prob of obs (rows) cond on states (col)
  matrix<Type> B1(4,4);
  B1(0,0) = Type(1.0)-pS;
  B1(0,1) = pS;
  B1(0,2) = Type(0.0);
  B1(0,3) = Type(0.0);
  B1(1,0) = Type(1.0)-pI;
  B1(1,1) = Type(0.0);
  B1(1,2) = pI;
  B1(1,3) = Type(0.0);
  B1(2,0) = Type(1.0)-pR;
  B1(2,1) = Type(0.0);
  B1(2,2) = Type(0.0);
  B1(2,3) = pR;
  B1(3,0) = Type(1.0);
  B1(3,1) = Type(0.0);
  B1(3,2) = Type(0.0);
  B1(3,3) = Type(0.0);
  
  matrix<Type> B2(4,5);
  B2(0,0) = Type(1.0);
  B2(0,1) = Type(0.0);
  B2(0,2) = Type(0.0);
  B2(0,3) = Type(0.0);
  B2(0,4) = Type(0.0);
  B2(1,0) = Type(0.0);
  B2(1,1) = deltaS;
  B2(1,2) = Type(0.0);
  B2(1,3) = Type(0.0);
  B2(1,4) = 1-deltaS;
  B2(2,0) = Type(0.0);
  B2(2,1) = Type(0.0);
  B2(2,2) = deltaI;
  B2(2,3) = Type(0.0);
  B2(2,4) = Type(1.0)-deltaI;
  B2(3,0) = Type(0.0);
  B2(3,1) = Type(0.0);
  B2(3,2) = Type(0.0);
  B2(3,3) = deltaR;
  B2(3,4) = Type(1.0)-deltaR;

  matrix<Type> BB(4, 5);
  BB = B1 * B2;
  matrix<Type> B(5, 4);
  B = BB.transpose();
  REPORT(B);
  
  // first encounter
  matrix<Type> BE1(4,4);
  BE1(0,0) = Type(0.0);
  BE1(0,1) = Type(1.0);
  BE1(0,2) = Type(0.0);
  BE1(0,3) = Type(0.0);
  BE1(1,0) = Type(0.0);
  BE1(1,1) = Type(0.0);
  BE1(1,2) = Type(1.0);
  BE1(1,3) = Type(0.0);
  BE1(2,0) = Type(0.0);
  BE1(2,1) = Type(0.0);
  BE1(2,2) = Type(0.0);
  BE1(2,3) = Type(1.0);
  BE1(3,0) = Type(1.0);
  BE1(3,1) = Type(0.0);
  BE1(3,2) = Type(0.0);
  BE1(3,3) = Type(0.0);
  
  matrix<Type> BE2(4,5);
  BE2(0,0) = Type(1.0);
  BE2(0,1) = Type(0.0);
  BE2(0,2) = Type(0.0);
  BE2(0,3) = Type(0.0);
  BE2(0,4) = Type(0.0);
  BE2(1,0) = Type(0.0);
  BE2(1,1) = deltaS;
  BE2(1,2) = Type(0.0);
  BE2(1,3) = Type(0.0);
  BE2(1,4) = Type(1.0)-deltaS;
  BE2(2,0) = Type(0.0);
  BE2(2,1) = Type(0.0);
  BE2(2,2) = deltaI;
  BE2(2,3) = Type(0.0);
  BE2(2,4) = Type(1.0)-deltaI;
  BE2(3,0) = Type(0.0);
  BE2(3,1) = Type(0.0);
  BE2(3,2) = Type(0.0);
  BE2(3,3) = deltaR;
  BE2(3,4) = Type(1.0)-deltaR;

  matrix<Type> BEE(4,5);
  BEE = BE1 * BE2;
  matrix<Type> BE(5,4);
  BE = BEE.transpose();
  REPORT(BE);
  
  // prob of states at t+1 given states at t
  matrix<Type> A1(4,4);
  A1(0,0) = phiS;
  A1(0,1) = Type(0.0);
  A1(0,2) = Type(0.0);
  A1(0,3) = Type(1.0)-phiS;
  A1(1,0) = Type(0.0);
  A1(1,1) = phiI;
  A1(1,2) = Type(0.0);
  A1(1,3) = Type(1.0)-phiI;
  A1(2,0) = Type(0.0);
  A1(2,1) = Type(0.0);
  A1(2,2) = phiR;
  A1(2,3) = Type(1.0)-phiR;
  A1(3,0) = Type(0.0);
  A1(3,1) = Type(0.0);
  A1(3,2) = Type(0.0);
  A1(3,3) = Type(1.0);

  matrix<Type> A2(4,4);
  A2(0,0) = Type(1.0)-betaSI;
  A2(0,1) = betaSI;
  A2(0,2) = Type(0.0);
  A2(0,3) = Type(0.0);
  A2(1,0) = Type(0.0);
  A2(1,1) = Type(0.0);
  A2(1,2) = Type(1.0);
  A2(1,3) = Type(0.0);
  A2(2,0) = Type(0.0);
  A2(2,1) = Type(0.0);
  A2(2,2) = Type(1.0);
  A2(2,3) = Type(0.0);
  A2(3,0) = Type(0.0);
  A2(3,1) = Type(0.0);
  A2(3,2) = Type(0.0);
  A2(3,3) = Type(1.0);

  matrix<Type> A(4,4);
  A = A1 * A2;
  REPORT(A);
  
  // init states
  vector<Type> PROP(4);
  PROP(0) = piS;
  PROP(1) = piI;
  PROP(2) = Type(1.0)-piS-piI;
  PROP(3) = Type(0.0);
  REPORT(PROP);
  
  // likelihood
  Type ll;
  Type nll;
  array<Type> ALPHA(4);
  for (int i = 0; i < nh; i++) {
    int ei = fc(i)-1;
    vector<int> evennt = ch.col(i);
    ALPHA = PROP * vector<Type>(BE.row(fs(i))); // element-wise vector product
    for (int j = ei+1; j < km; j++) {
      ALPHA = multvecmat(ALPHA,A) * vector<Type>(B.row(evennt(j))); // vector matrix product, then element-wise vector product
    }
    ll += log(sum(ALPHA));
  }
  nll = -ll;
  return nll;
}"
write(tmb_model, file = "sir_tmb.cpp")
```

Then load the model template:

```r
library(TMB)
compile("sir_tmb.cpp")
```

```
## [1] 0
```

```r
dyn.load(dynlib("sir_tmb"))
```

# Benchmarking

We do some benchmarking by hand. Note that we could have used a more formal approach using the `microbenchmark` package for example.
We will repeat the simulation and parameter estimation 5 times to evaluate the time it takes to run these steps in `TMB` vs. native `R`.
The input parameters are set up by default, see the `simul` function above. We consider 20 sampling occasions. The results are stored in 3 matrices that have as many rows as the number of simulations:   

* `tab_optim` contains the parameter estimates from the native `R` runs;   

* `tab_tmb` contains the parameter estimates from the native `TMB` runs;   

* `tab_time` contains the running times for `R` in the first column and for `TMB` in the second column.   



```r
MCiter <- 5
tab_optim <- matrix(0,MCiter,12)
tab_tmb <- matrix(0,MCiter,12)
tab_time <- matrix(0,MCiter,2) 

for(z in 1:MCiter){
	data <- simul(n.occasions=20)
	nh <- dim(data)[1]
	k <- dim(data)[2]
	km1 <- k-1
	eff <- rep(1,nh)
    fc <- NULL
    init.state <- NULL
    for (i in 1:nh){
      temp <- 1:k
      fc <- c(fc,min(which(data[i,]!=0)))
      init.state <- c(init.state,data[i,fc[i]])
    }
    binit <- runif(12,-1,0)
    data <- t(data)
    
    # native R optim-like optimisation
    deb <- Sys.time()
    tmp_optim <- optim(binit,devMULTIEVENT,NULL,hessian=T,data,eff,fc,init.state,nh,km1,method="BFGS")
    fin <- Sys.time()
    res_optim <- fin-deb
    x <- tmp_optim$par
    piS <- plogis(x[1])
    piI <- plogis(x[2])
    # 1-piS-piI = 1 / (1+exp(par[1])+exp(par[2]))
    phiS <- plogis(x[3])
    phiI <- plogis(x[4])
    phiR <- plogis(x[5])
    psiSI <- plogis(x[6])
    pS <- plogis(x[7])
    pI <- plogis(x[8])
    pR <- plogis(x[9])
    deltaS <- plogis(x[10])
    deltaI <- plogis(x[11])
    deltaR <- plogis(x[12])
    par_optim <- c(piS, piI, phiS, phiI, phiR, psiSI, pS, pI, pR, deltaS, deltaI, deltaR) 

    # TMB-like optimisation    
    deb <- Sys.time()
    f <- MakeADFun(data = list(ch = data, fc = fc, fs = init.state),parameters = list(b = binit),DLL = "sir_tmb")
    opt <- do.call("optim", f) # optimisation
    fin <- Sys.time()
    res_tmb <- fin-deb
    x <- opt$par
    piS <- plogis(x[1])
    piI <- plogis(x[2])
    # 1-piS-piI = 1 / (1+exp(par[1])+exp(par[2]))
    phiS <- plogis(x[3])
    phiI <- plogis(x[4])
    phiR <- plogis(x[5])
    psiSI <- plogis(x[6])
    pS <- plogis(x[7])
    pI <- plogis(x[8])
    pR <- plogis(x[9])
    deltaS <- plogis(x[10])
    deltaI <- plogis(x[11])
    deltaR <- plogis(x[12])
    par_tmb <- c(piS, piI, phiS, phiI, phiR, psiSI, pS, pI, pR, deltaS, deltaI, deltaR) 

    tab_optim[z,] <- par_optim
    tab_tmb[z,] <- par_tmb
    tab_time[z,1] <- res_optim
    tab_time[z,2] <- res_tmb
    
    }
```

```
## outer mgc:  1494.51 
## outer mgc:  260.9873 
## outer mgc:  2388.715 
## outer mgc:  249.3063 
## outer mgc:  240.9871 
## outer mgc:  156.0476 
## outer mgc:  196.6287 
## outer mgc:  173.0885 
## outer mgc:  83.93802 
## outer mgc:  44.59302 
## outer mgc:  54.32551 
## outer mgc:  62.23167 
## outer mgc:  58.03629 
## outer mgc:  49.95518 
## outer mgc:  56.03692 
## outer mgc:  39.08593 
## outer mgc:  77.00832 
## outer mgc:  97.07438 
## outer mgc:  114.6476 
## outer mgc:  167.8703 
## outer mgc:  101.3867 
## outer mgc:  94.33787 
## outer mgc:  82.68403 
## outer mgc:  106.2041 
## outer mgc:  83.30716 
## outer mgc:  65.91171 
## outer mgc:  65.00071 
## outer mgc:  41.3513 
## outer mgc:  38.83523 
## outer mgc:  28.71136 
## outer mgc:  25.51793 
## outer mgc:  15.34457 
## outer mgc:  19.51505 
## outer mgc:  22.78218 
## outer mgc:  23.40066 
## outer mgc:  13.7356 
## outer mgc:  15.93122 
## outer mgc:  14.12033 
## outer mgc:  4.03297 
## outer mgc:  1.944249 
## outer mgc:  1.261293 
## outer mgc:  0.404293 
## outer mgc:  0.2077742 
## outer mgc:  0.2964813 
## outer mgc:  0.1516369 
## outer mgc:  0.2876795 
## outer mgc:  0.2650137 
## outer mgc:  0.2082536 
## outer mgc:  0.1001712 
## outer mgc:  0.08147795 
## outer mgc:  0.12333 
## outer mgc:  1530.766 
## outer mgc:  261.3896 
## outer mgc:  260.7276 
## outer mgc:  270.4689 
## outer mgc:  291.9915 
## outer mgc:  312.7961 
## outer mgc:  380.5812 
## outer mgc:  177.673 
## outer mgc:  130.6905 
## outer mgc:  119.091 
## outer mgc:  106.0872 
## outer mgc:  86.0109 
## outer mgc:  72.18761 
## outer mgc:  78.15996 
## outer mgc:  75.7967 
## outer mgc:  36.99999 
## outer mgc:  43.75011 
## outer mgc:  105.3225 
## outer mgc:  63.76284 
## outer mgc:  60.27513 
## outer mgc:  62.50276 
## outer mgc:  79.75144 
## outer mgc:  72.69936 
## outer mgc:  35.37551 
## outer mgc:  34.09912 
## outer mgc:  48.99097 
## outer mgc:  28.36175 
## outer mgc:  13.94525 
## outer mgc:  5.774393 
## outer mgc:  9.971669 
## outer mgc:  6.781587 
## outer mgc:  5.86015 
## outer mgc:  5.673892 
## outer mgc:  5.803527 
## outer mgc:  4.507319 
## outer mgc:  4.209737 
## outer mgc:  2.904546 
## outer mgc:  1.002112 
## outer mgc:  0.1623214 
## outer mgc:  0.04144346 
## outer mgc:  0.1294262 
## outer mgc:  1.342533 
## outer mgc:  1.740284 
## outer mgc:  0.5586181 
## outer mgc:  0.3950054 
## outer mgc:  0.676709 
## outer mgc:  0.7760445 
## outer mgc:  0.5857662 
## outer mgc:  0.4157866 
## outer mgc:  0.3113676 
## outer mgc:  0.4413554 
## outer mgc:  0.6707272 
## outer mgc:  1.491758 
## outer mgc:  1.7234 
## outer mgc:  0.9675715 
## outer mgc:  1.328106 
## outer mgc:  0.9302313 
## outer mgc:  1.099889 
## outer mgc:  0.8010097 
## outer mgc:  2.167108 
## outer mgc:  2.16343 
## outer mgc:  1.732482 
## outer mgc:  2.212341 
## outer mgc:  2.223618 
## outer mgc:  2.391895 
## outer mgc:  0.6263815 
## outer mgc:  0.4294523 
## outer mgc:  0.2546599 
## outer mgc:  0.1581098 
## outer mgc:  0.02936748 
## outer mgc:  0.02794714 
## outer mgc:  0.008340797 
## outer mgc:  1307.673 
## outer mgc:  190.9755 
## outer mgc:  452.9882 
## outer mgc:  239.5127 
## outer mgc:  259.443 
## outer mgc:  231.5778 
## outer mgc:  172.6195 
## outer mgc:  116.7391 
## outer mgc:  83.0225 
## outer mgc:  47.50723 
## outer mgc:  48.73639 
## outer mgc:  21.4441 
## outer mgc:  24.23182 
## outer mgc:  38.84503 
## outer mgc:  16.05417 
## outer mgc:  12.2775 
## outer mgc:  10.46805 
## outer mgc:  2.219772 
## outer mgc:  5.520827 
## outer mgc:  7.221116 
## outer mgc:  4.504224 
## outer mgc:  2.242252 
## outer mgc:  1.743387 
## outer mgc:  1.195546 
## outer mgc:  0.9267127 
## outer mgc:  0.8882692 
## outer mgc:  1.017257 
## outer mgc:  0.5554739 
## outer mgc:  0.4739385 
## outer mgc:  0.4857425 
## outer mgc:  0.4809783 
## outer mgc:  0.2409444 
## outer mgc:  0.2169449 
## outer mgc:  0.1643512 
## outer mgc:  0.07596796 
## outer mgc:  0.1398887 
## outer mgc:  0.1195718 
## outer mgc:  0.0709888 
## outer mgc:  0.008896112 
## outer mgc:  0.008809542 
## outer mgc:  0.05021833 
## outer mgc:  0.08146955 
## outer mgc:  0.09754336 
## outer mgc:  0.08053354 
## outer mgc:  0.03157755 
## outer mgc:  0.006591604 
## outer mgc:  0.01668185 
## outer mgc:  0.01698102 
## outer mgc:  0.00819973 
## outer mgc:  1381.445 
## outer mgc:  239.1278 
## outer mgc:  241.5177 
## outer mgc:  220.0127 
## outer mgc:  176.9585 
## outer mgc:  174.8797 
## outer mgc:  156.0975 
## outer mgc:  151.6988 
## outer mgc:  91.22869 
## outer mgc:  97.85711 
## outer mgc:  76.96097 
## outer mgc:  45.61274 
## outer mgc:  54.70538 
## outer mgc:  46.00466 
## outer mgc:  32 
## outer mgc:  100.2498 
## outer mgc:  113.2193 
## outer mgc:  43.77251 
## outer mgc:  150.5963 
## outer mgc:  37.32339 
## outer mgc:  36.20879 
## outer mgc:  46.41986 
## outer mgc:  46.1904 
## outer mgc:  46.29101 
## outer mgc:  45.8399 
## outer mgc:  26.71226 
## outer mgc:  27.46837 
## outer mgc:  17.66209 
## outer mgc:  14.47132 
## outer mgc:  9.780992 
## outer mgc:  18.98239 
## outer mgc:  13.19483 
## outer mgc:  12.53809 
## outer mgc:  9.676407 
## outer mgc:  8.156408 
## outer mgc:  6.366084 
## outer mgc:  4.004398 
## outer mgc:  3.803764 
## outer mgc:  0.8559731 
## outer mgc:  0.878299 
## outer mgc:  0.2680634 
## outer mgc:  0.1246542 
## outer mgc:  0.4281125 
## outer mgc:  0.4814843 
## outer mgc:  0.3130686 
## outer mgc:  0.1181288 
## outer mgc:  0.07992162 
## outer mgc:  0.05851297 
## outer mgc:  0.02094862 
## outer mgc:  0.008959037 
## outer mgc:  1097.55 
## outer mgc:  268.6841 
## outer mgc:  268.5449 
## outer mgc:  258.3604 
## outer mgc:  228.5439 
## outer mgc:  221.5505 
## outer mgc:  263.2446 
## outer mgc:  252.8912 
## outer mgc:  248.5706 
## outer mgc:  112.8767 
## outer mgc:  101.591 
## outer mgc:  107.1722 
## outer mgc:  58.42133 
## outer mgc:  59.58718 
## outer mgc:  74.97324 
## outer mgc:  89.62676 
## outer mgc:  85.26705 
## outer mgc:  68.74189 
## outer mgc:  49.20698 
## outer mgc:  144.8106 
## outer mgc:  183.3379 
## outer mgc:  153.9409 
## outer mgc:  128.7477 
## outer mgc:  75.30231 
## outer mgc:  35.4622 
## outer mgc:  13.14304 
## outer mgc:  12.19207 
## outer mgc:  7.866475 
## outer mgc:  7.842383 
## outer mgc:  12.24108 
## outer mgc:  8.937254 
## outer mgc:  6.958744 
## outer mgc:  6.353576 
## outer mgc:  6.054076 
## outer mgc:  5.460202 
## outer mgc:  4.020603 
## outer mgc:  3.335044 
## outer mgc:  5.466193 
## outer mgc:  1.422917 
## outer mgc:  0.9003633 
## outer mgc:  0.6288708 
## outer mgc:  0.2731197 
## outer mgc:  0.3049809 
## outer mgc:  0.05112294
```

The estimates for $(\phi_S,\phi_I,\phi_R,\beta_{SI},p_S,p_I,p_R,\delta_S,\delta_I,\delta_R)$ we got with `R` are:

```r
apply(tab_optim,2,mean)[-c(1,2)]
```

```
##  [1] 0.96341206 0.52945051 0.84148944 0.57102704 0.95315620 0.91918038
##  [7] 0.94819202 0.07892201 0.09685499 0.10397751
```

while those obtained with `TMB` are:

```r
apply(tab_tmb,2,mean)[-c(1,2)]
```

```
##  [1] 0.96491367 0.52897752 0.84145416 0.57193338 0.95380304 0.91536888
##  [7] 0.94832532 0.07874303 0.09685467 0.10402125
```

which are to be compared with the values we used to simulate the data:

```r
c(0.96,0.53,0.84,0.60,0.95,0.95,0.95,0.1,0.1,0.1)
```

```
##  [1] 0.96 0.53 0.84 0.60 0.95 0.95 0.95 0.10 0.10 0.10
```


The time elapsed for `R` is 20.6571664 while 0.1811735 for `TMB`; in other words, `TMB` is 114.0187173 faster than `R`.

