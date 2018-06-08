
# If you have any problem with the function compile() from the TMB package (error status 127), 
# we recommend to follow the recommendations provided here: https://github.com/nwfsc-assess/geostatistical_delta-GLMM/wiki/Steps-to-install-TMB
# Alternatively running these 2 lines may work:

Sys.setenv(PATH = paste("C:/Rtools/bin", Sys.getenv("PATH"), sep=";"))
Sys.setenv(BINPREF = "C:/Rtools/mingw_$(WIN)/bin/")

#------------------------------------------------------------ Start SIR model


# Bias and MSE of SIR multi-event capture-recapture parameter estimates

## Data simulation

#We first build a function to simulate data from a SIR multi-event capture-recapture model:
  
simul <- function(
  n.occasions = 5, # number of occasions
  n.states = 4, # number of states: S, I and R plus dead
  n.obs = 4, # number of events: 1) detected and diagnosed as S, 2) as I, 3) as R, 4) detected and undiagnosed
  phiS = 0.90, # survival susceptible
  phiI = 0.50, # survival infected
  phiR = 0.90, # survival recovered
  betaSI = 0.90, # infection prob 
  gammaIR = 0.30, # recovery probability 
  pS = 0.50, # detection susceptible
  pI = 0.50, # detection infected
  pR = 0.50, # detection recovered
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
  parsim <- c(marked[1,1]/tot, marked[1,2]/tot, phiS, phiI, phiR, betaSI, gammaIR, pS, pI, pR, deltaS, deltaI, deltaR) 
  
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
      
        # PSI.STATE[,,i,t] <- matrix(c( 
        # phiS*(1-betaSI), phiS*(betaSI), 0, 1-phiS, 
        # 0, 0, phiI, 1-phiI, 
        # 0, 0, phiR, 1-phiR, 
        # 0, 0, 0, 1 ), nrow = n.states, byrow = TRUE) 
    
        PSI.STATE[,,i,t] <- matrix(c( 
        phiS*(1-betaSI), phiS*(betaSI), 0, 1-phiS, 
        0, phiI*(1-gammaIR),phiI*(gammaIR), 1-phiI, # now with recovery
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
  
  return(list(titi2, parsim))
}


## Likelihood function in `TMB`

#First, create the model template:

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
Type gammaIR = par(6);
Type pS = par(7);
Type pI = par(8);
Type pR = par(9);
Type deltaS = par(10);
Type deltaI = par(11);
Type deltaR = par(12);

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
A2(1,1) = Type(1.0)-gammaIR;
A2(1,2) = gammaIR;
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


#Then load the model template:
library(TMB)
compile("sir_tmb.cpp") 


dyn.load(dynlib("sir_tmb"))


## Monte Carlo simulations 

#First, let us set up the simulations scenario (homogeneous or heterogeneous; here homogeneous):
  
# Homogeneous assignement probability of infection states:
hetS <- hetI <- hetR <- 1 
# Heterogeneous assignement probability of infection states
#hetS<- 0.5 
#hetI<- 0.5
#hetR<- 0.5


#Now set up the values for input parameters (cf. Table 1)
  
n.occ <- 5 # Number of capture occasions 
phiS_par <- 0.9 # Survival probability of S
phiI_par <- 0.5 # Survival probability of I
phiR_par <- 0.9 # Survival probability of R
betaSI_par <- 0.9 # Infection probability
gammaIR_par <- 0.3 # Recovery probability
pS_par <- 0.5 # Detection probability of S
pI_par <- 0.5 # Detection probability of I
pR_par <- 0.5 # Detection probability of R


#We need to do a few things to get ready to run the Monte Carlo simulations:
  
# sequence of value of assignement probabilities to represent the uncertainty gradient
uncertainty <- seq(0.1, 0.8, 0.1)

# number of Monte Carlo iterations (here just a few for illustration; use 1000 if you would like serious results, see paper)
MCiter <- 1000

# vector of the true parameter value (with assignement parameter updated at each level of uncertainty )
truevalue <- NULL 

# vector of estimates averaged over the MCiter simulations
estimates <- NULL 

# vector of all estimates (bias)
all_estimates <- NULL

# table or results
tab <- array(data = 0, dim = c(length(uncertainty), MCiter, 13)) 


#Let us rool!
  
# increment index; to be used in the uncertainty loop to fill in the table that stores the averaged parameter estimates for each level of assignment probability

inc <- 0 

for(d in unique(uncertainty)){
  
  inc <- inc + 1
  tab_tmb <- array(data=0, dim= c( MCiter, 13))
  
  for(z in 1:MCiter){
    
    output_simu <-simul(
      n.occasions = n.occ, 
      phiS = phiS_par, 
      phiI = phiI_par, 
      phiR = phiR_par, 
      betaSI = betaSI_par, 
      gammaIR = gammaIR_par,
      pS = pS_par,
      pI = pI_par,  
      pR = pR_par,  
      deltaS = d*hetS, 
      deltaI = d*hetI, 
      deltaR = d*hetR)
    
    data <- output_simu[[1]]
    if (z == MCiter) parsim <- output_simu[[2]]
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
    } # i loop
    
    # optimisation    
    binit <- runif(13,-1,0)
    data <- t(data)
    f <- MakeADFun(data = list(ch = data, fc = fc, fs = init.state),parameters = list(b = binit),DLL = "sir_tmb")
    opt <- do.call("optim", f) # optimisation
    x <- opt$par
    piS <- plogis(x[1])
    piI <- plogis(x[2])
    phiS <- plogis(x[3])
    phiI <- plogis(x[4])
    phiR <- plogis(x[5])
    psiSI <- plogis(x[6])
    gammaIR <- plogis(x[7])
    pS <- plogis(x[8])
    pI <- plogis(x[9])
    pR <- plogis(x[10])
    deltaS <- plogis(x[11])
    deltaI <- plogis(x[12])
    deltaR <- plogis(x[13])
    par_tmb <- c(piS, piI, phiS, phiI, phiR, psiSI, gammaIR,pS, pI, pR, deltaS, deltaI, deltaR) 
    tab_tmb[z,] <- par_tmb
  } # z loop
  
  tab[inc,,] <- tab_tmb
  truevalue <- rbind(truevalue, parsim)
  estimates <- rbind(estimates, apply(tab_tmb,2,mean))
  
  all_estimates <-rbind(estimates)
} # d loop


#Get bias:

bias <- truevalue - estimates
colnames(bias) <-
  c("piS",
    "piI",
    "phiS",
    "phiI",
    "phiR",
    "psiSI",
    "gammaIR",
    "pS",
    "pI",
    "pR",
    "deltaS",
    "deltaI",
    "deltaR")
rownames(bias) <- c("90%", "80%", "70%", "60%", "50%", "40%", "30%", "20%")
round(bias,2)

#Get MSE:

MSE <- matrix(0,nrow=length(uncertainty), ncol=13)
for(p in 1:13){
  for (inc in 1:(length(uncertainty))){
    MSE[inc,p] <- sum(tab[inc,,p] - truevalue[inc,p])^2 / MCiter 
  }
}
colnames(MSE) <-
  c("piS",
    "piI",
    "phiS",
    "phiI",
    "phiR",
    "psiSI",
    "gammaIR",
    "pS",
    "pI",
    "pR",
    "deltaS",
    "deltaI",
    "deltaR")
rownames(MSE) <- c("90%", "80%", "70%", "60%", "50%", "40%", "30%", "20%")
round(MSE,2)


