## this script contains helper functions used by
## slice samplers for all models:
## 1) binomial
## 2) poisson
## 3) categorical

## the helper functions come from the following categories:
## 1) compute MH acceptance ratio
## 2) compute prior densities
## 3) slice sampling steps to update s and mu given the rest
## 4) gibbs sampling steps to update A and D given the rest 
## 5) initialize slice sampler after (A,D) has been initialized for a specific model
## 6) remove duplicates in the dictionary D and track the assignments


## the following functions are in the methods files:
## 1) initialize the slice sampler from a valid state
## 2) compute the likelihood l(Y | A, D)

## last updated: 4/20/2023
## ??: nianqiao@purdue.edu

## helper functions------

## Metropolis-Hastings----
## MH ratio, which is P(x = 1) / P(x = 0) 
get_mhratio <- function(logp1, logp0){
  if(is.infinite(logp1)) return(0);
  if(is.infinite(logp0)) return(1);

  maxlogp <- max(logp1, logp0);
  mhratio <- exp(logp1 - maxlogp) / (exp(logp1 - maxlogp) + exp(logp0 - maxlogp));
  return(mhratio);
}

logpriorA <- function(A, mu){
  ## length of mu is the same as column-size of A
  lpr <- sum(t(A) * log(mu));
  lpr <- lpr + sum(t(1-A) * log(1-mu));
  return(lpr);
}

logpriorMu <- function(mu){
  lpr <- length(mu) * log(alpha) + alpha * log(tail(mu,1)) - sum(log(mu))
  return(lpr);
}

logpriorD <- function(D){
  sumD <- sum(D);
  lpr <- sumD * log(rho) + (prod(dim(D)) - sumD) * log(1-rho)
  return(lpr);
}

get_m <- function(A){
  ## compute the colSum sum of A
  ## this is also the frequency/count of strains
  return(colSums(A));
}

get_kstar <- function(A){
  ## kstar is the maximal feature index with mu[k] > s. 
  ## which means a[n,k] = 0 for all k > kstar
  ## this is a function of A only
  return(tail(which(colSums(A)>0),1));
}

get_mustar <- function(A, mu){
  ## mustar is the mu value corresponding to the last active feature
  ## this is equivalent to mu[get_star(A)];
  return(mu[tail(which(colSums(A) != 0),1)]);
}

## for ARS
logf_newfeature <- function(x){
  ## log density used to generate mu for a new & inactive feature
  logf <- alpha * sum(sapply(1 : N, function(i) (1-x)**i / i));
  logf <- logf + (alpha - 1) * log(x) + N * log(1-x);
  return(logf);
}

logf_oldfeature <- function(x, m){ 
  ## log density used to update mu an existing feature
  ## m is the observed count of that feature appearing in A
  return((m-1) * log(x) + (N - m) * log(1-x));
}

sampleIndexLogWeights <- function(lw){
  ## random sampling from the categorical distribution
  ## returns an index according to log-weights;
  ## lw: logweights
  ## this function is used in gridsample
  lw <- lw - max(lw);
  return(sample.int(n = length(lw), size = 1, prob = exp(lw) / sum(exp(lw))));
}

gridsample <- function(lb, ub, logf, nsample = 100){
  ## sample from a grid, given lower bound, upper bound, and log density (unnormalized)
  ## this function cannot return values == lb or ub 
  ## we can also use adaptive rejection sampling to achieve the same goal
  xgrid <- seq(lb, ub, length.out = nsample); ## create an equally spaced grid between end points
  xgrid <- xgrid[-c(1,nsample)]; ## exclude end points
  lw <- logf(xgrid); ## compute log density at each point on the grid
  return(xgrid[sampleIndexLogWeights(lw)]);  ## sample according to log weights
}

refreshFeature <- function(state, k){
  ## for the kth feature, set allocation to zero
  ## and instantiate the feature from prior for the dictionary
  state$A[, k] <- 0;
  state$D[k, ] <- (runif(P) < rho);
  return(state);
}


## used for initialization 



removeDuplicates <- function(dd){
  ## dd is a binary matrix, each row is a strain, potentially there are duplicated rows
  metric <- 1 -  (dd %*% t(dd)  + (1-dd) %*% t(1-dd)) / dim(dd)[2];
  ## metric is the pairwise distance between strains
  ## if two strains are the same, then metric[k1,k2]=0
  assignments <- rep(NA, dim(dd)[1]);
  ## assignment tracks "who is who"
  counter <- 0;
  firstappear <- c();
  for(i in 1 : dim(dd)[1]){
    ## for every person, 
    ## if not inspected yet (check assignment)
    ## find its twins 
    if(is.na(assignments[i])){
      counter <- counter + 1;
      twins <- which(metric[i,] == 0);
      assignments[twins] <- counter;
      firstappear <- c(firstappear, i);
    }
  }
  return(list(assignments = assignments,
              D = dd[firstappear,]));
}


slice_init <- function(state){
  ## prepare state for slice sampler
  ## state: an list object with A, D, already initialized for the observation model
  ## state$A: a matrix of dimension N * K
  ## state$D: a matrix of dimension K * P
  ## state$mu: stick length 
  ## state$ktrunc: length of mu
  ## state$kstar: index for last active feature A[i,k]>0
  ## state$kplus: first inactive feature. kplus = kstar + 1
  ## state$logpost: log (unormalized) posterior density
  ## state$kmin: lowerbound of global K. This is (# of single infections + 1)
  ## state$mixed: index of mixed infection individuals


  ## initialize mu for slice breaking construction
  state$ktrunc <- dim(state$A)[2];
  mu <- rbeta(state$ktrunc + 1, alpha / (state$ktrunc + 1), 1)
  muSort <- sort(mu, decreasing = TRUE);
  state$mu <- muSort;
  ## store values of kstar and kplus
  state$kstar <- get_kstar(state$A);## maximal feature index for active features A[i,k]>0
  state$kplus <- state$kstar + 1; ## this index is larger than that of any active feature. 
  state$A <- cbind(state$A, 0);
  state$D <- rbind(state$D, (runif(P) < rho));
  state$ktrunc <- dim(state$A)[2];
  return(state);
}


## update S-----
## s is the slice variable
sliceUpdateS <- function(state){
  ## update the auxiliary variable according to the Uniform distribution
  s <- get_mustar(state$A, state$mu) * runif(1); ## 0 < s < mu[kstar]
  ## expand the length of mu (if needed)
  state$ktrunc <- dim(state$A)[2]; ## this is the old truncation level
  k <- state$ktrunc; ## k will change
  while(s < state$mu[k]){ ## if we need to break more sticks
    ## sample new values of mu, using a grid
    ## they correspond to 'empty' features and are smaller than existing mu values
    munext <- gridsample(0, state$mu[k], logf_newfeature); 
    state$mu <- c(state$mu, munext);
    k <- k+1;
  }
  rm(k);
  ## expand A and D to match kplus, if necessary
  if(state$ktrunc < length(state$mu)){## if we broke more sticks
    numNewColumn <- length(state$mu) - state$ktrunc;
    state$A <- cbind(state$A, matrix(0, nrow = N, ncol = numNewColumn));
    state$D <- rbind(state$D, matrix(0, ncol = P, nrow = numNewColumn));
    ## can improve the implementation later
    for(k in (state$ktrunc  + 1) : length(state$mu)){
      state <- refreshFeature(state, k); ## refresh new features
    }
  }
  state$ktrunc <- length(state$mu);
  state$kplus <- which(state$mu < s)[1]; ## kplus is the first index with mu[k]<s
  # if(is.na(state$kplus)) state$kplus <- state$ktrunc; ## this line is used during debug
  return(state);
}


## update A given the rest--------
sliceUpdateALocal <- function(state, i, k){
  ad0 <- state$A[i,-k] %*% state$D[-k, ];
  a0 <- sum(state$A[i,-k]);
  if(a0 == 0){
    state$A[i,k] <- 1;
    state$kstar <- get_kstar(state$A);
    return(state);
  }
  logp0 <- loglikelihood_vector(as.vector(ad0/a0), 
                               as.vector(y[i,]),
                                as.vector(r[i,]));
  logp1 <- loglikelihood_vector(as.vector((ad0 + state$D[k,]) / (a0 + 1)),
                                as.vector(y[i,]), 
                                as.vector(r[i,]));
  logp0 <- logp0 + log(1-state$mu[k]);
  logp1 <- logp1 + log(state$mu[k]);
  ## there are two cases where changing A[i,k] can change kstar
  if(k == state$kstar & state$A[i,k] == 1 & sum(state$A[-i,k]) == 0){
    logp1 <- logp1 - log(state$mu[k]);
    nextkstar <- tail(which(colSums(state$A) !=0),2)[1];
    logp0 <- logp0 - log(state$mu[nextkstar]);
    # cat("updating kstar, i= ", i, "k=", k,"\n");
  }else if(k > state$kstar){
    logp1 <- logp1 - log(state$mu[k]);
    logp0 <- logp0 - log(state$mu[state$kstar]);
    # cat("updating an inactive feature, i= ", i, "k=", k,"\n");
  }
  p1 <- get_mhratio(logp1, logp0);
  ## determine A[i,k] with coin flip
  u <- runif(1);
  state$A[i,k] = (u < p1);
  state$kstar <- get_kstar(state$A);
  return(state);
}

sliceUpdateA <- function(state){
  for(k in 1 : state$kplus){ ## update active features
    for(i in state$mixed){ ## update mixed infection hosts
      state <- sliceUpdateALocal(state, i, k);
    }
  }
  state$kstar <- get_kstar(state$A);
  return(state);
}

## update D given the rest-----
sliceUpdateDLocal <- function(state, k , p){
  ad0 <- state$A[,-k] %*% state$D[-k, p];
  an <- rowSums(state$A);## there is an error if any rowSum is 0
  logp1 <- log(rho); logp0 <- log(1-rho);
  ## compute p(y | d) under the observation model
  logp1 <- logp1 + loglikelihood_vector( (ad0 + state$A[,k]) / an, y[,p], r[ ,p]);
  logp0 <- logp0 + loglikelihood_vector(ad0 / an, y[,p], r[,p]);
  p1 <- get_mhratio(logp1, logp0);
  u <- runif(1);
  state$D[k,p] = (u < p1);
  return(state);
}

sliceUpdateD <- function(state){
  ## only need up update D upto k = kstar
  for(k in state$kmin : state$kstar){ ## update strains upto kstar, fix single infection dictionary
    for(p in 1 : P){
      state <- sliceUpdateDLocal(state, k,p);
    }
  }
  return(state);
}

## update mu-----
sliceUpdateMu <- function(state){
  m <- get_m(state$A);
  ## first Mu k = 1, mu[1] <= 1
  state$mu[1] <- gridsample(state$mu[2], 1, function(x) logf_oldfeature(x, m[1]));
  ## Other index, k = 2,3,...., Kplus - 1
  for(k in 2 : (state$kplus-1)){
    state$mu[k] <- gridsample(state$mu[k+1], state$mu[k-1], function(x) logf_oldfeature(x, m[k]));
  }
  ## kplus is the first inactive feature with mu[k] < s.
  ## update the tail (inactive features): kplus, kplus+1,...,ktrunc
  for(k in state$kplus : state$ktrunc){
    state$mu[k] <- gridsample(0, state$mu[k - 1], logf_newfeature);
  }
  return(state);
}


## slice sampling wrapper
sliceIter <- function(state){
  state <- sliceUpdateS(state);
  state <- sliceUpdateA(state);
  state <- sliceUpdateD(state);
  state <- sliceUpdateMu(state);
  state$loglik <- loglikelihood_matrix(state$A, state$D)
  state$logpost <- state$loglik + logpriorA(state$A, state$mu) + logpriorMu(state$mu) + logpriorD(state$D);
  return(state);
}
