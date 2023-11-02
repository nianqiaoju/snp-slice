## this script contains functions specific to the binomial model
## it contains functions that 
## 1) compute the likelihood p(Y | A, D), with matrix or vector inputs
## 2) initialize the MCMC chain
## 
## updated: 5/15/2023
## ???: nianqiao@purdue.edu

## compute loglikelihood---
loglikelihood_matrix <- function(A, D){
  proportions_ <- (A %*% D) / rowSums(A); 
  return(sum(dbinom(y, r , proportions_,TRUE), na.rm = T));
}

loglikelihood_vector <- function(propvec, yvec, rvec){
  return(sum(dbinom(yvec, rvec ,  propvec, TRUE), na.rm = T));
}


## initialize a valid state --------
## such that loglikelihood is not zero

resolve_exceptions <- function(state){
  ## modify (A,D) from random initialization
  ## without touching single infections
  ratios <- (state$A %*% state$D) / rowSums(state$A); 
  ## when y != r, cannot have ratios == 1 ;
  exceptions2 <- which(ratios == 1 & y != r, arr.ind = TRUE);
  state$A[exceptions2[,1],state$kmin] <- 1;
  state$D[state$kmin,exceptions2[,2]] <- 0;
  ratios <- (state$A %*% state$D) / rowSums(state$A); 
  ## when y > 0, we must have ratios > 0  
  exceptions1 <- which(y!=0 & ratios == 0, arr.ind = TRUE);
  state$A[exceptions1[,1], state$kmin + 1] <- 1;
  state$D[state$kmin + 1, exceptions1[,2]] <- 1;
  state$loglik <- loglikelihood_matrix(state$A, state$D);
  return(state);
}

snp_init <- function(thres = threshold){
  ratios <- y / r;
  ## can use ratios to build a dictionary
  cate <- (y/r > 0.5);
  if(any(is.na(cate))){
    cat("NA detected\n");
    cate[which(is.na(cate))] <- 0;  
  }
  ratios2dict <- removeDuplicates(cate);
  assignments <- ratios2dict$assignments;
  nstrain <- length(unique(assignments));
  state <- list();
  state$D <- ratios2dict$D;
  state$A <- matrix(0, nrow = N, ncol = nstrain);
  ## first find patients that have MOI = 1
  is_single <- apply(ratios, 1, function(x) all(pmin(x, 1 - x) < thres));
  nsingleppl <- sum(is_single, na.rm = TRUE);
  if(nsingleppl >= 1){
    ## if there are single infections
    which_single <- which(is_single);
    which_mixed <- c(which(!is_single), which(is.na(is_single)));
    nsinglestrain <- length(unique(assignments[which_single]));
    state$A[which_mixed,] <- 1;
    for(isingle in 1 : nsingleppl){
      state$A[which_single[isingle], assignments[which_single[isingle]]] <- 1;
    }    
    ord <- rep(NA, nstrain);
    ord[1 : nsinglestrain] <- unique(assignments[which_single]);
    ord[(nsinglestrain + 1) : nstrain] <- setdiff(unique(assignments), unique(assignments[which_single]));
    state$A <- state$A[, ord];
    state$D <- state$D[ord,];
    state$mixed <- which_mixed;
    state$kmin <- nsinglestrain + 1;
    ## reorder the strains so that single infection strains appear first
    ## single infections have only one strain
    ## resolve exceptions
    state$loglik <- loglikelihood_matrix(state$A, state$D);
    iter <- 1;
    while(is.infinite(state$loglik)){
      state <- resolve_exceptions(state);
      iter <- iter + 1;
      cat("iteration", iter, "of initialization\n");
      state$loglik <- loglikelihood_matrix(state$A, state$D);
    }
    ## mixed infections have all strains
  }else{
    state$mixed <- 1 : N;
    state$kmin <- 1;
    state$A[state$mixed,] <- 1;
    state$loglik <-  loglikelihood_matrix(state$A, state$D);
    ## resolve exceptions
    iter <- 1;
    while(is.infinite(state$loglik)){
      state <- resolve_exceptions(state);
      iter <- iter + 1;
      cat("iteration", iter, "of initialization\n");
    }
  }
  return(state);
}
