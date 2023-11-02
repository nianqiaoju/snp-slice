## this script contains functions specific to the categorical model
## it contains functions that 
## 1) compute the likelihood p(Y | A, D), with matrix or vector inputs
## 2) initialize the MCMC chain

## updated: 5/15/2023
## ???: nianqiao@purdue.edu

# loglikelihood function for this model
## prepare the log-likelihood table
logf0 <- log(c(1-e1,e1,0));
logf1 <- log(c(0,e1,1-e1));
logfmix <- log(c(e2/2,1-e2,e2/2));
llik_tab <- matrix(0, nrow = 3, ncol = 3);
llik_tab[1,] <- logf0;
llik_tab[2,] <- logfmix;
llik_tab[3,] <- logf1;


loglikelihood_vector <- function(propvec, yvec, rvec = NA){
  propvec <- cut(propvec, c(-Inf, 0,0.99,1.1), labels = c(0,0.5,1));
  yvec <- factor(yvec, levels = c(0,0.5,1));
  return(sum(table(propvec, yvec) * llik_tab, na.rm = T));
}


loglikelihood_matrix <- function(A, D){
  proportions_ <- (A %*% D) / rowSums(A); 
  return(sum(sapply(1 : P, function(p) loglikelihood_vector(proportions_[,p], y[,p]))));
}

## initialization for the categorical model
## cat model has different exceptions from bin and pois models
resolve_exceptions <- function(state){
  ratios <- (state$A %*% state$D) / rowSums(state$A); 
  ## when y == 0, cannot have pnp == 1 ;
  exceptions2 <- which(ratios == 1 & y ==0, arr.ind = TRUE);
  state$A[exceptions2[,1],state$kmin] <- 1;
  state$D[state$kmin,exceptions2[,2]] <- 0;
  ratios <- (state$A %*% state$D) / rowSums(state$A); 
  ## when y == 1, we cannot have pnp == 0;
  exceptions1 <- which(y == 1 & ratios == 0, arr.ind = TRUE);
  state$A[exceptions1[,1], state$kmin + 1] <- 1;
  state$D[state$kmin + 1 , exceptions1[,2]] <- 1;
  state$loglik <- loglikelihood_matrix(state$A, state$D);
  return(state);
}

snp_init <- function(thres = threshold){
  ## replace NA or 0.5 in y with 0 or 1 (randomly)
  cate <- y;
  for(i in 1 : N){
    for(j in 1 : P){
      if(is.na(y[i,j]) | y[i,j] == 0.5){
        cate[i,j] <- (runif(1) < 0.5);
      }
    }
  }  
  ratios2dict <- removeDuplicates(cate);
  assignments <- ratios2dict$assignments;
  nstrain <- length(unique(assignments));
  state <- list();
  state$D <- ratios2dict$D;
  state$A <- matrix(0, nrow = N, ncol = nstrain);
  ## first find patients that have MOI = 1
    ## find candidates for single infections
  is_single <- apply(y, 1, function(x) all(x == 0 | x == 1));
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

