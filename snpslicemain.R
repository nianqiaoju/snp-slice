## execution file for the slice sampler
## this file can 
## 1) read input data;
## 2) read command line tuning parameters
## 3) run the mcmc algorithm 
## 4) save output to txt;
## last update: 5/15/2023
## ???: nianqiao@purdue.edu

## set the correct working directory before using the code
# setwd()
#
## import data -----
rm(list = ls());


################# 
################# start of user determined values
################# 
## input file name----
prefix <- "scenario1";
inputdataname <- paste("inputdata/", prefix, sep = "");
## model----
# model <- 0 ## cat
# model <- 1 ## "pois"
# model <- 2 ## "bin";
model <- 3 ## negbin
## prior parameters
rho <- 0.5;
alpha <- 2.6;
## mcmc parameters
nmcmc <- 10000;
gap <- 20; ## terminate the chain if lpost does not change in > gap steps
burnin <- 0; ## this is the minimum number of MCMC iterations, default is 0.
rep <- 1;
## decide which information to store
store_mcmc <- FALSE;
store_sink <- FALSE;
store_map <- TRUE;

################# 
################# end of user determined values
################# 

## read values from command lind input
invisible(eval(parse(text=commandArgs(TRUE))));
if(is.na(burnin)) burnin <- floor(nmcmc / 2); ## default for burnin
## threshold to determine single infections
threshold <- 0.001;
## parameters for the observation model
e1 <- 0.05;
e2 <- 0.05;

set.seed(rep);
## data processing-----
model <- c("cat", "pois", "bin", "neg", "jpois")[model + 1];
if(model == "cat"){
  dfy <- read.delim(paste(prefix, "_cat.txt", sep = ""));
  # remove individual id
  dfy$host <- NULL;
  dfy$ind_name <- NULL;
  N <- dim(dfy)[1];
  P <- dim(dfy)[2];
  y <- as.matrix(dfy);
  r <- NA;
}else{
  ## for pois and bin models
  y <- as.matrix(read.delim(paste(prefix, "_read1.txt", sep = ""))[,-1]);
  r <- as.matrix(read.delim(paste(prefix, "_read0.txt", sep = ""))[,-1]) + y;
  N <- dim(y)[1];
  P <- dim(y)[2];
  rho <- sum(y, na.rm = T) / sum(r, na.rm = T);
}


# load helper functions for slice sampler-----
source(paste("source/snpslice", model, ".R", sep = ""));
source("source/snpslicecommon.R");


cat("running", model, "model on data", prefix, "\n");
cat("N =", N, "P=", P, "\n");
cat("initialization begin\n");
state <- snp_init(thres = threshold);

cat("initialization done, at", state$loglik, "\n");
cat("starting with",  sum(rowSums(state$A) ==1), "single infections\n")

state <- slice_init(state);
cat("plan to run", nmcmc, "iterations of MCMC with", burnin, "iterations for burn in, and gap =", gap, "\n");

map <- list();
lpostmax <- -Inf;
mapiter <- 0;
mapktrunc <- state$ktrunc;
if(store_sink){
  sink(file = paste("output/sink_", prefix,"_", model, "_nmcmc", nmcmc, "_gap", gap,"_rep", rep, ".txt", sep = ""), split = TRUE, type = "output");
}
for(iter in 1 : nmcmc){
	state <- sliceIter(state);
	cat(iter,  sum(colSums(state$A) > 0), state$kstar, dim(state$A)[2],  sum(rowSums(state$A) ==1), state$logpost, lpostmax, "\n");
	if(state$ktrunc > mapktrunc){## if ktrunc changes, restart map
	  map <- state;
	  mapiter <- iter;
	  mapktrunc <- state$ktrunc;
	  lpostmax <- state$logpost;
	}else if(state$logpost > lpostmax){## if still the same ktrunc, but lpost increases,
    map <- state;
    mapiter <- iter;
    lpostmax <- state$logpost;
  }
	if(iter > burnin & mapiter < iter - gap) break;
}
if(store_sink){
  sink();
}

## save the final sample (MAP estimator)
A <- map$A[,which(colSums(map$A) > 0)];
write.table(A, paste("output/", prefix, "_", model, "_A_nmcmc", nmcmc, "_gap", gap, "_rep", rep, ".txt", sep = ""),
            append = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE);
D <- map$D[which(colSums(map$A) > 0),];
write.table(D, paste("output/", prefix, "_", model, "_D_nmcmc", nmcmc, "_gap", gap, "_rep", rep, ".txt", sep = ""),
            append = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE);

## save MCMC state to RData file (for warm starts)
if(store_mcmc){
  filename <- paste("output/mcmc_", prefix, "_", model, "_nmcmc", nmcmc, "_gap", gap, "_rep", rep,  ".RData", sep = "");
  save(state, rho, alpha, e1, e2, model,r, y,  rep, prefix, file = filename);
}
