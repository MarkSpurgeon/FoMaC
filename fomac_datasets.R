############################# description ####################################

# fomac_datasets.R simulates datasets to be used to evaluate the performance 
# of F-test of Concordance and Magnitude (FoMaC) and comparable methods. Each 
# dataset contains [n.pheno] sets of [n.reps] phenotypes, each of which is
# associated with the interaction of [order] genotypes. The generative model
# is consistent for every phenotype in a set, giving rise to what we call
# "concordant epistasis". The associations are also designed to be largely 
# independent of marginal effects. A key is also generated to identify which 
# set of genotypes interact to produce a phenotype.

############################# packages ####################################

library(abind)
library(mvtnorm)

############################ functions ####################################

# Input: n, m
# Output: all binary vectors (can only take on value of 0 or 1) of length n that sum to m
comb.binary <- function(n, m) {
  ind <- combn(seq_len(n), m)
  ind <- t(ind) + (seq_len(ncol(ind)) - 1) * n
  res <- rep(0, nrow(ind) * n)
  res[ind] <- 1
  matrix(res, ncol = n, nrow = nrow(ind), byrow = TRUE)
}

# Input:  set of genotypes, coded additively [-1,0,1];
# Output: Design matrix of full linear epistasis model (intercept excluded)
additive2EpistasisDesign <- function(geno.a,rm.nzcols=T) {
  geno.d <- (geno.a == 0) - 0.5*(geno.a!=0)
  geno <- abind(geno.a,geno.d,along=3)
  binary <- NULL
  for (ORDER in 0:order) {binary <- rbind(binary,comb.binary(order,ORDER)+1)}
  X <- sapply(1:nrow(binary),function(x) apply(sapply(1:ncol(binary),function(y) geno[,y,binary[x,y]]),1,prod))
  if (rm.nzcols == T) {
    nzcols <- apply(X,2,function(x) var(x)>0)
    X <- X[, nzcols]
  }
  return(X)
}

############################ parameters ####################################

# total number of simulations
n.sims <- 50

#sample size
n <- 99

# number of phenotypes per set
n.reps <- 4

# number of phenotype sets
n.pheno <- 100

# number of genotypes
n.geno <- 15

# order of interaction to test for (i.e. number of interacting genotypes in a genomic association)
order <- 2

# number of possible combinations of genotypes of size [order]
n.tests <- choose(n.geno,order)

# noise: alter this to change signal/noise ratio
noise <- 75

######################### simulate datasets #################################

# genotypes are randomly sampled 3-category variables with equal samples in each category
geno.evendist <- rep(c(-1,0,1),n/3)

# initialized key to be used for each simulated phenotype
key.init <- matrix(0,n.tests,n.pheno)
rownames(key.init) <- apply(combn(n.geno,2),2,function(x) paste(x,collapse="/"))
colnames(key.init) <- paste("Phenotype",1:n.pheno)

# this object will hold all of the simulated data
sims <- vector(mode="list",length=n.sims)
names(sims) <- paste("Dataset",1:n.sims)

# set seed for easy recovery of randomly chosen genotype sets
set.seed(1)

for (s in 1:n.sims) {
  
  # each dataset entry contains a list of 3 entries: 
  # [[1]] genotypes, [[2]] list of phenotype sets, [[3]] key indicating condition true and false examples
  sims[[s]] <- vector(mode="list",length=3)
  names(sims[[s]]) <- c("Genotypes","Phenotypes","Key")
  
  # genotypes are 3-category variables with equal samples per category
  genotypes <- sapply(1:n.geno,function(x) sample(geno.evendist))
  sims[[s]][[1]] <- genotypes
  
  # entry for phenotype sets is a list since each set itself is a matrix 
  sims[[s]][[2]] <- vector(mode="list",length=n.pheno)
  names(sims[[s]][[2]]) <- paste("Phenotype",1:n.pheno)
  
  # reset key for this dataset to the initialized key
  key <- key.init
  
  for (p in 1:n.pheno) {
    
    # determine which genotype interaction will be used to generate phenotypes
    inter.idx <- sample(1:n.tests,1)
    
    # alter key to represent this new condition true example
    key[inter.idx,p] <- 1
    
    # define set of interacting genotypes using randomly chosen index
    geno.a <- genotypes[,as.numeric(strsplit(rownames(key)[inter.idx],split="/")[[1]])]
    
    # generate full linear epistasis model terms
    X <- additive2EpistasisDesign(geno.a)
    
    # signal is evenly distributed across parameters
    beta.stoch <- rep(1,ncol(X))
    sigma <- diag(noise,n.reps)
    sims[[s]][[2]][[p]] <- sapply(1:n.reps,function(x) X%*%beta.stoch) + rmvnorm(n,rep(0,n.reps),sigma)
    
  }
  
  # store the key for this dataset
  sims[[s]][[3]] <- key
  
}

setwd("C:/Users/Mark/Desktop/fomac_sim")

save(file="data_sim_small.RData",sims)

