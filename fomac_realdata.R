rm(list=ls())

args <- commandArgs(TRUE)
task.id <- as.numeric(args[1])
n.tasks <- as.numeric(args[2])
prefix <- args[3]


dir.create(sprintf("/zenodotus/dat01/mezeylab_scratch/mcs2008/%s/%s_l=%d_maf=%s",dir,dir,num.marginal,mafcutoffcode))

set.seed(task.id*3)

setwd("ftest")
setwd("F:/muther")

############################# packages ####################################

library(mvtnorm)
library(Matrix)
#library(MVB)
library(abind)
#library(lrgpr)
#library(gtools)
library(randomForest)

############################ functions ####################################

# Input: set of genotypes to be tested for interaction, set of phenotypes
# Output: FoMaC p-value, p-values from many other comparable methods
fomacFull <- function(genotype.mat,Y,expected.median=1) {
  
  n <- nrow(genotype.mat)
  order <- ncol(genotype.mat)
  n.reps <- ncol(Y)
  
  maf <- function(genotypes) {
    maf <- sum(genotypes + 1)/(2*length(genotypes))
    if (maf > 0.5) {maf <- 1-maf}
    return(maf)
  }
  
  geno.d <- apply(genotype.mat,2,function(x) as.integer(x == 1))
  geno.c <- apply(genotype.mat+1, 1, function (x) sum(x*3^(0:(order-1))))+1
  X <- diag(3^order)[geno.c, ]
  nzcols <- colSums(X)>0
  N <- sum(nzcols)-1
  X <- X[, nzcols][, -1]
  
  # perform linear regression
  itxx <- chol2inv(chol(crossprod(X)))
  #itxx <- solve(t(X)%*%X)
  beta.hat <- itxx%*%crossprod(X,Y)
  var.hat <- diag(crossprod(Y-X%*%beta.hat))/n
  var.beta.hat <- diag(itxx)%*%t(var.hat)
  beta.hat.adj <- beta.hat/sqrt(var.beta.hat)
  var.samp <- sum(apply(beta.hat.adj,1,function(x) sum((x - mean(x))^2)))
  beta.means <- sum(rowMeans(beta.hat.adj)^2)*n.reps
  fstat <- (beta.means/N)/(var.samp/(N*(n.reps-1)))
  pval.complementary <- pf(fstat/expected.median,df1=N,df2=N*(n.reps-1),lower.tail=F)
  
  #calculate lm complementary coding p values
  pheno <- Y[,1]
  datanull <- rep(1,n)
  datatst <- X
  L0 <- lm(pheno~.,data=data.frame(datanull))
  L1 <- lm(pheno~.,data=data.frame(datatst))
  pval.lm.comp.1 <-anova(L0,L1)$"Pr(>F)"[2]
  
  pheno <- Y[,2]
  datanull <- rep(1,n)
  datatst <- X
  L0 <- lm(pheno~.,data=data.frame(datanull))
  L1 <- lm(pheno~.,data=data.frame(datatst))
  pval.lm.comp.2 <-anova(L0,L1)$"Pr(>F)"[2]
  
  pheno <- Y[,3]
  datanull <- rep(1,n)
  datatst <- X
  L0 <- lm(pheno~.,data=data.frame(datanull))
  L1 <- lm(pheno~.,data=data.frame(datatst))
  pval.lm.comp.3 <-anova(L0,L1)$"Pr(>F)"[2]
  
  # combine p-values
  pval.lm.comp.fisher <- pchisq(-2*(log(pval.lm.comp.1)+log(pval.lm.comp.2)+log(pval.lm.comp.3)),df=2*n.reps,lower.tail=F)
  
  #calculate lm marginal p values
  pvals.marginal <- NULL
  for (geno.idx in 1:ncol(genotype.mat)) {
    datanull <- rep(1,n)
    datatst <- cbind(genotype.mat[,geno.idx],geno.d[,geno.idx])
    for (pheno.idx in 1:ncol(Y)) {
      pheno <- Y[,pheno.idx]
      L0 <- lm(pheno~.,data=data.frame(datanull))
      L1 <- lm(pheno~.,data=data.frame(datatst))
      pvals.marginal <- c(pvals.marginal,anova(L0,L1)$"Pr(>F)"[2])
      names(pvals.marginal[length(pvals.marginal)]) <- sprintf("geno %d, pheno %d",geno.idx,pheno.idx)
    }
  }
  
  ##############################################################################################################
  # calculate additive-by-additive-by-additive epistasis p values (PLINK, Fitzpatrick 2015, Lareau 2016)
  geno.a <- genotype.mat + 1
  X <- as.matrix(apply(geno.a,1,prod))
  #nzcols <- apply(X,2,function(x) var(x)>0)
  #X <- X[, nzcols]
  
  pheno <- Y[,1]
  datanull <- rep(1,n)
  datatst <- X
  L0 <- lm(pheno~.,data=data.frame(datanull))
  L1 <- lm(pheno~.,data=data.frame(datatst))
  pval.lm.epiadd.1 <-anova(L0,L1)$"Pr(>F)"[2]
  
  pheno <- Y[,2]
  datanull <- rep(1,n)
  datatst <- X
  L0 <- lm(pheno~.,data=data.frame(datanull))
  L1 <- lm(pheno~.,data=data.frame(datatst))
  pval.lm.epiadd.2 <-anova(L0,L1)$"Pr(>F)"[2]
  
  pheno <- Y[,3]
  datanull <- rep(1,n)
  datatst <- X
  L0 <- lm(pheno~.,data=data.frame(datanull))
  L1 <- lm(pheno~.,data=data.frame(datatst))
  pval.lm.epiadd.3 <-anova(L0,L1)$"Pr(>F)"[2]
  
  # combine p-values
  pval.lm.epiadd.fisher <- pchisq(-2*(log(pval.lm.epiadd.1)+log(pval.lm.epiadd.2)+log(pval.lm.epiadd.3)),df=2*n.reps,lower.tail=F)
  ######################################################################################################
  
  ######################################################################################################
  # calculate full lm epistasis model p values (additive coding [-1,0,1] dominant coding [-0.5,1,-0.5]) (Spurgeon - similar to Fish 2016)
  geno.a <- genotype.mat
  geno.d <- apply(geno.a,2,function(x) as.integer(x == 0)) - 0.5*apply(geno.a,2,function(x) as.integer(x != 0))
  geno <- abind(geno.a,geno.d,along=3)
  
  comb.binary <- function(n, m) {
    ind <- combn(seq_len(n), m)
    ind <- t(ind) + (seq_len(ncol(ind)) - 1) * n
    res <- rep(0, nrow(ind) * n)
    res[ind] <- 1
    matrix(res, ncol = n, nrow = nrow(ind), byrow = TRUE)
  }
  
  binary <- NULL
  for (ORDER in 0:order) {binary <- rbind(binary,comb.binary(order,ORDER)+1)}
  
  X <- sapply(1:nrow(binary),function(x) apply(sapply(1:ncol(binary),function(y) geno[,y,binary[x,y]]),1,prod))
  
  nzcols <- apply(X,2,function(x) var(x)>0)
  X <- X[, nzcols]
  N <- ncol(X)
  
  # calculate ftest p-value for spurgeon coding
  itxx <- chol2inv(chol(crossprod(X)))
  #itxx <- solve(t(X)%*%X)
  beta.hat <- itxx%*%crossprod(X,Y)
  var.hat <- diag(crossprod(Y-X%*%beta.hat))/n
  var.beta.hat <- diag(itxx)%*%t(var.hat)
  beta.hat.adj <- beta.hat/sqrt(var.beta.hat)
  var.samp <- sum(apply(beta.hat.adj,1,function(x) sum((x - mean(x))^2)))
  beta.means <- sum(rowMeans(beta.hat.adj)^2)*n.reps
  fstat <- (beta.means/N)/(var.samp/(N*(n.reps-1)))
  pval.spurgeon <- pf(fstat/expected.median,df1=N,df2=N*(n.reps-1),lower.tail=F)
  
  pheno <- Y[,1]
  datanull <- rep(1,n)
  datatst <- X
  L0 <- lm(pheno~.,data=data.frame(datanull))
  L1 <- lm(pheno~.,data=data.frame(datatst))
  pval.lm.fish.1 <-anova(L0,L1)$"Pr(>F)"[2]
  
  pheno <- Y[,2]
  datanull <- rep(1,n)
  datatst <- X
  L0 <- lm(pheno~.,data=data.frame(datanull))
  L1 <- lm(pheno~.,data=data.frame(datatst))
  pval.lm.fish.2 <-anova(L0,L1)$"Pr(>F)"[2]
  
  pheno <- Y[,3]
  datanull <- rep(1,n)
  datatst <- X
  L0 <- lm(pheno~.,data=data.frame(datanull))
  L1 <- lm(pheno~.,data=data.frame(datatst))
  pval.lm.fish.3 <-anova(L0,L1)$"Pr(>F)"[2]
  
  # combine p-values
  pval.lm.fish.fisher <- pchisq(-2*(log(pval.lm.fish.1)+log(pval.lm.fish.2)+log(pval.lm.fish.3)),df=2*n.reps,lower.tail=F)
  #############################################################################################
  
  # finally, a readout of the number of table entries and min number of samples in an entry
  output <- c(colnames(genotype.mat),
              pval.spurgeon,pval.complementary,
              pvals.marginal,
              pval.lm.comp.1,pval.lm.comp.2,pval.lm.comp.3,pval.lm.comp.fisher,
              pval.lm.epiadd.1,pval.lm.epiadd.2,pval.lm.epiadd.3,pval.lm.epiadd.fisher,
              pval.lm.fish.1,pval.lm.fish.2,pval.lm.fish.3,pval.lm.fish.fisher,
              #pval.lm.epistnd.1,pval.lm.epistnd.2,pval.epistnd.3,
              min(table(geno.c)),
              length(unique(geno.c)),
              expected.median,
              apply(genotype.mat,2,maf))
  
  # add geno.c to the output
  #output <- c(output,geno.c)
  #names(output) <- "Pheno ID","Median","Geno ID 1","Geno ID 2",
  
  return(output)
}

# Input: set of genotypes to be tested for interaction, set of phenotypes
# Output: FoMaC p-value
fomacFast <- function(genotype.mat,Y,expected.median=1,output="pval",model="additivedominant") {
  
  n <- nrow(genotype.mat)
  order <- ncol(genotype.mat)
  
  if (model == "complementary") {
    
    # transform set of genotypes into design matrix of binary variables
    geno.c <- apply(genotype.mat+1, 1, function (x) sum(x*3^(0:(order-1))))+1
    X <- diag(3^order)[geno.c, ]
    
    # remove terms with no samples
    nzcols <- colSums(X)>0
    X <- X[, nzcols][, -1]
    N <- ncol(X)
    
  } else if (model == "additivedominant") {
    
    geno.a <- genotype.mat
    
    # generate full linear epistasis model terms
    X <- additive2EpistasisDesign(geno.a)
    
    N <- ncol(X)
    
    # remove terms with zero variance
    nzcols <- apply(X,2,function(x) var(x)>0)
    X <- X[, nzcols]
    
    # tack on intercept term
    X <- cbind(rep(1,n),X)
  }
  
  # fit linear regressions & calculate statistics
  itxx <- chol2inv(chol(crossprod(X)))
  beta.hat <- itxx%*%crossprod(X,Y)
  var.hat <- diag(crossprod(Y-X%*%beta.hat))/n
  var.beta.hat <- diag(itxx)%*%t(var.hat)
  beta.hat.adj <- (beta.hat/sqrt(var.beta.hat))[-1,]
  var.samp <- sum(apply(beta.hat.adj,1,function(x) sum((x - mean(x))^2)))
  beta.means <- sum(rowMeans(beta.hat.adj)^2)*n.reps
  fstat <- (beta.means/N)/(var.samp/(N*(n.reps-1)))
  
  if (output == "pval") {pval <- pf(fstat/expected.median,df1=N,df2=N*(n.reps-1),lower.tail=F); return(pval)}
  if (output == "fstat") {return(list(fstat,N))}
}

# Input: one phenotype, all genotypes
# Output: interaction scores for every possible genotype pair
forestInteractions <- function(phenotype,genotypes,NTREE=10000,method="branch.interactions") {
  
  order <- 2
  n <- nrow(genotypes)
  
  if (method == "branch.interactions") {
    
    extractInteractions <- function(tree) {
      
      tree.bin <- rbind(tree[,-2],tree[,-1])
      
      # identify all candidate nodes (non-terminal and non-root)
      c.indices <- which(tree[,3] !=0)[-1]
      
      # follow paths back to root node for all candidate nodes
      interactions <- NULL
      for (c.idx in c.indices) {
        interacting.vars <- NULL
        c.var <- tree[c.idx,3]
        current.node <- c.idx
        while (current.node !=1) {
          parent.node <- as.numeric(rownames(tree.bin)[which(tree.bin[,1] ==current.node)])
          interacting.vars <- c(interacting.vars,tree[parent.node,3])
          current.node <- parent.node
        }
        interactions <- cbind(interactions,rbind(rep(c.var,length(interacting.vars)),interacting.vars))
      }
      
      return(t(interactions))
    }
    
    predictors <- data.frame(genotypes)
    for (m in 1:ncol(predictors)) {predictors[,m] <- as.factor(predictors[,m])}
    
    #start.time <- proc.time()[3]
    
    rfout <- randomForest(x=predictors,y=phenotype,ntree=NTREE,maxnodes=16)
    
    # extract interactions from all trees, merging all results
    inter.all <- NULL
    for (tree.idx in 1:NTREE) {
      inter.all <- rbind(inter.all,extractInteractions(getTree(rfout,k=tree.idx)))
    }
    
    # remove self interactions
    inter.all <- inter.all[apply(inter.all,1,function(x) x[1]!=x[2]),]
    
    # add on a phantom interaction for each possible interaction (prevents errors - will be removed shortly)
    inter.all <- rbind(inter.all,t(combn(1:ncol(genotypes),2)))
    
    # sum up all unordered counts (and remove phantom interaction that was added in previous step)
    counts <- table(apply(inter.all,1,function(x) paste(sort(x),collapse="/"))) - 1
    
    ordering <- apply(combn(1:ncol(genotypes),2),2,function(x) paste(x,collapse="/"))
    
    output <- sapply(ordering,function(x) counts[which(names(counts) == x)])
    names(output) <- ordering
    
    return(output)
  } else if (method == "variable.importance") {
    start.time <- proc.time()[3]
    
    comb <- combn(1:ncol(genotypes),2)
    
    predictors <- data.frame(rep(1,n))
    for (i in 1:ncol(comb)) {
      genotype.mat <- genotypes[,comb[,i]]
      geno.c <- apply(genotype.mat+1, 1, function (x) sum(x*3^(0:(order-1))))+1
      predictors <- data.frame(predictors,as.factor(geno.c))
    }
    predictors <- predictors[,-1]
    
    colnames(predictors) <- apply(comb,2,function(x) paste(x,collapse="/"))
    
    rfout <- randomForest(x=predictors,y=phenotype,ntree=100,maxnodes=2,mtry=ncol(predictors))
    
    print(proc.time()[3]-start.time)
    return(rfout$interactions)
  }
}

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

# Input: one genotype variable
# Output: minor allele frequency
maf <- function(genotypes) {
  maf <- sum(genotypes + 1)/(2*length(genotypes))
  if (maf > 0.5) {maf <- 1-maf}
  return(maf)
}

# Input: one genotype variable
# Output: minimum zygotic fraction
mzf <- function(genotypes,categories) {
  min(sum(genotypes == categories[1]),
      sum(genotypes == categories[2]),
      sum(genotypes == categories[3]))/length(genotypes)
}

# Input: p-values
# Output: plots a QQ plot of the -log10(p-values) vs theoretical expectation
qqUnif <- function(pvals,key=NULL,title=NULL,cicolor="blue",ylim=c(0,-log10(min(pvals)))) {
  
  COL <- rep(rgb(0,0,0,0),length(pvals))
  if (!is.null(key)) {
    for (i in which(key == 1)) {COL[i] <- "red"}
  }
  
  n.tests <- length(pvals)
  
  logp.exp <- -log10(1:n.tests/n.tests)
  
  # Confidence intervals ala
  # Casella & Berger, 2002, 
  # 2nd edition, pg 230, Duxbury
  cint.upper <- sapply(1:n.tests, function (x) { qbeta(0.90, x, n.tests - x + 1) })
  cint.lower <- sapply(1:n.tests, function (x) { qbeta(0.10, x, n.tests - x + 1) })
  
  # Make QQ Plot
  logp.obs <- -log10(sort(pvals))
  plot(logp.exp, logp.obs, type="n", pch=19, cex=0.5,
       xlab='-log10(expected p-values)', ylab='-log10(observed p-values)',
       xlim=c(0, max(logp.exp)), ylim=ylim, main=title)
  
  # Add x = y diagonal line
  abline(0, 1, col='gray', lty=2)
  
  # Add lines for confidence intervals
  lines(logp.exp, -log10(cint.upper), lty=2, col=cicolor)
  lines(logp.exp, -log10(cint.lower), lty=2, col=cicolor)
  
  points(logp.exp, logp.obs, pch=19, cex=0.5)
  
  points(logp.exp, logp.obs, pch=19, cex=0.5, col=COL)
}

############################# housekeeping ################################

task.id <- 1
n.tasks <- 100
prefix <- "muther_fatlclskin"

order <- 2
num.marginal <- 100
minmzf <- 0.2
minmzfcode <- "02"
cor.max <- 0.2

# data
load(sprintf("%s.RData",prefix))

# transform genotypes to [-1,0,1]
genotypes <- genotypes - 1

rm(covars)

n <- nrow(genotypes)
n.reps <- length(phenotypes)

# used to determine optimal minmzf
if (1 == 0) {
genotypes.temp <- genotypes[,which(zf > 0.2 & zf < 0.4)]
count <- 1
count.max <- 1000
store <- rep(NA,count.max)
while(count <= count.max) {
  cor.temp <- 1
  while(cor.temp > 0.2) {
    geno.a <- genotypes.temp[,sample(1:ncol(genotypes.temp),2)]
    cor.temp <- cor(geno.a[,1],geno.a[,2])
  }
  geno.c <- apply(geno.a, 1, function (x) sum(x*3^(0:(order-1))))+1
  store[count] <- min(table(geno.c))
  count <- count + 1
}
#rm(count,count.max,genotypes,temp,store,geno.a,cor.temp,geno.c)
}

# get rid of genotypes with mzf < minmzf
genotypes <- genotypes[,which(apply(genotypes,2,function(x) mzf(x,-1:1)) > minmzf)]

# determine which phenotypes this job should test
if (task.id <= (n.tasks-2)) {pheno.indices <- matrix(1:((n.tasks-2)*floor((ncol(phenotypes[[1]])-1)/(n.tasks-2))),floor((ncol(phenotypes[[1]])-1)/(n.tasks-2)),n.tasks-2)[,task.id]}
if (task.id == n.tasks-1) {pheno.indices <- ((n.tasks-2)*floor((ncol(phenotypes[[1]])-1)/(n.tasks-2))+1):(ncol(phenotypes[[1]])-1)}
if (task.id == n.tasks) {pheno.indices <- ncol(phenotypes[[1]])}

# detects unfinished analysis and picks up where it left off
if (1 == 1) {
  if (file.exists(sprintf("/zenodotus/dat01/mezeylab_scratch/mcs2008/%s/%s_l=%d_maf=%s/progress_task%d.RData",dir,dir,num.marginal,mafcutoffcode,task.id))) {
    load(sprintf("/zenodotus/dat01/mezeylab_scratch/mcs2008/%s/%s_l=%d_maf=%s/progress_task%d.RData",dir,dir,num.marginal,mafcutoffcode,task.id))
    pheno.indices <- pheno.indices[-which(pheno.indices %in% progress)]
  } else {
    progress <- NULL
  }
}



if (length(pheno.indices) > 0) {
  
  print(pheno.indices)
  results <- NULL
  pvals <- pvals.comp <- NULL
  scale.factor <- NULL
  n.hits <- n.tests <- rep(0,length(pheno.indices))
  names(n.hits) <- names(n.tests) <- colnames(phenotypes[[1]])[pheno.indices]
  start.time <- proc.time()[3]
  for (pheno.idx in pheno.indices) {
    
    ##################################################
    # Prepare phenotypes
    ##################################################
    
    pheno.idx <- pheno.indices[1]
    pheno.id <- colnames(phenotypes[[1]])[pheno.idx]
    
    # normalize variance
    Y <- sapply(1:n.reps,function(x) phenotypes[[x]][,pheno.idx])
    
    if (exists("covars")) {
      if (is.list(covars)) {
        Y <- sapply(1:n.reps,function(y) {
          Ytemp <- Y[,y]
          lm(Ytemp~.,data=data.frame(Ytemp,covars[[y]]))$residuals
        })
      } else {
        Y <- apply(Y,2,function(y) lm(y~.,data=data.frame(y,covars))$residuals)
      }
    }

    Y <- apply(Y,2,function(x) x/sd(x))
    colnames(Y) <- rep(colnames(phenotypes[[1]])[pheno.idx],ncol(Y))
    
    ##################################################
    # Find out how much the f-statistic will be inflated due to inter-phenotype correlations
    ##################################################
    
    n.tests <- 0
    n.tests.max <- 1000
    
    fomac.temp <- N.temp <- rep(NA,n.tests.max)
    geno.indices.temp <- matrix(NA,n.tests.max,order)
    
    while (n.tests < n.tests.max) {
      geno.idx <- sample(1:ncol(genotypes),order)
      genotypes.temp <- genotypes[,geno.idx]
      #print(sprintf("n.tests = %d     genotypes = %s/%s     geno.idx = %d/%d",
      #              n.tests,colnames(genotypes.temp)[1],colnames(genotypes.temp)[2],geno.idx[1],geno.idx[2]))
      #geno.c <- apply(genotypes.temp, 1, function (x) sum(x*c(1, 3)))+1
      #if (length(unique(geno.c)) == 9 & min(table(geno.c)) >= 4) {
      if (max(abs(cor(genotypes.temp)-diag(ncol(genotypes.temp)))) < cor.max) {
        n.tests <- n.tests + 1
        #getFstatFast.out <- getFstatFast(geno.c,Y,expected.median,output="fstat",model="additivedominant")
        fomacFast.out <- fomacFast(genotypes.temp,Y,expected.median,output="fstat",model="additivedominant")
        fomac.temp[n.tests] <- fomacFast.out[[1]]
        N.temp[n.tests] <- fomacFast.out[[2]]
        geno.indices.temp[n.tests,] <- geno.idx
      }
    }
    expected.median <- median(fomac.temp)/median(rf(1000000,mean(N.temp),mean(N.temp)*(n.reps-1)))
    #print(sprintf("expected median: %f",expected.median))
    
    ##################################################
    # Determine which genotype interactions to test
    ##################################################

    # marginal p-values
    if (1 == 1) {
      # separate datasets (marginal LRT)
      pvals.marginal <- matrix(NA,ncol(genotypes),length(phenotypes))
      n <- nrow(Y)
      
      for (data.idx in 1:length(phenotypes)) {
        pheno <- Y[,data.idx]
        print(sprintf("calculating marginal p-value: dataset %d",data.idx))
        pvals.marginal[,data.idx] <- as.vector(glmApply(pheno ~ SNP,features=genotypes)[[1]])
        #pvals.marginal[,data.idx] <- apply(genotypes,2,function(x) getPvalLinear(pheno,rep(1,n),cbind(rep(1,n),x,(x==1))))
        #pheno <- phenotypes[[data.idx]][,pheno.idx]
        #pvals.marginal[,data.idx] <- apply(genotypes,2,function(x) getPvalLinear(pheno,covars[[data.idx]],cbind(covars[[data.idx]],x,(x==1))))
      }
      
      pvals.marginal.vec <- as.vector(pvals.marginal)
      names(pvals.marginal.vec) <- rep(colnames(genotypes),ncol(pvals.marginal))
      pvals.marginal.vec <- sort(pvals.marginal.vec)
      
      genotype.ids <- names(pvals.marginal.vec)[1]
      count <- 1
      while(length(genotype.ids) < num.marginal) {
        count <- count + 1
        genotype.ids.temp <- c(genotype.ids,names(pvals.marginal.vec)[count])
        if (max(abs(cor(genotypes[,genotype.ids.temp])-diag(length(genotype.ids.temp)))) < cor.max) {genotype.ids <- genotype.ids.temp}
      }
      
      #n.unique <- sapply(1:length(pvals.marginal.vec),function(x) length(unique(names(pvals.marginal.vec[1:x]))))
      #genotype.ids <- unique(names(pvals.marginal.vec)[1:which(n.unique == 50)])
      genotype.indices <- sapply(genotype.ids,function(x) which(colnames(genotypes) == x))
      
    }
    
    # marginal FoMaC p-values
    if (1 == 0) {
      pvals.marginal <- apply(genotypes,2,function(x) fomacFast(as.matrix(x),Y,expected.median=1,output="pval",model="marginal"))
      genotype.indices <- order(pvals.marginal)[1:num.marginal]
      genotype.ids <- colnames(genotypes)[genotype.indices]
    }
    
    # random forest interaction scores
    if (1 == 0) {
      scores.forest <- rep(0,choose(ncol(genotypes),order))
      for (i in 1:n.reps) {
        scores.forest <- scores.forest + forestInteractions(Y[,i],genotypes,
                                                            NTREE=round(1000/n.reps),
                                                            method="branch.interactions")
      }
    }
    
    # randomly chosen
    if (1 == 0)
      genotype.ids.all <- sample(colnames(genotypes))
      genotype.ids <- genotype.ids.all[1]
      count <- 1
      while(length(genotype.ids) < num.marginal) {
        count <- count + 1
        genotype.ids.temp <- c(genotype.ids,genotype.ids.all[count])
        if (max(abs(cor(genotypes[,genotype.ids.temp])-diag(length(genotype.ids.temp)))) < cor.max) {genotype.ids <- genotype.ids.temp}
      }
      #n.unique <- sapply(1:length(pvals.marginal.vec),function(x) length(unique(names(pvals.marginal.vec[1:x]))))
      #genotype.ids <- unique(names(pvals.marginal.vec)[1:which(n.unique == 50)])
      genotype.indices <- sapply(genotype.ids,function(x) which(colnames(genotypes) == x))
    }
    
    #############################################
    geno.comb <- combn(genotype.indices,order)
    pvals.temp <- rep(NA,ncol(geno.comb))
    for (i in 1:length(pvals.temp)) {
      genotypes.temp <- genotypes[,geno.comb[,i]]
      pvals.temp[i] <- fomacFast(genotypes.temp,Y,expected.median,model="additivedominant")
    }
    names(pvals.temp) <- rep(pheno.id,length(pvals.temp))
    pvals <- c(pvals,pvals.temp)
    
    scale.factor.temp <- rep(expected.median,length(pvals.temp))
    names(scale.factor.temp) <- rep(pheno.id,length(pvals.temp))
    scale.factor <- c(scale.factor,scale.factor.temp)
    
    results.pheno <- NULL
    if (sum(pvals.temp < 0.05/sum(!is.na(pvals.temp)),na.rm=T) > 0) {
      which(pvals.temp < 0.05/sum(!is.na(pvals.temp)),arr.ind=T)
      
      hits <- as.matrix(geno.comb[,which(pvals.temp < 0.05/sum(!is.na(pvals.temp)))],order,sum(pvals.temp < 0.05/sum(!is.na(pvals.temp)),na.rm=T))
      
      for (hit.idx in 1:ncol(hits)) {
        genotypes.temp <- genotypes[,hits[,hit.idx]]
        
        # find annotation info for gene and snps
        distances <- rep(NA,ncol(genotypes.temp))
        for (peanut in 1:ncol(genotypes.temp)) {
          
          chr.idx.snp <- as.character(snp.info$chromosome[which(snp.info$id == colnames(genotypes.temp)[peanut])])
          
          if (identical(chr.idx,chr.idx.snp)) {
            distances[peanut] <- abs(snp.info$position[which(snp.info$id == colnames(genotypes.temp)[peanut])]-gene.info[which(gene.info[,which(colnames(gene.info)=="probe")]==colnames(phenotypes[[1]])[pheno.idx]),which(colnames(gene.info)=="start")])
          }
        }
        
        
        results.pheno <- rbind(results.pheno,c(pheno.id,getFstatFull(genotypes.temp,Y,expected.median),distances,0.05/sum(!is.na(pvals.temp))))
      }
      
      colnames(results.pheno) <- c("id:phenotype",paste("id:genotype",1:order,sep=""),"jointPval:FstatSpurgeon","jointPval:FstatComplementary",
                                   paste(paste("marginalPval:","geno",as.vector(sapply(1:order,function(x) rep(x,n.reps))),"data",sep=""),rep(1:n.reps,order),sep=""),
                                   paste("complementaryPval: data",1:n.reps,sep=""),"complementaryPvalCombined",
                                   paste("additiveepistasisPval:data",1:n.reps,sep=""),"additiveepistasisPvalCombined",
                                   paste("spurgeonepistasisPval:data",1:n.reps,sep=""),"spurgeonepistasisPvalCombined",
                                   "min samples","N","scale factor",paste("maf:geno",1:order,sep=""),
                                   paste("dist2TSS:geno",1:order,sep=""),"bonferroni")
      
      results.pheno <- results.pheno[order(as.numeric(results.pheno[,5])),]
      #results <- rbind(results,results.pheno)
      print(results.pheno)
      results <- rbind(results,results.pheno)
      #save(file=sprintf("hits/hits_%d_%d_%s_%d.RData",task.id,pheno.idx,pheno.id,sum(pvals.temp < 0.05/sum(!is.na(pvals.temp)))),results.pheno)
      #save(file=sprintf("pvals/pvals_%d_%d_%s.RData",task.id,pheno.idx,pheno.id),pvals.temp)
    }
    n.hits[which(names(n.hits) == pheno.idx)] <- sum(pvals.temp < 0.05/sum(!is.na(pvals.temp)),na.rm=T)
    n.tests[which(names(n.tests) == pheno.idx)] <- sum(!is.na(pvals.temp))
    
    
    progress <- c(progress,pheno.idx)
    save(file=sprintf("/zenodotus/dat01/mezeylab_scratch/mcs2008/%s/%s_l=%d_maf=%s/progress_task%d.RData",dir,dir,num.marginal,mafcutoffcode,task.id),progress)
    print(sprintf("Phenotype %d/%d completed; %f hours elapsed",which(pheno.indices == pheno.idx),length(pheno.indices),(proc.time()[3]-start.time)/60^2))
    save(file=sprintf("/zenodotus/dat01/mezeylab_scratch/mcs2008/%s/%s_l=%d_maf=%s/%s_%s_results_%d.RData",dir,dir,num.marginal,mafcutoffcode,dir,prefix,task.id),results,pvals,scale.factor,n.hits,n.tests)
  }
  
} else {
  print("task already complete")
}
