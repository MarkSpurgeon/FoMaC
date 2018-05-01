# Version for array job on cluster. Differences between this and local script:
# (1) working directory
# (2) d is set to task.id
# 

########### command line arguments (unnecessary for local jobs) ###########
args <- commandArgs(TRUE)
task.id <- as.numeric(args[1])

############################# packages ####################################

library(abind)
library(aod)
library(randomForest)
library(flux)

############################ functions ####################################

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

# Input: inverse(t(X)*X), parameter estimates, residuals
# Output: p-value
waldPval <- function(ITXX,b,residuals,terms=NULL) {
  n <- length(residuals)
  if(is.null(terms)) {terms <- 2:length(b)}
  k <- length(terms)
  
  Sigma <- ((1/(n-k))*sum(residuals^2))*ITXX
  pval <- exp(pchisq(wald.test(Sigma,b,Terms=terms)$result[[1]][1],df=k,lower.tail=F,log.p=T))
  return(pval)
}

# Input: time in seconds
# Output: time in hours:minutes:seconds
formatTime <- function(time) {
  hours <- floor(elapsed.time/60^2)
  minutes <- floor((elapsed.time - hours*60^2)/60)
  seconds <- floor(time - minutes*60 - hours*60^2)
  return(paste(paste(hours,"h",sep=""),paste(minutes,"m",sep=""),paste(seconds,"s",sep=""),sep=" : "))
}

# Input: condition labels, scores (high scores representing higher confidence)
# Output: receiver operating characteristic (TPR & FPR at all thresholds)
simple_roc <- function(labels, scores){
  labels <- labels[order(scores, decreasing=TRUE)]
  data.frame(TPR=cumsum(labels)/sum(labels),FPR=cumsum(!labels)/sum(!labels))
}

# Input: condition labels, scores (high scores representing higher confidence)
# Output: precision and recall at all thresholds
simple_pr <- function(labels, scores) {
  labels <- labels[order(scores, decreasing=TRUE)]
  data.frame(precision=cumsum(labels)/1:length(labels),recall=cumsum(labels)/sum(labels))
}

# Input: merged FPR and TPR from multiple roc curves, sorted according to FPR
# Output: single roc curve, with TPR values averaged at each unique FPR
curveMean <- function(dv,iv) {
  breakpoints <- cbind(c(1,1+which(diff(iv)!=0)),c(c(1,1+which(diff(iv)!=0))[-1]-1,length(iv)))
  dv.mean <- apply(breakpoints,1,function(x) mean(dv[x[1]:x[2]]))
  return(dv.mean,cbind(unique(iv)))
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

########################### load data ######################################

#setwd("scratch/fomac_sim/round9")

load("data_sim_small.RData")

############################ parameters ####################################

# dataset index
d <- task.id

# order of interaction to test for (i.e. number of interacting genotypes in a genomic association)
order <- 2

# sample size
n <- nrow(datasets[[d]][[1]])

# number of phenotypes per set
n.reps <- ncol(datasets[[d]][[2]][[1]])

# which phenotype sets to analyze (default all sets)
pheno.indices <- 1:length(datasets[[d]][[2]])

# number of phenotype sets
n.pheno <- length(pheno.indices)

# number of genotypes
n.geno <- ncol(datasets[[d]][[1]])

# all possible genotype combinations of size [order]
comb <- combn(1:n.geno,order)

# number of association tests per phenotype set
n.tests <- ncol(comb)

# number of different methods being tested
n.methods <- 7

####################### test for epistasis ####################################
# preallocate storage for all results of association tests
pvals <-  array(NA,dim=c(n.tests,n.pheno,n.methods))
rownames(pvals) <- apply(combn(n.geno,2),2,function(x) paste(x,collapse="/"))

# set seed because of stochastic components, which are:
# (1) correction estimation for FoMac, (2) random forest interactions
set.seed(1)

start.time <- proc.time()[3]
for (p in pheno.indices) {
  
  # select pth set of phenotypes
  Y <- datasets[[d]][[2]][[p]]
  genotypes <- datasets[[d]][[1]]
  
  # preallocate method-specific storage
  pvals.fomac <- pvals.mvlm <- pvals.lm.comb <- pvals.plink.comb <- rep(NA,n.tests)
  pvals.lm <- pvals.plink <- matrix(NA,n.tests,n.reps)
  
  # calculate expected median (FoMaC correction heuristic) for current set of phenotypes
  count.em <- 1
  count.em.max <- 100
  fstats.temp <- N.temp <- rep(NA,count.em.max)
  while (count.em <= count.em.max) {
    
    # choose set of random genotypes with same distribution as sample genotypes
    geno.evendist <- rep(c(-1,0,1),n/3)
    
    # additive coding
    geno.a.em <- sapply(1:order,function(x) sample(geno.evendist))
    
    # generate full linear epistasis model terms
    X <- additive2EpistasisDesign(geno.a.em)
    
    # number of linear model terms being tested for significance
    N <- ncol(X)
    
    # tack on intercept term
    X <- cbind(rep(1,n),X)
    
    ### fit linear model & calculate statistics ###
    # calculate (t(X)*X)^-1
    itxx <- chol2inv(chol(crossprod(X)))
    
    # parameter estimates
    beta.hat <- itxx%*%crossprod(X,Y)
    
    # variance estimate
    var.hat <- diag(crossprod(Y-X%*%beta.hat))/n
    
    # estimate of parameter estimate variances
    var.beta.hat <- diag(itxx)%*%t(var.hat)
    
    # scale parameter estimates to transform into standard normal distributions
    beta.hat.adj <- (beta.hat/sqrt(var.beta.hat))[-1,]
    
    # fstat numerator - chi-square distributed with df=N
    var.samp <- sum(apply(beta.hat.adj,1,function(x) sum((x - mean(x))^2)))
    
    # fstat denominator - chi-square distributed with df=N*(n.reps-1)
    beta.means <- sum(rowMeans(beta.hat.adj)^2)*n.reps
    
    # fstat - ratio of chi-square distributions scaled by their dfs
    fstat <- (beta.means/N)/(var.samp/(N*(n.reps-1)))
    
    fstats.temp[count.em] <- fstat
    N.temp[count.em] <- N
    
    count.em <- count.em + 1
  }
  # expected median is ratio of empirically-derived / theoretical expectation
  em <- median(fstats.temp)/median(rf(100000,mean(N.temp),mean(N.temp)*(n.reps-1)))
  
  # test all possible genotype sets for association with phenotype
  for (i in 1:n.tests) {
    geno.a <- genotypes[,comb[,i]]
    
    # generate full linear epistasis model terms
    X <- additive2EpistasisDesign(geno.a)
    
    # number of linear model terms being tested for significance
    N <- ncol(X)
    
    # tack on intercept term
    X <- cbind(rep(1,n),X)
    
    # calculate (t(X)*X)^-1
    itxx <- chol2inv(chol(crossprod(X)))
    
    # parameter estimates
    beta.hat <- itxx%*%crossprod(X,Y)
    
    # variance estimate
    var.hat <- diag(crossprod(Y-X%*%beta.hat))/n
    
    # estimate of variances of parameter estimates
    var.beta.hat <- diag(itxx)%*%t(var.hat)
    
    # scale parameter estimates to transform into standard normal distributions
    beta.hat.adj <- (beta.hat/sqrt(var.beta.hat))[-1,]
    
    # numerator of fstat - chi-square distributed with df=N
    var.samp <- sum(apply(beta.hat.adj,1,function(x) sum((x - mean(x))^2)))
    
    # denominator of fstat - chi-square distributed with df=N*(n.reps-1)
    beta.means <- sum(rowMeans(beta.hat.adj)^2)*n.reps
    
    # fstat - ratio of chi-square distributions scaled by their dfs
    fstat <- (beta.means/N)/(var.samp/(N*(n.reps-1)))
    
    # p-value produced by testing against upper tail of f-distribution
    pvals.fomac[i] <- pf(fstat/em,df1=N,df2=N*(n.reps-1),lower.tail=F)
    
    # lm pvals for each phenotype in set
    resid <- Y - X%*%beta.hat
    pvals.lm[i,] <- sapply(1:n.reps,function(x) waldPval(itxx,beta.hat[,x],resid[,x]))
    
    # combine pvals from separate phenotypes into one p-value
    pvals.lm.comb[i] <- pchisq(-2*sum(log(pvals.lm[i,])),df=2*n.reps,lower.tail=F)
    
    # multivariate regression
    pvals.mvlm[i] <- summary(manova(Y ~ X))$stats[1,6]
    
    # PLINK p-values (not defined for epistasis order > 2)
    if (order == 2) {
      # design matrix contains intercept, additive genotype codings, and product of additive codings
      X.plink <- cbind(rep(1,n),geno.a,apply(geno.a+1,1,prod))
      
      # 
      itxx.plink <- chol2inv(chol(crossprod(X.plink)))
      beta.hat.plink <- itxx.plink%*%crossprod(X.plink,Y)
      resid.plink <- Y - X.plink%*%beta.hat.plink
      
      # PLINK pvals for each phenotype in set
      pvals.plink[i,] <- sapply(1:n.reps,function(x) waldPval(itxx.plink,beta.hat.plink[,x],resid.plink[,x],terms=4))
    }
    
    pvals.plink.comb[i] <- pchisq(-2*sum(log(pvals.plink[i,])),df=2*n.reps,lower.tail=F)
    
  }
  
  # run random forest interactions, summing scores from all phenotypes in set
  scores.forest <- rep(0,n.tests)
  for (i in 1:n.reps) {
    scores.forest <- scores.forest + forestInteractions(Y[,i],genotypes,
                                                        NTREE=round(500/n.reps),
                                                        method="branch.interactions")
  }
  
  # store p-values/scores
  pvals[,p,1] <- pvals.fomac
  pvals[,p,2] <- pvals.mvlm
  pvals[,p,3] <- pvals.lm[,1]
  pvals[,p,4] <- pvals.lm.comb
  pvals[,p,5] <- pvals.plink[,1]
  pvals[,p,6] <- pvals.plink.comb
  pvals[,p,7] <- -scores.forest
  
  # display progress
  elapsed.time <- proc.time()[3] - start.time
  elapsed.time.fmt <- formatTime(elapsed.time)
  print(sprintf("Phenotype %d/%d Completed; %s Elapsed",p,n.pheno,elapsed.time.fmt))
}

save(file=sprintf("fomac_sim_pvals_dataset%d.RData",d),pvals)

########################## Performance assessment ################################

# produce ROCs and calculate AUC for all methods
roc.d <- pr.d <- vector(mode="list",length=n.methods)
auc.roc.d <- auc.pr.d <- rep(NA,n.methods)
key.temp <- as.vector(datasets[[d]][[3]])
for (midx in 1:n.methods) {
  scores.temp <- -(as.vector(pvals[,,midx]))
  roc.temp <- simple_roc(key.temp,scores.temp)
  roc.d[[midx]] <- roc.temp
  auc.roc.d[midx] <- auc(roc.temp[,2],roc.temp[,1])
  
  pr.temp <- simple_pr(key.temp,scores.temp)
  pr.d[[midx]] <- pr.temp
  auc.pr.d[[midx]] <- auc(pr.temp[,2],pr.temp[,1])
}

save(file=sprintf("fomac_sim_performance_dataset%d.RData",d),roc.d,auc.roc.d,pr.d,auc.pr.d)

if (task.id == 1) {
  
  n.datasets <- length(datasets)
  while (length(grep("fomac_sim_performance_dataset",list.files())) < n.datasets) {
    Sys.sleep(10)
  }
  
  roc.merge <- pr.merge <- vector(mode="list",length=n.methods)
  auc.roc.all <- auc.pr.all <- NULL
  for (d in 1:n.datasets) {
    load(sprintf("fomac_sim_performance_dataset%d.RData",d))
    for (midx in 1:n.methods) {
      roc.merge[[midx]] <- rbind(roc.merge[[midx]],roc.d[[midx]])
      auc.roc.all <- rbind(auc.roc.all,auc.roc.d)
      
      pr.merge[[midx]] <- rbind(pr.merge[[midx]],pr.d[[midx]])
      auc.pr.all <- rbind(auc.pr.all,auc.pr.d)
    }
  }
  
  roc.means <- pr.means
  for (midx in 1:n.methods) {
    roc.merge[[midx]] <- roc.merge[[midx]][order(roc.merge[[midx]][,2]),]
    roc.means[[midx]] <- curveMean(roc.merge[[midx]][,1],roc.merge[[midx]][,2])
    
    pr.merge[[midx]] <- pr.merge[[midx]][order(pr.merge[[midx]][,2]),]
    pr.means[[midx]] <- curveMean(pr.merge[[midx]][,1],pr.merge[[midx]][,2])
  }
  
  save(file="fomac_sim_performance_ALL.RData",roc.means,auc.roc.all,pr.means,auc.pr.all)

  
  col.key <- c("darkgreen","purple","red","yellow","orange","blue","cyan")
  leg.key <- c("FoMaC","MVar LM","UVar LM","Cmbn LM","PLINK","Cmbn PLINK","RF")
  
  png("performance_roc+auc.png",height=2000,width=1300,pointsize=40)
  layout(matrix(1:2,2,1))
  
  par(font.lab=2)
  method.indices <- c(1,2,3,4,5,6,7)
  plot(0,0,type="n",xlim=c(0,1),ylim=c(0,1),cex.lab=1.3,xlab="FPR",ylab="TPR")
  grid(NULL,NULL)
  lines(0:1,0:1,col="gray",lty="dotted")
  #plot(0,0,type="n",xlim=c(0,0.1),ylim=c(0,0.6),xlab="False Positive Rate",ylab="True Positive Rate")
  par(font=2)
  legend("bottomright",legend=leg.key[method.indices],col=col.key[method.indices],cex=0.75,lty="solid",lwd=5)
  par(font=1)
  for (midx in method.indices) {
    lines(roc.means[[midx]][,2],roc.means[[midx]][,1],col=col.key[midx])
  }
  
  auc.all.subset <- auc.roc.all[,method.indices]
  YLIM <- range(auc.all.subset)
  plot(0,0,xlim=c(0.5,length(method.indices)+0.5),ylim=YLIM,type="n",
       cex.lab=1.3,ylab="ROC AUC",xlab=NA,xaxt="n")
  axis(side=1,at=1:length(method.indices),tick=T,labels=F)
  grid(NULL,NULL)
  x <- boxplot(auc.all.subset,col=col.key[method.indices],outline=F,add=T,xaxt="n")
  #axis(side=2,labels=)
  # need to manually adjust y here so that it sits right below the x axis
  text(cex=1, x=1:length(method.indices)+0.2, y=YLIM[1]-.05, leg.key[method.indices],
       xpd=TRUE, srt=45, pos=2)
  par(font.lab=1)
  
  dev.off()
  
  png("performance_pr+auc.png",height=2000,width=1300,pointsize=40)
  layout(matrix(1:2,2,1))
  par(font.lab=2)
  method.indices <- c(1,2,3,4,5,6,7)
  YLIM <- c(0,max(sapply(1:length(method.indices),function(x) max(pr.means[[x]][,1]))))
  plot(0,0,type="n",xlim=c(0,1),ylim=YLIM,cex.lab=1.3,xlab="Recall",ylab="Precision")
  grid(NULL,NULL)
  #plot(0,0,type="n",xlim=c(0,0.1),ylim=c(0,0.6),xlab="False Positive Rate",ylab="True Positive Rate")
  par(font=2)
  legend("topright",legend=leg.key[method.indices],col=col.key[method.indices],cex=0.75,lty="solid",lwd=5)
  par(font=1)
  for (midx in method.indices) {
    lines(pr.means[[midx]][,2],pr.means[[midx]][,1],col=col.key[midx])
  }
  
  auc.all.subset <- auc.pr.all[,method.indices]
  YLIM <- range(auc.all.subset)
  plot(0,0,xlim=c(0.5,length(method.indices)+0.5),ylim=YLIM,type="n",
       cex.lab=1.3,ylab="PR AUC",xlab=NA,xaxt="n")
  axis(side=1,at=1:length(method.indices),tick=T,labels=F)
  grid(NULL,NULL)
  x <- boxplot(auc.all.subset,col=col.key[method.indices],outline=F,add=T,xaxt="n")
  #axis(side=2,labels=)
  # need to manually adjust y here so that it sits right below the x axis
  text(cex=1, x=1:length(method.indices)+0.2, y=YLIM[1]-0.001, leg.key[method.indices],
       xpd=TRUE, srt=45, pos=2)
  par(font.lab=1)
  
  dev.off()

  
  png(sprintf("QQplots_dataset%d",1))
  layout(matrix(1:6,2,3,byrow=F))
  key <- as.vector(datasets[[d]][[3]])
  for (midx in 1:6) {
  qqUnif(as.vector(pvals[,,midx]),title=leg.key[midx])
  }
  
}
