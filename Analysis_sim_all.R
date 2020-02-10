# Author: Camilo Broc (camilo.broc@univ-pau.fr)
# Description: Application of sgPLS, ASSET and metaSKAT on Simulated data
# -----

# This code simulates data following the article. The part "Simulate data" sets the parameters of the simulation.
# In the following "structured sgPLS" stands for the joint sgPLS method.

#####
### Clean phase
#####

rm(list=ls())
graphics.off()
options(stringsAsFactors = FALSE)
options(na.strings = c("", "NA", " "))
gc(reset = T)
today = format(Sys.Date(),"%Y%m%d")
invisible(lapply(paste0('package:', names(sessionInfo()$otherPkgs)), detach, character.only=TRUE, unload=TRUE))
cat("\014")

# path <-dirname(rstudioapi::getActiveDocumentContext()$path)
# setwd(path)

#####
### library statements
#####

Sys.time()-> time1

library(mixOmics)
library(sgspls)
library(Rcpp)
library(MASS)
library(magrittr)
library(tidyverse)
library(MetaSKAT)
library(ASSET)
require(testit)
require(BhGLM)
assign("last.warning", NULL, envir = baseenv())
#####
### sourcing files
#####

path_function <- "./Source_code"

source(paste0(path_function,"/" ,"sgPLS_structure_v2.R"))
sourceCpp(paste0(path_function,"/" ,"sgsPLS.cpp"))

source(paste0(path_function,"/" ,"coef.R"))
source(paste0(path_function,"/" ,"internal.R"))
source(paste0(path_function,"/" ,"internalPLS.R"))
source(paste0(path_function,"/" ,"sgspls.R"))
source(paste0(path_function,"/" ,"predict.R"))
source(paste0(path_function,"/" ,"perf.R"))
source(paste0(path_function,"/" ,"tuning.R"))

source(paste0(path_function,"/" ,"Simulation2.R"))
#####
### Simulate data
#####

seed.start <- 1
seed.end<- 1
nseed <- seed.end -seed.start+1
res_template <- matrix(0,ncol=4,nrow=2)
colnames(res_template) <- c(                                "TP",
                                                            "FP",
                                                            "FN",
                                                            "TN")
rownames(res_template) <- c("SNP","gene")

res_tot_sgPLS <- res_tot_ASSET<- res_tot_metaSKAT <-res_tot_jsgPLS <- res_template
# <- res_tot_sgPLS.true <- res_tot_jsgPLS.true 
a1 <- Sys.time()  
for (iseed in seed.start:seed.end){
  
# Generation of simulated data
N = 200 # Number of individuals
P = 500 # Number of variables
K = 2 # Number of studies
ngrp = 25 # Number of group
group = parallel::splitIndices(P, ngrp) # groups as a list of indices

set.seed(iseed)
groupSparsity = 0.8# percentage of group being sparse
intraGroupSparsity = 0.5 # percentage of marker being sparse within non-sparse group
corr = c(0.5, 0) # Intragroup correlation, between group correlation
MAF = 0.3 # Minor Allele Frequency
typeEffect = "fixed" # Or "variable"
# typeEffect = "variable" # Or "variable"
maxEffect = exp(1) # maximum odd ratio for marker non sparse. If "fixed", effect is maxEffect
verbose = TRUE # Make it chatty

# Number of true groups
grpNonNull = 1:(floor(length(group)*(1-groupSparsity))+1)

nonNull = lapply(1:length(group), function(x) { # for each group
  len = length(group[[x]]) # Save length
  if (x %in% grpNonNull) { # If group is in non-null group
    markerNonNull = 1:floor(len*(1-intraGroupSparsity)) # First markers taken
    markerNonNull = ifelse(1:len %in% markerNonNull, 1, 0) # convert into a vector of 0/1
    return(markerNonNull)
  } else { # If group is 0
    return(rep(0, len))
  }
}) %>% do.call(c, .) 

nonNull <- unlist(lapply(nonNull,function(x){sample(c(-1,1),1)*x}))
pAssoc <- replicate(2, nonNull)


simulation = sim_sparse_group_pleio(N=N, 
                                    P=P, 
                                    K=K, 
                                    group=group, 
                                    pAssoc = pAssoc,
                                    groupSparsity = groupSparsity,
                                    intraGroupSparsity = intraGroupSparsity,
                                    corr = corr, 
                                    MAF = MAF,
                                    typeEffect = typeEffect,
                                    minEffect = minEffect,
                                    maxEffect = maxEffect,
                                    verbose = verbose)

x1 = as.matrix(simulation$Study1$X)
y1 = as.matrix(simulation$Study1$Y)
x2 = as.matrix(simulation$Study2$X)
y2 = as.matrix(simulation$Study2$Y)
simulation$simBeta
X1 <- x1
Y1 <- unmap(y1)
X2 <- x2
Y2 <- unmap(y2)

#####
### Apply sgPLS : retreat select and compare
#####

######
### Preparation
######

# iseed <- 1
n_sample1 <- N
n_sample2 <- N
ncomp <- 1
q <- 2
p <- P
is.tune <- TRUE
group_seq <-c(20)
grid.tune <- c(0.01,0.05,0.1,0.3,0.5,0.7,0.9,0.95,0.99)

group_seq <- c(1:10)
group_seq <- c(1:10)
grid.tune <- c(0.1,0.5,0.9)

list_group_var <- unlist(lapply(1:ngrp,function(x){rep(x,P/ngrp)}))

set.seed(iseed)

sample1 <-1:dim(x1)[1]
sample2 <-1:dim(x2)[1]
da <- TRUE


params_run <- list(
  iseed = iseed,
  group_seq = group_seq,
  sample1 = sample1,
  sample2 = sample2,
  da = da,
  is.tune = is.tune,
  grid.tune = grid.tune,
  ncomp= ncomp
)



#######

run <- function(params_run){
  
  res <- NULL
  
  ### sgPLS
  
  iseed = params_run$iseed
  group_seq = params_run$group_seq
  sample1 = params_run$sample1
  sample2 = params_run$sample2
  da = params_run$da
  is.tune = params_run$is.tune
  grid.tune = params_run$grid.tune
  ncomp <- params_run$ncomp
  
  X <- rbind(X1[sample1,],X2[sample2,])
  Y <- rbind(Y1[sample1,],Y2[sample2,])
  
  l_param <- list(p = dim(X)[2],
                  q = dim(Y)[2],
                  n = dim(X)[1],
                  s = 2,
                  group_seq = group_seq,
                  list_group_var = list_group_var,
                  list_group_obs = c(rep(1,length(sample1)),rep(2,length(sample2))),
                  da = da,
                  ncomp = ncomp,
                  iseed = iseed
  )
  
    t1<- Sys.time()
    gc()
    system.time(res.sgpls.12 <- my.perf.spls(X = X, Y = Y, l.param = l_param,method="sgpls", is.tune = is.tune, grid.tune = grid.tune))
    res$res.sgpls.12 <- res.sgpls.12 
    t2<- Sys.time()
    print(paste0("Time sgPLS ",t2-t1))
  
    ### joint sgPLS
  
    l_param <- list(p = dim(X)[2],
                    q = dim(Y)[2],
                    n = dim(X)[1],
                    s = 2,
                    group_seq = group_seq,
                    # group_seq1 = group_seq1[4],
                    #  grid.tune1 = grid.tune1[4],
                    list_group_var = list_group_var,
                    list_group_obs = c(rep(1,length(sample1)),rep(2,length(sample2))),
                    da = da,
                    ncomp = ncomp,
                    iseed = iseed
    )
    
    t1<- Sys.time()
    print(gc())
    system.time(res.sgplss <- my.perf.spls(X = X, Y = Y, l.param = l_param, method = "sgpls_structure", is.tune = is.tune, grid.tune = grid.tune))
    res$res.sgplss <- res.sgplss 
    t2<- Sys.time()
    print(paste0("Time structured sgPLS ",t2-t1))
  
  return(res)
  
}

Sys.time()-> time1
n_rep <- 1

  print(iseed)
  a <- Sys.time()
  system.time(temp <-run(params_run))
  b <- Sys.time()
  b-a
  print(b-a)
  
  res.sgpls.12.true <- temp$res.sgpls.12.true
  res.sgplss.true <- temp$res.sgplss.true  
  res.sgpls.12 <- temp$res.sgpls.12
  res.sgplss <- temp$res.sgplss
  
  res.sgpls.12.true$model.tune$parameters <- NULL
  res.sgplss.true$model.tune$parameters <- NULL
  res.sgpls.12$model.tune$parameters <- NULL
  res.sgplss$model.tune$parameters <- NULL
  
  res.sgpls.12.true$model$parameters <- NULL
  res.sgplss.true$model$parameters <- NULL  
  res.sgpls.12$model$parameters <- NULL
  res.sgplss$model$parameters <- NULL
 
  res_template <- matrix(0,ncol=4,nrow=2)
  colnames(res_template) <- c(                                "TP",
                                                           "FP",
                                                           "FN",
                                                           "TN")
  rownames(res_template) <- c("SNP","gene")
  true_select_SNP <- which(pAssoc[,1]!=0)
  others_SNP <- setdiff(1:P,true_select_SNP)
  
  true_select_gene <- grpNonNull 
  others_gene <- setdiff(1:length(group),true_select_gene)
  
  selected_SNP_sgPLS <- res.sgpls.12$selected.C
  selected_SNP <- selected_SNP_sgPLS
  selected_gene <- selected_SNP_sgPLS
  
  tabled <- function(true_SNP, true_gene, s_SNP, s_gene, t){
    t0 <- t
    t0[1,1] <- length(which(s_SNP %in% true_SNP))
    t0[1,2] <- length(which(s_SNP %in% setdiff(1:P, true_SNP))) 
    t0[1,3] <- length(which(setdiff(1:P,s_SNP) %in% true_SNP))
    t0[1,4] <- length(which(setdiff(1:P,s_SNP) %in% setdiff(1:P, true_SNP)))
    
    t0[2,1] <- length(which(s_gene %in% true_gene))
    t0[2,2] <- length(which(s_gene %in% setdiff(1:ngrp, true_gene))) 
    t0[2,3] <- length(which(setdiff(1:ngrp,s_gene) %in% true_gene))
    t0[2,4] <- length(which(setdiff(1:ngrp,s_gene) %in% setdiff(1:P, true_gene)))
    
    return(t0)
    
  }
 selected <- res.sgpls.12$selected.C
  res_sgPLS <- tabled(true_select_SNP, true_select_gene, selected, unique(floor(( selected/(P/ngrp) ) ))+1,res_template)
 selected <- res.sgplss$selected.C
  res_jsgPLS <- tabled(true_select_SNP, true_select_gene, selected, unique(floor(( selected/(P/ngrp) ) ))+1,res_template)

  
#####
### Apply ASSET : retreat select and compare
#####

  library(ASSET)
  t1<- Sys.time()
  t0 = Sys.time()
  tmp = apply(x1, 2, function(x) {
    x = matrix(x, ncol = 1)
    fit = glm(y1 ~ x, family=binomial())
    co = summary(fit)$coefficients
    return(c(co["x", "Estimate"], co["x", "Std. Error"], co["x", "Pr(>|z|)"]))
  })
  tmp2 = t(tmp) %>% as_tibble(rownames = "snp") %>% 
    rename(beta=V1, sigma=V2, pval=V3) %>%
    mutate(z = beta / sigma)
  Sys.time() - t0 # 26 secs OMG
summary1 = tmp2
t0 = Sys.time()
tmp = apply(x2, 2, function(x) {
  x = matrix(x, ncol = 1)
  fit = glm(y2 ~ x, family=binomial())
  co = summary(fit)$coefficients
  return(c(co["x", "Estimate"], co["x", "Std. Error"], co["x", "Pr(>|z|)"]))
})

tmp2 = t(tmp) %>% as_tibble(rownames = "snp") %>% 
  rename(beta=V1, sigma=V2, pval=V3) %>%
  mutate(z = beta / sigma)
Sys.time() - t0 # 26 secs OMG
summary2 = tmp2
case1 <- length(which(y1==1))
control1 <- length(which(y1==0))
case2 <- length(which(y2==1))
control2 <- length(which(y2==0))

case_1 = rep(case1,P)
control_1 = rep(control1,P)
case_2 = rep(case2,P)
control_2 = rep(control2,P)

beta = as.matrix(data.frame(data1 = summary1$beta, data2 = summary2$beta, 
                            row.names = summary1$snp))
sigma = as.matrix(data.frame(data1 = summary1$sigma, data2 = summary2$sigma, 
                             row.names = summary1$snp))
case = as.matrix(data.frame(data1 = case_1, data2 = case_2, 
                             row.names = summary1$snp))
control = as.matrix(data.frame(data1 = control_1, data2 = control_2, 
                                row.names = summary1$snp))
Study = c("Study1","Study2")

SNPs <- rownames(beta)

t0 = Sys.time()
res_ASSET = h.traits(SNPs, Study, beta, sigma, case, control, meta=TRUE)
Sys.time() - t0 # 20 secs

res_2sides = h.summary(res_ASSET)$Subset.2sided

t2<- Sys.time()
print(paste0("Time sgPLS ASSET ",t2-t1))
selected <- which(res_2sides$Pvalue <0.05)
res_ASSET <- tabled(true_select_SNP, true_select_gene, selected, NULL,res_template)

#####
### Apply metaSKAT : retreat select and compare
#####
t1<- Sys.time()

y.list = list(as.vector(y1), as.vector(y2))
cov.list = list(rep(0,N),rep(0,N))
X.list = list()
for (ig in 1:length(group)) {

  X.list[[ig]] = rbind(x1[,group[[ig]]],x2[,group[[ig]]])
}


t0 = Sys.time()
null = Meta_Null_Model(y.list,cov.list, n.cohort=2, out_type="D")
res = list()

for (ig in 1:length(group)) {
  pval_SKAT = MetaSKAT_wZ(X.list[[ig]], null, is.separate = TRUE)$p.value # Het-Meta-SKAT
  pval_burden = MetaSKAT_wZ(X.list[[ig]], null, r.corr=1, is.separate = TRUE)$p.value # Het-Meta-Burden
  pval_SKATO = MetaSKAT_wZ(X.list[[ig]], null, is.separate = TRUE, method="optimal")$p.value # Het-Meta-SKATO
  res_g = c(pval_SKAT, pval_burden, pval_SKATO)
  names(res_g) = c("pval_SKAT", "pval_burden", "pval_SKATO")
  res[[ig]] = res_g
}
res_SKAT_gene = rbind.data.frame(res) %>% t() %>% as_tibble(rownames = NA) 
Sys.time() - t0 # 2.1 mins

selected<- which(res_SKAT_gene$pval_SKATO < 0.05)
res_metaSKAT <- tabled(true_select_SNP, true_select_gene, NULL, selected,res_template)

t2<- Sys.time()
print(paste0("Time sgPLS mSKAT ",t2-t1))

#######
#######
#######

res_tot_sgPLS <- res_tot_sgPLS + 1/nseed*res_sgPLS
res_tot_jsgPLS <- res_tot_jsgPLS + 1/nseed*res_jsgPLS
res_tot_ASSET <- res_tot_ASSET + 1/nseed*res_ASSET
res_tot_metaSKAT <- res_tot_metaSKAT + 1/nseed*res_metaSKAT

}
b1<- Sys.time()
print("time")
print(b1-a1)

res_tot <- rbind(
  res_tot_sgPLS,
  res_tot_jsgPLS,
  res_tot_ASSET,
  res_tot_metaSKAT)
res <- NULL
res$res_tot <- res_tot
save(res,file="./Results/res.RData")
Sys.time()-> time2
print(time2-time1)

  res$res_tot