##%######################################################%##
#                                                          #
####                   Help functions                   ####
#                                                          #
##%######################################################%##

# Return a matrix P X K with 0 or 1. 0 for non association, 1 for association 
# with the marker P and the study K
simulate_association = function(K, group, pAssoc, groupSparsity, intraGroupSparsity) {
  
  grpNonNull = sample(1:length(group), floor(length(group)*(1-groupSparsity)))
  grpNonNull = 1:(floor(length(group)*(1-groupSparsity)))
  
  # Randomly pick markers according groupEffect and intraGroupSparsity proba
  nonNull = lapply(1:length(group), function(x) { # for each group
    len = length(group[[x]]) # Save length
    if (x %in% grpNonNull) { # If group is in non-null group
      markerNonNull = 1:floor(len*(1-intraGroupSparsity)) # randomly pick markers
      markerNonNull = ifelse(1:len %in% markerNonNull, 1, 0) # convert into a vector of 0/1
      return(markerNonNull)
    } else { # If group is 0
      return(rep(0, len))
    }
  }) %>% do.call(c, .) 
  # We got a vector telling if the marker is associated or not with at least 1 phenotype
  
  # Loop over these associated markers to randomly choose how many phenotypes the marker is associated with
  
  out <- list()
  temp = lapply(nonNull, function(x) {
    nAssoc = sample(1:K, 1, prob=pAssoc) # Randomly pick number of pheno associated
    studyPicked = sample(1:K, nAssoc) # Randomly pick which pheno are associated
    studyAssoc = rep(0, K) # Initiate vector
    studyAssoc[studyPicked] = 1 # Attribute 1 to pheno associated
    return(x*studyAssoc) # multiply by x to set non-associated marker to 0
  }) %>% do.call(rbind, .)
  out[["res"]] <- temp
  out[["grpNonNull"]] <- grpNonNull
  out[["nonNull"]] <- nonNull
  return(out)
}


# Take the output of simulate_association
# Output a matrix of the same dimension with simulated coefficients
# typeEffect can be fixed or variable
# If effect is fixed, the effect is maxEffect
# If effect is variable, the effect is randomly picked between minEffect and maxEffect
# simulate_beta = function(assoc, typeEffect, minEffect=2, maxEffect=5) {
#   nSim = ncol(assoc) * nrow(assoc)
#   random = runif(nSim) # Generate random number
#   # In 50% of case, OR is picked between 1.05 and 1.25
#   # In the other 50% of case, OR is picked between 1/1.25 and 1/1.05
#   OR = vector("numeric", nSim)
#   if (typeEffect=="variable") {
#     OR[random > 0.5] = runif(sum(random > 0.5), minEffect, maxEffect)
#     OR[random <= 0.5] = runif(sum(random <= 0.5), 1/maxEffect, 1/minEffect) 
#   } else if (typeEffect=="fixed") {
#     OR[random > 0.5] = maxEffect
#     OR[random <= 0.5] = 1/maxEffect
#   }
#   beta = log(OR)
#   # Multiply by assoc to set non-associated marker to 0
#   beta = matrix(assoc*beta, ncol=ncol(assoc), nrow=nrow(assoc)) 
#   return(beta)
# }

simulate_beta = function(assoc, typeEffect, minEffect=2, maxEffect=5) {
  # nSim = ncol(assoc) * nrow(assoc)
  nSim = nrow(assoc)
  random = runif(nSim) # Generate random number
  # In 50% of case, OR is picked between 1.05 and 1.25
  # In the other 50% of case, OR is picked between 1/1.25 and 1/1.05
  OR = vector("numeric", nSim)
  if (typeEffect=="variable") {
    OR[random > 0.5] = runif(sum(random > 0.5), minEffect, maxEffect)
    OR[random <= 0.5] = runif(sum(random <= 0.5), 1/maxEffect, 1/minEffect) 
  } else if (typeEffect=="fixed") {
    OR[random > 0.5] = maxEffect
    OR[random <= 0.5] = 1/maxEffect
  }
  beta = log(OR)
  beta = replicate(ncol(assoc),beta)
  #print("dim")
  #print(dim(beta))
  #print(dim(assoc))
  # Multiply by assoc to set non-associated marker to 0
  beta = matrix(assoc*beta, ncol=ncol(assoc), nrow=nrow(assoc)) 
  return(beta)
}
##%######################################################%##
#                                                          #
####                   Main function                    ####
#                                                          #
##%######################################################%##



# N = 50 # Number of individuals
# P = 101 # Number of variables
# K = 2 # Number of studies
# ngrp = 5 # Number of group
# group = parallel::splitIndices(P, ngrp) # groups as a list of indices
# pAssoc = rep(1/K, K) # probability for marker with effect to be associated with n study
# groupSparsity = 0.4 # percentage of group being sparse
# intraGroupSparsity = 0.5 # percentage of marker being sparse within non-sparse group
# corr = c(0.5, 0) # Intragroup correlation, between group correlation
# MAF = 0.3 # Minor Allele Frequency
# typeEffect = "fixed" # Or "variable"
# minEffect = 2 # minimum odd ratio for marker non sparse
# maxEffect = 5 # maximum odd ratio for marker non sparse. If "fixed", effect is maxEffect
# verbose = TRUE # Make it chatty

sim_sparse_group_pleio =  function(N, 
                                   P, 
                                   K, 
                                   group, 
                                   pAssoc = rep(1/K, K),
                                   groupSparsity = 0.4,
                                   intraGroupSparsity = 0.5,
                                   corr = c(0.5, 0), 
                                   MAF = 0.3,
                                   typeEffect = "fixed",
                                   minEffect = 2,
                                   maxEffect = 5,
                                   verbose=FALSE) {
    #print(pAssoc)
  # Check inputs
  # assert(K==length(pAssoc))
  # print(pAssoc)
  assert(groupSparsity>=0 & groupSparsity<=1)
  assert(intraGroupSparsity>=0 & intraGroupSparsity<=1)
  assert(length(unlist(group)) == P)
  
  t0 = Sys.time()
  
  # Simulate association and beta for each study
  if (verbose) cat("Simulating associations\n")
  # print(pAssoc)
  # simAssoc = simulate_association(K, group, pAssoc, groupSparsity, intraGroupSparsity)
  simAssoc = pAssoc
  # print(simAssoc)
  if (verbose) cat("Simulating coefficients\n")
  print(dim(pAssoc))
  simBeta = simulate_beta(simAssoc, typeEffect, minEffect, maxEffect)
  # simBeta = simulate_beta(pAssoc, typeEffect, minEffect, maxEffect)
  
  if (verbose) cat("Starting simulating X et Y for each study\n")
  out = list()
  for (i in 1:K) {
    if (verbose) cat("  Study",i,"\n")
    beta = simBeta[,i]
    if (verbose) cat("    X\n")
    simX = sim.x(N, P, corr=corr, p=MAF, group=group, genotype = 1:P)
    if (verbose) cat("    Y\n")
    simY = sim.y(simX, coefs=beta)$y.ordinal
    simStudy = list(X=simX, Y=simY, coefs=beta)
    out[[paste0("Study",i)]] = simStudy
  }
  
  t1 = Sys.time()
  minutes <- round(difftime(t1, t0, units = "min"), 3)
  out[["simAssoc"]] <- simAssoc 
  out[["simBeta"]] <- simBeta 
  if (verbose) cat("Computation took",minutes,"minutes\n")
  return(out)
}



##%######################################################%##
#                                                          #
####                 Plotting Functions                 ####
#                                                          #
##%######################################################%##



plot_coefs = function(coefs, group=NULL) {
  require(testit)
  require(tidyverse)
  
  # test input
  assert("X is tabular", any(class(coefs) %in% c("matrix", "data.frame")))
  # Test group here
  
  K = ncol(coefs) # K studies
  if (K==1) {
    colnames(coefs) = "useless"
  }
  Kn = colnames(coefs) # names of studies
  Vn = rownames(coefs) # names of variables
  coefs2 = coefs %>% as_tibble(rownames = "Vn") %>% 
    tidyr::gather(., Kn, key="Study", value="coef") %>%
    mutate(Vn=factor(Vn, levels=unique(Vn)))
  
  p = ggplot(coefs2, aes(colour=Study)) + 
    geom_hline(yintercept = 0, colour = gray(1/2), lty = 2) +
    geom_linerange(aes(x=Vn, ymin=coef, ymax=0), lwd=1,
                   position=position_dodge(width = 1/2)) +
    coord_flip() + theme_classic() + xlab("") +
    ggtitle("Coefficient plot")
  if (K==1) {
    p = p + theme(legend.position = "none") 
  } else {
    p = p
  }
  p
}


plot_genotype = function(X) {
  require(testit)
  require(tidyverse)
  
  # check input
  assert("X is tabular", any(class(X) %in% c("matrix", "data.frame")))
  assert("X only contains 0,1,2", all(unique(unlist(X)) %in% c(0,1,2)))
  
  Xmelted = X %>% as_tibble() %>% mutate(id=1:nrow(X)) %>%
    gather(key="Variable", value="genotype", -id) %>%
    mutate(Variable = factor(Variable, levels=unique(Variable)))
  
  
  ggplot(Xmelted) +
    geom_tile(aes(x=Variable, y=id, fill=factor(genotype)), colour="white") +
    scale_fill_manual(values = c("white", "sienna1", "red4")) +
    labs(x = "", y = "", fill="Genotype") + 
    scale_x_discrete(expand = c(0, 0), position = "bottom") +
    scale_y_discrete(expand = c(0, 0)) +
    theme(axis.text.x = element_text(size=10, angle = 90, hjust = 1, vjust = 0.5))
}
