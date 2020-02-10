# Author: Camilo Broc (camilo.broc@univ-pau.fr)
# Description: Performs sPLS, MINT, sgPLS, sPLS for structured data and sgPLS for structured data and evaluate their performance.

# Library statements ---------------------------


# library(sgPLS)
# library(sgspls)

# Sources ---------------------------

# source(file = "auxiliar.R")
# source(file = "group_index.R")


my.perf.spls <- function(X,Y,l.param,method,is.tune,grid.tune){
  
  # Compute the methods and evaluate their performance (MSEP, TPR, TD) on synthetic data. 
  # Those data are generated following the methodology of the article. 
  
  ### Args
  # keepX : For sPLS, MINT and sPLS for structured data, refers to the number of variables kept. For the others refers to the number of groups kept.
  # n.train :  total number of observations in the training set
  # n.test : total number of observations in the test set
  # l.param : list of parameters for the generation of the synthetic data
  # seed :  seed to reproduce the results
  # method : "spls", "MINT", "sgpls", "spls_structure", "sgpls_structure"
  # method.generation : whether "pb5" for first simulation case and "pb6" for second simualtion case. 
  
  ### Out
  # X,Y : standardized matrices
  # mean.X, sd.X : (s x dim(X)[1]) matrices of resp. mean and standard deviation by which X is standardized. each row correspond to a different observation set
  # mean.Y, sd.Y : same for Y
  
  ###
  # ncomp <- 2
  
  ###
  # l.param <- l_param
  p <- l.param$p
  s <- l.param$s
  q <- l.param$q
  n <- l.param$n
  da <- l.param$da
  group_seq<- l.param$group_seq
  ncomp <- l.param$ncomp
  
  # D <- l.param$D
  
  list.group.var <- l.param$list_group_var
  list.group.obs <- l.param$list_group_obs
  
  
  if (method == "sgpls"){
    
    #####
    ### sgpls with Matthew's Package
    #####
    
    #############################
    
    groupX <- subgroupX <- list.group.var
    
    
    if(is.tune){
      sparsities <- cbind(grid.tune,
                          rep(0,length(grid.tune)),
                          1-grid.tune)
      
      model.tune = list(X = X, Y = Y,
                        groupX=groupX, subgroupX=subgroupX,da = da,
                        group_obs = list.group.obs, structured = FALSE,
                        n_study = max(list.group.obs),ncomp = 0)
      for (icomp in 1:ncomp){
        a <- Sys.time()
        model.tune<- tune_sgspls(pls_obj = model.tune,
                                 sparsities = sparsities,
                                 group_seq = group_seq,
                                 folds= 5,
                                 scale_resp = FALSE,
                                 progressBar = FALSE,
                                 setseed = l.param$iseed)
        b <- Sys.time()
        print("time sgspls cv")
        print(b-a)

      }
      
      keepX = model.tune$parameters$keepX
      indiv_sparsity_x =  model.tune$parameters$indiv_sparsity_x
      
    } else {
      # print("sgPLS_str 2")
      model.tune <- NULL
      keepX = group_seq
      indiv_sparsity_x = grid.tune
    }
    
    a <- Sys.time()
    model<- sgspls(X = X, Y = Y,
                   ncomp = ncomp,
                   mode = "regression",
                   # keepX = model.tune$parameters$keepX,
                   keepX =  keepX,
                   groupX = groupX,
                   subgroupX = subgroupX, 
                   da = da, group_obs = list.group.obs, structured = FALSE,n_study = max(list.group.obs),
                   indiv_sparsity_x =  indiv_sparsity_x,
                   # indiv_sparsity_x =  c(0.8,0,0.2),
                   subgroup_sparsity_x = rep(0,ncomp))
    b <- Sys.time()
    print("time sgspls simple")
    print(b-a)
    
    #############################
    
  } else if( (method == "spls_structure")||(method == "sgpls_structure")){
    
    
 
    
    # shifting X and Y
    
    X_shift <- shift(X,s,group_obs = list.group.obs)
    Y_shift <- shift(Y,s,group_obs = list.group.obs)
    
    if (method == "sgpls_structure"){
      
      #####
      ### sgspls structured with my code
      #####
      
      #############################
      
      groupX<- rep(list.group.var,s) 
      subgroupX<-rep(1:p,s) 
      if(is.tune){
        sparsities <- cbind(grid.tune,
                            1-grid.tune,
                            rep(0,length(grid.tune)))
        model.tune <- list(X = X_shift, Y = Y_shift,
                        groupX=groupX, subgroupX=subgroupX, da = da,
                        group_obs = list.group.obs,
                        structured = TRUE,
                        n_study = max(list.group.obs),ncomp =0)
        for (icomp in 1:ncomp){
          a <- Sys.time()
          model.tune<- tune_sgspls(pls_obj =  model.tune,
                                   sparsities = sparsities,
                                   group_seq = group_seq,
                                   folds= 5,
                                   scale_resp = FALSE,
                                   progressBar = TRUE,
                                   setseed = l.param$iseed)
          b <- Sys.time()
          print("time sgspls structure cv")
          print(b-a)
        }
        
        
        
        keepX = model.tune$parameters$keepX
        subgroup_sparsity_x =  model.tune$parameters$subgroup_sparsity_x
        
        #############################
        
      } else {
        model.tune <- NULL
        keepX = group_seq
        subgroup_sparsity_x = grid.tune
      }

      a <- Sys.time()
      model<- sgspls(X = X_shift, Y = Y_shift,
                     ncomp = ncomp,
                     mode = "regression",
                     # keepX = model.tune$parameters$keepX,
                     keepX = keepX,
                     groupX = groupX,
                     subgroupX = subgroupX,
                     da = da, group_obs = list.group.obs, structured = TRUE, n_study = max(list.group.obs),
                     indiv_sparsity_x =  rep(0,ncomp),
                     subgroup_sparsity_x = subgroup_sparsity_x)

      b <- Sys.time()
      print("time_model")
      print(b-a)
    }
  } else {
    model <- spls(X,
                  Y,
                  ncomp=1,
                  mode = "regression"
                  ,keepX=keepX
    )
  }
  
  
  if( (method == "spls_structure")||(method == "sgpls_structure")){
    # The variables must be reordered
    temp.p <- p
    weights.var.0 <- model$weights$X
    weights.var <- rep(0,p)
    for (i in 1:p){
      weights.var[i] <-  sqrt(sum( (weights.var.0[(0:(s-1))*p +i])^2))
    }
    selected.C <- which(weights.var !=0)
  } else if(method == "sgpls"){
    
    weights.var  <- model$weights$X
    selected.C <- which(weights.var !=0)
  } else {
    
    weights.var <- model$weights$X
    selected.C <- which(weights.var !=0)
  }
  
  selected.C
  
  return(list(
    selected.C = selected.C,
    weights.var = weights.var,
    model = model
    ,model.tune = model.tune
  ))
  
}

