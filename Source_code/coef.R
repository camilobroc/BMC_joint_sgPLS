get_loading <- 
  function(A, scores){
    tt <- scale(scores, scale =  apply(scores, 2, crossprod))
    return(crossprod(A, tt))
    
  }

get_adj_weights <- 
  function(x_loads, x_weights, nresp){
    
    ncomp <- ncol(x_loads)
    if(ncomp == 1) {return(x_weights)} 
    else {
      PW <- crossprod(x_loads, x_weights)
      if(nresp == 1){
        # Using trick from the pls package to speed up calculations
        PWinv <- diag(ncomp)
        bidiag <- - PW[row(PW) == col(PW)-1]
        for (a in 1:(ncomp - 1)){ 
        }
      } else {
        PWinv <- backsolve(PW, diag(ncomp))
      }
      return(x_weights%*%PWinv)
    }
  }

#' Compute regression coefficients, loadings and 
#' 
#' Function to extract the coefficients of the PLS models.
#'
#' @param object a sgspls object, parameters will be extracted from the object.
#' @param type  a vector of the type of coefficients to be extracted. Can be 
#' the regression \code{coefficients}, the \code{loadings}, the \code{adjusted_weights} or
#' any subset of these.
#' @param comps a vector of the components that should be returned.
#' @param ... not used currently.
#' 
#' @export
#' @return \code{coef} returns a list with the required measures.
#' 
#' @references Liquet Benoit, Lafaye de Micheaux, Boris Hejblum, Rodolphe
#'   Thiebaut. A group and Sparse Group Partial Least Square approach applied in
#'   Genomics context. \emph{Submitted}.
#' 
#' @examples
#'  set.seed(1)
#'  n = 50; p = 500; 
#'  size.groups = 30; size.subgroups = 5
#'  groupX <- ceiling(1:p / size.groups)
#'  subgroupX <- ceiling(1:p / size.subgroups)
#'  
#'  X = matrix(rnorm(n * p), ncol = p, nrow = n)
#'  
#'  beta <- rep(0,p)
#'  bSG <- -2:2; b0 <- rep(0,length(bSG))
#'  betaG <- c(bSG, b0, bSG, b0, bSG, b0)
#'  beta[1:size.groups] <- betaG
#'  
#'  y = X %*% beta + 0.1*rnorm(n)
#'  
#' model <- sgspls(X, y, ncomp = 3, mode = "regression", keepX = 1,
#'                 groupX = groupX, subgroupX = subgroupX,
#'                 indiv_sparsity_x = 0.8, subgroup_sparsity_x = 0.15)
#'  
#'  # get the regression coefficients
#'  model_coef <- coef(model, type = "coefficients", comps = 3)
#'  
#'  # check fit
#'  cbind(beta, model_coef$B)

coef.sgspls <- function(object, type = c("coefficients", "loadings", "adjusted_weights"), comps = NULL, ...){
    # Get data matrices with correct scale
    X <- object$parameters$X
    Y <- object$parameters$Y
    # scale.x : TRUE/FALSE
    # scale.y : TRUE/FALSE
    x_scale <- NULL
    y_scale <- NULL
    
    n_study <- object$parameters$n_study 
    group_obs <- object$parameters$group_obs
    if(n_study > 1){
      if (object$parameters$structured){
        X_s <- X
        Y_s <- Y
        
        p <- dim(X)[2]/n_study
        q <- dim(Y)[2]/n_study
        
        xmeans <- matrix(0, ncol = n_study * p, nrow = n_study)
        ymeans <- matrix(0, ncol = n_study * q, nrow = n_study)
        for (iobs in 1:n_study){
          X_study <- X[which(group_obs == iobs),1:p + (iobs-1)*p]
          Y_study <- Y[which(group_obs == iobs),1:q + (iobs-1)*q]
          
          X_s[which(group_obs == iobs),1:p + (iobs-1)*p] <- scale(X_study, center = T, scale =  object$parameters$scale.x)
          Y_s[which(group_obs == iobs),1:q + (iobs-1)*q] <- scale(Y_study, center = T, scale =  object$parameters$scale.y)
          
          xmeans[iobs, 1:p + (iobs-1)*p] <- colMeans(X_study)
          ymeans[iobs, 1:q + (iobs-1)*q] <- colMeans(Y_study)
          if(object$parameters$scale.x){
            x_scale <- cbind(rbind(x_scale), rbind(apply(X_study, 2, sd)))
          } else {
            x_scale <- cbind(rbind(x_scale), rbind(rep(1, p)))
          }
          if(object$parameters$scale.y){
            y_scale <- cbind(rbind(y_scale), rbind(apply(Y_study, 2, sd)))
          } else {
            y_scale <- cbind(rbind(y_scale), rbind(rep(1, q)))
          }
        }
      } else {
        
     
        p <- dim(X)[2]
        q <- dim(Y)[2]
       
        X_s <- X
        Y_s <- Y
        
        
        xmeans <- matrix(0, ncol = n_study * p, nrow = n_study)
        ymeans <- matrix(0, ncol = n_study * q, nrow = n_study)
        
        for (iobs in 1:n_study){
          X_study <- X[which(group_obs == iobs),]
          Y_study <- Y[which(group_obs == iobs),]
          
          X_s[which(group_obs == iobs),] <- scale(X_study, center = T, scale = object$parameters$scale.x)
          Y_s[which(group_obs == iobs),] <- scale(Y_study, center = T, scale = object$parameters$scale.y)
          
          xmeans[iobs, 1:p + (iobs-1)*p] <- colMeans(X_study)
          ymeans[iobs, 1:q + (iobs-1)*q] <- colMeans(Y_study)
          if(object$parameters$scale.x){
            x_scale <- cbind(rbind(x_scale), rbind(apply(X_study, 2, sd)))
          } else {
            x_scale <- cbind(rbind(x_scale), rbind(rep(1, p)))
          }
          if(object$parameters$scale.y){
            y_scale <- cbind(rbind(y_scale), rbind(apply(Y_study, 2, sd)))
          } else {
            y_scale <- cbind(rbind(y_scale), rbind(rep(1, q)))
          }
        }
      }
    } else {
      
      x_scale <- rep(1, ncol(X))
      y_scale <- rep(1, ncol(Y))
      X_s <- scale(X, T, object$parameters$scale.x)
      if(object$parameters$scale.x){
        x_scale <- apply(X, 2, sd)
      }
      Y_s <- scale(Y, T, object$parameters$scale.y)
      if(object$parameters$scale.y){
        y_scale <- apply(Y, 2, sd)
      }
      xmeans <- rbind(colMeans(X))
      ymeans <- rbind(colMeans(Y))
    }
    
    ########
    
    nresp <- ncol(Y)
    npred <- ncol(X)
    # Get ncomp and mode
    ncomp <- object$parameters$ncomp
    mode <- object$parameters$mode
    
    if(is.null(comps)) comps <- 1:ncomp
    
    # Get scores
    x_scores <- object$scores$X
    y_scores <- if(mode == "regression") object$scores$X else object$scores$Y
    
    # Get weights
    x_weights <- object$weights$X
    y_weights <- object$weights$Y
    
    res <- list()
    x_loads <- y_loads <- x_adjusted_weights <- B <- NULL
    
    if( "loadings" %in% type ){
      x_loads <- get_loading(X_s, x_scores)
      y_loads <- get_loading(Y_s, y_scores)
      
      res$x_loads <- x_loads[ , comps , drop = F]
      res$y_loads <- y_loads[ , comps , drop = F]
      
    }
    
    if( "adjusted_weights" %in% type ){
      x_loads <- if(is.null(x_loads)) get_loading(X_s, x_scores) else x_loads
      
      if(ncomp == 1){
        x_adjusted_weights <- x_weights[,, drop = F]
      } else {
        x_adjusted_weights <- get_adj_weights(x_loads, x_weights, nresp)
      }
      res$x_adjusted_weights <- x_adjusted_weights[, comps, drop = F]
    }
    
    if( "coefficients" %in% type & mode == "regression" ){
      
     
      x_loads <- if(is.null(x_loads)) get_loading(X_s, x_scores) else x_loads
      y_loads <- if(is.null(y_loads)) get_loading(Y_s, x_scores) else y_loads
      x_adjusted_weights <- if(is.null(x_adjusted_weights)) get_adj_weights(x_loads, x_weights, nresp) else x_adjusted_weights

      if (n_study > 1){
        if (object$parameters$structured){
          B_hat <- array(0, dim = c(npred, nresp, ncomp))
          B0 <- array(0, dim = c(n_study, nresp , ncomp))
        } else {
          B_hat <- array(0, dim = c(npred * n_study, nresp * n_study, ncomp))
          B0 <- array(0, dim = c(n_study, nresp * n_study, ncomp))
        }
      } else {
        B_hat <- array(0, dim = c(npred, nresp, ncomp))
        B0 <- array(0, dim = c(1, nresp, ncomp))
      }
      #####
      for( a in 1:ncomp ){
        B <- tcrossprod(x_adjusted_weights[,1:a, drop = F], y_loads[,1:a, drop = F])
        
        if ((!object$parameters$structured) && (n_study > 1) ){
          temp_B <- matrix(0, nrow = npred *n_study, ncol = nresp*n_study )
          for (iobs in 1:n_study){
            temp_B[1:p + (iobs-1)*p, 1:q + (iobs-1)*q] <- B
          }
          B <- temp_B
        }
        if ((object$parameters$structured) && (n_study > 1) ){
          for (iobs in 1:n_study){
            B[1:p + (iobs-1)*p, -(1:q + (iobs-1)*q)] <- 0
          }
        }
        B <- apply( B, 2, function(b) scale(t(b),scale = x_scale, center = F) )
        B_hat[ , , a] <- scale(B, scale = 1/y_scale, center = F) 
        
        B0[, , a] <- ymeans - xmeans%*%B_hat[ , , a]
        
      }
      
      res$B <- B_hat[,,comps, drop = F]
      res$B0 <- B0
      #####
    }
    class(res) = c("sgspls_coef")
    return(res)
  }
