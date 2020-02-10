#' Predict Method for sgspls
#'
#' Predicted values based on sparse group subgroup PLS. New responses are predicted using a fitted model and a new matrix of observations.
#'
#' @param object Object of class inheriting from \code{"sgspls"}. 
#' @param newdata Data matrix in which to look for for explanatory variables to be used for prediction.
#' @param ... Not currently used.
#' 
#' @export
#' @return \code{perf} returns a list that contains the following performance measures: 
#' \code{predict} function produces predicted values, obtained by evaluating the sparse group subgroup PLS. 
#' The prediction values are calculated based on the regression coefficients of \code{object$Y} onto \code{object$variates$X}.
#' 

predict.sgspls <-
  function(object, newdata, group_obs_test,  ...)  {
    
    # predict function for sgspls
    newdata <-  as.matrix(newdata)
    nobs <- nrow(newdata)
    p <- ncol(newdata)
    nresp <- ncol(object$parameters$Y)
    npred <- ncol(object$parameters$X)
    ncomp <- object$parameters$ncomp

    #-- validation des arguments --#
    if (missing(newdata)){
      stop("No new data available.")
    }
    
    if (length(dim(newdata)) == 0) {
      if (length(newdata) != p)
        stop("'newdata' must be a numeric matrix with ncol = ", p,
             " or a vector of length = ", p, ".")
      dim(newdata) = c(1, p)
    }

    B <- array(0, dim = c(npred, nresp, ncomp))
    B_coef <- coef(object, type = "coefficients")
    B <- B_coef$B
    B0 <- B_coef$B0
    
    pred <- array(dim = c(nobs, nresp, ncomp))
    
    design <- matrix(1, nrow = dim(newdata)[1], ncol = object$parameters$n_study)
    
    if(object$parameters$n_study > 1){
      for (is in 1:object$parameters$n_study){
        design[-which(group_obs_test == is),is] <- 0
      }
      if(!object$parameters$structured){
        newdata <- shift(newdata,object$parameters$n_study,group_obs_test)
        pred <- array(dim = c(nobs, nresp * object$parameters$n_study, ncomp))
      }
    } 
    
    # C'est B et B0 qui pèchent parce que un  x_scores est pas bien fait en amont
    for( i in 1:ncomp ){
      pred[,,i] <- newdata%*%B[ , , i] + design %*% B0[,,i]
    }
    pred0 <- NULL
    if(object$parameters$da){
      temp <- pred
      pred0 <- pred
      for( i in 1:ncomp ){
        temp[,,i] = t(apply(pred[,,i], 1, FUN = function(x){as.numeric(x==max(x))}))
      }
      pred <- temp
    }
    
    if(object$parameters$n_study > 1 && !(object$parameters$structured)){
      max_group_obs <- object$parameters$n_study 
      temp <- array(dim = c(nobs, nresp , ncomp))
      for( i in 1:ncomp ){
        temp[,,i] <- unshift(pred[,,i],object$parameters$n_study,group_obs_test)
      }
      pred <- temp
    }
    return(list(pred=pred, pred0 = pred0, B = B, B0 = B0))
  }

