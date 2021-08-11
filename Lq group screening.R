## MLqE of coefficients of regression using each group of variables.
MLqE.est <- function(X, Y, q, eps1){
  ##X: the group of variables
  ##Y: the response
  ##q: the distortion parameter
  ##eps1: the iteration coverage criterion
  beta0 <- solve(t(X)%*%X)%*%t(X)%*%Y
  sigma0 <- crossprod(Y-X%*%beta0)/length(Y)
  t <- 1
  beta_old <- beta0
  sigma_old <- sigma0
  repeat
  {
    omega_hat <- NULL
    for(i in 1:length(Y)){
      omega_hat[i] <- (1/sqrt(2*pi*sigma_old)*exp(-1/(2*sigma_old)*(Y[i]-X[i,]%*%beta_old)^2))^(1-q)
    }
    OMEGA_new <- diag(omega_hat)
    beta_new <- solve(t(X)%*%OMEGA_new%*%X)%*%t(X)%*%OMEGA_new%*%Y
    sigma_new <- sum(omega_hat*(Y-X%*%beta_new)^2)/sum(omega_hat)
    if ((crossprod(beta_new - beta_old) <= eps1))
      break
    t <- t + 1
    beta_old <- beta_new
    sigma_old <- sigma_new
  }
  return(list(t=t, beta_hat = beta_new, sigma_hat = sigma_new, OMEGA_hat = OMEGA_new))
}

## Group screening of LqG
blosc.MLqE <- function(data, p, m, n, q){
  ##data: the list of predictor X=X and reponse Y=Y
  ##p: the dimension of predictors
  ##m: the number of groups
  ##n: sample size
  ##q: the distortion parameter
  pp = rep(p/m, m); psum <- cumsum(c(0, pp))
  Y <- data$Y; X <- data$X
  Y <- Y - mean(Y); X <- apply(X, 2, function(c) c-mean(c))
  X <- apply(X, 2, scale)
  eps1 = 1e-6
  beta.group <- NULL
  for(l in 1:m){
    Xb = X[ ,(psum[l]+1):psum[l+1]]
    betab <- MLqE.est(Xb, Y, q, eps1)$beta_hat
    beta_mean <- sqrt(crossprod(betab))/ncol(Xb)
    beta.group <- c(beta.group, beta_mean)
  }
  return(beta.group)
}

## Group screening of Lq1
blosc.marg.MLqE <- function(data, p, m, n, q){
  pp = rep(p/m, m); psum <- cumsum(c(0, pp))
  Y <- data$Y; X <- data$X
  Y <- Y - mean(Y); X <- apply(X, 2, function(c) c-mean(c))
  X <- apply(X, 2, scale)
  eps1 = 1e-6
  beta.marginal <- NULL
  for(l in 1:p){
    Xb = as.matrix(X[ ,l])
    betab <- MLqE.est(Xb, Y, q, eps1)$beta_hat
    beta.marginal <- c(beta.marginal, betab)
  }
  beta.group <- NULL
  for(l in 1:m){
    beta_mean <- sqrt(crossprod(beta.marginal[(psum[l]+1):psum[l+1]]))/(p/m)
    beta.group <- c(beta.group, beta_mean)
  }
  return(beta.group)
}