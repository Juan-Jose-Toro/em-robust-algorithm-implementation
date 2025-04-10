
# Helper function
dmvnorm <- function(Mu, Sigma, x) {
  # computes f(x; Mu, Sigma)
  #
  # Args: 
  #   Mu: mean 
  #   Sigma: covariance matrix
  #   x: observed data
  #
  # Returns:
  #   f(x; Mu, Sigma), density function of N(Mu, Sigma) evaled at X
  
  # d variate GMM
  d <- length(Mu)
  # compute f
  1/sqrt((2*pi)^(d)*det(Sigma))*exp(-1/2*t(x - Mu)%*%solve(Sigma)%*%(x - Mu)) 
}

EM <- function(X, C, Z, tol=1e-10, m_iter=1e3) {
  # Estimates centers of Gaussian Mixture models:
  #
  # Args:
  #   X: Data, d by n matrix, n obs, d variate
  #   C: Number of clusters
  #   Z: Labels (suspected clusters) each X_i belongs to
  #   tol: tolerance for convergence
  #   m_iter: max iterations
  #
  # Returns:
  #   The estimated K centers from data X
  
  n = dim(X)[2] # Num of observations
  d = dim(X)[1] # Data dimensions
  
  # Initialize Z matrix
  Z_update <- matrix(0, n, C)
  for (i in 1:C) {
    Z_update[which(Z == i),i] <- 1
  }
  Z = Z_update
  # Initialize estimated contributions, mu, and Sigma for each suspected cluster
  alpha = NULL
  mu = NULL
  Sigma = NULL  
  
  ### for snap shots, Examples:
  ###   Density estimation, 
  ###   Keep mu and Sigmas for Ellipsoid plotting
  # Snap <- list() 
  ###
  
  for (i in 1:m_iter) {
    # update alpha 
    alpha = colSums(Z)/n
    
    # update mu
    mu <- lapply(1:C, function(k, X, Z) {
      X %*% Z[,k]/sum(Z[,k])
    }, X, Z)
    
    # update Sigma
    Sigma = lapply(1:C, function(k, X, mu, Z_update,d) {
      # compute residuals
      residuals <- sweep(X, MARGIN = 1, STATS = mu[[k]]) 
      matrix(as.matrix(apply(residuals, 2, function(x) {
        x %*% t(x)
      }))%*%Z_update[,k],d,d)/sum(Z_update[,k])
    }, X, mu, Z_update,d)
    
    # update Z 
    for (j in 1:n) { # eq. (4)
      for (k in 1:C) {
        Z_update[j,k] = alpha[k]*dmvnorm(mu[[k]], Sigma[[k]], X[,j])
      }
      Z_update[j,] <- Z_update[j,]/sum(Z_update[j,])
    }
    
    if (norm(Z_update - Z, "2") < tol) {
      # collects mu, sigma, alpha, and final iteration number 
      OUT <- list(mu=mu, Sigma=Sigma, alpha=alpha,iter=i)
      return(OUT)
    }
    Z <- Z_update
  }
  message(paste("failed to converge in", m_iter,"iterations."))
  return(NULL)
}

# EM Robust algorithm
EM_Robust <- function(X, nC, tol = 1e-10, m_iter=1e4) {
  # Estimates centers of Gaussian Mixture models without manual initialization:
  #
  # Args:
  #   X: Data, d by n matrix, n obs, d variate
  #   tol: tolerance for convergence
  #   m_iter: max iterations
  #   nC: number of clusters, should be <= n
  #
  # Returns:
  #   The estimated K centers from data X, variances, and contributions 
  
  n = dim(X)[2] # number of observations
  d = dim(X)[1] # dimension of data
  nC <- n
  # init beta:
  beta = 1
  # init C, C idx is shifted by 1, C[[iter_t+1]] is C[[t]] 
  C <- list(nC)
  # init alpha as 1/c(initial)
  alpha = rep(1/nC, nC)
  # init mu with every observations in X
  mu_update <- lapply(1:n, function(k) { X[,k] })
  # init Sigma
  tmp = ceiling(sqrt(n))
  d_m = .Machine$double.xmax
  
  outliers <- c()
  for (i in 1:(n-1)) {
    for (j in (i+1):n){
      dif <- norm(X[,i] - X[,j], "2")^2
      if (dif > 0 & dif < d_m) {
        d_m = dif
      }
    }
  }
  Sigma <- lapply(1:n, function(k, X, tmp, d, d_m) {
    # Use eq 28 to avoid issue with matrix being close to singular
    g <- 0.0001
    (1-g)*norm(X[,k] - X[,tmp], "2")^2 * diag(d) + g*d_m*diag(d)
  }, X, tmp, d, d_m)
  # Initialize Z
  Z <- as.data.frame(matrix(0, n, n))
  for (j in 1:n) {
    for (k in 1:n) {
      Z[j,k] = alpha[k]*dmvnorm(mu_update[[k]], Sigma[[k]], X[,j])
    }
    Z[j,] <- Z[j,]/sum(Z[j,])
  }
  # Compute first iteration of mu
  mu_update <- lapply(1:C[[1]], function(k, X, Z) {
    X %*% Z[,k]/sum(Z[,k])
  }, X, Z)
  mu <- mu_update
  
  # Loop: 
  for (i in 1:m_iter) {
    # update of alpha
    a_em = colSums(Z)/n 
    alpha_update = sapply(1:C[[i]], function(k, alpha, a_em) {
      a_em[k] + alpha[k]*(log(alpha[k]) - t(alpha)%*%(log(alpha)))
    }, alpha, a_em)
    
    # update beta
    tmp1 <- (1-max(a_em))/(-max(alpha)*t(alpha)%*%alpha)
    eta <- min(1, 0.5^floor(d/2 - 1))
    tmp2 <- 0
    for (k in 1:C[[i]]) { # eq 24
      tmp2<-tmp2+exp(-eta*n*abs(alpha_update[k] - alpha[k]))/C[[i]]
    }
    # set to perma 0 when beta is set to 0 (stable C)
    beta <- min(tmp1, tmp2, beta)
    
    # discard all clusters with contribution < 1/n
    keep <- which(alpha_update >= 1/n)
    
    alpha_update <- alpha_update[keep]
    # update mu (discard clusters)
    mu_update <- mu_update[keep]
    # update number of clusters
    C[[i+1]] <- length(alpha_update)
    # update alpha, eq 15
    alpha <- alpha_update/sum(alpha_update)
    # update Z, eq 16
    Z <- as.matrix(Z[,keep])
    for (j in 1:n) {
      Z[j,] <- Z[j,]/sum(Z[j,])
    }
    # check stability of C
    if (i >= 60) {
      if (C[[i-59]] == C[[i+1]]) {
        beta <- 0
      }
    }
    
    outlier <- which(is.na(Z[,1]))
    if (length(outlier) > 0) {
      Z <- Z[-outlier,]
      outliers <- cbind(outliers, X[,outlier])
      X <- X[,-outlier]
      n <- n - length(outlier)
    }
    
    # update Sigma
    Sigma <- lapply(1:C[[i+1]], function(k, X, mu_update, Z, d, d_m) {
      gamma <- 0.0001
      # compute all X_i - Mu
      x_diff_mu <- sweep(X, MARGIN = 1, STATS = mu_update[[k]]) 
      # eq 26
      
      vcovs <- matrix(as.matrix(apply(x_diff_mu, 2, function(x) {
        x %*% t(x)
      }))%*%Z[,k], d, d)/sum(Z[,k])
      # eq 28
      (1 - gamma)*vcovs + gamma*d_m*diag(d)
    }, X, mu_update, Z,d,d_m)
    
    # update Z
    for (j in 1:n) {
      for (k in 1:C[[i+1]]) {
        Z[j,k] = alpha[k]*dmvnorm(mu_update[[k]], Sigma[[k]], X[,j])
      }
      Z[j,] <- Z[j,]/sum(Z[j,])
    }
    # update mu
    mu_update <- lapply(1:C[[i+1]], function(k, X, Z) {
      X %*% Z[,k]/sum(Z[,k])
    }, X, Z)
    
    dif <- as.matrix(sapply(mu_update, cbind) - sapply(mu[keep], cbind))
    m_dif <- max(apply(dif, 2, function(d) { norm(d, "2") }))
    if (m_dif < tol) {
      # check against tolerance
      OUT <- list(mu = mu_update, 
                  Sigma = Sigma, 
                  alpha = alpha,
                  iter = i,
                  outliers = outliers) 
      # return a named list of elements of interest
      return(OUT)
    }
    mu <- mu_update
  }
  # If failed to converge
  message(paste("failed to converge in", m_iter,"iterations."))
  OUT <- list(mu = mu_update, 
              Sigma = Sigma, 
              alpha = alpha,
              iter = i,
              outliers = outliers) 
  # return a named list of elements of interest
  return(OUT)
}

# helper functions
library(shape)
draw_ellipse <- function(EM_res, i) {
  colour <- c("steelblue", "aquamarine", 
              "darkorchid", "brown1", 
              "deeppink3", "orange")
  center <- EM_res$mu[[i]]
  tmp <- sqrt(eigen(EM_res$Sigma[[i]])$values*qchisq(p=0.95, df = 2))
  
  lines(getellipse(rx = tmp[2], ry = tmp[1], mid = center), col = colour[i],
        lwd = 2)
  points(x = center[1], y = center[2], pch = 19, col = "red")
}


