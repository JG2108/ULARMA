SEfromHessian <- function(a, hessian=FALSE, silent=FALSE) {
  namesp <- colnames(a)
  mathessian <- a
  
  mathessian <- ifelse(mathessian==-Inf, -1E9, mathessian)
  mathessian <- ifelse(mathessian==+Inf, 1E9, mathessian)
  
  sigma  <- try(solve(mathessian), silent=TRUE)
  # Add all: 22/4/2020; any !
  if (inherits(sigma, "try-error")) {
    if (!silent) warning("Error in Hessian matrix inversion")
    mathessianx <- try(as.matrix(getFromNamespace("nearPD", ns="Matrix")(mathessian)$mat), silent=TRUE)
    # 29/1/2021
    if (inherits(mathessianx, "try-error")) {
      if (!silent) warning("Error in estimation of the Nearest Positive Definite Matrix. Calculates the Moore-Penrose generalized inverse. Use result with caution.")
      sigma  <- try(ginv(mathessian), silent=TRUE)
    } else {
      if (!silent) warning("Calculates the Nearest Positive Definite Matrix. Use result with caution.")
      sigma  <- try(solve(mathessianx), silent=TRUE)
    }
  } 
  
  # Add all: 22/4/2020
  if (!inherits(sigma, "try-error")) {
    
    if (all(diag(sigma)>=0)) {
      # méthode classique
      res <- sqrt(diag(sigma))
      
    } else {
      # Autre essai
      s. <- svd(sigma)
      R <- t(s.$v %*% (t(s.$u) * sqrt(pmax(s.$d, 0))))
      res <- structure((matrix(rep(1, nrow(R)), nrow = 1, byrow = TRUE) %*% R)[1, ], .Names=colnames(mathessian))
      
      # Si j'ai des SE négatif... pas bon signe
      if (any(res<0)) {
        d <- diag(as.matrix(getFromNamespace("nearPD", ns="Matrix")(sigma)$mat))
        names(d) <- colnames(mathessian)
        res <- ifelse(d<0, NA, sqrt(d))
      }
      if (any(is.na(res))) {
        
        
        a <- sigma
        
        n = dim(a)[1];
        root = matrix(0,n,n);
        
        for (i in 1:n){
          sum = 0;
          if (i>1){
            sum = sum(root[i,1:(i-1)]^2);
          }
          
          x = a[i,i] - sum;
          
          if (x<0){
            x = 0;
          }
          
          root[i,i] = sqrt(x);
          
          if (i < n){
            for (j in (i+1):n){
              
              if (root[i,i] == 0){
                x=0;
              }
              else{
                sum = 0;
                if (i>1) {
                  sum = root[i,1:(i-1)] %*% t(t(root[j,1:(i-1)]))
                }
                x = (a[i,j] - sum)/root[i,i];
              }
              
              root[j,i] = x;
            }
          }
        }
        colnames(root) <- rownames(root) <- colnames(mathessian)
        
        pseudoV <- root %*% t(root)
        
        d <- diag(pseudoV)
        
        if (any(d != 0) & all(d >= 0)) {
          res <- sqrt(d)
          
          if (!silent) warning("Estimates using pseudo-variance based on Gill & King (2004)")
          
        } else {
          if (!silent) warning("Approximation of Cholesky matrix based on Rebonato and Jackel (2000)")
          if (!silent) warning("Estimates using pseudo-variance based on Gill & King (2004)")
          
          # The paper by Rebonato and Jackel, “The most general methodology for creating a valid correlation matrix for risk management and option pricing purposes”, Journal of Risk, Vol 2, No 2, 2000, presents a methodology to create a positive definite matrix out of a non-positive definite matrix.
          # FP Brissette, M Khalili, R Leconte, Journal of Hydrology, 2007, “Efficient stochastic generation of multi-site synthetic precipitation data”
          # GA Baigorria, JW Jones, Journal of Climate, 2010, “GiST: A stochastic model for generating spatially and temporally correlated daily rainfall data”
          # M Mhanna and W Bauwens, International Journal of Climatology, 2012, “A stochastic space-time model for the generation of daily rainfall in the Gaza Strip”
          # fix the correl matrix
          newMat <- a
          cholError <- TRUE
          iter <- 0
          while (cholError) {
            
            iter <- iter + 1
            # cat("iteration ", iter, "\n")
            
            # replace -ve eigen values with small +ve number
            newEig <- eigen(newMat)
            newEig2 <- ifelse(newEig$values < 0, 0, newEig$values)
            
            # create modified matrix eqn 5 from Brissette et al 2007, inv = transp for
            # eig vectors
            newMat <- newEig$vectors %*% diag(newEig2) %*% t(newEig$vectors)
            
            # normalize modified matrix eqn 6 from Brissette et al 2007
            newMat <- newMat/sqrt(diag(newMat) %*% t(diag(newMat)))
            
            # try chol again
            cholStatus <- try(u <- chol(newMat), silent = TRUE)
            cholError <- ifelse(inherits(cholStatus, "try-error"), TRUE, FALSE)
          }
          root <- cholStatus
          
          colnames(root) <- rownames(root) <- colnames(mathessian)
          
          pseudoV <- root %*% t(root)
          
          res <- sqrt(diag(pseudoV))
        }
      }
    }
  }
  
  
  SEInf <- namesp[!namesp %in% names(res)]
  
  res <- c(res, structure(rep(+Inf, length(SEInf)), .Names=SEInf))
  
  if (hessian) {
    return(list(SE=res, hessian=mathessian))
  } else {
    return(res)
  }
  
}