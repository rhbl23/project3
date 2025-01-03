gendat_select <- function(seed, n){
  
  desired_SNR = 5
  
  beta = c(3, 2, 1.2, 0.5, 0.2, 1.5, rep(0, 14))
  
  # Create correlation matrix for first 6 variables
  cor_matrix <- matrix(0.9, 6, 6)
  diag(cor_matrix) <- 1
  
  # Cholesky decomposition
  chol_matrix <- chol(cor_matrix)
  
  # Generate correlated variables for first 6 predictors
  X_corr <- matrix(rnorm(n * 6), nrow = n, ncol = 6) %*% chol_matrix
  
  X_indep <- matrix(rnorm(n * 14), nrow = n, ncol = 14)
  
  X <- cbind(X_corr,X_indep)
  colnames(X) <- paste0("x", 1:20)
  
  # Generate response with intercept
  y_signal <- 0.2 + X %*% beta
  y_signal <- scale(y_signal, center = TRUE, scale = FALSE)  # Centers y_signal (mean zero)
  
  var_signal <- var(as.vector(y_signal))
  var_noise <- var_signal / desired_SNR
  sigma <- sqrt(var_noise)
  noise <- rnorm(n, mean = 0, sd = sigma)
  y <- y_signal + noise + 0.2  # Adds back the intercept
  
  # Create and return final dataset
  dataset <- data.frame(y = as.vector(y), X)
  
  return(dataset)
}


nn <- 100  # sample size
replicate <- 50
s <- 4     # starting seed

library(glmnet)
# Create an array for storing results
nnr <- length(nn) * replicate # total number of runs
res_hi <- array(data=NA, dim=c(nnr, 24), dimnames=list(paste(1:nnr), c(
  "nn", "sam", "intercept", "b1", "b2", "b3", "b4", "b5", "b6", "b7", "b8", "b9", "b10",
  "b11", "b12", "b13", "b14", "b15", "b16", "b17", "b18", "b19", "b20","sigma2"))) 

counter <- 1
for (i in 1:length(nn)) {
  for (j in s:(s + replicate - 1)) {
    res_hi[counter, 1] <- nn[i]  # Store the sample size
    res_hi[counter, 2] <- j      # Store the seed number
    
    # Generate the dataset with covariates x1 to x20
    dd <- gendat_select(j, nn[i])
    
    # Prepare the data for Hi Lasso
    x_matrix <- as.matrix(dd[, -1])  # Covariates (x1 to x20)
    y_vector <- dd$y                 # Outcome variable
    
    # Fit Hi Lasso model
    hi_lasso <- HiLasso(q1 = 10, q2 = 10, random_state = 123, L = 30)
    hi_lasso <- fit.HiLasso(hi_lasso, x_matrix, y_vector)
    
    # Extract the coefficients for intercept and x1 to x20
    res_hi[counter,3] <- hi_lasso$intercept_
    res_hi[counter, 4:23] <- hi_lasso$coef_[1:20]  # 1 intercept + 20 coefficients
    
    
    # Calculate and store the estimated residual variance (sigma^2)
    hl_coefficients <- as.numeric(hi_lasso$coefficients[-1])
    # Convert intercept_ from vector to scalar
    intercept <- as.numeric(hi_lasso$intercept_)[1]
    y_pred <- x_matrix %*% hl_coefficients + intercept
    res_hi[counter, 24] <- var(y_vector - y_pred)
    
    # Increment the counter
    if (counter < nrow(res_hi)) counter <- counter + 1
  }
}
compute_aggregated_metrics(res_hi)