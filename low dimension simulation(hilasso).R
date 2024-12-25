gendat_select <- function(seed, n) {
  # Set seed for reproducibility
  set.seed(seed)
  
  # Generate covariates
  x1 <- rnorm(n, 0, 1)
  x2 <- rnorm(n, 0, 1)
  x3 <- rnorm(n, 0, 1)
  x4 <- rnorm(n, 0, 1)
  x5 <- rnorm(n, 0, 1)
  x6 <- rnorm(n, 0, 1)
  # Noise variables (x7 to x10)
  x7 <- rnorm(n, 0, 1)
  x8 <- rnorm(n, 0, 1)
  x9 <- rnorm(n, 0, 1)
  x10 <- rnorm(n, 0, 1)
  x11 <- rnorm(n, 0, 1)
  x12 <- rnorm(n, 0, 1)
  x13 <- rnorm(n, 0, 1)
  x14 <- rnorm(n, 0, 1)
  x15 <- rnorm(n, 0, 1)
  x16 <- rnorm(n, 0, 1)
  x17 <- rnorm(n, 0, 1)
  x18 <- rnorm(n, 0, 1)
  x19 <- rnorm(n, 0, 1)
  x20 <- rnorm(n, 0, 1)
  
  # Generate error term
  epsilon <- rnorm(n, 0, 1)
  
  # Linear model to generate y (x1 to x6 contribute to y, while x7 to x20 are noise)
  y <- 0.2 + 3*x1 + 2*x2 + 1.2*x3 + 0.5*x4 + 0.2*x5 + 1.5*x6+epsilon
  
  # Return a data frame with y and covariates
  data <- data.frame(y, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10,
                     x11, x12, x13, x14, x15, x16, x17, x18, x19, x20)
  
  return(data)
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
    res_hi[counter, 4:23] <- hi_lasso$coefficients[1:20]  # 1 intercept + 20 coefficients
    
    
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
$sensitivity
[1] 0.6766667

$specificity
[1] 0.9314286

$MCC
[1] 0.6426929

$g_mean
[1] 0.7938934

$f1
[1] 0.7368421
