############################################################################
# Data generation for variable selection
########################################################################
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
  # Noise variables (x7 to x20)
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
res_lasso <- array(data=NA, dim=c(nnr, 24), dimnames=list(paste(1:nnr), c(
  "nn", "sam", "intercept", "b1", "b2", "b3", "b4", "b5", "b6", "b7", "b8", "b9", "b10",
  "b11", "b12", "b13", "b14", "b15", "b16", "b17", "b18", "b19", "b20","sigma2"))) 

counter <- 1
for (i in 1:length(nn)) {
  for (j in s:(s + replicate - 1)) {
    res_lasso[counter, 1] <- nn[i]  # Store the sample size
    res_lasso[counter, 2] <- j      # Store the seed number
    
    # Generate the dataset with covariates x1 to x20
    dd <- gendat_select(j, nn[i])
    
    # Prepare the data for glmnet
    x_matrix <- as.matrix(dd[, -1])  # Covariates (x1 to x20)
    y_vector <- dd$y                 # Outcome variable
    
    # Fit LASSO using cross-validation to select lambda
    cv_fit <- cv.glmnet(x_matrix, y_vector, alpha = 1)  # alpha=1 for LASSO
    
    
    # Get the coefficients at the lambda that gives minimum cross-validated error
    best_lambda <- cv_fit$lambda.min
    lasso_coef <- coef(cv_fit, s = best_lambda)
    
    # Extract the coefficients for intercept and x1 to x20
    res_lasso[counter, 3:23] <- as.vector(lasso_coef)[1:21]  # 1 intercept + 20 coefficients
    
    # Calculate and store the estimated residual variance (sigma^2)
    y_pred <- predict(cv_fit, s = best_lambda, newx = x_matrix)
    res_lasso[counter, 24] <- var(y_vector - y_pred)
    
    # Increment the counter
    if (counter < nrow(res_lasso)) counter <- counter + 1
  }
}

#######################################################################
#Performance metrics
#########################################################################
# Function to compute MCC
compute_mcc <- function(tp, tn, fp, fn) {
  numerator <- (tp * tn) - (fp * fn)
  denominator <- sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
  
  if (denominator == 0) {
    return(0)  # Prevent division by zero
  } else {
    return(numerator / denominator)
  }
}

f1_score <- function(tp, fp, fn) {
  
  tpr <- tp / (tp + fn)
  ppv <- tp/ (tp + fp)
  
  # Calculate F1 Score
  f1_score <- 2*(ppv*tpr)/(ppv+tpr)
  
  return(f1_score)
}


# Function to compute aggregated sensitivity, specificity, MCC, and g-mean
compute_aggregated_metrics <- function(res_lasso) {
  nnr <- nrow(res_lasso)
  
  # True variables (x1 to x6 are relevant, x7 to x20 are noise)
  true_relevant <- c(1, 1, 1, 1, 1, 1, rep(0,14))  # True relevance of variables x1 to x20
  
  # Initialize counters for TP, TN, FP, and FN across all simulations
  total_tp <- 0
  total_tn <- 0
  total_fp <- 0
  total_fn <- 0
  
  for (counter in 1:nnr) {
    # Extract the coefficients (excluding intercept)
    selected_vars <- as.numeric(res_lasso[counter, 4:23] != 0)  # Non-zero coefficients indicate selection
    
    # Calculate True Positives (TP), True Negatives (TN), False Positives (FP), and False Negatives (FN)
    tp <- sum(selected_vars[1:6] == 1)  # Correctly selected relevant variables
    tn <- sum(selected_vars[7:20] == 0) # Correctly excluded noise variables
    fp <- sum(selected_vars[7:20] == 1) # Incorrectly selected noise variables
    fn <- sum(selected_vars[1:6] == 0)  # Missed relevant variables
    
    # Accumulate totals for TP, TN, FP, and FN
    total_tp <- total_tp + tp
    total_tn <- total_tn + tn
    total_fp <- total_fp + fp
    total_fn <- total_fn + fn
  }
  
  # Compute aggregated sensitivity, specificity, MCC, and g-mean
  sensitivity <- total_tp / (total_tp + total_fn)
  specificity <- total_tn / (total_tn + total_fp)
  mcc <- compute_mcc(total_tp, total_tn, total_fp, total_fn)
  g_mean <- sqrt(sensitivity * specificity)
  f1 <- f1_score(total_tp, total_fp, total_fn)
  
  # Return the aggregated metrics as a named list
  return(list(
    sensitivity = sensitivity,
    specificity = specificity,
    MCC = mcc,
    g_mean = g_mean,
    f1=f1
  ))
}

compute_aggregated_metrics(res_lasso)

#result
#$sensitivity
#[1] 0.9666667

#$specificity
#[1] 0.5657143

#$MCC
#[1] 0.4967935

#$g_mean
#[1] 0.7394979

#$f1
#[1] 0.6487696

