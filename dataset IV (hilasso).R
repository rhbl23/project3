gendat_select <- function(seed, n){
  
  desired_SNR = 2
  
  beta = c(3, 2, 1.2, 0.5, 0.2, 1.5, 5, 0.3, 2.5, 4, rep(0, 90))
  
  # Create correlation matrix for first 10 variables
  # Create correlation matrix for the first 10 variables
  cor_matrix <- matrix(0, 10, 10)
  diag(cor_matrix) <- 1  # Diagonal elements are 1, indicating each variable is perfectly correlated with itself
  
  # Set the correlation coefficients
  cor_matrix[1, 2] <- cor_matrix[2, 1] <- 0.8  # High correlation between x1 and x2
  cor_matrix[1, 3] <- cor_matrix[3, 1] <- 0.6  # Moderate correlation between x1 and x3
  cor_matrix[2, 3] <- cor_matrix[3, 2] <- 0.5  # Moderate correlation between x2 and x3
  cor_matrix[4, 5] <- cor_matrix[5, 4] <- 0.7  # Moderate correlation between x4 and x5
  cor_matrix[6, 7] <- cor_matrix[7, 6] <- 0.6  # Moderate correlation between x6 and x7
  cor_matrix[8, 9] <- cor_matrix[9, 8] <- 0.5  # Moderate correlation between x8 and x9
  cor_matrix[9, 10] <- cor_matrix[10, 9] <- 0.7  # Moderate correlation between x9 and x10
  
  # Cholesky decomposition
  chol_matrix <- chol(cor_matrix)
  
  # Generate correlated variables for first 6 predictors
  X_corr <- matrix(rnorm(n * 10), nrow = n, ncol = 10) %*% chol_matrix
  
  X_indep <- matrix(rnorm(n * 90), nrow = n, ncol = 90)
  
  X <- cbind(X_corr,X_indep)
  colnames(X) <- paste0("x", 1:100)
  
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

nn <- 50  # sample size
replicate <- 50
s <- 4     # starting seed
nnr <- length(nn) * replicate # total number of runs
res_hi <- array(data=NA, dim=c(nnr, 105), dimnames=list(paste(1:nnr), c(
  "nn", "sam", 'intercept', paste0("b", 1:100), "sigma2","mse")))


counter <- 1
for (i in 1:length(nn)) {
  for (j in s:(s + replicate - 1)) {
    res_hi[counter, 1] <- nn[i]  # Store the sample size
    res_hi[counter, 2] <- j      # Store the seed number
    
    # Generate the dataset with covariates x1 to x100
    dd <- gendat_select(j, nn[i])
    
    # Prepare the data for Hi Lasso
    x_matrix <- as.matrix(dd[, -1])  # Covariates (x1 to x100)
    y_vector <- dd$y                 # Outcome variable
    
    # Fit Hi Lasso model
    hi_lasso <- HiLasso(q1 = 50, q2 = 50, random_state = 123, L = 30)
    hi_lasso <- fit.HiLasso(hi_lasso, x_matrix, y_vector)
    
    # Extract the coefficients for intercept and x1 to x100
    res_hi[counter,3] <- hi_lasso$intercept_
    res_hi[counter, 4:103] <- hi_lasso$coef_[1:100]  # 1 intercept + 100 coefficients
    
    # Calculate and store the estimated residual variance (sigma^2)
    hl_coefficients <- as.numeric(hi_lasso$coefficients[-1])
    # Convert intercept_ from vector to scalar
    intercept <- as.numeric(hi_lasso$intercept_)[1]
    y_pred <- x_matrix %*% hl_coefficients + intercept
    res_hi[counter, 104] <- var(y_vector - y_pred)
    
    # Calculate and store Mean Squared Error (MSE)
    mse <- mean((y_vector - y_pred)^2)
    res_hi[counter, 105] <- mse
    
    # Increment the counter
    if (counter < nrow(res_hi)) counter <- counter + 1
  }
}

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
  
  # True variables (x1 to x10 are relevant, x11 to x100 are noise)
  true_relevant <- c(rep(1,10), rep(0,90))  # True relevance of variables x1 to x100
  
  # Initialize counters for TP, TN, FP, and FN across all simulations
  total_tp <- 0
  total_tn <- 0
  total_fp <- 0
  total_fn <- 0
  
  print("res_random:")
  print(res_lasso)
  
  for (counter in 1:nnr) {
    # Extract the coefficients (excluding intercept)
    selected_vars <- as.numeric(res_lasso[counter, 4:103] != 0)  # Non-zero coefficients indicate selection
    
    # Calculate True Positives (TP), True Negatives (TN), False Positives (FP), and False Negatives (FN)
    tp <- sum(selected_vars[1:10] == 1)  # Correctly selected relevant variables
    tn <- sum(selected_vars[11:100] == 0) # Correctly excluded noise variables
    fp <- sum(selected_vars[11:100] == 1) # Incorrectly selected noise variables
    fn <- sum(selected_vars[1:10] == 0)  # Missed relevant variables
    
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
  print("g_mean:")
  print(g_mean)
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

compute_aggregated_metrics(res_hi)
