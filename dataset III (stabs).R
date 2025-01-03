gendat_select <- function(seed, n){
  
  desired_SNR = 2
  
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
nnr <- length(nn) * replicate # total number of runs
res_stabs <- array(data=NA, dim=c(nnr, 22), dimnames=list(paste(1:nnr), c(
  "nn", "sam", "p1", "p2", "p3", "p4", "p5", "p6", "p7", "p8", "p9", "p10",
  "p11", "p12", "p13", "p14", "p15", "p16", "p17", "p18", "p19", "p20")))

counter <- 1
for (i in 1:length(nn)) {
  for (j in s:(s + replicate - 1)) {
    res_stabs[counter, 1] <- nn[i]  # Store the sample size
    res_stabs[counter, 2] <- j      # Store the seed number
    
    # Generate the dataset with covariates x1 to x20
    dd <- gendat_select(j, nn[i])
    
    # Prepare the data for Stability Selection
    x_matrix <- as.matrix(dd[, -1])  # Covariates (x1 to x20)
    y_vector <- dd$y                 # Outcome variable
    
    # Fit Stability Selection Lasso model
    stabs_lasso <- stabsel(x_matrix, y_vector, fitfun = glmnet.lasso, cutoff = 0.75,PFER = 1)
    
    # Extract the probabilities for intercept and x1 to x20
    res_stabs[counter, 3:22] <- stabs_lasso$max[1:20] 
    
    # Increment the counter
    if (counter < nrow(res_stabs)) counter <- counter + 1
  }
}


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
compute_aggregated_metrics_stabs <- function(res_stabs) {
  nnr <- nrow(res_stabs)
  
  # True variables (x1 to x5 are relevant, x6 to x10 are noise)
  true_relevant <- c(1, 1, 1, 1, 1, 1, rep(0,14))  # True relevance of variables x1 to x10
  
  # Initialize counters for TP, TN, FP, and FN across all simulations
  total_tp <- 0
  total_tn <- 0
  total_fp <- 0
  total_fn <- 0
  
  for (counter in 1:nnr) {
    # Extract the coefficients (excluding intercept)
    selected_vars <- as.numeric(res_stabs[counter, 3:22] >= 0.75)  #selection probability greater than 0.75 indicate selection
    
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

compute_aggregated_metrics_stabs(res_stabs)
$sensitivity
[1] 0.2966667

$specificity
[1] 1

$MCC
[1] 0.4774459

$g_mean
[1] 0.5446712

$f1
[1] 0.4575835
