gendat_select <- function(seed, n) {
  # Set seed for reproducibility
 set.seed(seed)
  l <- list()
  for (i in 1:20) {
    l[[paste0("x", i)]] <- rnorm(n, 0, 1)
  }
  
  # Generate error term
  epsilon <- rnorm(n, 0, 1)
  
  # Linear model to generate y (x1 to x6 contribute to y, while x7 to x20 are noise)
  y <- 0.2 + 3*l$x1 + 2*l$x2 + 1.2*l$x3 + 0.5*l$x4 + 0.2*l$x5 + 1.5*l$x6+epsilon
  
  # Return a data frame with y and covariates
  data <- data.frame(y, l)
  
  return(data)
}

nn <- 100  # sample size
replicate <- 50
s <- 4     # starting seed

# Create an array for storing results
nnr <- length(nn) * replicate # total number of runs
res_random <- array(data=NA, dim=c(nnr, 24), dimnames=list(paste(1:nnr), c(
  "nn", "sam", "intercept", "b1", "b2", "b3", "b4", "b5", "b6", "b7", "b8", "b9", "b10",
  "b11", "b12", "b13", "b14", "b15", "b16", "b17", "b18", "b19", "b20","sigma2"))) 

counter <- 1
for (i in 1:length(nn)) {
  for (j in s:(s + replicate - 1)) {
    res_random[counter, 1] <- nn[i]  # Store the sample size
    res_random[counter, 2] <- j      # Store the seed number
    
    # Generate the dataset with covariates x1 to x20
    dd <- gendat_select(j, nn[i])
    
    # Prepare the data for Random Lasso
    x_matrix <- as.matrix(dd[, -1])  # Covariates (x1 to x20)
    y_vector <- dd$y                 # Outcome variable
    
    # Fit Random Lasso model
    random_lasso <- RandomLasso(q1 = 10, q2 = 10, B = 'auto', 
                                random_state = 123, L = 30, threshold = 0.01)
    random_lasso <- fit.RandomLasso(random_lasso, x_matrix, y_vector)
    
    # Extract the coefficients for intercept and x1 to x20
    res_random[counter,3] <- random_lasso$intercept_
    res_random[counter, 4:23] <- random_lasso$coef_[1:20]  # 1 intercept + 20 coefficients
    
    # Calculate and store the estimated residual variance (sigma^2)
    y_pred <- x_matrix %*% random_lasso$coefficients[-1] + random_lasso$intercept_
    res_random[counter, 24] <- var(y_vector - y_pred)
    
    # Increment the counter
    if (counter < nrow(res_random)) counter <- counter + 1
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

compute_aggregated_metrics(res_random)
$sensitivity
[1] 0.98

$specificity
[1] 0.4728571

$MCC
[1] 0.439035

$g_mean
[1] 0.6807349

$f1
[1] 0.6105919
