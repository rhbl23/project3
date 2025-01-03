library(glmnet)

gendat_select <- function(seed, n) {
  # Set seed for reproducibility
  set.seed(seed)
  desired_SNR = 2
  
  l <- list()
  for (i in 1:100) {
    l[[paste0("x", i)]] <- rnorm(n, 0, 1)
  }
  
  # Linear model to generate y (x1 to x10 contribute to y, while x11 to x100 are noise)
  y_signal <- 0.2 + 3 * l$x1 + 2 * l$x2 + 1.2 * l$x3 + 0.5 * l$x4 + 0.2 * l$x5 +
    1.5 * l$x6 + 5 * l$x7 + 0.3 * l$x8 + 2.5 * l$x9 + 4 * l$x10 
  
  y_signal <- scale(y_signal, center = TRUE, scale = FALSE) 
  var_signal <- var(as.vector(y_signal))
  var_noise <- var_signal / desired_SNR
  sigma <- sqrt(var_noise)
  noise <- rnorm(n, mean = 0, sd = sigma)
  y <- y_signal + noise + 0.2  # Adds back the intercept
  data <- data.frame(y, l)
  
  return(data)
}

nn <- 50  # sample size
replicate <- 50
s <- 4     # starting seed
nnr <- length(nn) * replicate # total number of runs
res_hi <- array(data=NA, dim=c(nnr, 104), dimnames=list(paste(1:nnr), c(
  "nn", "sam", 'intercept', paste0("b", 1:100), "sigma2")))


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
$sensitivity
[1] 0.52

$specificity
[1] 0.8795556

$MCC
[1] 0.326633

$g_mean
[1] 0.6762905

$f1
[1] 0.3993856
