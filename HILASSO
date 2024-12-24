rm(list = ls())

# Load necessary packages
library(glmnet)
library(progress)

# Standardization function
standardization <- function(X, y) {
  X_mean <- colMeans(X)
  X_sd <- apply(X, 2, sd)
  X_sc <- scale(X, center = X_mean, scale = X_sd)
  
  y_mean <- mean(y)
  y_sc <- y - y_mean
  
  list(X_sc = X_sc, y_sc = y_sc, X_sd = X_sd, y_mean = y_mean)
}

# Constructor function for Hi-LASSO class
HiLasso <- function(q1 = "auto", q2 = "auto", L = 30, alpha = 0.05, 
                    logistic = FALSE, random_state = NULL, n_jobs = 1) {
  hl <- list(
    q1 = q1,
    q2 = q2,
    L = L,
    alpha = alpha,
    logistic = logistic,
    random_state = random_state,
    n_jobs = n_jobs,
    n = NULL,
    p = NULL,
    X = NULL,
    y = NULL,
    sample_weight = NULL,
    select_prob = NULL,
    coef_ = NULL,
    intercept_ = NULL,
    coefficients = NULL
  )
  class(hl) <- "HiLasso"
  return(hl)
}

# Hi-LASSO fit method
fit.HiLasso <- function(hl, X, y, sample_weight = NULL) {
  # Initialize data
  hl$n <- nrow(X)
  hl$p <- ncol(X)
  hl$X <- as.matrix(X)
  hl$y <- as.numeric(y)
  
  # Set q1 and q2
  hl$q1 <- ifelse(hl$q1 == "auto", hl$n, hl$q1)
  hl$q2 <- ifelse(hl$q2 == "auto", hl$n, hl$q2)
  
  # Cap q1 and q2 at p to avoid sampling errors
  hl$q1 <- min(hl$q1, hl$p)
  hl$q2 <- min(hl$q2, hl$p)
  
  # Set sample weights
  hl$sample_weight <- if (is.null(sample_weight)) {
    rep(1, hl$n)
  } else {
    as.numeric(sample_weight)
  }
  
  # Step 1: Compute importance scores
  cat("Procedure 1: Computing importance scores...\n")
  beta1_res <- bootstrapping_h(hl, mode = "procedure1")
  beta1 <- beta1_res[-1, ]
  beta1_mean <- rowMeans(abs(beta1), na.rm = TRUE)
  
  hl$select_prob <- ifelse(beta1_mean == 0, 1e-10, beta1_mean) / sum(beta1_mean)
  
  # Step 2: Compute coefficients and select variables
  cat("Procedure 2: Computing coefficients...\n")
  beta2_res <- bootstrapping_h(hl, mode = "procedure2")
  beta2 <- beta2_res[-1, ]
  beta2_mean <- rowMeans(beta2, na.rm = TRUE)
  
  # ---------------
  
  b2 <- beta2_res  # Extract betas
  
  # Compute mean coefficients across bootstraps (including intercept)
  b2_mean <- suppressWarnings(apply(b2, 1, function(x) mean(x, na.rm = TRUE)))
  b2_se <- suppressWarnings(apply(b2, 1, function(x) sd(x, na.rm = TRUE)))
  
  # Degrees of freedom for OLS (adjusted for intercept)
  n <- nrow(X)
  df <- n - length(b2_mean)
  
  # Compute t-statistics
  t_stats <- b2_mean / b2_se
  
  # Compute p-values based on t-distribution
  p_values <- 2 * pt(-abs(t_stats), df = df)
  
  # Compute confidence intervals
  alpha <- 0.05 # 95% Confidence intervel
  t_value <- qt(1 - alpha / 2, df = df)
  hl$lower_bound <- b2_mean - t_value * b2_se
  hl$upper_bound <- b2_mean + t_value * b2_se
  hl$Std_Error <- b2_se
  
  # ---------------
  
  # Compute p-values
  hl$p_values_ <- compute_p_values(hl, beta2)
  hl$p_values_with_intercept <- p_values
  
  hl$coef_ <- ifelse(hl$p_values_ < hl$alpha, beta2_mean, 0)
  
  hl$intercept_ <- mean(hl$y) - colMeans(hl$X) %*% hl$coef_

  
  # Assign names to coefficients
  names(hl$coef_) <- paste0("X", 1:hl$p)
  
  # Combine intercept and coefficients
  hl$coefficients <- c(Intercept = hl$intercept_, hl$coef_)
  return(hl)
}

# Bootstrapping method
bootstrapping_h <- function(hl, mode) {
  q <- ifelse(mode == "procedure1", hl$q1, hl$q2)
  method <- ifelse(mode == "procedure1", "ElasticNet", "AdaptiveLASSO")
  B <- floor(hl$L * hl$p / q)
  
  # Initialize progress bar
  pb <- progress_bar$new(
    format = "  [:bar] :percent eta: :eta",
    total = B, clear = FALSE, width = 60
  )
  
  # Initialize coefficient matrix
  betas <- matrix(NA, nrow = hl$p + 1, ncol = B)
  
  # Perform bootstrapping
  for (b in 1:B) {
    betas[, b] <- estimate_coef_h(hl, b, q, method)
    pb$tick()
  }
  
  return(betas)
}

# Coefficient estimation method
estimate_coef_h <- function(hl, bootstrap_number, q, method) {
  beta <- rep(NA, hl$p + 1)
  
  # Set random seed
  if (!is.null(hl$random_state)) {
    set.seed(hl$random_state + bootstrap_number)
  } else {
    set.seed(NULL)
  }
  
  
  # Sample indices for observations and variables
  bst_sample_idx <- sample(1:hl$n, size = hl$n, replace = TRUE)
  bst_predictor_idx <- sample(1:hl$p, size = q, replace = FALSE, prob = hl$select_prob)
  
  # Subset X and y
  X_bst <- hl$X[bst_sample_idx, bst_predictor_idx, drop = FALSE]
  y_bst <- hl$y[bst_sample_idx]
  
  # Standardization
  std <- standardization(X_bst, y_bst)
  X_sc <- std$X_sc
  y_sc <- std$y_sc
  X_sd <- std$X_sd
  
  # Select model
  if (method == "ElasticNet") {
    fit <- cv.glmnet(X_sc, y_sc, alpha = 0.5, weights = hl$sample_weight[bst_sample_idx], 
                     standardize = FALSE, nfolds = 5,)
    coef_all <- coef(fit, s = "lambda.min")
    coef_fit <- as.numeric(coef_all)[-1]
    
    coef_all_vec <- as.vector(coef(fit, s = 'lambda.min'))
    beta[1] <- as.numeric(coef_all)[1]  # Intercept
  } else if (method == "AdaptiveLASSO") {
    penalty_weights <- 1 / (hl$select_prob[bst_predictor_idx] * 100)
    fit <- cv.glmnet(X_sc, y_sc, alpha = 1, penalty.factor = penalty_weights, 
                     weights = hl$sample_weight[bst_sample_idx],
                     standardize = FALSE)
    coef_all <- coef(fit, s = "lambda.min")
    coef_fit <- as.numeric(coef_all)[-1]
    
    coef_all_vec <- as.vector(coef(fit, s = 'lambda.min'))
    beta[1] <- as.numeric(coef_all)[1]  # Intercept
  }
  
  # Assign coefficients back to beta vector
  beta[bst_predictor_idx + 1] <- coef_fit / X_sd

  return(beta)
}

# Method for calculating p-values
compute_p_values <- function(hl, betas) {
  not_null <- !is.na(betas)
  d_j <- rowSums(betas != 0 & not_null)
  pi <- sum(d_j) / sum(not_null)
  B <- ncol(betas)
  
  p_values <- pbinom(d_j - 1, size = B, prob = pi, lower.tail = FALSE)
  return(p_values)
}
  
# Print method for HPLasso in table format
print.HiLasso <- function(hl, ...) {
  cat("HiLasso Model\n")
  cat("Parameters:\n")
  cat("  q1 =", hl$q1, "\n")
  cat("  q2 =", hl$q2, "\n")
  cat("  random_state =", ifelse(is.null(hl$random_state), "NULL", hl$random_state), "\n")
  cat("  L =", hl$L, "\n\n")
  
  if (!is.null(hl$coef_)) {
    cat("Model has been fitted.\n")
    cat("Number of samples (n):", hl$n, "\n")
    cat("Number of predictors (p):", hl$p, "\n\n")
    
    # Prepare the table for coefficients and p-values
    coefficients <- c(hl$intercept_, hl$coef_)
    names(coefficients) <- c("Intercept", paste0("Beta_", 1:hl$p))
    
    p_values <- c(hl$p_values_with_intercept[1], hl$p_values_)
    names(p_values) <- c("Intercept", paste0("Beta_", 1:hl$p))
    
    lower_bound <- hl$lower_bound
    upper_bound <- hl$upper_bound
    
    
    # Create a data frame
    table <- data.frame(
      Variable = names(coefficients),
      Coefficient = coefficients,
      Std_Error = hl$Std_Error,
      Lower_bound = lower_bound,
      Upper_bound = upper_bound,
      P_Value = p_values,
      stringsAsFactors = FALSE
    )
    
    # Print the table
    print(table, row.names = FALSE)
  } else {
    cat("Model has not been fitted yet.\n")
  }
}

# Example usage

set.seed(123)   # For reproducibility

# Step 1: Generate a synthetic dataset
n <- 100  # Number of observations
p <- 20   # Number of predictors
X <- matrix(rnorm(n * p), n, p) # Predictor matrix with n rows and p columns

# Define the true beta vector with 6 non-zero coefficients and 14 zeros
beta <- c(3, 2, 1.2, 0.5, 0.2, 1.5, rep(0, p - 6)) # 6 non-zero coefficients + 14 zeros = 20
inter <- 0.2
y <- inter + X %*% beta + rnorm(n)  # Response with noise

hi_lasso <- HiLasso(q1 = 10, q2 = 10,
                    random_state = 123, L = 30)

# Fit the model
hi_lasso <- fit.HiLasso(hi_lasso, X, y, sample_weight = NULL)

# Print the model summary
print(hi_lasso)

