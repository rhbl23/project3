
############################################################################
# Python Code translated to R for RLASSO
##########################################################################
rm(list=ls())

library(glmnet)
library(progress)

# Utility function for standardization
standardization <- function(X, y) {
  X_mean <- colMeans(X)
  X_sd <- apply(X, 2, sd)
  X_sc <- scale(X, center = X_mean, scale = X_sd)
  
  y_mean <- mean(y)
  y_sc <- y - y_mean
  
  list(X_sc = X_sc, y_sc = y_sc, X_sd = X_sd, y_mean = y_mean)
}

# Constructor for RandomLasso S3 class
RandomLasso <- function(q1 = 'auto', q2 = 'auto', B = 'auto', 
                        random_state = NULL, L = 30, threshold = 0.1) {
  rl <- list(
    q1 = q1,
    q2 = q2,
    B = B,
    random_state = random_state,
    L = L,
    threshold = threshold,  # Threshold for setting coefficients to zero
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
  class(rl) <- "RandomLasso"
  return(rl)
}

# Fit method for RandomLasso
fit.RandomLasso <- function(rl, X, y, sample_weight = NULL) {
  # Initialize number of samples and predictors
  rl$n <- nrow(X)
  rl$p <- ncol(X)
  rl$X <- as.matrix(X)
  rl$y <- as.numeric(y)
  
  # Set q1 and q2
  rl$q1 <- ifelse(rl$q1 == 'auto', rl$n, rl$q1)
  rl$q2 <- ifelse(rl$q2 == 'auto', rl$n, rl$q2)
  
  
  # Cap q1 and q2 at p to avoid sampling errors
  rl$q1 <- min(rl$q1, rl$p)
  rl$q2 <- min(rl$q2, rl$p)
  
  # Set B
  rl$B <- ifelse(rl$B == 'auto', floor(rl$L * rl$p / rl$q1), rl$B)
  
  # Set sample weights
  rl$sample_weight <- if (is.null(sample_weight)) {
    rep(1, rl$n)
  } else {
    as.numeric(sample_weight)
  }
  
  rl$select_prob <- NULL
  
  # Store the mean of y for intercept
  # rl$intercept_ <- mean(rl$y)
  rl$intercept_ <- NULL
  
  # Procedure 1: Bootstrapping for Importance Scores
  cat("Procedure 1: Bootstrapping for Importance Scores\n")
  proc1_result <- bootstrapping(rl, mode = "procedure1")
  beta1 <- proc1_result$betas
  beta1_mean <- rowMeans(abs(beta1))
  rl$importance_score <- ifelse(beta1_mean == 0, 1e-10, beta1_mean)
  
  # Compute selection probability
  rl$select_prob <- rl$importance_score / sum(rl$importance_score)
  
  # Procedure 2: Bootstrapping for Coefficients
  cat("Procedure 2: Bootstrapping for Coefficients\n")
  proc2_result <- bootstrapping(rl, mode = "procedure2")
  beta2 <- proc2_result$betas
  intercepts <- proc2_result$intercepts
  rl$intercept_ <- mean(intercepts)  # Use the mean of tmp_intercept values
  
  beta2_mean <- rowMeans(beta2)
  rl$coef_ <- ifelse(abs(beta2_mean) > rl$threshold, beta2_mean, 0)
  
  # Assign names to coefficients
  names(rl$coef_) <- paste0("X", 1:rl$p)
  # Combine intercept and coefficients
  rl$coefficients <- c(Intercept = rl$intercept_, rl$coef_)
  
  return(rl)
}

# Bootstrapping function
bootstrapping <- function(rl, mode) {
  if (mode == 'procedure1') {
    q <- rl$q1
    method <- 'ElasticNet'
  } else {
    q <- rl$q2
    method <- 'AdaptiveLASSO'
  }
  
  # Initialize progress bar
  pb <- progress_bar$new(
    format = "  [:bar] :percent eta: :eta",
    total = rl$B, clear = FALSE, width=60
  )
  
  # Define the coefficient matrix
  betas <- matrix(0, nrow = rl$p, ncol = rl$B)
  rownames(betas) <- paste0("X", 1:rl$p)
  
  # Store intercepts
  intercepts <- numeric(rl$B)
  
  # Sequential bootstrapping
  for (b in 1:rl$B) {
    result <- estimate_coef(rl, bootstrap_number = b, q = q, method = method)
    betas[, b] <- result$beta
    intercepts[b] <- result$tmp_intercept  # Store tmp_intercept for each iteration
    pb$tick()
  }
  
  list(betas = betas, intercepts = intercepts)
}

# Function to estimate coefficients for each bootstrap sample
estimate_coef <- function(rl, bootstrap_number, q, method) {
  # Initialize beta vector
  beta <- rep(0, rl$p)
  
  # Set random seed for reproducibility
  if (!is.null(rl$random_state)) {
    set.seed(rl$random_state + bootstrap_number)
  } else {
    set.seed(NULL)
  }
  
  # Generate bootstrap sample indices
  bst_sample_idx <- sample(1:rl$n, size = rl$n, replace = FALSE)
  
  # Generate predictor indices based on select_prob
  bst_predictor_idx <- sample(1:rl$p, size = q, replace = FALSE, prob = rl$select_prob)
  
  # Subset X and y
  X_bst <- rl$X[bst_sample_idx, bst_predictor_idx, drop = FALSE]
  y_bst <- rl$y[bst_sample_idx]
  
  # Standardize X and y
  std <- standardization(X_bst, y_bst)
  X_sc <- std$X_sc
  y_sc <- std$y_sc
  X_sd <- std$X_sd
  
  # Fit the model
  if (method == 'ElasticNet') {
    # Elastic Net with alpha=0.5
    fit <- tryCatch({
      glmnet(X_sc, y_sc, alpha = 0.5, weights = rl$sample_weight[bst_sample_idx],
             standardize = FALSE, intercept = TRUE)
    }, error = function(e) { NULL })
    
    if (!is.null(fit)) {
      # Select lambda that minimizes cross-validated error
      cv_fit <- tryCatch({
        cv.glmnet(X_sc, y_sc, alpha = 0.5, weights = rl$sample_weight[bst_sample_idx],
                  standardize = FALSE, intercept = TRUE)
      }, error = function(e) { NULL })
      
      if (!is.null(cv_fit)) {
        coefs <- as.numeric(coef(cv_fit, s = "lambda.min"))
        tmp_intercept <- coefs[1]
        coef_fit <- coefs[-1]
      } else {
        coef_fit <- rep(0, length(bst_predictor_idx))
      }
    } else {
      coef_fit <- rep(0, length(bst_predictor_idx))
    }
    
  } else if (method == 'AdaptiveLASSO') {
    # Adaptive Lasso with penalty.factor
    weights_adaptive <- 1 / rl$importance_score[bst_predictor_idx]
    fit <- tryCatch({
      glmnet(X_sc, y_sc, alpha = 1, penalty.factor = weights_adaptive, 
             weights = rl$sample_weight[bst_sample_idx],
             standardize = FALSE, intercept = FALSE)
    }, error = function(e) { NULL })
    
    if (!is.null(fit)) {
      # Select lambda that minimizes cross-validated error
      cv_fit <- tryCatch({
        cv.glmnet(X_sc, y_sc, alpha = 1, penalty.factor = weights_adaptive, 
                  weights = rl$sample_weight[bst_sample_idx],
                  standardize = FALSE, intercept = FALSE)
      }, error = function(e) { NULL })
      
      if (!is.null(cv_fit)) {
        coefs <- as.numeric(coef(cv_fit, s = "lambda.min"))
        tmp_intercept <- coefs[1]
        coef_fit <- coefs[-1]
      } else {
        coef_fit <- rep(0, length(bst_predictor_idx))
      }
    } else {
      coef_fit <- rep(0, length(bst_predictor_idx))
    }
  }
  
  # Assign coefficients back to the beta vector
  if (length(coef_fit) != length(bst_predictor_idx)) {
    # In case of fitting failure, skip assigning
    return(list(beta = beta, tmp_intercept = tmp_intercept))
  }

  beta[bst_predictor_idx] <- coef_fit / X_sd
  return(list(beta = beta, tmp_intercept = tmp_intercept))
}

# Print method for RandomLasso
print.RandomLasso <- function(rl, ...) {
  cat("RandomLasso Model\n")
  cat("Parameters:\n")
  cat("  q1 =", rl$q1, "\n")
  cat("  q2 =", rl$q2, "\n")
  cat("  B =", rl$B, "\n")
  cat("  random_state =", ifelse(is.null(rl$random_state), "NULL", rl$random_state), "\n")
  cat("  L =", rl$L, "\n\n")
  
  if (!is.null(rl$coef_)) {
    cat("Model has been fitted.\n")
    cat("Number of samples (n):", rl$n, "\n")
    cat("Number of predictors (p):", rl$p, "\n\n")
    
    cat("Intercept and Coefficients:\n")
    print(rl$coefficients)
  } else {
    cat("Model has not been fitted yet.\n")
  }
}

# Example usage
# Uncomment and run the following lines to test the RandomLasso implementation

set.seed(123)   # For reproducibility

# Step 1: Generate a synthetic dataset
n <- 50  # Number of observations
p <- 100   # Number of predictors
X <- matrix(rnorm(n * p), n, p) # Predictor matrix with n rows and p columns

# Define the true beta vector with 6 non-zero coefficients and 14 zeros
beta <- c(3, 2, 1.2, 0.5, 0.2, 1.5, 5, 0.3, 2.5, 4, rep(0, p - 10)) # 6 non-zero coefficients + 14 zeros = 20
inter <- 0.2
y <- inter + X %*% beta + rnorm(n)  # Response with noise
threshold <- 0.1 # Threshold for selecting predictors

# Initialize the RandomLasso object with the threshold
random_lasso <- RandomLasso(q1 = 50, q2 = 50, B = 'auto', 
                            random_state = 123, L = 50, threshold = threshold)
start_time <- Sys.time()
# Fit the model
random_lasso <- fit.RandomLasso(random_lasso, X, y, sample_weight = NULL)
end_time <- Sys.time()
# Print the model summary
print(random_lasso)
processing_time <- end_time - start_time

# Access the coefficients
coefficients <- random_lasso$coefficients
print(coefficients)

