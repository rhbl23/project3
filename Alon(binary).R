library(datamicroarray)
# more data in https://github.com/ramhiser/datamicroarray
data('alon', package = 'datamicroarray')
X <- alon$x
y <- factor(alon$y)


#lasso
cv_lasso <- cv.glmnet(X, y, alpha = 1, family = "binomial")
best_lambda <- cv_lasso$lambda.min
lasso_coefficients <- coef(cv_lasso, s = "lambda.min")
coef_matrix <- as.matrix(lasso_coefficients)

coef_lasso <- coef_lasso[coef_lasso$Gene != "(Intercept)", ]

# Add a column for the absolute value of coefficients
coef_lasso$AbsCoefficient <- abs(coef_lasso$Coefficient)
coef_lasso <- coef_lasso[order(coef_lasso$AbsCoefficient, decreasing = TRUE), ]
top_genes_lasso <- head(coef_lasso, 10)
print(top_genes_lasso)


#stabs
library(stabs)
stabs_lasso <- stabsel(
  x = X,
  y = y,
  fitfun = glmnet.lasso,
  args.fitfun = list(family = "binomial"), # Logistic regression
  cutoff = 0.75,  # Selection frequency threshold
  PFER = 1        # Per-family error rate
)

stabs_prob <- as.matrix(stabs_lasso$max)

prob <- data.frame(
  Gene = rownames(stabs_prob),
  prbability = as.numeric(stabs_prob)
)
prob<- prob[order(prob$prbability, decreasing = TRUE), ]
top_genes_stabs <- head(prob, 10)
print(top_genes_stabs)


#randomlasso
random_lasso <- RandomLasso(q1 = 1000, q2 = 1000, B = 'auto', 
                            random_state = 123, L = 50, threshold = 1/500, family="binomial")
random_lasso <- fit.RandomLasso(random_lasso, X, y, sample_weight = NULL)
random_coefficients <- as.matrix(random_lasso$coefficients)
rownames(random_coefficients) <- rownames(lasso_coefficients)

coef_random <- data.frame(
  Gene = rownames(random_coefficients),
  Coefficient = as.numeric(random_coefficients)
)
coef_random$AbsCoefficient <- abs(coef_random$Coefficient)
coef_random <- coef_random[order(coef_random$AbsCoefficient, decreasing = TRUE), ]
top_genes_random <- head(coef_random, 10)
print(top_genes_random)


#hilasso
hi_lasso <- HiLasso(q1 = 1500, q2 = 1500,
                    random_state = 123, L = 30, family="binomial")
hi_lasso <- fit.HiLasso(hi_lasso, X, y, sample_weight = NULL)
hi_coefficients <- as.matrix(hi_lasso$coefficients)

coef_hi <- data.frame(
  Gene = rownames(lasso_coefficients),
  Coefficient = as.numeric(hi_coefficients)
)

coef_hi <- coef_hi [coef_hi $Gene != "(Intercept)", ]
coef_hi$AbsCoefficient <- abs(coef_hi$Coefficient)
coef_hi<- coef_hi[order(coef_hi$AbsCoefficient, decreasing = TRUE), ]
top_genes_hi <- head(coef_hi, 10)
print(top_genes_hi)



