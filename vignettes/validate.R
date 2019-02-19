## ---- echo = FALSE-------------------------------------------------------
library(beset)
suppressPackageStartupMessages(library(ggplot2))

## ------------------------------------------------------------------------
lin_mod <- lm(Fertility ~ ., data = swiss)
summary(lin_mod)

## ------------------------------------------------------------------------
cv_results <- validate(lin_mod)
cv_results

## ------------------------------------------------------------------------
null_model <- lm(Fertility ~ 1, data = swiss)
mce <- as.numeric(-logLik(null_model) / nobs(null_model))

## ------------------------------------------------------------------------
log_mod <- glm(tumor ~ ., data = prostate, family = "binomial")
summary(log_mod)
validate(log_mod)

## ------------------------------------------------------------------------
validate(log_mod, n_folds = 5, n_reps = 5)

## ------------------------------------------------------------------------
validate(log_mod, n_folds = nobs(log_mod), n_reps = 1)

## ------------------------------------------------------------------------
cv_results <- validate(lin_mod)

## ------------------------------------------------------------------------
cv_results$predictions

## ------------------------------------------------------------------------
cv_results$fold_assignments

