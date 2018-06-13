## ---- echo = FALSE-------------------------------------------------------
suppressPackageStartupMessages(library(beset))
suppressPackageStartupMessages(library(ggplot2))

## ------------------------------------------------------------------------
library(ElemStatLearn)
data(prostate)
prostate$gleason <- factor(prostate$gleason)
qplot(gleason, lpsa, data = prostate, geom = "boxplot") + geom_jitter()

## ------------------------------------------------------------------------
prostate$gleason <- forcats::fct_collapse(
  prostate$gleason, gleason_6 = "6", gleason_above_6 = c("7", "8", "9"))
qplot(gleason, lpsa, data = prostate, geom = "boxplot") + geom_jitter()

## ------------------------------------------------------------------------
lin_mod <- lm(lpsa ~ gleason, data = prostate)
summary(lin_mod)

## ------------------------------------------------------------------------
validate(lin_mod)

## ------------------------------------------------------------------------
null_model <- lm(lpsa ~ 1, data = prostate)
-logLik(null_model) / nobs(null_model)

## ------------------------------------------------------------------------
log_mod <- glm(gleason ~ lpsa, data = prostate, family = "binomial")
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

