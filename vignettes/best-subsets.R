## ---- echo = FALSE-------------------------------------------------------
suppressPackageStartupMessages(library(beset))
suppressPackageStartupMessages(library(ggplot2))

## ------------------------------------------------------------------------
data("swiss", package = "datasets")
train_data <- swiss
set.seed(42)
train_data <- cbind(train_data, 
                    matrix(replicate(5, rnorm(nrow(train_data))), ncol = 5))
names(train_data)[7:11] <- paste0("noise", names(train_data)[7:11])

## ---- eval = FALSE-------------------------------------------------------
#  mod <- beset_lm(Fertility ~ ., train_data)

## ---- eval = FALSE-------------------------------------------------------
#  mod <- beset_lm(Fertility ~ ., train_data, test_data = NULL, p_max = 10,
#                  n_cores = 2, n_folds = 10, n_repeats = 10, seed = 42)

## ---- eval = FALSE-------------------------------------------------------
#  mod <- beset_glm(Fertility ~ ., train_data)

## ------------------------------------------------------------------------
mod <- beset_glm(Fertility ~ ., train_data, test_data = NULL, p_max = 10,
                 family = "gaussian", link = "identity", n_cores = 2, 
                 n_folds = 10, n_repeats = 10, seed = 42)

## ---- fig.height=4, fig.width=5------------------------------------------
plot(mod)

## ---- fig.height=4, fig.width=5------------------------------------------
plot(mod, "r2")

## ---- eval = FALSE-------------------------------------------------------
#  summary(mod)

## ------------------------------------------------------------------------
summary(mod, metric = "mse", n_pred = NULL, oneSE = TRUE, n_cores = 2)

## ------------------------------------------------------------------------
summary(mod, oneSE = FALSE)

## ------------------------------------------------------------------------
summary(mod, n_pred = 3)

## ------------------------------------------------------------------------
summary(mod, metric = "r2")

## ------------------------------------------------------------------------
data("Boston", package = "MASS")

## ------------------------------------------------------------------------
Boston.1 <- partition(data = Boston, y = medv, p = 0.1, seed = 42)
Boston.3 <- partition(data = Boston, y = medv, p = 0.3, seed = 42)
Boston.5 <- partition(data = Boston, y = medv, p = 0.5, seed = 42)

## ------------------------------------------------------------------------
mod.1 <- beset_lm(medv ~ ., Boston.1$train, Boston.1$test)
mod.3 <- beset_lm(medv ~ ., Boston.3$train, Boston.3$test)
mod.5 <- beset_lm(medv ~ ., Boston.5$train, Boston.5$test)

## ---- eval = FALSE-------------------------------------------------------
#  mod.5 <- beset_lm(medv ~ ., Boston.5$train, Boston.5$test, p_max = 13)

## ---- fig.height=4, fig.width=5------------------------------------------
plot(mod.1) + ggtitle("Train sample = 52, Test sample = 454")
plot(mod.3) + ggtitle("Train sample = 154, Test sample = 352")
plot(mod.5) + ggtitle("Train sample = 254, Test sample = 252")

## ------------------------------------------------------------------------
summary(mod.1)

## ------------------------------------------------------------------------
summary(mod.5)

## ------------------------------------------------------------------------
data("biopsy", package = "MASS")
biopsy <- biopsy[,-1]

## ------------------------------------------------------------------------
biopsy <- partition(data = biopsy, y = class)

## ------------------------------------------------------------------------
mod <- beset_glm(class ~ ., biopsy$train, biopsy$test, family = "binomial")

## ------------------------------------------------------------------------
summary(mod)

## ---- eval = FALSE-------------------------------------------------------
#  # Incorrect syntax
#  beset_glm(class ~ ., biopsy$train, biopsy$test,
#            family = binomial("probit"))

## ---- eval = FALSE-------------------------------------------------------
#  # Correct syntax
#  beset_glm(class ~ ., biopsy$train, biopsy$test,
#            family = "binomial", link = "probit")

## ---- fig.height=4, fig.width=5------------------------------------------
plot(mod)

## ---- fig.height=4, fig.width=5------------------------------------------
plot(mod) + ylab("Log-loss")

## ---- fig.height=4, fig.width=5------------------------------------------
plot(mod, metric = "r2")

## ------------------------------------------------------------------------
data("bioChemists", package = "pscl")
qplot(art, data = bioChemists, binwidth = 1) + theme_bw()

## ------------------------------------------------------------------------
bioChemists <- partition(bioChemists, art)
mod <- beset_glm(art ~ ., bioChemists$train, bioChemists$test,
                 family = "poisson")
summary(mod)

## ------------------------------------------------------------------------
mod <- beset_glm(art ~ ., bioChemists$train, bioChemists$test,
                 family = "negbin")
summary(mod)

## ---- fig.height=4, fig.width=5------------------------------------------
plot(mod, metric = "r2") + 
  ggtitle("Cross-validated subsets for
negative binomial regression")

## ------------------------------------------------------------------------
mod <- beset_zeroinfl(art ~ ., bioChemists$train, bioChemists$test)

## ---- eval = FALSE-------------------------------------------------------
#  mod <- beset_zeroinfl(art ~ ., bioChemists$train, bioChemists$test,
#                        family = "poisson", link = "logit", p_count_max = 10,
#                        p_zero_max = 10, n_cores = 2, n_folds = 10,
#                        n_repeats = 10, seed = 42)

## ------------------------------------------------------------------------
summary(mod)

## ------------------------------------------------------------------------
summary(mod, metric = "aic")

## ------------------------------------------------------------------------
summary(mod, oneSE = FALSE)

## ------------------------------------------------------------------------
summary(mod, n_count_pred = 2, n_zero_pred = 2)

## ---- fig.height=4, fig.width=5------------------------------------------
plot(mod, type = "train")

## ---- fig.height=4, fig.width=5------------------------------------------
plot(mod, type = "train", metric = "mse")

## ---- fig.height=4, fig.width=5------------------------------------------
plot(mod, type = "cv", metric = "r2")

## ---- fig.height=4, fig.width=5------------------------------------------
plot(mod, type = "cv", metric = "r2", se = FALSE)

## ------------------------------------------------------------------------
mod <- beset_zeroinfl(art ~ ., bioChemists$train, bioChemists$test,
                      family = "negbin", p_count_max = 3, p_zero_max = 2)
summary(mod)

## ---- fig.height=4, fig.width=5------------------------------------------
plot(mod, type = "cv", metric = "r2")

