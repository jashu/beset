## ---- echo = FALSE, message=FALSE---------------------------------------------
library(beset)
suppressPackageStartupMessages(library(tidyverse))

## -----------------------------------------------------------------------------
set.seed(42)
train_data <- cbind(swiss, 
                    matrix(replicate(5, rnorm(nrow(swiss))), ncol = 5))
names(train_data)[7:11] <- paste0("noise", names(train_data)[7:11])

## -----------------------------------------------------------------------------
mod <- beset_lm(Fertility ~ ., train_data)

## ---- fig.height=4, fig.width=5-----------------------------------------------
plot(mod)

## ---- fig.height=4, fig.width=5-----------------------------------------------
plot(mod, "rsq")

## -----------------------------------------------------------------------------
summary(mod)

## ---- eval = FALSE------------------------------------------------------------
#  summary(mod, n_pred = NULL, metric = "mse", oneSE = TRUE)

## -----------------------------------------------------------------------------
summary(mod, oneSE = FALSE)

## -----------------------------------------------------------------------------
summary(mod, metric = "aic")

## -----------------------------------------------------------------------------
summary(mod, metric = "mae")

## ---- fig.height=4, fig.width=5-----------------------------------------------
plot(mod)

## -----------------------------------------------------------------------------
summary(mod, n_pred = 4)

## -----------------------------------------------------------------------------
data <- partition(train_data, y = "Fertility", seed = 42, frac = .75)

## ---- fig.height=4, fig.width=5-----------------------------------------------
mod <- beset_lm(Fertility ~ ., data = data)
plot(mod)

## ---- fig.height=4, fig.width=5-----------------------------------------------
mod <- beset_lm(Fertility ~ ., data = train_data, nest_cv = TRUE, p_max = 5)
plot(mod)

## ---- fig.height=4, fig.width=5-----------------------------------------------
plot(mod, se = FALSE)

## -----------------------------------------------------------------------------
summary(mod)

## -----------------------------------------------------------------------------
summary(mod) %>% print(standardize = FALSE)

## -----------------------------------------------------------------------------
summary(mod, oneSE = FALSE)

## -----------------------------------------------------------------------------
summary(mod) %>% print(metric = "mse")

## -----------------------------------------------------------------------------
validate(mod, metric = "mse", oneSE = TRUE)

## -----------------------------------------------------------------------------
mod <- beset_lm(Fertility ~ Agriculture + Examination + Education +
                  Catholic + Infant.Mortality, train_data)
summary(mod)

## ---- eval = FALSE------------------------------------------------------------
#  ### DO NOT RUN!
#  mod <- beset_lm(Fertility ~ Agriculture * Examination * Education *
#                    Catholic * Infant.Mortality, train_data)

## -----------------------------------------------------------------------------
mod <- beset_lm(Fertility ~ Education * Catholic * Infant.Mortality,
                train_data)

## -----------------------------------------------------------------------------
summary(mod)

## -----------------------------------------------------------------------------
mod <- beset_lm(Fertility ~ Education + Catholic + Infant.Mortality +
                  Education:Catholic + Education:Infant.Mortality +
                  Catholic:Infant.Mortality, train_data, nest_cv = TRUE,
                force_in = c("Education", "Catholic", "Infant.Mortality"))
summary(mod)

## -----------------------------------------------------------------------------
mod <- beset_lm(Fertility ~ ., train_data, n_folds = 5, n_reps = 5)
summary(mod)

## ---- fig.height=4, fig.width=5-----------------------------------------------
plot(mod)

## -----------------------------------------------------------------------------
data("prostate")
summary(prostate)

## -----------------------------------------------------------------------------
mod <- beset_glm(tumor ~ ., data = prostate, family = "binomial")

## -----------------------------------------------------------------------------
summary(mod)

## ---- eval = FALSE------------------------------------------------------------
#  # Incorrect syntax
#  beset_glm(tumor ~ ., data = prostate, family = binomial("probit"))

## ---- eval = FALSE------------------------------------------------------------
#  # Correct syntax
#  beset_glm(tumor ~ ., data = prostate, family = "binomial",
#            link = "probit") %>% summary()

## ---- fig.height=4, fig.width=5-----------------------------------------------
plot(mod)

## ---- fig.height=4, fig.width=5-----------------------------------------------
plot(mod) + ylab("Log-loss")

## ---- fig.height=4, fig.width=5-----------------------------------------------
plot(mod, metric = "auc")

## -----------------------------------------------------------------------------
summary(mod, metric = "auc")

## -----------------------------------------------------------------------------
data("adolescents")
qplot(x = dep, bins = 10, data = na.omit(adolescents)) + theme_classic()

## -----------------------------------------------------------------------------
mod_poisson <- beset_glm(dep ~ ., data = adolescents, family = "poisson")

## ---- fig.height=4, fig.width=5-----------------------------------------------
plot(mod_poisson)

## -----------------------------------------------------------------------------
poisson_summary <- summary(mod_poisson, oneSE = FALSE)
poisson_summary

## -----------------------------------------------------------------------------
mod_negbin <- beset_glm(dep ~ ., data = adolescents, family = "negbin")

## ---- fig.height=4, fig.width=5-----------------------------------------------
plot(mod_negbin) 

## -----------------------------------------------------------------------------
negbin_summary <- summary(mod_negbin, oneSE = FALSE)
negbin_summary

