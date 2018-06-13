## ---- echo = FALSE-------------------------------------------------------
library(beset)
library(tidyverse)
library(ElemStatLearn)
data(prostate)
summary(prostate)

## ------------------------------------------------------------------------
train <- prostate %>% filter(train) %>% select(-train)
test <- prostate %>% filter(!train) %>% select(-train)
prostate <- prostate %>% select(-train)
data <- data_partition(train, test, y = "lpsa")

## ------------------------------------------------------------------------
mod <- beset_elnet(lpsa ~ ., data)

## ---- fig.height=4, fig.width=5------------------------------------------
plot(mod)

## ------------------------------------------------------------------------
summary(mod, oneSE = FALSE)

## ------------------------------------------------------------------------
summary(mod)

## ------------------------------------------------------------------------
summary(mod, alpha = 0.01)

## ------------------------------------------------------------------------
range(train$lpsa)
range(test$lpsa)

## ------------------------------------------------------------------------
var(train$lpsa)
var(test$lpsa)

## ---- fig.height=4, fig.width=5------------------------------------------
plot(mod, "rsq")

## ------------------------------------------------------------------------
set.seed(42)
noise <- sample_frac(prostate) %>% select(-lpsa)
names(noise) <- paste("noise", 1:ncol(noise), sep = "")
prostate <- bind_cols(prostate, noise)

## ---- fig.height=4, fig.width=5------------------------------------------
mod <- beset_elnet(lpsa ~ ., prostate, nest_cv = TRUE)
plot(mod, "rsq")

## ------------------------------------------------------------------------
summary(mod)

## ------------------------------------------------------------------------
summary(mod, oneSE = FALSE)

## ------------------------------------------------------------------------
summary(mod) %>% print(metric = "mse")

## ------------------------------------------------------------------------
validate(mod, metric = "mce", oneSE = TRUE, alpha = NULL, lambda = NULL)

## ------------------------------------------------------------------------
prostate <- read_csv("https://raw.github.com/0xdata/h2o/master/smalldata/logreg/prostate.csv") %>% 
  select(-ID) %>%
  mutate(RACE = factor(RACE, levels = c(1,2), labels = c("white", "black")),
         DCAPS = factor(DCAPS, labels = c("no", "yes")))
summary(prostate)

## ------------------------------------------------------------------------
mod <- beset_elnet(CAPSULE ~ ., data = prostate, family = "binomial",
                   nest_cv = TRUE)
summary(mod)

## ---- fig.height=4, fig.width=5------------------------------------------
plot(mod)

## ---- fig.height=4, fig.width=5------------------------------------------
plot(mod) + ylab("Log-loss")

## ---- fig.height=4, fig.width=5------------------------------------------
plot(mod, metric = "auc")

## ------------------------------------------------------------------------
summary(mod, metric = "auc") %>% print(metric = "auc")

