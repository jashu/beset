library(beset)
context("create_folds")

y <- c(rep(0, 10), rep(1:5, each = 2))
seeds <- 1:1000
fold_ids <- map(seeds, ~ beset:::create_folds(y, 10, 1, seed = .x))
has_empty_fold <- map_lgl(fold_ids, ~ any(map_lgl(.x, is_empty)))

test_that("there are no empty folds without any hold-out cases", {
 expect_false(any(has_empty_fold))
})

y <- pscl::bioChemists$art
fold_ids <- map(seeds, ~ beset:::create_folds(y, 10, 1, seed = .x)) %>%
  reduce(c)
n_0s <- map_int(fold_ids, ~ sum(y[.x] == 0))
expected_n_0s <- sum(y == 0) / 10
test_that("all folds contain approximately equal numbers of 0s",{
  expect_true(all(between(n_0s, expected_n_0s - 1, expected_n_0s + 1)))
})
