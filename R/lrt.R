#' @export
#'
lrt <- function(object1, object2){
  resid_df_1 <- object1$df.residual
  resid_df_2 <- object2$df.residual
  resid_dev_1 <- -2 * logLik(object1)
  resid_dev_2 <- -2 * logLik(object2)
  ii <- order(c(resid_df_1, resid_df_2), decreasing = TRUE)
  resid_df <- c(resid_df_1, resid_df_2)[ii]
  resid_dev <- c(resid_dev_1, resid_dev_2)[ii]
  df <- c(NA, resid_df[1] - resid_df[2])
  dev <- c(NA, resid_dev[1] - resid_dev[2])
  p <- c(NA, pchisq(dev[2], df = df[2], lower.tail=F))
  structure(data.frame(`Resid. Df` = resid_df,
                       `Resid. Dev` = resid_dev,
                       `Df` = df,
                       `Deviance` = dev,
                       `Pr(>Chi)` = p,
                       check.names = FALSE),
            class = c("anova", "data.frame"))
}
