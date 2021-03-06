% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{prostate}
\alias{prostate}
\title{Prostate Cancer Study}
\format{
A data frame with 380 rows and 8 variables:
\describe{
  \item{tumor}{tumor penetration of prostatic capsule, 0 = no penetration}
  \item{age}{in years}
  \item{race}{white or black}
  \item{dpros}{results of the digital rectal exam, 4 levels}
  \item{dcaps}{detection of capsular involvement in rectal exam (yes, no)}
  \item{psa}{antigen (mg/ml)}
  \item{vol}{tumor volume obtained from ultrasound (cm3)}
  \item{gleason}{total gleason score (0 - 10)}
}
}
\source{
\url{https://raw.github.com/0xdata/h2o/master/smalldata/logreg/prostate.csv}
}
\usage{
prostate
}
\description{
Baseline exam results on prostate cancer patients from Dr. Donn Young at The
Ohio State University Comprehensive Cancer Center.
}
\references{
Hosmer and Lemeshow (2000) Applied Logistic Regression: Second
Edition.
}
\keyword{datasets}
