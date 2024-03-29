% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ict.hausman.test.R
\name{ict.hausman.test}
\alias{ict.hausman.test}
\title{Hausman Specification Test for Two List Experiment Regression Fit Objects}
\usage{
ict.hausman.test(ml, nls, abs = FALSE, psd = FALSE)
}
\arguments{
\item{ml}{Maximum likelihood model fit, for example from ictreg(method = "ml")}

\item{nls}{NLS model fit, for example from ictreg(method = "nls")}

\item{abs}{Flag to override error when Hausman statistic is negative, which may indicate misspecification. Set to \code{FALSE} to override.}

\item{psd}{Flag to override error when variance-covariance difference is non-positive semidefinite, suggesting misspecification.  Set to \code{TRUE} to override.}
}
\value{
List containing Hausman statistic, p-value, and the degrees of freedom of the test.
}
\description{
Hausman Specification Test for Two List Experiment Regression Fit Objects
}
