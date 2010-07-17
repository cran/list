\name{ictreg} 
\alias{ictreg} 
\alias{ict}
\alias{list}
\title{Item Count Technique} 
\description{ 
Function to conduct multivariate regression analyses of survey data with the item count technique, also known as the list experiment and the unmatched count technique.
} 
\usage{
ictreg(formula, data, treat = "treat", J, method = "nls", 
       overdispersed = FALSE, constrained = TRUE, maxIter = 5000, verbose = FALSE, ...)
} 

\arguments{ 
\item{formula}{An object of class "formula": a symbolic description of the model to be fitted.} 
\item{data}{A data frame containing the variables in the model}
\item{treat}{Name of binary treatment indicator.}
\item{J}{Number of non-sensitive (control) survey items.}
\item{method}{ Method for regression, either \code{ml} for the Maximum Likelihood (ML) estimation with the Expectation-Maximization algorithm or \code{nls} for the Non-linear Least Squares (NLS) estimation with the two-step procedure.}
\item{overdispersed}{Indicator for the presence of overdispersion. If \code{TRUE}, the beta-binomial model is used in the EM algorithm, if \code{FALSE} the binomial model is used. Not relevant for the NLS method.}
\item{constrained}{a logical value indicating whether the control group parameters are constrained to be equal.  Not relevant for the NLS method}
\item{maxIter}{Maximum number of iterations for the Expectation-Maximization algorithm of the ML estimation.  The default is 5000.} 
\item{verbose}{ a logical value indicating whether model diagnostics are printed out during fitting.}
\item{...}{further arguments to be passed to NLS regression commands.}
}

\details{This function allows the user to perform regression analysis on data from the item count technique, also known as the list experiment and the unmatched count technique. 

Two methods are implemented in this function, the Maximum Liklihoood (ML) estimation for the Expectation-Maximization (EM) algorithm and the nonlinear least squares (NLS) estimation with the two-step procedure both proposed in Imai (2010). The constrained and unconstrained ML models presented in Imai (2010) are available through the \code{constrained} option, and the user can specify if overdispersion is present in the data using the \code{overdispersed} option to control whether a beta-binomial or binomial model is used in the EM algorithm to model the item counts.}


\value{\code{ictreg} returns an object of class "ictreg".  The function \code{summary} is used to obtain a table of the results.  The object \code{ictreg} is a list that contains the following components.  Some of these elements are not available depending on which method is used (\code{nls} or \code{ml}) and whether the constrained model is used (\code{constrained = TRUE}).

\item{par.treat}{point estimate for effect of covariate on item count fitted on treatment group}
\item{se.treat}{standard error for estimate of effect of covariate on item count fitted on treatment group}
\item{par.control}{point estimate for effect of covariate on item count fitted on control group}
\item{se.control}{standard error for estimate of effect of covariate on item count fitted on control group}
\item{coef.names}{variable names as defined in the data frame}
\item{method}{call of the method used}
\item{overdispersed}{call indicating whether data is overdispersed}
\item{constrained}{call indicating whether the constrained model is used}

For the unconstrained model, the \code{par.control} and \code{se.control} output is replaced by:
\item{par.control.phi0}{point estimate for effect of covariate on item count fitted on treatment group}
\item{se.control.phi0}{standard error for estimate of effect of covariate on item count fitted on treatment group}
\item{par.control.phi1}{point estimate for effect of covariate on item count fitted on treatment group}
\item{se.control.phi1}{standard error for estimate of effect of covariate on item count fitted on treatment group}
}

\references{Imai, Kosuke. (2010) ``Statistical Inference for the Item
  Count Technique.''  Technical Report, Department of Politics,
  Princeton University. available at \url{http://imai.princeton.edu/research/list.html}}

\author{Graeme Blair, Princeton University,
\email{gblair@princeton.edu} and Kosuke Imai, Princeton University,
\email{kimai@princeton.edu} 
}

\examples{

data(race)

# Fit two-step non-linear least squares regression
# Replicates Table 1 Columns 1-2 Imai (2010); note that age is divided by 10

nls.results <- ictreg(y ~ south + age + male + college, data = race, treat = "treat", J=3, method = "nls")

summary(nls.results)

\dontrun{

# Fit EM algorithm ML model with constraint
# Replicates Table 1 Columns 3-4, Imai (2010); note that age is divided by 10

ml.constrained.results <- ictreg(y ~ south + age + male + college, data = race, treat = "treat", J=3, method = "ml", overdispersed = FALSE, constrained = TRUE)

summary(ml.constrained.results)

# Fit EM algorithm ML model with no constraint
# Replicates Table 1 Columns 5-8, Imai (2010); not that age is divided by 10

ml.unconstrained.results <- ictreg(y ~ south + age + male + college, data = race, treat = "treat", J=3, method = "ml", overdispersed = FALSE, constrained = FALSE)

summary(ml.unconstrained.results)

}
}
 
\seealso{
\code{\link{predict.ictreg}} for fitted values
}

\keyword{models}
\keyword{regression}