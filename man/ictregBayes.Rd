\name{ictregBayes} 
\alias{ictregBayes}
\title{Item Count Technique} 
\description{ 
Function to conduct multivariate regression analyses of survey data with the item count technique, also known as the list experiment and the unmatched count technique.
} 
\usage{
ictregBayes(formula, data = parent.frame(), treat = "treat", J, constrained = TRUE,
            fit.start = "lm", n.draws = 10000, burnin = 5000, 
            thin = 100, delta.start, psi.start, delta.mu0, psi.mu0, delta.A0, psi.A0, 
            delta.tune, psi.tune, verbose = FALSE, ...)
} 

\arguments{ 
\item{formula}{An object of class "formula": a symbolic description of the model to be fitted.} 
\item{data}{A data frame containing the variables in the model}
\item{treat}{Name of treatment indicator as a string. The treatment indicator can be a numerical or factor variable, and must contain only two conditions, treatment and control (which must be indicated be zero if numeric).}
\item{J}{Number of non-sensitive (control) survey items. This will be set automatically to the maximum value of the outcome variable in the treatment group if no input is sent by the user.}
\item{constrained}{A logical value indicating whether the control group parameters are constrained to be equal}
\item{fit.start}{Fit method for starting values. The options are \code{lm}, \code{glm}, \code{nls}, and \code{ml}, which use OLS, logistic regression, non-linear least squares, and maximum likelihood estimation to generate starting values, respectively. The default is \code{lm}.}
\item{n.draws}{Number of MCMC iterations after the burnin.}
\item{burnin}{The number of initial MCMC iterations that are discarded.}
\item{thin}{The interval of thinning, in which every other (\code{thin} = 1) or more iterations are discarded in the output object}
\item{delta.start}{Optional starting values for the sensitive item fit. This should be a vector with the length of the number of covariates. The default runs an \code{ictreg} fit with the method set by the \code{fit.start} option.}
\item{psi.start}{Optional starting values for the control items fit. This should be a vector of length the number of covariates. The default runs an \code{ictreg} fit with the method set by the \code{fit.start} option.}
\item{delta.mu0}{Optional vector of prior means for the sensitive item fit parameters, a vector of length the number of covariates.}
\item{psi.mu0}{Optional vector of prior means for the control item fit parameters, a vector of length the number of covariates.}
\item{delta.A0}{Optional matrix of prior precisions for the sensitive item fit parameters, a matrix of dimension the number of covariates.}
\item{psi.A0}{Optional matrix of prior precisions for the control items fit parameters, a matrix of dimension the number of covariates.}
\item{delta.tune}{A required vector of tuning parameters for the Metropolis algorithm for the sensitive item fit. This must be set and refined by the user until the acceptance ratios are approximately .4 (reported in the output).}
\item{psi.tune}{A required vector of tuning parameters for the Metropolis algorithm for the control item fit. This must be set and refined by the user until the acceptance ratios are approximately .4 (reported in the output).}
\item{verbose}{A logical value indicating whether model diagnostics are printed out during fitting.}
\item{...}{further arguments to be passed to NLS regression commands.}
}

\details{This function allows the user to perform regression analysis on data from the item count technique, also known as the list experiment and the unmatched count technique using a Bayesian MCMC algorithm.

In the first version (4.0), only the standard design list experiment regression is estimated with a single sensitive item.

Convergence is at times difficult to achieve, so we recommend running multiple chains from overdispersed starting values by, for example, running an MLE or linear model using the ictreg() function, and then generating a set of overdispersed starting values using those estimates and their estimated variance-covariance matrix. An example is provided below.

}

\value{\code{ictregBayes} returns an object of class "ictregBayes".  The function \code{summary} is used to obtain a table of the results, using the \code{coda} package. Two attributes are also included, the data ("x"), the call ("call"), which can be extracted using the command, e.g., attr(ictregBayes.object, "x").

}

\references{
Blair, Graeme and Kosuke Imai. (2010) ``Statistical Analysis of List Experiments." 
  Technical Report, Department of Politics, Princeton University. 
  available at \url{http://imai.princeton.edu/research/listP.html}

Imai, Kosuke. (2011) ``Multivariate Regression Analysis for the Item
  Count Technique.''  Journal of the American Statistical Association,
  Forthcoming. available at
  \url{http://imai.princeton.edu/research/list.html}}

\author{Graeme Blair, Princeton University,
\email{gblair@princeton.edu} and Kosuke Imai, Princeton University,
\email{kimai@princeton.edu} 
}

\examples{

data(race)

\dontrun{

## Multiple chain MCMC list experiment regression
## starts with overdispersed MLE starting values

mle.estimates <- ictreg(y ~ male + college + age + south, data = race)

draws <- mvrnorm(n = 4, mu = coef(mle.estimates), Sigma = vcov(mle.estimates) * 9)

bayesDraws <- list()
for (i in 1:4)
  bayesDraws[[i]] <- ictregBayes(y ~ male + college + age + south, data = race, delta.start = draws[i, 1:5], psi.start = draws[i, 6:10], burnin = 10000, n.draws = 100000, delta.tune = diag(.002, 5), psi.tune = diag(.00025, 5))

coda.object <- as.mcmc.list(bayesDraws)


}
}
 
\seealso{
\code{\link{predict.ictreg}} for fitted values
}

\keyword{models}
\keyword{regression}