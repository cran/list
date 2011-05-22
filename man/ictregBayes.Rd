\name{ictregBayes} 
\alias{ictregBayes}
\title{Item Count Technique} 
\description{ 
Function to conduct multivariate regression analyses of survey data with the item count technique, also known as the list experiment and the unmatched count technique.
} 
\usage{
ictregBayes(formula, data = parent.frame(), treat = "treat", J, 
            constrained.single = c("full","none","intercept"), 
	    constrained.multi = TRUE, fit.start = "lm", 
	    n.draws = 10000, burnin = 5000, thin = 0, delta.start, 
	    psi.start, delta.mu0, psi.mu0, delta.A0, psi.A0, 
            delta.tune, psi.tune, verbose = TRUE, ...)
} 

\arguments{ 
\item{formula}{An object of class "formula": a symbolic description of the model to be fitted.} 
\item{data}{A data frame containing the variables in the model}
\item{treat}{Name of treatment indicator as a string. For single sensitive item models, this refers to a binary indicator, and for multiple sensitive item models it refers to a multi-valued variable with zero representing the control condition. This can be an integer (with 0 for the control group) or a factor (with "control" for the control group).}
\item{J}{Number of non-sensitive (control) survey items. This will be set automatically to the maximum value of the outcome variable in the treatment group if no input is sent by the user.}
\item{constrained.single}{A logical value indicating whether the control group parameters are constrained to be equal in the single sensitive item design.}
\item{constrained.multi}{A logical value indicating whether the non-sensitive item count is included as a predictor in the sensitive item fits for the multiple sensitive item design.}
\item{fit.start}{Fit method for starting values. The options are \code{lm}, \code{glm}, \code{nls}, and \code{ml}, which use OLS, logistic regression, non-linear least squares, and maximum likelihood estimation to generate starting values, respectively. The default is \code{lm}.}
\item{n.draws}{Number of MCMC iterations after the burnin.}
\item{burnin}{The number of initial MCMC iterations that are discarded.}
\item{thin}{The interval of thinning, in which every other (\code{thin} = 1) or more iterations are discarded in the output object}
\item{delta.start}{Optional starting values for the sensitive item fit. This should be a vector with the length of the number of covariates for the single sensitive item design, and either a vector or a list with a vector of starting values for each of the sensitive items. The default runs an \code{ictreg} fit with the method set by the \code{fit.start} option.}
\item{psi.start}{Optional starting values for the control items fit. This should be a vector of length the number of covariates for the constrained models. The default runs an \code{ictreg} fit with the method set by the \code{fit.start} option.}
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

For the single sensitive item design, the model can constrain all control parameters to be equal (\code{constrained = "full"}), or just the intercept (\code{constrained = "intercept"}) or all the control fit parameters can be allowed to vary across the potential sensitive item values (\code{constrained = "none"}).

In the first version (4.0), only the standard design list experiment regression is estimated with a single sensitive item.

Convergence is at times difficult to achieve, so we recommend running multiple chains from overdispersed starting values by, for example, running an MLE or linear model using the ictreg() function, and then generating a set of overdispersed starting values using those estimates and their estimated variance-covariance matrix. An example is provided below.

}

\value{\code{ictregBayes} returns an object of class "ictregBayes".  The function \code{summary} is used to obtain a table of the results, using the \code{coda} package. Two attributes are also included, the data ("x"), the call ("call"), which can be extracted using the command, e.g., attr(ictregBayes.object, "x").

\item{mcmc}{an object of class "mcmc" that can be analyzed using the \code{coda} package.}
\item{x}{the design matrix}
\item{multi}{a logical value indicating whether the data included multiple sensitive items.}
\item{constrained}{a logical or character value indicating whether the control group parameters are constrained to be equal in the single sensitive item design, and whether the non-sensitive item count is included as a predictor in the sensitive item fits for the multiple sensitive item design.}
\item{delta.start}{Optional starting values for the sensitive item fit. This should be a vector with the length of the number of covariates. The default runs an \code{ictreg} fit with the method set by the \code{fit.start} option.}
\item{psi.start}{Optional starting values for the control items fit. This should be a vector of length the number of covariates. The default runs an \code{ictreg} fit with the method set by the \code{fit.start} option.}
\item{delta.mu0}{Optional vector of prior means for the sensitive item fit parameters, a vector of length the number of covariates.}
\item{psi.mu0}{Optional vector of prior means for the control item fit parameters, a vector of length the number of covariates.}
\item{delta.A0}{Optional matrix of prior precisions for the sensitive item fit parameters, a matrix of dimension the number of covariates.}
\item{psi.A0}{Optional matrix of prior precisions for the control items fit parameters, a matrix of dimension the number of covariates.}
\item{delta.tune}{A required vector of tuning parameters for the Metropolis algorithm for the sensitive item fit. This must be set and refined by the user until the acceptance ratios are approximately .4 (reported in the output).}
\item{psi.tune}{A required vector of tuning parameters for the Metropolis algorithm for the control item fit. This must be set and refined by the user until the acceptance ratios are approximately .4 (reported in the output).}
\item{call}{the matched call}

If the data includes multiple sensitive items, the following object is also included:
\item{treat.values}{a vector of the values used in the \code{treat} vector for the sensitive items, either character or numeric depending on the class of \code{treat}. Does not include the value for the control status}
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

## Standard single sensitive-item design

mle.estimates <- ictreg(y ~ male + college + age + south, data = race)

draws <- mvrnorm(n = 3, mu = coef(mle.estimates), 
  Sigma = vcov(mle.estimates) * 9)

bayesDraws.1 <- ictregBayes(y ~ male + college + age + south, data = race, 
  delta.start = draws[1, 1:5], psi.start = draws[1, 6:10], burnin = 10000, 
  n.draws = 100000, delta.tune = diag(.002, 5), psi.tune = diag(.00025, 5))

bayesDraws.2 <- ictregBayes(y ~ male + college + age + south, data = race, 
  delta.start = draws[2, 1:5], psi.start = draws[2, 6:10], burnin = 10000, 
  n.draws = 100000, delta.tune = diag(.002, 5), psi.tune = diag(.00025, 5))

bayesDraws.3 <- ictregBayes(y ~ male + college + age + south, data = race, 
  delta.start = draws[3, 1:5], psi.start = draws[3, 6:10], burnin = 10000, 
  n.draws = 100000, delta.tune = diag(.002, 5), psi.tune = diag(.00025, 5))

bayesDraws <- as.list(bayesDraws.1, bayesDraws.2, bayesDraws.3)

summary(bayesDraws)

## Multiple sensitive item design

mle.estimates.multi <- ictreg(y ~ male + college + age + south, data = multi)

draws <- mvrnorm(n = 3, mu = coef(mle.estimates.multi), 
  Sigma = vcov(mle.estimates.multi) * 9)

bayesMultiDraws.1 <- ictregBayes(y ~ male + college + age + south, data = multi, 
  delta.start = list(draws[1, 6:10], draws[1, 11:15]), psi.start = draws[1, 1:5], 
  burnin = 10000, n.draws = 50000, delta.tune = diag(.002, 5), 
  psi.tune = diag(.001, 5))

bayesMultiDraws.2 <- ictregBayes(y ~ male + college + age + south, data = multi, 
  delta.start = list(draws[2, 6:10], draws[2, 11:15]), psi.start = draws[2, 1:5], 
  burnin = 10000, n.draws = 50000, delta.tune = diag(.002, 5), 
  psi.tune = diag(.001, 5))

bayesMultiDraws.3 <- ictregBayes(y ~ male + college + age + south, data = multi, 
  delta.start = list(draws[3, 6:10], draws[3, 11:15]), psi.start = draws[3, 1:5], 
  burnin = 10000, n.draws = 50000, delta.tune = diag(.002, 5), 
  psi.tune = diag(.001, 5))

bayesMultiDraws <- as.list(bayesMultiDraws.1, bayesMultiDraws.2, bayesMultiDraws.3)

summary(bayesMultiDraws)

}
}
 
\seealso{
\code{\link{predict.ictreg}} for fitted values
}

\keyword{models}
\keyword{regression}