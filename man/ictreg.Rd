\name{ictreg} 
\alias{ictreg} 
\alias{ict}
\alias{list}
\title{Item Count Technique} 
\description{ 
Function to conduct multivariate regression analyses of survey data with the item count technique, also known as the list experiment and the unmatched count technique.
} 
\usage{
ictreg(formula, data = parent.frame(), treat = "treat", J, method = "ml", 
       overdispersed = FALSE, constrained = TRUE, floor = FALSE, ceiling = FALSE,
       ceiling.fit = "glm", floor.fit = "glm", fit.start = "nls", 
       fit.nonsensitive = "nls", multi.condition = "none", maxIter = 5000, 
       verbose = FALSE, ...)
} 

\arguments{ 
\item{formula}{An object of class "formula": a symbolic description of the model to be fitted.} 
\item{data}{A data frame containing the variables in the model}
\item{treat}{Name of binary treatment indicator.}
\item{J}{Number of non-sensitive (control) survey items.}
\item{method}{ Method for regression, either \code{ml} for the Maximum Likelihood (ML) estimation with the Expectation-Maximization algorithm; \code{lm} for linear model estimation; or \code{nls} for the Non-linear Least Squares (NLS) estimation with the two-step procedure.}
\item{overdispersed}{Indicator for the presence of overdispersion. If \code{TRUE}, the beta-binomial model is used in the EM algorithm, if \code{FALSE} the binomial model is used. Not relevant for the \code{NLS} or \code{lm} methods.}
\item{constrained}{A logical value indicating whether the control group parameters are constrained to be equal.  Not relevant for the \code{NLS} or \code{lm} methods}
\item{floor}{A logical value indicating whether the floor liar model should be used to adjust for the possible presence of respondents dishonestly reporting a negative preference for the sensitive item among those who hold negative views of all the non-sensitive items.}
\item{ceiling}{A logical value indicating whether the ceiling liar model should be used to adjust for the possible presence of respondents dishonestly reporting a negative preference for the sensitive item among those who hold affirmative views of all the non-sensitive items.}
\item{ceiling.fit}{Fit method for the M step in the EM algorithm used to fit the ceiling liar model. \code{glm} uses standard logistic regression, while \code{bayesglm} uses logistic regression with a weakly informative prior over the parameters.}
\item{floor.fit}{Fit method for the M step in the EM algorithm used to fit the floor liar model. \code{glm} uses standard logistic regression, while \code{bayesglm} uses logistic regression with a weakly informative prior over the parameters.}
\item{fit.start}{Fit method for starting values for standard design \code{ml} model. The options are \code{lm}, \code{glm}, and \code{nls}, which use OLS, logistic regression, and non-linear least squares to generate starting values, respectively. The default is \code{nls}.}
\item{fit.nonsensitive}{Fit method for the non-sensitive item fit for the \code{nls} method and the starting values for the \code{ml} method for the \code{modified} design. Options are \code{glm} and \code{nls}, and the default is \code{nls}.}
\item{multi.condition}{For the multiple sensitive item design, covariates representing the estimated count of affirmative responses for each respondent can be included directly as a level variable by choosing \code{level}, or as indicator variables for each value but one by choosing \code{indicators}. The default is \code{none}.}
\item{maxIter}{Maximum number of iterations for the Expectation-Maximization algorithm of the ML estimation.  The default is 5000.} 
\item{verbose}{ a logical value indicating whether model diagnostics are printed out during fitting.}
\item{...}{further arguments to be passed to NLS regression commands.}
}

\details{This function allows the user to perform regression analysis on data from the item count technique, also known as the list experiment and the unmatched count technique.

Three list experiment designs are accepted by this function: the standard design; the multiple sensitive item standard design; and the modified design proposed by Corstange (2009). 

For the standard design, three methods are implemented in this function: the linear model; the Maximum Likelihood (ML) estimation for the Expectation-Maximization (EM) algorithm; the nonlinear least squares (NLS) estimation with the two-step procedure both proposed in Imai (2010); and the Maximum Likelihood (ML) estimator in the presence of two types of dishonest responses, "ceiling" and "floor" liars. The ceiling model, floor model, or both, as described in Blair and Imai (2010) can be activated by using the \code{ceiling} and \code{floor} options. The constrained and unconstrained ML models presented in Imai (2010) are available through the \code{constrained} option, and the user can specify if overdispersion is present in the data for the no liars models using the \code{overdispersed} option to control whether a beta-binomial or binomial model is used in the EM algorithm to model the item counts.

The modified design and the multiple sensitive item design are automatically detected by the function, and only the binomial model without overdispersion is available. }

\value{\code{ictreg} returns an object of class "ictreg".  The function \code{summary} is used to obtain a table of the results.  The object \code{ictreg} is a list that contains the following components.  Some of these elements are not available depending on which method is used (\code{lm}, \code{nls} or \code{ml}), which design is used (\code{standard}, \code{modified}), whether multiple sensitive items are include (\code{multi}), and whether the constrained model is used (\code{constrained = TRUE}).

\item{par.treat}{point estimate for effect of covariate on item count fitted on treatment group}
\item{se.treat}{standard error for estimate of effect of covariate on item count fitted on treatment group}
\item{par.control}{point estimate for effect of covariate on item count fitted on control group}
\item{se.control}{standard error for estimate of effect of covariate on item count fitted on control group}
\item{coef.names}{variable names as defined in the data frame}
\item{design}{call indicating whether the \code{standard} design as proposed in Imai (2010) or thee \code{modified} design as proposed in Corstange (2009) is used}
\item{method}{call of the method used}
\item{overdispersed}{call indicating whether data is overdispersed}
\item{constrained}{call indicating whether the constrained model is used}
\item{boundary}{call indicating whether the floor/ceiling boundary models are used}
\item{multi}{indicator for whether multiple sensitive items were included in the data frame}
\item{call}{the matched call}
\item{data}{the \code{data} argument}
\item{x}{the design matrix}
\item{y}{the response vector}
\item{treat}{the vector indicating treatment status}

For the floor/ceiling models, several additional output objects are included:
\item{ceiling}{call indicating whether the assumption of no ceiling liars is relaxed, and ceiling parameters are estimated}
\item{par.ceiling}{point estimate for effect of covariate on whether respondents who answered affirmatively to all non-sensitive items and hold a true affirmative opinion toward the sensitive item lied and reported a negative response to the sensitive item }
\item{se.ceiling}{standard error for estimate for effect of covariate on whether respondents who answered affirmatively to all non-sensitive items and hold a true affirmative opinion toward the sensitive item lied and reported a negative response to the sensitive item}
\item{floor}{call indicating whether the assumption of no floor liars is relaxed, and floor parameters are estimated}
\item{par.ceiling}{point estimate for effect of covariate on whether respondents who answered negatively to all non-sensitive items and hold a true affirmative opinion toward the sensitive item lied and reported a negative response to the sensitive item }
\item{se.ceiling}{standard error for estimate for effect of covariate on whether respondents who answered negatively to all non-sensitive items and hold a true affirmative opinion toward the sensitive item lied and reported a negative response to the sensitive item}

For the multiple sensitive item design, the \code{par.treat} and \code{se.treat} vectors are returned as lists of vectors, one for each sensitive item. An additional output object is also included:
\item{treat.values}{a vector of the values used in the \code{treat} vector for the sensitive items. does not include the value for the control status}


For the unconstrained model, the \code{par.control} and \code{se.control} output is replaced by:
\item{par.control.phi0}{point estimate for effect of covariate on item count fitted on treatment group}
\item{se.control.phi0}{standard error for estimate of effect of covariate on item count fitted on treatment group}
\item{par.control.phi1}{point estimate for effect of covariate on item count fitted on treatment group}
\item{se.control.phi1}{standard error for estimate of effect of covariate on item count fitted on treatment group}

Depending upon the estimator requested by the user, model fit statistics are also included:
\item{llik}{the log likelihood of the model, if \code{ml} is used}
\item{resid.se}{the residual standard error, if \code{nls} or \code{lm} are used. This will be a scalar if the standard design was used, and a vector if the multiple sensitive item design was used}
\item{resid.df}{the residual degrees of freedom, if \code{nls} or \code{lm} are used. This will be a scalar if the standard design was used, and a vector if the multiple sensitive item design was used}

}

\references{
Blair, Graeme and Kosuke Imai. (2010) ``Statistical Analysis of List Experiments." 
  Technical Report, Department of Politics, Princeton University. 
  available at \url{http://imai.princeton.edu/research/listP.html}

Imai, Kosuke. (2010) ``Statistical Inference for the Item
  Count Technique.''  Technical Report, Department of Politics,
  Princeton University. available at \url{http://imai.princeton.edu/research/list.html}}

\author{Graeme Blair, Princeton University,
\email{gblair@princeton.edu} and Kosuke Imai, Princeton University,
\email{kimai@princeton.edu} 
}

\examples{

data(race)

# Fit linear regression
# Replicates Table 1 Columns 1-2 Imai (2010); note that age is divided by 10

lm.results <- ictreg(y ~ south + age + male + college, data = race, 
	       	      treat = "treat", J=3, method = "lm")

summary(lm.results)

# Fit two-step non-linear least squares regression
# Replicates Table 1 Columns 3-4 Imai (2010); note that age is divided by 10

nls.results <- ictreg(y ~ south + age + male + college, data = race, 
	       	      treat = "treat", J=3, method = "nls")

summary(nls.results)

\dontrun{

# Fit EM algorithm ML model with constraint
# Replicates Table 1 Columns 5-6, Imai (2010); note that age is divided by 10

ml.constrained.results <- ictreg(y ~ south + age + male + college, data = race, 
		       	  	 treat = "treat", J=3, method = "ml", 
				 overdispersed = FALSE, constrained = TRUE)

summary(ml.constrained.results)

# Fit EM algorithm ML model with no constraint
# Replicates Table 1 Columns 7-10, Imai (2010); note that age is divided by 10

ml.unconstrained.results <- ictreg(y ~ south + age + male + college, data = race, 
			    	   treat = "treat", J=3, method = "ml", 
				   overdispersed = FALSE, constrained = FALSE)

summary(ml.unconstrained.results)

# Fit EM algorithm ML model for multiple sensitive items
# Replicates Table 3 in Blair and Imai (2010)

multi.results <- ictreg(y ~ male + college + age + south, treat = "treat", J = 3, data = multi, method = "ml", multi.condition = "level")

summary(multi.results)

# Fit standard design ML model
# Replicates Table 7 Columns 1-2 in Blair and Imai (2010)

noboundary.results <- ictreg(y ~ age + college + male + south, treat = "treat", J = 3, data = affirm, method = "ml", overdispersed = FALSE)

summary(noboundary.results)

# Fit standard design ML model with ceiling effects alone
# Replicates Table 7 Columns 3-4 in Blair and Imai (2010)

ceiling.results <- ictreg(y ~ age + college + male + south, treat = "treat", J = 3, data = affirm, method = "ml", overdispersed = FALSE, ceiling = TRUE, ceiling.fit = "bayesglm")

summary(ceiling.results)

# Fit standard design ML model with floor effects alone
# Replicates Table 7 Columns 5-6 in Blair and Imai (2010)

floor.results <- ictreg(y ~ age + college + male + south, treat = "treat", J = 3, data = affirm, method = "ml", overdispersed = FALSE, floor = TRUE, floor.fit = "bayesglm")

summary(floor.results)

# Fit standard design ML model with floor and ceiling effects
# Replicates Table 7 Columns 7-8 in Blair and Imai (2010)

both.results <- ictreg(y ~ age + college + male + south, treat = "treat", J = 3, data = affirm, method = "ml", overdispersed = FALSE, floor = TRUE, ceiling = TRUE, floor.fit = "bayesglm", ceiling.fit = "bayesglm")

summary(both.results)

}
}
 
\seealso{
\code{\link{predict.ictreg}} for fitted values
}

\keyword{models}
\keyword{regression}