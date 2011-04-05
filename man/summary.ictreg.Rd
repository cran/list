\name{summary.ictreg} 
\alias{summary.ictreg} 
\title{Summary Method for the Item Count Technique} 
\description{ 
Function to summarize results from list experiment regression based on the ictreg() function, and to produce proportions of liars estimates.
} 
\usage{

  \method{summary}{ictreg}(object, boundary.proportions = FALSE, ...) 

} 

\arguments{ 
\item{object}{Object of class inheriting from "ictreg"}
\item{boundary.proportions}{A switch indicating whether, for models with ceiling effects, floor effects, or both (indicated by the \code{floor = TRUE}, \code{ceiling = TRUE} options in ictreg), the conditional probability of lying and the population proportions of liars are calculated.}
\item{...}{further arguments to be passed to or from other methods.}
}

\details{
\code{predict.ictreg} produces a summary of the results from an \code{ictreg} object. It displays the coefficients, standard errors, and fit statistics for any model from \code{ictreg}.

\code{predict.ictreg} also produces estimates of the conditional probability of lying and of the population proportion of liars for boundary models from \code{ictreg()} if \code{ceiling = TRUE} or \code{floor = TRUE}. 

The conditional probability of lying for the ceiling model is the probability that a respondent with true affirmative views of all the sensitive and non-sensitive items lies and responds negatively to the sensitive item. The conditional probability for the floor model is the probability that a respondent lies to conceal her true affirmative views of the sensitive item when she also holds true negative views of all the non-sensitive items. In both cases, the respondent may believe her privacy is not protected, so may conceal her true affirmative views of the sensitive item. 
}

\seealso{\code{\link{ictreg}} for model fitting}

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
# Fit standard design ML model with floor and ceiling effects
# Replicates Table 7 Columns 7-8 in Blair and Imai (2010)

floor.results <- ictreg(y ~ age + college + male + south, 
   treat = "treat", J = 3, data = affirm, method = "ml",
   fit.start = "nls", floor = TRUE, floor.fit = "bayesglm")

# Summarize fit object and generate conditional probability 
# of floor liars the population proportion of floor liars,
# both with standard errors.
# Replicates Table 7 Columns 7-8 last row in Blair and Imai (2010)


summary(floor.results, boundary.proportions = TRUE)
}

}
 
\keyword{models}
\keyword{regression}