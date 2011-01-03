\name{predict.ictreg} 
\alias{predict.ictreg} 
\title{Predict Method for Item Count Technique} 
\description{ 
Function to calculate predictions and uncertainties of predictions from estimates from multivariate regression analysis of survey data with the item count technique.
} 
\usage{

  \method{predict}{ictreg}(object, newdata, se.fit = FALSE, interval = c("none","prediction"), level = .95, avg = FALSE, ...) 

} 

\arguments{ 
\item{object}{Object of class inheriting from "ictreg"}
\item{newdata}{An optional data frame in which to look for variables with which to predict. If omitted, the fitted values are used.}
\item{se.fit}{A switch indicating if standard errors are required.}
\item{interval}{Type of interval calculation.}
\item{level}{Tolerance level for prediction (tolerance) intervals.}
\item{avg}{A switch indicating if the mean prediction and associated statistics across all obserations in the dataframe will be returned instead of predictions for each observation.}
\item{...}{further arguments to be passed to or from other methods.}
}

\details{\code{predict.ictreg} produces predicted values, obtained by evaluating the regression function in the frame newdata (which defaults to \code{model.frame(object)}. If the logical \code{se.fit} is \code{TRUE}, standard errors of the predictions are calculated. Setting \code{interval} specifies computation of prediction (tolerance) intervals at the specified level or no intervals.

If \code{avg} is set to \code{TRUE}, the mean prediction across all observations in the dataset will be calculated, and if the \code{se.fit} option is set to \code{TRUE} a standard error for this mean estimate will be provided. The \code{interval} option will output prediction intervals instead of only the point estimate if set to \code{TRUE}.}

\value{\code{predict.ictreg} produces a vector of predictions or a matrix of predictions and bounds with column names fit, lwr, and upr if interval is set. If se.fit is TRUE, a list with the following components is returned:

\item{fit}{vector or matrix as above}
\item{se.fit}{standard error of prediction}
}

\seealso{\code{\link{ictreg}} for model fitting}

\references{Imai, Kosuke. (2010) ``Statistical Inference for the Item
  Count Technique.''  Technical Report, Department of Politics,
  Princeton University. available at \url{http://imai.princeton.edu/research/list.html}}

\author{Graeme Blair, Princeton University,
\email{gblair@princeton.edu} and Kosuke Imai, Princeton University,
\email{kimai@princeton.edu} 
}

\examples{
data(race)

race.south <- race.nonsouth <- race

race.south[, "south"] <- 1
race.nonsouth[, "south"] <- 0

\dontrun{

# Fit EM algorithm ML model with constraint with no covariates

ml.results.south.nocov <- ictreg(y ~ 1, 
   data = race[race$south == 1, ], method = "ml", treat = "treat", 
   J = 3, overdispersed = FALSE, constrained = TRUE)
ml.results.nonsouth.nocov <- ictreg(y ~ 1, 
   data = race[race$south == 0, ], method = "ml", treat = "treat", 
   J = 3, overdispersed = FALSE, constrained = TRUE)

# Calculate average predictions for respondents in the South 
# and the the North of the US for the MLE no covariates 
# model, replicating the estimates presented in Figure 1, 
# Imai (2010)

avg.pred.south.nocov <- predict(ml.results.south.nocov,
   newdata = as.data.frame(matrix(1, 1, 1)), se.fit = TRUE, 
   avg = TRUE)
avg.pred.nonsouth.nocov <- predict(ml.results.nonsouth.nocov,
   newdata = as.data.frame(matrix(1, 1, 1)), se.fit = TRUE, 
   avg = TRUE)

# Fit linear regression

lm.results <- ictreg(y ~ south + age + male + college, 
   data = race, treat = "treat", J=3, method = "lm")

# Calculate average predictions for respondents in the 
# South and the the North of the US for the lm model, 
# replicating the estimates presented in Figure 1, Imai (2010)

avg.pred.south.lm <- predict(lm.results, newdata = race.south, 
   se.fit = TRUE, avg = TRUE)
avg.pred.nonsouth.lm <- predict(lm.results, newdata = race.nonsouth, 
   se.fit = TRUE, avg = TRUE)

# Fit two-step non-linear least squares regression

nls.results <- ictreg(y ~ south + age + male + college, 
   data = race, treat = "treat", J=3, method = "nls")

# Calculate average predictions for respondents in the South 
# and the the North of the US for the NLS model, replicating
# the estimates presented in Figure 1, Imai (2010)

avg.pred.south.nls <- predict(nls.results, newdata = race.south, 
   se.fit = TRUE, avg = TRUE)
avg.pred.nonsouth.nls <- predict(nls.results, newdata = race.nonsouth, 
   se.fit = TRUE, avg = TRUE)

# Fit EM algorithm ML model with constraint

ml.constrained.results <- ictreg(y ~ south + age + male + college, 
   data = race, treat = "treat", J=3, method = "ml", 
   overdispersed = FALSE, constrained = TRUE)

# Calculate average predictions for respondents in the South 
# and the the North of the US for the MLE model, replicating the 
# estimates presented in Figure 1, Imai (2010)

avg.pred.south.mle <- predict(ml.constrained.results, 
   newdata = race.south, se.fit = TRUE, avg = TRUE)
avg.pred.nonsouth.mle <- predict(ml.constrained.results, 
   newdata = race.nonsouth, se.fit = TRUE, avg = TRUE)

}

}
 
\keyword{models}
\keyword{regression}