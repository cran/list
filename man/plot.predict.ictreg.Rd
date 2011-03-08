\name{plot.predict.ictreg} 
\alias{plot.predict.ictreg} 
\title{Plot Method for the Item Count Technique} 
\description{ 
Function to plot predictions and confidence intervals of predictions from estimates from multivariate regression analysis of survey data with the item count technique.
} 
\usage{

  \method{plot}{predict.ictreg}(x, labels = NA, axes.ict = TRUE,
                                xlim = NULL, ylim = NULL, xlab = NULL, 
				ylab = "Estimated Proportion",
                                axes = F, pch = 19, xvec = NULL, ...) 

} 

\arguments{ 
\item{x}{object of class inheriting from "predict.ictreg". Either a single \code{predict} object from an \code{ictreg()} model fit, or a set of predict objects combined with \code{c()}.}
\item{labels}{a vector of labels for each prediction, plotted at the x axis.}
\item{axes.ict}{a switch indicating if custom plot axes are to be used with the user-provided estimate \code{labels}.}
\item{xlim}{a title for the y axis.}
\item{ylim}{a title for the y axis.}
\item{xlab}{a title for the x axis.}
\item{ylab}{a title for the y axis.}
\item{axes}{an indicator for whether default plot axes are included.}
\item{pch}{either an integer specifying a symbol or a single character to be used as the default in plotting points.}
\item{xvec}{a vector of x values at which the proportions will be printed.}
\item{...}{other graphical parameters for the \code{plot()} command are accepted.}
}

\details{\code{plot.predict.ictreg} produces plots with estimated population proportions of respondents answering the sensitive item in a list experiment in the affirmative, with confidence intervals. 

The function accepts a set of \code{predict.ictreg} objects calculated in the following manner:

\code{predict(ictreg.object, avg = TRUE, interval = "confidence")}

For each average prediction, a point estimate and its confidence interval is plotted at equally spaced intervals. The x location of the points can be manipulated with the \code{xvec} option.

Either a single predict object can be plotted, or a group of them combined with \code{c(predict.object1, predict.object2)}.
}

\seealso{\code{\link{ictreg}} for model fitting and \code{\link{predict.ictreg}} for predictions based on the model fits.}


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
race.south <- race.nonsouth <- race
race.south[, "south"] <- 1
race.nonsouth[, "south"] <- 0

\dontrun{

# Fit EM algorithm ML model with constraint
ml.constrained.results <- ictreg(y ~ south + age + male + college, 
   data = race, treat = "treat", J=3, method = "ml", 
   overdispersed = FALSE, constrained = TRUE)

# Calculate average predictions for respondents in the South 
# and the the North of the US for the MLE model, replicating the 
# estimates presented in Figure 1, Imai (2011)
avg.pred.south.mle <- predict(ml.constrained.results, 
   newdata = race.south, avg = TRUE, interval = "confidence")
avg.pred.nonsouth.mle <- predict(ml.constrained.results, 
   newdata = race.nonsouth, avg = TRUE, interval = "confidence")

# A plot of a single estimate and its confidence interval
plot(avg.pred.south.mle, labels = "South")
# A  plot of the two estimates and their confidence intervals
plot(c(avg.pred.south.mle, avg.pred.nonsouth.mle), labels = c("South", "Non-South"))

}
}
 
\keyword{models}
\keyword{regression}