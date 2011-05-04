
ictregBayes <- function(formula, data = parent.frame(), treat="treat", J, constrained = TRUE, fit.start = "lm", n.draws = 10000, burnin = 5000, thin = 100, delta.start, psi.start, delta.mu0, psi.mu0, delta.A0, psi.A0, delta.tune, psi.tune, verbose = FALSE, ...){

  ictreg.call <- match.call()
  
  require(magic)
  
  # set up data frame, with support for standard and modified responses
  mf <- match.call(expand.dots = FALSE)
  
  # make all other call elements null in mf <- NULL in next line
  mf$treat <- mf$J <- mf$constrained <- mf$verbose <- mf$n.draws <- mf$burnin <- mf$thin <- mf$delta.start <- mf$psi.start <- mf$delta.mu0 <- mf$psi.mu0 <- mf$delta.A0 <- mf$psi.A0 <- mf$delta.tune <- mf$psi.tune <- mf$fit.start <- NULL
  mf[[1]] <- as.name("model.frame")
  mf$na.action <- 'na.pass'
  mf <- eval.parent(mf)
  
  # define design, response data frames
  x.all <- model.matrix.default(attr(mf, "terms"), mf)
  y.all <- model.response(mf)
 
  # list-wise missing deletion
  na.x <- apply(is.na(x.all), 1, sum)
  na.y <- is.na(y.all)
  
  ## treatment indicator for subsetting the dataframe
  t <- data[na.x==0 & na.y==0, paste(treat)]

  if(class(t) == "factor") {
    
    levels(t) <- tolower(levels(t))
    
    if (length(which(levels(t) == "control")) == 1) {
      t <- relevel(t, ref = "control")
    } else {
      warning("Note: using the first level as the control condition, but it is not labeled 'control'.")
    }
        
    condition.labels <- levels(t)
    t <- as.numeric(t) - 1
    treatment.labels <- condition.labels[2:length(condition.labels)]
    
  } else {
    condition.labels <- as.character(paste("Sensitive Item",sort(unique(t))))
    condition.labels[which(condition.labels == "Sensitive Item 0")] <- "control"
    treatment.labels <- condition.labels[condition.labels != "control"]
  }
  
  ## list wise delete
  y.all <- y.all[na.x==0 & na.y==0]
  x.all <- x.all[na.x==0 & na.y==0, , drop = FALSE]
  
  ## set up data objects for y and x for each group from user input
  x.treatment <- x.all[t != 0, , drop = FALSE]
  y.treatment <- subset(y.all, t != 0)
  x.control <- x.all[t == 0 , , drop = FALSE]
  y.control <- subset(y.all, t==0)

  ## J is defined by user for the standard design
  if(missing("J")) {
    J <- max(y.treatment) - 1
  }

  condition.values <- sort(unique(t))
  treatment.values <- 1:length(condition.values[condition.values!=0])

  if(length(treatment.values) > 1) {
    multi <- TRUE
    stop("The multi design is not yet supported. Please use a single sensitive item dataset.")
  } else {
    multi <- FALSE
  }
  
  n <- nrow(x.treatment) + nrow(x.control)

  coef.names <- colnames(x.all)

  nPar <- ncol(x.all)

  intercept.only <- ncol(x.all)==1 & sum(x.all[,1]==1) == n

  if (intercept.only == TRUE) {
    x.vars <- "1"
  } else {
    x.vars <- coef.names[-1]
  }
  
  logistic <- function(x) exp(x)/(1+exp(x))
  logit <- function(x) return(log(x)-log(1-x))

  if (missing("delta.mu0")) {
    if (multi == FALSE)
      delta.mu0 <- rep(0, nPar)
    else
      delta.mu0 <- rep(0, nPar + (1-constrained))
  }
  if (missing("psi.mu0")) {
    psi.mu0 <- rep(0, nPar)
  }
  if (missing("delta.A0")) {
    if (multi == FALSE)
      delta.A0 <- .01 * diag(nPar)
    else
      delta.A0 <- .01 * diag(nPar + (1-constrained))
  }
  if (missing("psi.A0")) {
    if (multi == FALSE)
      psi.A0 <- .01 * diag(nPar)
    else
      psi.A0 <- 0.01*diag(nPar)
  }
  if (missing("delta.tune")) {
    stop("The Metropolis tuning input object delta.tune is required.")

    ##if (multi == FALSE)
    ##  delta.tune <- diag(.002, nPar)
    ##else
    ##  delta.tune <- diag(.002, nPar + (1-constrained))
  }
  if (missing("psi.tune")) {
    stop("The Metropolis tuning input object psi.tune is required.") 

    ##if (multi == FALSE)
    ##  psi.tune <- diag(.0005, nPar)
    ##else
    ##  psi.tune <- diag(.0002, nPar)
  }
    
  if (multi == FALSE) {

    ictreg.fit <- ictreg(formula, data = data, treat = treat, J = J, method = fit.start, constrained = constrained)

    if (missing(delta.start))
      delta.start <- ictreg.fit$par.treat
    if (missing(psi.start))
      psi.start <- ictreg.fit$par.control

    ##ictregBayes.fit(reg$y, reg$treat, reg$x, 3, FALSE, 20000, 10000, 0, TRUE, coef(reg)[1:3], coef(reg)[4:6],
    ##                    rep(0, 3), rep(0, 3), 0.01*diag(3), 0.01*diag(3), diag(0.002, 3), diag(0.0005, 3))
    ictregBayes.fit <- ictregBayes.fit(Y = y.all, treat = t, X = x.all, J = J, constrained = constrained,
                                       n.draws = n.draws, burnin = burnin, thin = thin, verbose = verbose,
                                       delta.start = delta.start, psi.start = psi.start,
                                       delta.mu0 = delta.mu0, psi.mu0 = psi.mu0, delta.A0 = delta.A0, psi.A0 = psi.A0,
                                       delta.tune = delta.tune, psi.tune = psi.tune, ...)

    if (constrained == TRUE)
      colnames(ictregBayes.fit) <- c(paste("sensitive.", coef.names,sep=""),
                                     paste("control.", coef.names,sep=""), "acceptance.sensitive", "acceptance.control")
    else
      colnames(ictregBayes.fit) <- c(paste("sensitive.", coef.names,sep=""),
                                     paste("control.psi0.", coef.names,sep=""),
                                     paste("control.psi1.", coef.names,sep=""),
                                     "acceptance.sensitive", "acceptance.psi0", "acceptance.psi1")

    bayes.fit <- ictregBayes.fit
  } else {
    ##ictregBayesMulti.fit(reg3$y, reg3$treat, reg3$x, 3, FALSE, 10000, 5000, 0, TRUE,
    ##                         list("1" = coef(reg3)[7:13], "2" = coef(reg3)[14:20]), coef(reg2)[1:6],
    ##                         rep(0, 7), rep(0, 6), 0.01*diag(7), 0.01*diag(6), diag(0.002, 7), diag(0.0002, 6))

    if (constrained == T)
      ictreg.fit <- ictreg(formula, data = data, treat = treat, J = J, method = "ml", multi.condition = "none")
    else
      ictreg.fit <- ictreg(formula, data = data, treat = treat, J = J, method = "ml", multi.condition = "level")

    par.treat <- list()
    for (m in 1:length(treatment.values))
      par.treat[[as.character(m)]] <- ictreg.fit$par.treat[[m]]

    if (missing(delta.start))
      delta.start <- par.treat
    if (missing(psi.start))
      psi.start <- ictreg.fit$par.control
    
    ictregBayesMulti.fit <- ictregBayesMulti.fit(Y = y.all, treat = t, X = x.all, J = J, constrained = constrained,
                                                 n.draws = n.draws, burnin = burnin, thin = thin,  verbose = verbose,
                                                 delta.start = delta.start, psi.start = psi.start,
                                                 delta.mu0 = delta.mu0, psi.mu0 = psi.mu0, delta.A0 = delta.A0,
                                                 psi.A0 = psi.A0, delta.tune = delta.tune, psi.tune = psi.tune, ...)

    ##res2 <- ictregBayesMulti.fit(reg2$y, reg2$treat, reg2$x, 3, TRUE, 10000, 5000, 0, TRUE,
    ##                         list("1" = coef(reg2)[7:12], "2" = coef(reg2)[13:18]), coef(reg2)[1:6],
    ##                         rep(0, 6), rep(0, 6), 0.01*diag(6), 0.01*diag(6), diag(0.002, 6), diag(0.0002, 6))
        
    colnamesTemp <- convergeNames <- c()
    for (m in 1:length(treatment.labels)) {
      if (constrained == TRUE)
        colnamesTemp <- c(colnamesTemp, paste("sensitive.", treatment.labels[m], ".", coef.names, sep = ""))
      else
        colnamesTemp <- c(colnamesTemp, paste("sensitive.", treatment.labels[m], ".", c(coef.names, "y_i(0)"), sep = ""))

      convergeNames <- c(convergeNames, paste("acceptance.sensitive.", treatment.labels[m], sep = ""))
    }

    colnames(ictregBayesMulti.fit) <- c(colnamesTemp,
                                        paste("control.", coef.names,sep=""), convergeNames, "acceptance.control")

    bayes.fit <- ictregBayesMulti.fit
    
  }
  
  rownames(bayes.fit) <- seq(from = n.draws - burnin + 1 + thin + 1, to = n.draws - burnin + thin + 1 + floor((n.draws - burnin)/(thin + 1)) * (thin + 1), by = thin + 1)                            

  coda.object <- mcmc(data=bayes.fit,
                    start = 1,
                    thin = 1,
                    end = nrow(bayes.fit)
                    )

  if (multi == FALSE)
    class(coda.object) <- c("ictregBayes", "mcmc")
  else {
    class(coda.object) <- c("ictregBayesMulti", "mcmc")
    attr(coda.object, "treat.values") <- treatment.labels
  }

  attr(coda.object, "x") <- x.all
  attr(coda.object, "call") <- match.call()
  
  return(coda.object)
  
}

coef.ictregBayes <- function(object, ...) {

  if (is.null(attr(object, "call")$constrained) == FALSE) {
    if (attr(object, "call")$constrained == TRUE) {
      constrained <- TRUE
    } else {
      constrained <- FALSE
    }
  } else {
    constrained <- TRUE
  }
  
  if (constrained == TRUE)
    apply(object[,1:(ncol(object)-2)], 2, mean)
  else
    apply(object[,1:(ncol(object)-3)], 2, mean)
  
}

print.ictregBayes <- function(x, ...) {
  
  cat("\nItem Count Technique Bayesian Regression \n\nCall: ")
  
  dput(attr(x, "call"))

  cat("\nCoefficient estimates\n")

  print(coef(x))

  cat("\n")
  
  invisible(x)
  
}

vcov.ictregBayes <- function(object, ...) {

  if (is.null(attr(object, "call")$constrained) == FALSE) {
    if (attr(object, "call")$constrained == TRUE) {
      constrained <- TRUE
    } else {
      constrained <- FALSE
    }
  } else {
    constrained <- TRUE
  }
 
  if (constrained == TRUE)
    cov(object[, 1:(ncol(object)-2)])
  else
    cov(object[, 1:(ncol(object)-3)])
  
}

ictregBayes.fit <- function(Y, treat, X, J, constrained, n.draws, burnin, thin, verbose,
                            delta.start, psi.start, delta.mu0, psi.mu0, delta.A0, psi.A0,
                            delta.tune, psi.tune) {

  n <- length(Y)
  k <- ncol(X)
  n.par <- 2*(k + 1)
  keep <- thin + 1
  
  if (!constrained) {
    if (is.list(psi.start)) {
      psi.start <- c(psi.start$psi0, psi.start$psi1)
    } else {
      psi.start <- rep(psi.start, 2)
    }      
    if (is.list(psi.mu0)) {
      psi.mu0 <- c(psi.mu0$psi0, psi.mu0$psi1)
    } else {
      psi.mu0 <- rep(psi.mu0, 2)
    }
    if (is.list(psi.A0)) {
      psi.A0 <- c(as.double(psi.A0$psi0), as.double(psi.A0$psi1))
    } else {
      psi.A0 <- rep(as.double(psi.A0), 2)
    }
    if (is.list(psi.tune)) {
      psi.tune <- c(as.double(psi.tune$psi0), as.double(psi.tune$psi1))
    } else {
      psi.tune <- rep(as.double(psi.tune), 2)
    }
    n.par <- 3*(k + 1)
  }

  res <- .C("ictregBinom", as.integer(Y), as.integer(J), as.integer(n),
            as.integer(n.draws), as.integer(treat), as.double(X), 
            as.double(delta.start), as.double(psi.start), as.integer(k),
            as.double(delta.mu0), as.double(psi.mu0), as.double(delta.A0),
            as.double(psi.A0), as.double(delta.tune), as.double(psi.tune),
            as.integer(!constrained), as.integer(burnin), as.integer(keep),
            as.integer(verbose),
            allresults = double(n.par*floor((n.draws - burnin)/keep)),
	    PACKAGE = "list")$allresults

  res <- matrix(res, byrow = TRUE, ncol = n.par)

  class(res) <- "ictregBayes"
  return(res)
            
}

predict.ictregBayes <- function(object, newdata, newdata.diff, direct.glm, se.fit = FALSE,
                                interval = c("none","confidence"), level = .95, ...){
  
  n.draws <- nrow(object)
  
  if(missing(interval)) interval <- "none"
  
  nPar <- sum(substr(colnames(object), 1, 8)=="control.")
  ##if (sum(substr(colnames(object), 1, 10)=="sensitive.1.") > 1)
  ##  nPar <- sum(substr(colnames(object), 1, 10)=="sensitive.1.")

  diff <- missing(newdata.diff) == FALSE
  direct <- missing(direct.glm) == FALSE
  if (diff == TRUE & direct == TRUE)
    stop("The difference and direct comparison options cannot be used simultaneously. Try just one.")
  
  logistic <- function(object) exp(object)/(1+exp(object))
  
  if(missing(newdata)) {
    xvar <- attr(object, "x")
  } else { 
    if(nrow(newdata)==0)
      stop("No data in the provided data frame.")
    xvar <- model.matrix(as.formula(paste("~", c(attr(object, "call")$formula[[3]]))), newdata)
  }

  draws.list <- object[ , 1:nPar]
  
  if (direct == TRUE) {
    beta.direct <- coef(direct.glm)
    vcov.direct <- vcov(direct.glm)
    
    xvar.direct <- model.matrix(as.formula(paste("~", c(object$call$formula[[3]]))),
                                direct.glm$data)

    if (ncol(xvar.direct) != nPar)
      stop("Different number of covariates in direct and list regressions.")

    draws.direct <- mvrnorm(n = n.draws, mu = beta.direct, Sigma = vcov.direct)

    pred.direct.mean <- pred.direct.diff.mean <- rep(NA, n.draws)

  }

  if (diff == TRUE) {
    xvar.diff <- model.matrix(as.formula(paste("~", c(attr(object, "call")$formula[[3]]))),
                              newdata.diff)
    pred.diff.mean <- pred.diff.diff.mean <- rep(NA, n.draws)
  }
  
  pred.list.mean <- rep(NA, n.draws)
  
  for (d in 1:n.draws) {
        
    pred.list <- logistic(as.matrix(xvar) %*% as.matrix(draws.list[d, ]))
    pred.list.mean[d] <- mean(pred.list)

    if (diff == TRUE) {
      pred.diff <- logistic(as.matrix(xvar.diff) %*% as.matrix(draws.list[d, ]))
      pred.diff.mean[d] <- mean(pred.list)
      pred.diff.diff.mean[d] <- pred.list.mean[d] - pred.diff.mean[d]
    }
    
    if (direct == TRUE) {
      pred.direct <- logistic(as.matrix(xvar.direct) %*% as.matrix(draws.direct[d,]))
      pred.direct.mean[d] <- mean(pred.direct)
      pred.direct.diff.mean[d] <- pred.list.mean[d] - pred.direct.mean[d]
    }
    
  }

  critical.value <- qt(1-(1-level)/2, df = nrow(xvar))
  
  est.list <- mean(pred.list.mean)
  se.list <- sd(pred.list.mean)
  ci.upper.list <- est.list + critical.value * se.list
  ci.lower.list <- est.list - critical.value * se.list

  if (direct == TRUE) {
    est.direct <- mean(pred.direct.mean)
    se.direct <- sd(pred.direct.mean)
    est.direct.diff <- mean(pred.direct.diff.mean)
    se.direct.diff <- sd(pred.direct.diff.mean)
    ci.upper.direct <- est.direct + critical.value * se.direct
    ci.lower.direct <- est.direct - critical.value * se.direct
    ci.upper.direct.diff <- est.direct.diff + critical.value * se.direct.diff
    ci.lower.direct.diff <- est.direct.diff - critical.value * se.direct.diff

    return.object <- list(fit = rbind(est.list, est.direct, est.direct.diff))
      
    if ( interval == "confidence") {
      return.object <- as.data.frame(rbind(c(est.list, ci.lower.list, ci.upper.list),
                                      c(est.direct, ci.lower.direct, ci.upper.direct),
                                      c(est.direct.diff, ci.lower.direct.diff, ci.upper.direct.diff)))
      colnames(return.object) <- c("fit","lwr","upr")
      rownames(return.object) <- c("List", "Direct", "Difference (list - direct)")
    }
      
    if (se.fit == T)
      return.object$se.fit <- c(se.list, se.direct, se.direct.diff)

    attr(return.object, "concat") <- TRUE
    
  }

  if (diff == TRUE) {
    est.diff <- mean(pred.diff.mean)
    se.diff <- sd(pred.diff.mean)
    est.diff.diff <- mean(pred.diff.diff.mean)
    se.diff.diff <- sd(pred.diff.diff.mean)
    ci.upper.diff <- est.diff + critical.value * se.diff
    ci.lower.diff <- est.diff - critical.value * se.diff
    ci.upper.diff.diff <- est.diff.diff + critical.value * se.diff.diff
    ci.lower.diff.diff <- est.diff.diff - critical.value * se.diff.diff

    return.object <- list(fit = rbind(est.list, est.diff, est.diff.diff))
      
    if ( interval == "confidence") {
      return.object <- as.data.frame(rbind(c(est.list, ci.lower.list, ci.upper.list),
                                      c(est.diff, ci.lower.diff, ci.upper.diff),
                                      c(est.diff.diff, ci.lower.diff.diff, ci.upper.diff.diff)))
      colnames(return.object) <- c("fit","lwr","upr")
      rownames(return.object) <- c("Dataset 1", "Dataset 2", "Difference (1 - 2)")
    }
      
    if (se.fit == T)
      return.object$se.fit <- c(se.list, se.diff, se.diff.diff)

    attr(return.object, "concat") <- TRUE
    
  }

  if (diff == FALSE & direct == FALSE) {

    return.object <- list(fit = est.list)
      
    if ( interval == "confidence") {
      return.object <- as.data.frame(rbind(c(est.list, ci.lower.list, ci.upper.list)))
      names(return.object) <- c("fit","lwr","upr")
    }
      
    if (se.fit == T)
      return.object$se.fit <- c(se.list)
  }

  class(return.object) <- "predict.ictreg"
  
  return.object
  
}


