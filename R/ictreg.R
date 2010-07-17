
ictreg <- function(formula, data, treat="treat", J, method = "nls", overdispersed = FALSE, constrained = TRUE, maxIter = 5000, verbose = FALSE, ...){

	ictreg.call <- match.call()
 	
	require(magic)
				
	n <- nrow(data)
	
	if(missing("J")) stop("J must be defined.")

	###
	###  Two-step Method of Moments estimator for the standard design
	###

	# treatment indicator for subsetting the dataframe
	t <- data[,paste(treat)]
	
	# set up data frame based on the formula sent
	data.all <- model.frame(formula=formula, data = data)
	
	# set up two subsets of data for treatment and control groups
	data.treatment <- model.frame(formula=formula, data=subset(data, t==1))
	data.control <- model.frame(formula=formula, data=subset(data, t==0))

	# set up the data for all points for the one-step estimator
	x.all <- model.matrix(attr(data.all, "terms"), data = data.all)
	y.all <- model.response(data.all)
	
	# set up data objects for y and x for each group from user input
	x.treatment <- model.matrix(attr(data.treatment, "terms"), data=data.treatment)
	y.treatment <- model.response(data.treatment)
	x.control <- model.matrix(attr(data.control, "terms"), data=data.control)
	y.control <- model.response(data.control)
			
	logit <- function(x) exp(x)/(1+exp(x))
	inv.logit <- function(x) return(log(x)-log(1-x))
	
	# ## 
	# Run two-step estimator
	
	## fit to the control group
	fit.glm <- glm(cbind(y.control, J-y.control) ~ x.control - 1, family = binomial(logit))
	coef.glm <- coef(fit.glm)
	names(coef.glm) <- paste("beta", 1:length(coef.glm), sep = "")
	
	if(is.null(ictreg.call$control)) fit.control <- nls( as.formula(paste("I(y.control/J) ~ logit(x.control %*% c(", paste(paste("beta", 1:length(coef.glm),sep=""), collapse= ","), "))")) , start = coef.glm, control = nls.control(maxiter=maxIter, warnOnly=TRUE), ... = ...)
	else fit.control <- nls( as.formula(paste("I(y.control/J) ~ logit(x.control %*% c(", paste(paste("beta", 1:length(coef.glm),sep=""), collapse= ","), "))")) , start = coef.glm, ... = ...)
	
	fit.control.coef <- summary(fit.control)$parameters[,1]
	
	## calculate the adjusted outcome
	y.treatment.pred <- y.treatment - logit(x.treatment %*% fit.control.coef)*J
	
	## fit to the treated
	y.treatment.start <- ifelse(y.treatment.pred > 0.5, 1, 0)
	fit.glm <- glm(y.treatment.start ~ x.treatment - 1, family = binomial(logit))
	coef.glm <- coef(fit.glm)
	names(coef.glm) <- paste("beta", 1:length(coef.glm), sep = "")
	 
	if(is.null(ictreg.call$control)) fit.treat <- nls( as.formula(paste("y.treatment.pred ~ logit(x.treatment %*% c(", paste(paste("beta", 1:length(coef.glm),sep=""), collapse= ","), "))")) , start = coef.glm, control = nls.control(maxiter = maxIter, warnOnly = TRUE), ... = ...)
	else fit.treat <- nls( as.formula(paste("y.treatment.pred ~ logit(x.treatment %*% c(", paste(paste("beta", 1:length(coef.glm),sep=""), collapse= ","), "))")) , start = coef.glm, ... = ...)
	
	vcov.twostep <- function(treat.fit, control.fit, J, y, treat, x) {

	  n <- length(y)
	  y1 <- y[treat == 1]
	  y0 <- y[treat == 0]
	  x1 <- x[treat == 1,]
	  x0 <- x[treat == 0,]
	  delta <- coef(treat.fit)
	  gamma <- coef(control.fit)
	
	  m1 <- c((y1 - J*logit(x1 %*% gamma) - logit(x1 %*% delta))*logit(x1 %*% delta)/(1+exp(x1 %*% delta))) * x1
	  m0 <- c((y0 - J*logit(x0 %*% gamma))*J*logit(x0 %*% gamma)/(1+exp(x0 %*% gamma))) * x0
	  Em1 <- t(m1) %*% m1 / n
	  Em0 <- t(m0) %*% m0 / n
	  F <- adiag(Em1, Em0)
	  Gtmp <- c(logit(x1 %*% delta)/(1 + exp(x1 %*% delta))) * x1
	  G1 <- -t(Gtmp) %*% Gtmp / n  
	  Gtmp <- c(sqrt(J*logit(x1 %*% delta)*logit(x1 %*% gamma)/((1 + exp(x1 %*% delta))*(1 + exp(x1 %*% gamma))))) * x1
	  G2 <- -t(Gtmp) %*% Gtmp / n
	  Gtmp <- c(J*logit(x0 %*% gamma)/(1 + exp(x0 %*% gamma))) * x0
	  G3 <- -t(Gtmp) %*% Gtmp / n
	  invG1 <- solve(G1)
	  invG3 <- solve(G3)
	  invG <- rbind(cbind(invG1, - invG1 %*% G2 %*% invG3), cbind(matrix(0, ncol = ncol(G1), nrow = nrow(G3)), invG3))
	
	  return(invG %*% F %*% t(invG) / n)
	  
	}
	
	vcov.nls <- vcov.twostep(fit.treat, fit.control, J, y.all, t, x.all)
	se.twostep <- sqrt(diag(vcov.nls))
	
	coef.control.start <- coef(fit.control)
	coef.treat.start <- coef(fit.treat)
				
	if(method=="ml") {
	
		# ##
		# Run EM estimator if requested
		
		#require(gamlss.dist)
		#require(VGAM)
		
		###
		### observed-data log-likelihood for beta binomial
		###
		obs.llik <- function(par, J, y, treat, x, const = FALSE) {

		  k <- ncol(x)
		  if (const) {
			coef.h0 <- coef.h1 <- par[1:k]
			rho.h0 <- rho.h1 <- logit(par[(k+1)])
			coef.g <- par[(k+2):(2*k+1)]
		  } else {
			coef.h0 <- par[1:k]
			rho.h0 <- logit(par[(k+1)])
			coef.h1 <- par[(k+2):(2*k+1)]
			rho.h1 <- logit(par[(2*k+2)])
			coef.g <- par[(2*k+3):(3*k+2)]
		  }
			
		  h0X <- logit(x %*% coef.h0)
		  h1X <- logit(x %*% coef.h1)
		  gX <- logit(x %*% coef.g)
		
		  ind10 <- ((treat == 1) & (y == 0))
		  ind1J1 <- ((treat == 1) & (y == (J+1)))
		  ind1y <- ((treat == 1) & (y > 0) & (y < (J+1)))
		
		  if (sum(ind10) > 0) {
			p10 <- sum(log(1-gX[ind10]) + dBB(x = 0, bd = J, mu = h0X[ind10], sigma = rho.h0/(1-rho.h0), log = TRUE))
		  } else {
			p10 <- 0
		  }
		  if (sum(ind1J1) > 0) {
			p1J1 <- sum(log(gX[ind1J1]) + dBB(J, bd = J, mu = h1X[ind1J1], sigma = rho.h1/(1-rho.h1), log = TRUE))
		  } else {
			p1J1 <- 0
		  }
		  if (sum(ind1y) > 0) {
			p1y <- sum(log(gX[ind1y]*dBB(y[ind1y]-1, bd = J, mu = h1X[ind1y], sigma = rho.h1/(1-rho.h1), log = FALSE) +
						   (1-gX[ind1y])*dBB(y[ind1y], bd = J, mu = h0X[ind1y], sigma = rho.h0/(1-rho.h0), log = FALSE)))
		  } else {
			p1y <- 0
		  }
		  if (sum(treat == 0) > 0) {
			p0y <- sum(log(gX[!treat]*dBB(y[!treat], bd = J, mu = h1X[!treat], sigma = rho.h1/(1-rho.h1), log = FALSE) +
						   (1-gX[!treat])*dBB(y[!treat], bd = J, mu = h0X[!treat], sigma = rho.h0/(1-rho.h0), log = FALSE)))
		  } else {
			p0y <- 0
		  }
		  
		  return(p10+p1J1+p1y+p0y)
		}
		###
		###  Observed data log-likelihood for binomial
		###
		obs.llik.binom <- function(par, J, y, treat, x, const = FALSE) {

		  k <- ncol(x)
		  if (const) {
			coef.h0 <- coef.h1 <- par[1:k]
			coef.g <- par[(k+1):(2*k)]
		  } else {
			coef.h0 <- par[1:k]
			coef.h1 <- par[(k+1):(2*k)]
			coef.g <- par[(2*k+1):(3*k)]
		  }
		  
		  h0X <- logit(x %*% coef.h0)
		  h1X <- logit(x %*% coef.h1)
		  gX <- logit(x %*% coef.g)
		
		  ind10 <- ((treat == 1) & (y == 0))
		  ind1J1 <- ((treat == 1) & (y == (J+1)))
		  ind1y <- ((treat == 1) & (y > 0) & (y < (J+1)))
		
		  if (sum(ind10) > 0) {
			p10 <- sum(log(1-gX[ind10]) + dbinom(x = 0, size = J, prob = h0X[ind10], log = TRUE))
		  } else {
			p10 <- 0
		  }
		  if (sum(ind1J1) > 0) {
			p1J1 <- sum(log(gX[ind1J1]) + dbinom(J, size = J, prob = h1X[ind1J1], log = TRUE))
		  } else {
			p1J1 <- 0
		  }
		  if (sum(ind1y) > 0) {
			p1y <- sum(log(gX[ind1y]*dbinom(y[ind1y]-1, size = J, prob = h1X[ind1y], log = FALSE) +
						   (1-gX[ind1y])*dbinom(y[ind1y], size = J, prob = h0X[ind1y], log = FALSE)))
		  } else {
			p1y <- 0
		  }
		  if (sum(treat == 0) > 0) {
			p0y <- sum(log(gX[!treat]*dbinom(y[!treat], size = J, prob = h1X[!treat], log = FALSE) +
						   (1-gX[!treat])*dbinom(y[!treat], size = J, prob = h0X[!treat], log = FALSE)))
		  } else {
			p0y <- 0
		  }
		
		  return(p10+p1J1+p1y+p0y)
		}
		
		###
		### EM algorithm
		###
		
		## Estep for beta-binomial
		Estep <- function(par, J, y, treat, x, const = FALSE) {

		  k <- ncol(x)
		  if (const) {
			coef.h0 <- coef.h1 <- par[1:k]
			rho.h0 <- rho.h1 <- logit(par[(k+1)])
			coef.g <- par[(k+2):(2*k+1)]
		  } else {
			coef.h0 <- par[1:k]
			rho.h0 <- logit(par[(k+1)])
			coef.h1 <- par[(k+2):(2*k+1)]
			rho.h1 <- logit(par[(2*k+2)])
			coef.g <- par[(2*k+3):(3*k+2)]
		  }
		  
		  h0X <- logit(x %*% coef.h0)
		  h1X <- logit(x %*% coef.h1)
		  gX <- logit(x %*% coef.g)
		
		  ind <- !((treat == 1) & ((y == 0) | (y == (J+1))))
		  w <- rep(NA, length(y))
		  w[ind] <- gX[ind]*dBB((y-treat)[ind], bd = J, mu = h1X[ind], sigma = rho.h1/(1-rho.h1), log = FALSE) /
			(gX[ind]*dBB((y-treat)[ind], bd = J, mu = h1X[ind], sigma = rho.h1/(1-rho.h1), log = FALSE) +
			 (1-gX[ind])*dBB(y[ind], bd = J, mu = h0X[ind], sigma = rho.h0/(1-rho.h0), log = FALSE))
		
		  w[(treat == 1) & (y == 0)] <- 0
		  w[(treat == 1) & (y == (J+1))] <- 1
		  
		  return(w)
		}

		
		## Estep for binomial
		Estep.binom <- function(par, J, y, treat, x, const = FALSE) {

		  k <- ncol(x)
		  if (const) {
			coef.h0 <- coef.h1 <- par[1:k]
			coef.g <- par[(k+1):(2*k)]
		  } else {
			coef.h0 <- par[1:k]
			coef.h1 <- par[(k+1):(2*k)]
			coef.g <- par[(2*k+1):(3*k)]
		  }
		  
		  h0X <- logit(x %*% coef.h0)
		  h1X <- logit(x %*% coef.h1)
		  gX <- logit(x %*% coef.g)
		
		  ind <- !((treat == 1) & ((y == 0) | (y == (J+1))))
		  w <- rep(NA, length(y))
		  w[ind] <- gX[ind]*dbinom((y-treat)[ind], size = J, prob = h1X[ind], log = FALSE) /
			(gX[ind]*dbinom((y-treat)[ind], size = J, prob = h1X[ind], log = FALSE) +
			 (1-gX[ind])*dbinom(y[ind], size = J, prob = h0X[ind], log = FALSE))
		  w[(treat == 1) & (y == 0)] <- 0
		  w[(treat == 1) & (y == (J+1))] <- 1
		  
		  return(w)
		}


		

		## Mstep 1: weighted MLE for logistic regression
		wlogit.fit <- function(y, treat, x, w, par = NULL) {

		  yrep <- rep(c(1,0), each = length(y))
		  xrep <- rbind(x, x)
		  wrep <- c(w, 1-w)
		  return(glm(cbind(yrep, 1-yrep) ~ xrep - 1, weights = wrep, family = binomial(logit), start = par))
		  		  
		}
		
		## Mstep 2: weighted MLE for beta-binomial regression
		wbb.fit <- function(J, y, treat, x, w, par0 = NULL, par1 = NULL) {

		  Not0 <- ((treat == 1) & (y == (J+1)))
		  y0 <- y[!Not0]
		  x0 <- x[!Not0,]
		  w0 <- 1-w[!Not0]
		  fit0 <- vglm(cbind(y0, J-y0) ~ x0, betabinomial, weights = w0, coefstart = par0)
		
		  Not1 <- ((treat == 1) & (y == 0))
		  y1 <- y
		  y1[treat == 1] <- y1[treat == 1] - 1
		  y1 <- y1[!Not1]
		  x1 <- x[!Not1]
		  w1 <- w[!Not1]
		  fit1 <- vglm(cbind(y1, J-y1) ~ x1, betabinomial, weights = w1, coefstart = par1)
		
		  return(list(fit1 = fit1, fit0 = fit0))
		}
		
		## Mstep2: weighted MLE for binomial regression
		wbinom.fit <- function(J, y, treat, x, w, par0 = NULL, par1 = NULL) {

		  Not0 <- ((treat == 1) & (y == (J+1)))
		  y0 <- y[!Not0]
		  x0 <- x[!Not0,]
		  w0 <- 1-w[!Not0]
		  fit0 <- glm(cbind(y0, J-y0) ~ x0, family = binomial(logit), weights = w0, start = par0)
		
		  Not1 <- ((treat == 1) & (y == 0))
		  y1 <- y
		  y1[treat == 1] <- y1[treat == 1] - 1
		  y1 <- y1[!Not1]
		  x1 <- x[!Not1,]
		  w1 <- w[!Not1]
		  fit1 <- glm(cbind(y1, J-y1) ~ x1, family = binomial(logit), weights = w1, start = par1)
		
		  return(list(fit1 = fit1, fit0 = fit0))
		}
		
		###
		### Running the EM algorithm
		###
		
		nPar <- ncol(x.all)
		
		if(constrained==F) {
		
			# Run unconstrained model
								
			if (overdispersed==T) {
			  par <- c(rep(c(coef.control.start, 0.5), 2), coef.treat.start)
			} else {
			  par <- c(rep(coef.control.start, 2), coef.treat.start)
			}
					
			## max number of iterations
			pllik <- -Inf
			
			if (overdispersed==T) {
				llik <- obs.llik(par, J = J, y = y.all, treat = t, x = x.all)
			} else {
				llik <- obs.llik.binom(par, J = J, y = y.all, treat = t, x = x.all)
			}
								
			Not0 <- (t & (y.all == (J+1)))
			Not1 <- (t & (y.all == 0))
			
			counter <- 0
			while (((llik - pllik) > 10^(-8)) & (counter < maxIter)) {
						
				if(overdispersed==T) {
				
					w <- Estep(par, J, y.all, t, x.all)
					
					lfit <- wlogit.fit(y.all, t, x.all, w, par = par[(nPar*2+3):length(par)])
								
					fit0 <- vglm(cbind(y.all[!Not0], J-y.all[!Not0]) ~ x.all[!Not0,] - 1, betabinomial, weights = 1-w[!Not0], coefstart = par[c(1,(nPar+1),2:(nPar))])
					
					y1 <- y.all
					
					y1[t==1] <- y1[t==1]-1
					
					fit1 <- vglm(cbind(y1[!Not1], J-y1[!Not1]) ~ x.all[!Not1,] - 1, betabinomial, weights = w[!Not1], coefstart = par[c((nPar+2),(2*nPar+2),(nPar+3):(2*nPar+1))])
					
					par <- c(coef(fit0)[c(1,3:(nPar+1),2)], coef(fit1)[c(1,3:(nPar+1),2)], coef(lfit))
					
				} else {
	
					w <- Estep.binom(par, J, y.all, t, x.all)
					
					lfit <- wlogit.fit(y.all, t, x.all, w, par = par[(nPar*2+1):length(par)])
								
					fit0 <- glm(cbind(y.all[!Not0], J-y.all[!Not0]) ~ x.all[!Not0,] - 1, family = binomial(logit), weights = 1-w[!Not0], start = par[1:(nPar)])
					
					y1 <- y.all
					
					y1[t==1] <- y1[t==1]-1
					
					fit1 <- glm(cbind(y1[!Not1], J-y1[!Not1]) ~ x.all[!Not1,] - 1, family = binomial(logit), weights = w[!Not1], start = par[(nPar+1):(2*nPar)])
					
					par <- c(coef(fit0), coef(fit1), coef(lfit))
				
				}
							
				pllik <- llik
				
				if(verbose==T) cat(paste(counter, round(llik, 4), "\n"))			
				
				if (overdispersed==T) {
					llik <- obs.llik(par, J = J, y = y.all, treat = t, x = x.all)
				} else {
					llik <- obs.llik.binom(par, J = J, y = y.all, treat = t, x = x.all)
				}
				
				counter <- counter + 1
							
				if (llik < pllik) {
					stop("log-likelihood is not monotonically increasing.")
				}
				
				if(counter == (maxIter-1)) warning("number of iterations exceeded maximum in ML")
			
			}
					
			## getting  standard errors
			
			if (overdispersed==T) {
			  MLEfit <- optim(par, obs.llik, method = "BFGS", J = J, y = y.all, treat = t, x = x.all, hessian = TRUE, control = list(maxit = 0))
			  vcov.mle <- solve(-MLEfit$hessian)
			  se.mle <- sqrt(diag(vcov.mle))
			} else {
			  MLEfit <- optim(par, obs.llik.binom, method = "BFGS", J = J, y = y.all, treat = t, x = x.all, hessian = TRUE, control = list(maxit = 0))
			  vcov.mle <- solve(-MLEfit$hessian)
			  se.mle <- sqrt(diag(vcov.mle))
			}
			
		} else { # end of constrained model
		
			# constrained model
			
			if (overdispersed==T) {
			  par <- c(coef.control.start, 0.5, coef.treat.start)
			} else {
			  par <- c(coef.control.start, coef.treat.start)
			}
			
			pllik.const <- -Inf
			
			if (overdispersed==T) {
			  llik.const <- obs.llik(par, J = J, y = y.all, treat = t, x = x.all, const = TRUE)
			} else {
			  llik.const <- obs.llik.binom(par, J = J, y = y.all, treat = t, x = x.all, const = TRUE)
			}
			  			
			counter <- 0
			while (((llik.const - pllik.const) > 10^(-8)) & (counter < maxIter)) {
			
			  if (overdispersed==T) {
			  
				w <- Estep(par, J, y.all, t, x.all, const = TRUE)
				lfit <- wlogit.fit(y.all, t, x.all, w, par = par[(nPar+2):(nPar*2+1)])
				
			  } else {
			  
				w <- Estep.binom(par, J, y.all, t, x.all, const = TRUE)
				lfit <- wlogit.fit(y.all, t, x.all, w, par = par[(nPar+1):(nPar*2)])
				
			  }
			  
			  
			  y.var <- as.character(formula)[[2]]
			  
			  x.vars <- strsplit(gsub(" ", "",as.character(formula)[[3]]), split="+", fixed=T)[[1]]
			  			  
			  data.all$w <- w
			  
			  dtmp <- rbind(cbind(data.all, t), cbind(data.all, t)[t==1, ])
			 
			  dtmp[((dtmp$t == 1) & ((1:nrow(dtmp)) <= n)), paste(y.var)] <-
				dtmp[((dtmp$t == 1) & ((1:nrow(dtmp)) <= n)), paste(y.var)] - 1
				
			  dtmp$w[((dtmp$t == 1) & ((1:nrow(dtmp)) > n))] <-
				1 - dtmp$w[((dtmp$t == 1) & ((1:nrow(dtmp)) > n))]
				
			  dtmp$w[dtmp$t == 0] <- 1
			  
			  dtmp <- dtmp[dtmp$w > 0, ]
			  			  
			  if (overdispersed==T) {
			  			  
				fit <- vglm(as.formula(paste("cbind(", y.var, ", J-", y.var, ") ~ ", paste(x.vars, collapse=" + "))), betabinomial, weights = dtmp$w, coefstart = par[c(1,(nPar+1),2:(nPar))], data = dtmp)
				
				par <- c(coef(fit)[c(1,3:(nPar+1),2)], coef(lfit))
				
			  } else {
	
				fit <- glm(as.formula(paste("cbind(", y.var, ", J-", y.var, ") ~ ", paste(x.vars, collapse=" + "))), family = binomial(logit), weights = dtmp$w, start = par[1:(nPar)], data = dtmp)
				
				par <- c(coef(fit), coef(lfit))
				
			  }
			  
			  pllik.const <- llik.const
			  
			  if(verbose==T) cat(paste(counter, round(llik.const, 4), "\n"))
			  
			  if (overdispersed==T) {
				llik.const <- obs.llik(par, J = J, y = y.all, treat = t, x = x.all, const = TRUE)
			  } else {
				llik.const <- obs.llik.binom(par, J = J, y = y.all, treat = t, x = x.all, const = TRUE)
			  }
			  
			  counter <- counter + 1
			  if (llik.const < pllik.const) {
				stop("log-likelihood is not monotonically increasing.")
			  }
			  
			  if(counter == (maxIter-1)) warning("number of iterations exceeded maximum in ML")
			  
			}
			
			if (overdispersed==T) {
			
			  MLEfit <- optim(par, obs.llik, method = "BFGS", J = J, y = y.all, treat = t, x = x.all, const = TRUE, hessian = TRUE, control = list(maxit = 0))
			  
			  vcov.mle <- solve(-MLEfit$hessian)
			  se.mle <- sqrt(diag(vcov.mle))
			  
			} else {
			
			  MLEfit <- optim(par, obs.llik.binom, method = "BFGS", J = J, y = y.all, treat = t,  x = x.all, const = TRUE, hessian = TRUE, control = list(maxit = 0))
			  
			  vcov.mle <- solve(-MLEfit$hessian)
			  se.mle <- sqrt(diag(vcov.mle))
			  
			}

		} # end of constrained model

	} # end of em algorithm
	
	# ## 
	# Set up return object
	
	coef.names <- c("(Intercept)", strsplit(gsub(" ", "",as.character(formula)[[3]]), split="+", fixed=T)[[1]])
	
	df <- data.all
	df$treat <- t
	
	if (method == "nls") {
	
		par.treat <- coef(fit.treat)
		se.treat <- se.twostep[1:(length(par.treat))]
		
		par.control <- coef(fit.control)
		se.control <- se.twostep[(length(par.treat)+1):(length(se.twostep))]
		
		names(par.treat) <- names(se.treat) <- names(par.control) <- names(se.control) <- coef.names
		
		sum.fit.treat <- summary(fit.treat)
		
		resid.se <- sum.fit.treat$sigma
		resid.df <- sum.fit.treat$df[2]
	
		return.object <- list(par.treat=par.treat, se.treat=se.treat, par.control=par.control, se.control=se.control, vcov=vcov.nls, resid.se=resid.se, resid.df=resid.df, coef.names=coef.names, method = method, overdispersed=overdispersed, call = match.call(), df = df)

	}
	
	if (method == "ml"){
	

		if (constrained == T) {

			par.treat <- MLEfit$par[(nPar+1):(nPar*2)]
			par.control <- MLEfit$par[1:(nPar)]
		
			if (overdispersed == T){
				se.treat <- se.mle[(nPar+2):(nPar*2+1)]
				se.control <- se.mle[1:(nPar)]
			} else {
				se.treat <- se.mle[(nPar+1):(nPar*2)]
				se.control <- se.mle[1:(nPar)]
			}
			
			names(par.treat) <- names(se.treat) <- names(par.control) <- names(se.control) <- coef.names
			
			return.object <- list(par.treat=par.treat, se.treat=se.treat, par.control=par.control, se.control=se.control, vcov=vcov.mle, llik=llik.const, coef.names=coef.names, method = method, overdispersed=overdispersed, constrained=constrained, call = match.call(), df = df)

		} else { 
		
			par.treat <- MLEfit$par[(nPar*2+1):(nPar*3)]
			par.control.phi0 <- MLEfit$par[1:(nPar)]
			par.control.phi1 <- MLEfit$par[(nPar+1):(nPar*2)]
			
			if (overdispersed == T){
				se.treat <- se.mle[(nPar*2+3):(nPar*3 + 2)]
				se.control.phi0 <- se.mle[1:(nPar)]
				se.control.phi1 <- se.mle[(nPar*2+2):(nPar*2+1)]
			} else {
				se.treat <- se.mle[(nPar*2+1):(nPar*3)]
				se.control.phi0 <- se.mle[1:(nPar)]
				se.control.phi1 <- se.mle[(nPar+1):(nPar*2)]
			}
			
			names(par.treat) <- names(se.treat) <- names(par.control.phi0) <- names(se.control.phi0) <- names(par.control.phi1) <- names(se.control.phi1) <- coef.names

			return.object <- list(par.treat=par.treat, se.treat=se.treat, par.control.phi0=par.control.phi0, se.control.phi0=se.control.phi0, par.control.phi1=par.control.phi1, se.control.phi1=se.control.phi1, vcov=vcov.mle, llik=llik, coef.names=coef.names, method = method, overdispersed=overdispersed, constrained=constrained, call = match.call(), df=df)

		}
	
	}
	
	
	class(return.object) <- "ictreg"

	return.object
	
}


print.ictreg <- function(x, ...){
	
	print(unlist(x))
	
	invisible(x)
	
}

predict.ictreg <- function(object, newdata, se.fit = FALSE, interval = c("none","prediction"), level = .95, avg = FALSE, ...){

	if(missing(interval)) interval <- "none"

	nPar <- length(object$coef.names)
	
	# dummy value for constraint
	if(object$method=="nls") object$constrained <- F
	
	logit <- function(object) exp(object)/(1+exp(object))
	
	# extract only treatment coef's, var/covar
	beta <- coef(object)[1:nPar]
	var.beta <- vcov(object)[1:nPar, 1:nPar]
	
	df <- model.frame(formula=as.formula(object$call$formula), data = object$df)
    if(missing(newdata)) { xvar <- model.matrix(attr(df, "terms"), data = df) }
	else { 
   		if(nrow(newdata)==0) stop("No data in the provided data frame.")
   		xvar <- model.matrix(as.formula(paste("~",as.character(object$call$formula)[[3]])), newdata)
	}
	
	n <- nrow(xvar)

	pix <- logit(xvar %*% beta)
 	v <- pix/(1+exp(xvar %*% beta))
 	
 	if(avg==FALSE){
		var.pred <- rep(NA, n)
		for (i in 1:n) {
			var.pred[i] <- v[i] * v[i] * (xvar[i,, drop = FALSE] %*% var.beta %*% t(xvar[i,, drop = FALSE]))
		}
		
		#var.pred.matrix <- v * v %*% (xvar %*% var.beta %*% t(xvar))
	} else {
		
		pix <- mean(pix)
		
		var.pred <- 0
		for (i in 1:n) {
			for (j in 1:n) {
				var.pred <- var.pred + v[i] * v[j] * (xvar[i,, drop = FALSE] %*% var.beta %*% t(xvar[j,, drop = FALSE]))
			}
		}
			
	}
	
	se <- c(sqrt(var.pred)/n)
	
	return.object <- list(fit=as.vector(pix))
	
	if(interval=="prediction"){
	
		critical.value <- qt(1-(1-level)/2, df = n)
		ci.upper <- pix + critical.value * se
		ci.lower <- pix - critical.value * se
		
		return.object <- list(fit=data.frame(fit = pix, lwr = ci.lower, upr = ci.upper))
		
	}
	
	if(se.fit==T) return.object$se.fit <- as.vector(se)

	return(return.object)

}

coef.ictreg <- function(object, ...){

	nPar <- length(object$coef.names)

	if(object$method=="nls"){
		coef <- c(object$par.treat,object$par.control)
		names(coef) <- c(paste("treat.",object$coef.names,sep=""), paste("control.",object$coef.names,sep=""))
	} else if(object$method=="ml") {
		if(object$constrained==F) {
			coef <- c(object$par.treat,object$par.control.phi0,object$par.control.phi1)
			names(coef) <- c(paste("treat.",object$coef.names,sep=""), paste("control.phi0.",object$coef.names,sep=""),paste("control.phi1.",object$coef.names,sep=""))
		} else {
			coef <- c(object$par.treat,object$par.control)
			names(coef) <- c(paste("treat.",object$coef.names,sep=""), paste("control.",object$coef.names,sep=""))
		}
	}
		
	return(coef)

}

vcov.ictreg <- function(object, ...){

	vcov <- object$vcov
	
	nPar <- length(object$coef.names)
	
	# dummy value for constraint
	if(object$method=="nls") object$constrained <- F
	
	#if(object$method=="nls") vcov <- rbind( cbind( vcov[(nPar+1):(nPar*2), (nPar+1):(nPar*2)], vcov[(nPar+1):(nPar*2), 1:nPar]  ), cbind(vcov[1:nPar, (nPar+1):(nPar*2)] ,vcov[1:nPar, 1:nPar])  )
		
	if(object$method=="ml" & object$constrained==T)	vcov <- rbind( cbind( vcov[(nPar+1):(nPar*2), (nPar+1):(nPar*2)], vcov[(nPar+1):(nPar*2), 1:nPar]  ), cbind(vcov[1:nPar, (nPar+1):(nPar*2)] ,vcov[1:nPar, 1:nPar])  )
	
	else if (object$method=="ml" & object$constrained==F) vcov <- rbind( cbind(vcov[(nPar*2+1):(nPar*3),(nPar*2+1):(nPar*3)], vcov[(nPar*2+1):(nPar*3),1:nPar], vcov[(nPar*2+1):(nPar*3), (nPar+1):(nPar*2)] ), cbind(vcov[1:nPar, (nPar*2+1):(nPar*3)] , vcov[1:nPar,1:nPar] , vcov[1:nPar, (nPar+1):(nPar*2)]), cbind(vcov[(nPar+1):(nPar*2), (nPar*2+1):(nPar*3)] , vcov[(nPar+1):(nPar*2), 1:nPar] , vcov[(nPar+1):(nPar*2),(nPar+1):(nPar*2)]) )

	if(object$method=="nls"){
		rownames(vcov) <- colnames(vcov)<- c(paste("treat.",object$coef.names,sep=""), paste("control.",object$coef.names,sep=""))
	} else if(object$method=="ml") {
		if(object$constrained==F) {
			rownames(vcov) <- colnames(vcov) <- c(paste("treat.",object$coef.names,sep=""), paste("control.phi0.",object$coef.names,sep=""),paste("control.phi1.",object$coef.names,sep=""))
		} else {
			rownames(vcov) <- colnames(vcov) <- c(paste("treat.",object$coef.names,sep=""), paste("control.",object$coef.names,sep=""))
		}
	}
	
	return(vcov)

}

summary.ictreg <- function(object, ...)
	structure(object, class = c("summary.ictreg", class(object)))
 
print.summary.ictreg <- function(x, ...){

	cat("\nItem Count Technique Regression \n\nCall: ")
	
	dput(x$call)
	
	cat("\n")
		
	if(x$method=="nls"){
		cat(rep(" ", max(nchar(x$coef.names))+7), "Treatment", rep(" ", 15), "Control \n", sep="")
		tb <- matrix(NA, ncol = 4, nrow = length(x$par.treat))
		colnames(tb) <- rep(c("Est.", "S.E."), 2)
		rownames(tb) <- x$coef.names
		tb[,1] <- x$par.treat
		tb[,2] <- x$se.treat
		tb[,3] <- x$par.control
		tb[,4] <- x$se.control
		
		summ.stat <- paste("Residual standard error:",signif(x$resid.se, digits=6),"with",x$resid.df,"degrees of freedom")
		
	} else if(x$method=="ml") {
	
		if(x$constrained==F) {
			cat(rep(" ", max(nchar(x$coef.names))+7), "Treatment", rep(" ", 7), "Control (phi0)", rep(" ",7), "Control (phi1)\n", sep="")
			tb <- matrix(NA, ncol = 6, nrow = length(x$par.treat))
			colnames(tb) <- rep(c("Est.", "S.E."), 3)
			rownames(tb) <- x$coef.names
			tb[,1] <- x$par.treat
			tb[,2] <- x$se.treat
			tb[,3] <- x$par.control.phi0
			tb[,4] <- x$se.control.phi0
			tb[,5] <- x$par.control.phi1
			tb[,6] <- x$se.control.phi1
		} else {
			cat(rep(" ", max(nchar(x$coef.names))+7), "Treatment", rep(" ", 15), "Control \n", sep="")
			tb <- matrix(NA, ncol = 4, nrow = length(x$par.treat))
			colnames(tb) <- rep(c("Est.", "S.E."), 2)
			rownames(tb) <- x$coef.names
			tb[,1] <- x$par.treat
			tb[,2] <- x$se.treat
			tb[,3] <- x$par.control
			tb[,4] <- x$se.control
		}
		
		summ.stat <- paste("Log-likelihood:", x$llik)
		
	}
	
	print(as.matrix(tb))
			
	cat("\n",summ.stat,"\n\n",sep="")
	
	invisible(x)
	
}
