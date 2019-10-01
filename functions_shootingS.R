
###############################################################
# The shooting S-estimator                                    #
# A regression method for cellwise contamination              #
#                                                             #
# Reference: Oellerer V., Alfons A., Croux C.                 #
#            'The shooting S-estimator for robust regression' #
#                                                             #
# Author: Viktoria Oellerer                                   #
#         KU Leuven                                           #
###############################################################


require("mvtnorm")
require('robustHD')
shooting <- function(x, y, tol = 10^(-2), maxIteration = 100, k=3.420, method='biweight') {
  # shooting S-estimator computed either for biweight loss or skipped Huber loss
  # INPUT:  x: design matrix (dimension nxp)
  #         y: response variable (length n)
  #         tol: numerical tolerance for convergence
  #         maxIteration: maximum number of iterations in coordinate descent loop
  #         k: tuning parameter of rho-function
  #                 (default to achieve breakdown point 20% - biweight: 3.420, skipped Huber: 2.176922)
  #         method: name of method for which shooting algorithm should be used
  #                 ('biweight': shooting S-estimator using biweight loss function,
  #                  'skHuber': shooting S-estimator using skipped Huber loss function)
  # OUTPUT: coef: coefficient estimates of the shooting S-estimator (including intercept -> length p+1)
  #         resid: residuals of the shooting S-estimator (length n)
  #         scaleVar: vector of scale estimates of residuals (one result per variable -> length p)
  #         weights: weight matrix (dimension nxp)
  #         nonconv: boolean indicating if algorithm converged (nonconv=0) or not converged (nonconv=1)


  ## initializations
  x <- as.matrix(x)
  n <- nrow(x)
  p <- ncol(x)
  if (method == 'biweight') {
    rho <- match.fun('rho_bi')
    weights <- match.fun('w_bi')
    delta <- (1 - 3/k^2 + 5/k^4 - k^2/3) *pnorm(k) + (4/(3*k) - k/3 - 5/k^3)*dnorm(k) - 1/2 + 3/(2*k^2) - 5/(2*k^4) + k^2/3
  } else if(method == 'skHuber') {
    rho <- match.fun('rho_skh')
    weights <- match.fun('w_skh')
    delta <- (1-k^2)*pnorm(k) - k*dnorm(k)+k^2-1/2
  }

  # initial estimate of coefficients
  fitini <- lmrob.hubx(x,y)
  betaEst <- fitini$coeforig[-1]
  intercept <- rep.int(fitini$coeforig[1], p)
  scaleVar <- rep.int(fitini$scale, p)
  scaleVarmat <- rbind(rep.int(Inf, p), scaleVar)


  # compute Huberized x
  xshub <- apply(x,2,function(z) pmax(median(z)-2*mad(z), pmin(z,median(z)+2*mad(z))))
  xtilde <- xshub

  # construct error terms y - x_1*beta_1 -.... - x_p*beta_p
  ytilde <- y - xtilde%*%betaEst

  # initialize weight matrix and xhat matrix (the imputed x's)
  wt <- xhat <- matrix(NA, ncol=p, nrow=n)

  # draw a new sequence in which loop has to be done
  sequencevector <- sample(1:p)

  ## loop iteratively through all variables until convergence of estimate (coordinate descent loop)
  r <- 0
  while(sum(abs(scaleVar-scaleVarmat[nrow(scaleVarmat)-1,])) > tol*mad(y) && r < maxIteration) {
    r <- r+1

    for(j in sequencevector) {
      # Regression step
      ytilde <- ytilde + xtilde[, j] * betaEst[j]
      resid <- drop(ytilde - x[,j] * betaEst[j] - median(ytilde - x[,j] * betaEst[j]))
      wt[,j] <- sapply(resid/scaleVar[j], weights, k=k)
      fit <- univ_est(x[, j, drop=FALSE], ytilde, intercept=TRUE, k = k, rho = rho, weights = weights, delta = delta,
                      wt = wt[,j], tol= 10^(-6)*mad(y), prevresid=resid)
      # regression fit for variable j (I-steps done inside univ_est)
      # update coefficient and intercept
      betaEst[j] <- fit$coef[2]
      intercept[j] <- fit$coef[1]
      # update scale, weights and xtilde
      scaleVar[j] <- fit$s
      wt[,j] <- fit$wt
      xhat[,j] <- median(x[,j])
      if (abs(betaEst[j]) > 10^(-4)*mad(y)/mad(x[,j])) xhat[,j] <- (ytilde-intercept[j])/betaEst[j]
      lambda <- (wt[,j]>weights(3,k))*1
      xtilde[,j] <- xhat[,j] + (x[,j] - xhat[,j])*lambda
      # update response
      ytilde <- ytilde - xtilde[, j] * betaEst[j]
    }
    scaleVarmat <- rbind(scaleVarmat, scaleVar) # save new scale estimates
  }

  # give warning if maximum number of iterations reached:
  if(r==maxIteration) {
    warning(paste('maximum number of iterations reached in shooting()'))
    nonconv <- 1
  } else {
    nonconv <- 0
  }

  # recompute intercept
  alpha <- median(ytilde)

  # return estimates
  list('coef' = c(alpha, betaEst), 'resid' = resid, 'scaleVar' = scaleVar, 'weights' = wt, 'nonconv' = nonconv)
}

univ_est <- function(x, y, k, weights, rho, wt, delta, maxIteration = 500, intercept=TRUE, tol, prevresiduals) {
  # computation of one step of coordinate descent loop (i.e. all coefficients except beta_j are held fixed)
  # INPUT:  x: current variable j (n)
  #         y: response variable (n)
  #         k: tuning parameter for rho-function
  #         weights: name of weight-function to use (w_bi, w_skh)
  #         rho: name of rho-function to use (rho_bi, rho_skh)
  #         wt: vector of initial weights to use in weighted least squares regression (n)
  #         delta: constant of consistency of M-scale (E[Z]=delta)
  #         maxIteration: maximum number of iterations in iteratively reweighted least squares (= maximum number of I-steps)
  #         intercept: boolean indicating if intercept should be added to regression (needs to be TRUE, if data is not standardized)
  #         tol: numerical tolerance for convergence
  #         prevresiduals: current residuals (y-x%*%betaEst-intercept)
  # OUTPUT: fit: last regression fit in iteratevely reweighted least squares (IRLS) algorithm
  #         s: scale estimate of residuals in last regression of IRLS
  #         wt: estimated weights in last regression of IRLS

  # initializations
  n <- length(x)
  rho <- match.fun(rho)
  weights <- match.fun(weights)

  # add column for intercept if requested
  x <- if(isTRUE(intercept)) cbind(rep.int(1, n), x) else as.matrix(x)

  # IRLS algorithm
  residuals <- rep(-Inf,n)
  r <- 0
  while(max(abs(prevresiduals - residuals)) > tol && r < maxIteration) {
    r <- r+1
    prevresiduals <- residuals
    # fit weighted least squares and extract coefficients
    fit <- lm.wfit(x, y, wt)
    # extract residuals and compute residual scale
    residuals <- fit$residuals
    if (r==1){
      s <- scale_iter(residuals, kp = delta, cc = k, rho = rho)
    } else {
      s <- scale_iter(residuals, kp = delta, cc = k, initial.sc = s, rho = rho)
    }
    wt <- sapply(residuals/s, weights, k = k) # update weights
  }

  # give warning if maximum number of iterations reached
  if(r==maxIteration) {
    warning('maximum number of iterations reached in univ_est()')
  }

  # return results
  fit <- c(fit, list('s'=s, 'wt'=wt))
  fit
}

scale_iter <- function(u, kp, cc, initial.sc=median(abs(u))/.6745, rho, max.it = 5000, tol = 1e-6) {
  # computation of M-scale estimate (see Salibian-Barrera and Yohai 2006 "A fast algorithm for S-Regression Estimates")
  # INPUT:  u: current residuals for which M-scale should be computed
  #         kp: constant for consistency of M-scale (E[Z]=kp)
  #         cc: tuning parameter for rho-function
  #         initial.sc: initial scale estimate to be updated
  #         rho: name of rho-function to use (rho_bi, rho_skh)
  #         max.it: maximum number of fixed point iterations
  #         tol: numerical tolerance for convergence
  # OUTPUT: M-scale estimate of residuals 'u'

  # initializations
  rho <- match.fun(rho)
  sc <- initial.sc
  i <- 0
  err <- 1
  # perform fixed point iteration
  while( ( (i <- i+1) < max.it ) && (err > tol) ) {
    sc2 <- sqrt( sc^2 * mean(sapply(u/sc,rho, k = cc)) / kp  )
    err <- abs(sc2/sc - 1)
    sc <- sc2
  }

  # give warning if maximum number of iterations reached
  if(i == max.it-1) warning('Maximum number of iterations reached in scale_iter()')

  # return M-scale estimate
  return(sc)
}

lmrob.hubx <- function(x,y, chub=2) {
  # computation of MM-estimate of Huberized data
  # INPUT: x: design matrix (dimension nxp)
  #        y: response variable (length n)
  #        chub: tuning parameter for Huberizing data
  # OUTPUT: rr: estimated coefficients (including intercept)

  # standardize and Huberize data
  xx <- robStandardize(x) # standardization
  mx <- attr(xx, "center")
  sx <- attr(xx, "scale")
  xh <- data.frame(apply(xx,2,function(x) pmin(pmax(x,-chub),chub))) # Huberization
  rr <- lmrob(y~.,data=cbind(y,xh), setting="KS2011", k.max = 100000,
              max.it = 5000, maxit.scale = 1000000) # MM-estimation

  # backtransform
  coeforig <- rr$coef[-1]/sx
  coeforigalpha <- rr$coef[1] - sum(coeforig * mx)

  # return estimates
  rr$coeforig <- c(coeforigalpha,coeforig)
  rr
}

w_bi <- function(z, k) {
  # biweight weight function for a single input value
  # INPUT:  z: numerical value
  #         k: tuning parameter
  # OUTPUT: value of biweight weight function for input value z
  if (abs(z) <= k) {
    (1-(z/k)^2)^2
  } else {
    0
  }
}

rho_bi <- function(z, k) {
  # biweight rho function for a single input value
  # INPUT:  z: numerical value
  #         k: tuning parameter
  # OUTPUT: value of biweight rho function for input value z
  if (abs(z) <= k) {
    (1-(1-(z/k)^2)^3)*k^2/6
  } else {
    k^2/6
  }
}

w_skh <- function(z, k) {
  # skipped Huber weight function for a single input value
  # INPUT:  z: numerical value
  #         k: tuning parameter
  # OUTPUT: value of skipped Huber weight function for input value z
  if (abs(z) <= k) {
    1
  } else {
    0
  }
}

rho_skh <- function(z, k) {
  # skipped Huber rho function for a single input value
  # INPUT:  z: numerical value
  #         k: tuning parameter
  # OUTPUT: value of skipped Huber rho function for input value z
  if (abs(z) <= k) {
    1/2*z^2
  } else {
    k^2/2
  }
}
