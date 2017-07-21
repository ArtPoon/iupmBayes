# examine the effect of partitioning positive cells

iupm.ll <- function(rate, positives, cell.counts) {
  # log likelihood for single-hit Poisson model
  pr <- exp(-1e-6*rate*cell.counts)
  pr[pr==0] <- 1e-320
  sum( positives*log(1-pr) + (1-positives)*log(pr) )
}

iupm <- function(positives, n.cells, lower=1e-3, upper=1e3, hessian=FALSE) {
  # Implementation of IUPMStats calculator
  # 
  # @arg positives: a binary vector indicating if i-th well was positive
  # @arg n.cells: an integer vector for number of cells in i-th well
  # @arg lower: IUPM lower bound for minimization
  # @arg upper: IUPM upper bound for minimization
  # @arg hessian: optim() optional return Hessian matrix
  
  if (!is.vector(n.cells) | !is.vector(positives)) {
    stop("all arguments must be vectors\n")
  }
  n <- length(n.cells)  # number of wells
  if (length(positives)!=n) {
    stop("n.cells must have same length as positives\n")
  }
  if (any(n.cells<0) | any(positives<0)) {
    stop("n.cells and positives must be non-negative\n")
  }
  
  # are all wells negative?
  if (all(positives==0)) {
    # return posterior median estimate
    res <- list(par=log(2)/sum(n.cells) * 1e6, lo.95=NA, hi.95=NA)
    return(res)
  }
  
  # reverse sign for objective function minimization
  obj.f <- function(rate) -iupm.ll(rate, positives, n.cells)
  res <- optim(par=1, fn=obj.f, method='Brent', lower=lower, upper=upper, hessian=hessian)

  # # calculate 95% confidence interval
  # l1 <- iupm.ll(res$par, positives, n.cells)
  # obj.f <- function(x) exp(iupm.ll(x, positives, n.cells) - l1) - 0.147
  # res$lo.95 <- NA
  # res$hi.95 <- NA
  # try({
  #   res$lo.95 <- uniroot(obj.f, c(1e-3, res$par))$root
  #   res$hi.95 <- uniroot(obj.f, c(res$par, 1e2))$root
  # })
  
  # asymptotic normal confidence interval
  lmle <- log(res$par*1e-6)
  term1 <- exp(-n.cells*res$par*1e-6) / (exp(-n.cells*res$par*1e-6)-1)
  d2 <- sum(n.cells^2 * positives * term1 * (1-term1))
  shift <- 1.96*1/sqrt(-d2)
  
  return(res)
}


# example, IUPM = 8.148 (8.147718 {1.862907, 35.635330})
run.example <- function() {
  cell.counts <- rep(c(1e6, 2e5, 4e4, 8000, 1600, 320), times=c(8, 2, 2, 2, 2, 2))
  positives <- c(rep(1,8), rep(1,2), rep(0,2), rep(0,2), rep(0,2), rep(0,2))
  iupm(positives, cell.counts)
}



# Is it possible to use this variant information in another way?
# For example, observing two variants implies that at least two cells are infected.
# The likelihood of the positive outcome changes from: 
#   1 - exp(-rN) 
# where r=IUPM and N is the number of cells in the well, to 
#   1 - exp(-rN) - rN exp(-rN)  = 1 - ppois(1,8)
# where ppois() is the Poisson cumulative distribution function.
# In other words, we've removed the case where only one cell is infected.

lcup <- function(k,L) {
  # log cumulative distribution function of Poisson
  # e.g., lcup(1,8) = log(1-ppois(1,8))
  if (k==0) {
    return(-L)
  } else {
    pgamma(L, floor(k), log=T) - lgamma(k)
  }
}

iupm.ll2 <- function(rate, variants, cells) {
  # log-likelihood for multi-hit Poisson model
  n <- length(variants)
  L <- 1e-6*rate*cells
  #sum(sapply(1:n, function(i) lcup(variants[i], L[i])))
  sum(sapply(1:n, function(i) {
    k <- variants[i]
    ifelse( k==0, -L[i], log(1-ppois(k-1, L[i])) )
  }))
}

iupm.2 <- function(n.variants, n.cells, upper=1e3, hessian=FALSE) {
  # n.variants:  a vector of number of variants observed per well
  # n.cells: a vector of number of cells per well
  if (length(n.variants) != length(n.cells)) {
    stop("Length of n.variants must match length of n.cells")
  }
  
  obj.f <- function(rate) { -iupm.ll2(rate, n.variants, n.cells) }
  res <- optim(par=1, fn=obj.f, method='Brent', lower=0, upper=upper, hessian=hessian)
  
  # calculate 95% confidence interval
  l1 <- iupm.ll2(res$par, n.variants, n.cells)
  obj.f <- function(x) exp(iupm.ll2(x, positives, n.cells) - l1) - 0.147
  try({
    res$lo.95 <- uniroot(obj.f, c(1e-3, res$par))$root
    res$hi.95 <- uniroot(obj.f, c(res$par, 1e2))$root
  })
  return(res)
}





# TODO: write function to coerce data frame into this list format
run.example2 <- function() {
  variants <- list('1e6'=c(1,2,1,1,1,1,2,2), '2e5'=c(1,1), '4e4'=c(0,0), '8000'=c(0,0), '1600'=c(0,0), '320'=c(0,0))
  iupm.2(variants)
}



