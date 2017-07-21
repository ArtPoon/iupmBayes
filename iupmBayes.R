require(MCMCpack)  # for ddirichlet()


# likelihood function
bern <- function(well, rate, probs) {
  # @arg well: binary vector indicating presence/absence of variants in well
  # @arg rate: Poisson rate, IUPM * N (number of cells in well)
  # @arg probs: a vector of probabilities for each variant
  
  # some sanity checks
  if (!is.vector(well) | any(well<0)) {
    stop("well argument must be a vector of non-negative values")
  }
  if (!is.vector(probs) | any(probs<0)) {
    stop("probs argument must be a vector of non-negative values")
  }
  if (length(well) != length(probs)) {
    stop("well and probs arguments must have the same length")
  }
  if (sum(probs)!=1) {
    probs <- probs/sum(probs)
  }
  pr <- exp(-rate*probs)
  prod(ifelse(well, 1-pr, pr))
}

ll <- function(data, params) {
  # ll:  log-likelihood
  # @arg data: S3 object containing:
  #   wells: matrix of K binary vectors (rows) for presence/absence
  #          of M variants
  #   cells: vector of length K, the number of cells in the k-th well
  # @arg params: S3 object containing:
  #   iupm: IUPM
  #   probs: vector of probabilities of length M, representing the relative
  #          proportion of cells infected by each variant
  #          The first value is fixed to 1.  All other values can vary 
  #          from -Inf to Inf, and is transformed to a positive value by 
  #          exp() and rescaled so that all probs sum to 1.
  p.vec <- c(1, params$probs)
  p.vec <- p.vec/ sum(p.vec)
  sum(sapply(1:nrow(data$wells), function(i) {
    wells <- data$wells[i,]
    cells <- data$cells[i]
    log(bern(wells, params$iupm/1e6*cells, p.vec))
  }))
}


iupm3.ml <- function(variants, cells, params) {
  # maximum-likelihood estimation
  data <- list(wells=variants, cells=cells)
  init.p <- c(params$iupm, params$probs)
  n <- length(init.p)
  obj.f <- function(p) {
    pv <- list(iupm=p[1], probs=c(p[2:(n-1)], 1-sum(p[2:(n-1)])))
    -ll(data, pv)
  }
  
  optim(par=init.p, fn=obj.f, lower=rep(1e-5,n), 
        upper=c(1e3, rep(1,n)), method='L-BFGS-B')
}


lprior <- function(params, hyper, use.gamma=FALSE) {
  # lprior:  log prior
  # @arg params: S3 object as in ll()
  # @arg hyper: S3 object containing:
  #   shape: shape hyperparameter for gamma prior (shape=1 gives exponential)
  #   rate: rate hyperparameter for gamma prior
  #   alpha: vector of length M, hyperparameters for Dirichlet prior
  pv <- c(1, params$probs)
  log(ddirichlet(pv/sum(pv), hyper$alpha)) + ifelse(use.gamma,
      dgamma(params$iupm, shape=hyper$shape, rate=hyper$rate, log=T),
      0)
}


propose <- function(params, sd=0.1, delta=0.005) {
  # shift IUPM by small amount, reflecting at zero bound
  step <- rnorm(1, 0, sd)
  params$iupm <- abs(params$iupm + step)
  # modify probability vector
  n.var <- length(params$probs)  # note first variant parameter is constrained to 1
  params$probs <- abs(params$probs + runif(n.var, -delta, delta))
  return(params)
}



# Metropolis-Hastings sampling
mh <- function(data, params, hyper, max.steps, logfile=NA, skip.console=100, skip.log=100, overwrite=FALSE) {
  res <- c()
  m <- length(params$probs) + 1  # number of variants
  
  labels <- c('step', 'posterior', 'iupm', paste0('prob', 1:m))
  if (!is.na(logfile)) {
    # prevent overwriting logfiles
    if (file.exists(logfile) & !overwrite)
      stop("To overwrite logfiles, set overwrite to TRUE.")
    write(labels, ncolumns=3+m, sep='\t', file=logfile)
  }
  
  lpost <- ll(data, params) + lprior(params, hyper)
  for (step in 0:max.steps) {
    next.par <- propose(params)
    next.lpost <- ll(data, next.par) + lprior(next.par, hyper)
    ratio <- next.lpost - lpost
    if (ratio >= 0 | runif(1) < exp(ratio)) {
      lpost <- next.lpost
      params <- next.par
    }
    
    # logging
    if (step %% skip.console == 0) {
      cat(step, '\t', lpost, '\t', params$iupm, '\t', params$probs, '\n')
    }
    if (step %% skip.log == 0) {
      pv <- c(1, params$probs)
      pv <- pv/sum(pv)
      row <- c(step, lpost, params$iupm, pv)
      if (!is.na(logfile)) {
        write(row, ncolumns=3+m, sep='\t', file=logfile, append=T)
      }
      res <- rbind(res, row)
    }
  }
  res <- as.data.frame(res)
  names(res) <- labels
  return(res)
}


run.example <- function() {
  # example
  params <- list(iupm=1.0, probs=c(1,1))
  hyper <- list(rate=1, shape=2, alpha=c(1,2,3))
  
  data <- list(wells=matrix(c(1,0,0, 1,1,0, 1,0,1, 1,0,0), ncol=3, byrow=T), cells=rep(1e5, 4))
  
  mh(data, params, hyper, 1e5, 'test.log', skip.console=100, skip.log=100, overwrite=T)
  trace <- read.table('test.log', sep='\t', header=T)
}


parse.data <- function(path, sep) {
  # Read tabular data from text file in expected format (see README.md).
  #
  # @arg path:  Absolute or relative path to data file
  # @arg sep:   Delimiting character in tabular data
  # Return:  a list keyed by region, containing data lists
  data <- read.table(path, header=T, sep=sep, comment.char='')
  
  result <- list()
  by.region <- split(data, data$Region)
  for (temp in by.region) {
    region <- unique(temp$Region)
    temp2 <- temp[order(temp$Well.number, temp$Variant..), ]
    
    wells <- split(temp2$Presence.of.variant, temp2$Well.number)
    wells <- t(sapply(wells, unlist))  # convert into matrix
    
    cells <- lapply(split(temp2$Cells.plated, temp2$Well.number), unique)
    cells <- unlist(cells)
    
    if (!all(row.names(wells) == names(cells))) {
      stop("Failed to parse well and cell data")
    }
    eval(parse(text=paste("result$", region, "<-list('data'=list('wells'=wells, 'cells'=cells))")))
  }
  return(result)
}



simulate.data <- function(iupm, probs, cells) {
  # iupm: rate per 10E6
  # probs: vector of probabilities (relative abundance of each variant)
  # cells: vector of number of cells per well
  n <- length(cells)
  n.var <- length(probs)  # number of variants to simulate
  #res <- matrix(NA, nrow=n, ncol=n.var)
  
  # number of infected cells per well
  n.inf <- rpois(n, iupm/1e6 * cells)
  res <- t(sapply(n.inf, function(m) {
    z <- sample(1:n.var, size=m, replace=TRUE, prob=probs)
    as.integer(is.element(1:n.var, z))
  }))
  
  # censor absent variants
  as.matrix(res[,apply(res, 2, sum)>0])
}


iupm.ngs <- function(variants) {
  # from SK Lee et al. (2017) JAIDS 74(2): 221-228
  # number of wells expressing i-th variant
  # note this estimator assumes that all wells have the same number of cells
  positives <- apply(variants, 2, sum)
  mle <- -log(1-positives/nrow(variants))
  return(sum(mle))
}

