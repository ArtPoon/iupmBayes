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
  # @arg data: a nested list containing:
  #     region: genomic region sequenced
  #       wells: list of K binary vectors (rows) for presence/absence
  #              of M variants, indexed by well label
  #       cells: vector of length K, the number of cells in the k-th well,
  #              indexed by well label
  # @arg params: S3 object containing:
  #   iupm: IUPM
  #   probs: vector of probabilities of length M, representing the relative
  #          proportion of cells infected by each variant
  #          The first value is fixed to 1.  All other values can vary 
  #          from -Inf to Inf, and is transformed to a positive value by 
  #          exp() and rescaled so that all probs sum to 1.

  res <- sapply(1:length(data), function(i) {
    region <- names(data)[i]
    rdat <- data[[region]]
    
    # unpack variant frequency vector
    var.freq <- params$freqs[[region]])
    
    # for each well, apply Bernoulli distribution
    sum(sapply(1:nrow(rdat$wells), function(j) {
      wells <- rdat$wells[j,]
      cells <- rdat$cells[j]
      log(bern(wells, params$iupm/1e6*cells, var.freq))
    }))
  })
  
  sum(res)
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
  #   alpha: a list of vectors of length M, hyperparameters for Dirichlet prior,
  #          keyed by region name
  res <- 0
  if (use.gamma) {
    dgamma(params$iupm, shape=hyper$shape, rate=hyper$rate, log=T)
  }
  
  for (i in 1:length(params$freqs)) {
    region <- names(params$freqs)[i]
    pv <- params$freqs[region]
    alpha <- hyper$alpha[[region]]
    if (length(alpha) != length(pv)) {
      stop(paste("Error: length of alpha hyperparameter != variant frequency vector for", region))
    }
    res <- res + log(ddirichlet(pv, alpha))
  }
  return(res)
}


propose <- function(params, sd=0.1, delta=0.005) {
  # shift IUPM by small amount, reflecting at zero bound
  step <- rnorm(1, 0, sd)
  params$iupm <- abs(params$iupm + step)
  
  # modify probability vectors
  for (i in 1:length(params$freqs)) {
    region <- names(params$freqs)[i]
    n.var <- length(params$freqs[[region]])
    temp <- abs(params$freqs[[region]] + runif(n.var, -delta, delta))
    params$freqs[[region]] <- temp / sum(temp)  # normalize
  }
  return(params)
}



mh <- function(data, iupm0, hyper, max.steps, logfile=NA, skip.console=100, 
               skip.log=100, overwrite=FALSE) {
  # Metropolis-Hastings sampling
  #
  # @param data: S3 object returned by parse.data()
  # @param iupm0: initial IUPM value
  # @param hyper: nested S3 object, see lprior() above
  # @param max.steps:  number of iterations to run chain sample
  # @param logfile:  path to write log
  # @param skip.console:  number of iterations between print-outs
  # @param skip.log:  number of iterations between chain samples to logfile
  # @param overwrite:  if TRUE, replace contents of logfile
  res <- c()
  m <- sapply(data, function(sub) ncol(sub$wells))  # number of variants
  labels <- c('step', 'posterior', 'iupm')
  
  # initialize model parameters
  params <- list(iupm=iupm0, freqs=list())
  regions <- c()
  for (i in 1:length(m)) {
    region <- names(m)[i]
    regions <- c(regions, region)
    nvar <- m[i]
    
    # we reduce the number of parameters by one because the frequency
    # vector is constrained to sum to one
    params[["freqs"]][[region]] <- rep(1, nvar-1)
    labels <- c(labels, paste(region, 1:nvar, sep='.'))
  }
  
  # prepare output file
  if (!is.na(logfile)) {
    # prevent overwriting logfiles
    if (file.exists(logfile) & !overwrite)
      stop("To overwrite logfiles, set overwrite to TRUE.")
    write(labels, ncolumns=length(labels), sep='\t', file=logfile)
  }
  
  # calculate initial posterior probability
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
      cat(step, '\t', lpost, '\t', params$iupm)
      for (region in regions) {
        cat('\t', params$freqs[[region]])
      }
      cat('\n')  # EOL
    }
    if (step %% skip.log == 0) {
      row <- c(step, lpost, params$iupm)
      
      # concatenate frequency vectors
      for (region in regions) {
        row <- c(row, params$freqs[[region]])
      }
      
      if (!is.na(logfile)) {
        write(row, ncolumns=length(row), sep='\t', file=logfile, append=T)
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
  params <- list(iupm=1.0, probs=c(1,1))  # initial values
  hyper <- list(rate=1, shape=2, alpha=c(1,2,3))  # hyperparameters for priors
  
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
  data$Participant <- as.character(data$Participant)
  
  result <- list()  # prepare output container
  
  by.subject <- split(data, data$Participant)
  for (part in by.subject) {
    participant <- unique(part$Participant)
    result[[participant]] <- list()
    
    by.region <- split(part, as.character(part$Region))
    for (temp in by.region) {
      region <- as.character(unique(temp$Region))
      
      temp$Well.number <- as.character(temp$Well.number)
      temp2 <- temp[order(temp$Well.number, temp$Variant), ]
      
      wells <- split(temp2$Presence.of.Variant, temp2$Well.number)
      wells <- t(sapply(wells, unlist))  # convert into matrix
      if (nrow(wells) == 1) {
        wells <- t(wells)
      }
      
      cells <- lapply(split(temp2$Cells.plated, temp2$Well.number), unique)
      cells <- unlist(cells)
      
      if (!setequal(row.names(wells), names(cells))) {
        print(wells)
        print(cells)
        stop("Failed to parse well and cell data")
      }
      
      result[[participant]][[region]] <- list('wells'=wells, 'cells'=cells)
      #eval(parse(text=paste("result$", region, "<-list('data'=list('wells'=wells, 'cells'=cells))")))
    }
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

