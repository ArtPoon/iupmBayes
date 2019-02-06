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
  if (any(probs<0)) {
    stop("probs argument must be a vector of non-negative values")
  }
  if (length(well) != length(probs)) {
    stop("well and probs arguments must have the same length")
  }
  if (sum(probs)!=1) {
    probs <- probs/sum(probs)
  }
  pr <- exp(-rate*probs)  # probability of count zero
  prod(ifelse(well, 1-pr, pr))
}


ll <- function(data, params) {
  # ll:  log-likelihood
  # @arg data: a nested list containing:
  #   region: genomic region sequenced
  #   wells: binary matrix with K rows for replicate wells with row 
  #          names and M columns for variants
  #   cells: named vector of length K, the number of cells in the k-th well;
  #          names must match row names of :wells: matrix
  # @arg params: S3 object containing:
  #   iupm: IUPM
  #   probs: vector of probabilities of length M, representing the relative
  #          proportion of cells infected by each variant
  res <- sapply(1:length(data), function(i) {
    region <- names(data)[i]
    rdat <- data[[region]]
    var.freq <- params$freqs[[region]]
    
    # for each well, apply Bernoulli distribution
    sum(sapply(1:nrow(rdat$wells), function(j) {
      wells <- rdat$wells[j,]
      cells <- rdat$cells[j]
      log(bern(wells, params$iupm/1e6*cells, var.freq))
    }))
  })
  
  sum(res)
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
    res <- dgamma(params$iupm, shape=hyper$shape, rate=hyper$rate, log=T)
  }
  
  for (i in 1:length(params$freqs)) {
    region <- names(params$freqs)[i]
    pv <- params$freqs[[region]]
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


make.hyper <- function(data, shape=NA, rate=NA) {
  # convenience function for specifying hyperparameters
  # alpha set to rep(1, times=<number of variants>)
  hyper <- list(shape=shape, rate=rate, alpha=list())
  if (is.na(shape) || is.na(rate)) {
    # a broad prior centered (mode) on 1 with 95% interval 0.12-18.4
    hyper['shape'] <- 1.
    hyper['rate'] <- 0.2
  }
  for (i in 1:length(data)) {
    region <- names(data)[i]
    part <- data[[region]]
    n.var <- ncol(part$wells)
    hyper[['alpha']][[region]] <- rep(1, times=n.var)
  }
  return(hyper)
}


mh <- function(data, params=list(), hyper=list(), max.steps=1e5, logfile=NA, skip.console=1000, 
               skip.log=1000, overwrite=FALSE) {
  # Metropolis-Hastings sampling
  #
  # @param data: S3 object returned by parse.data() for one participant
  # @param iupm0: initial IUPM value
  # @param hyper: nested S3 object, see lprior() above
  # @param max.steps:  number of iterations to run chain sample
  # @param logfile:  path to write log
  # @param skip.console:  number of iterations between print-outs
  # @param skip.log:  number of iterations between chain samples to logfile
  # @param overwrite:  if TRUE, replace contents of logfile
  res <- c()
  # number of variants
  m <- lapply(data, function(sub) ncol(sub$wells))
  regions <- names(m)
  labels <- c('step', 'posterior', 'iupm')
  
  # check params
  if (!is.list(params)) {
    stop('params must be list()')
  }
  
  # set default model parameters
  if (length(params)==0) {
    params <- list(iupm=5, freqs=list())
    for (i in 1:length(m)) {
      region <- names(m)[i]
      nvar <- as.integer(m[i])
      params[["freqs"]][[region]] <- rep(1/nvar, nvar)
      labels <- c(labels, paste(region, 1:nvar, sep='.'))
    }
  } else {
    if (!all(is.element(c('iupm', 'freqs'), names(params)))) {
      stop('params list must have iupm and freqs entries')
    }
  }
  
  # use default hyperparameters if not set
  if (length(hyper)==0) {
    hyper <- list(alpha=list())
    for (region in regions) {
      hyper[['alpha']][[region]] <- rep(1, times=as.integer(m[region]))
    }
    use.gamma <- FALSE
  } else {
    if (!all(is.element(c('shape', 'rate', 'alpha'), names(hyper)))) {
      stop('hyper list must have shape, rate and alpha entries')
    }
    use.gamma <- !is.na(hyper$rate) & !is.na(hyper$shape)
  }
  
  # prepare output file
  if (!is.na(logfile)) {
    # prevent overwriting logfiles
    if (file.exists(logfile) & !overwrite)
      stop("To overwrite logfiles, set overwrite to TRUE.")
    write(labels, ncolumns=length(labels), sep='\t', file=logfile)
  }
  
  # calculate initial posterior probability
  lpost <- ll(data, params) + lprior(params, hyper, use.gamma)
  
  for (step in 0:max.steps) {
    next.par <- propose(params)
    next.lpost <- ll(data, next.par) + lprior(next.par, hyper, use.gamma)
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



parse.data <- function(path, sep) {
  # Read tabular data from text file in expected format (see README.md).
  # @arg path:  Absolute or relative path to data file
  # @arg sep:   Delimiting character in tabular data
  # Return:  a nested list keyed by subject and region, containing data as follows:
  #   $Participant1
  #   $Participant1$region1
  #   $Participant1$region1$wells
  #       [,1] [,2] [,3]
  #   1A     1    0    1
  #   1B     1    0    1
  #   1C     0    1    0
  #   $Participant$region1$cells
  #        1A      1B      1C
  #   1000000 1000000 1000000
  data <- read.table(path, header=T, sep=sep, comment.char='')
  
  # remove extra characters from column names
  #names(data) <- gsub('^\.+(.+)\.+$', "\\1", names(data))
  
  data$Participant <- as.character(data$Participant)
  
  result <- list()  # prepare output container
  
  by.subject <- split(data, data$Participant)
  for (part in by.subject) {
    participant <- as.character(unique(part$Participant))
    result[[participant]] <- list()
    
    by.region <- split(part, as.character(part$Region))
    for (temp in by.region) {
      region <- as.character(unique(temp$Region))
      
      temp$Well.number <- as.character(temp$Well.number)
      temp$Variant <- as.integer(temp$Variant)
      
      temp2 <- temp[order(temp$Well.number, temp$Variant), ]
      
      wells <- split(temp2$Presence.of.Variant, temp2$Well.number)
      wells <- t(sapply(wells, unlist))  # convert into matrix
      if (nrow(wells) == 1 && is.null(row.names(wells))) {
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
  list(wells=as.matrix(res[,apply(res, 2, sum)>0]), cells=cells)
}



run.example <- function(iupm0, iupm, shape=1, rate=1, seed=1) {
  set.seed(seed)
  data1 <- simulate.data(iupm=iupm, probs=c(.5,.2,.1,.1,.1), cells=rep(1e6,8))
  data2 <- simulate.data(iupm=iupm, probs=c(.3,.1,rep(.05,12)), cells=rep(1e6,6))
  joint.data <- list(region1=data1, region2=data2)
  
  # initialize parameters to random values
  nvar1 <- ncol(data1$wells)
  nvar2 <- ncol(data2$wells)
  params <- list(iupm=iupm0, freqs=list(
    region1=rdirichlet(1, rep(1, nvar1))[1,],
    region2=rdirichlet(1, rep(1, nvar2))[1,]
  ))
  
  hyper <- list(shape=shape, rate=rate, alpha=list(
    region1=rep(1, nvar1), region2=rep(1, nvar2)
  ))
  mh(joint.data, params, hyper=hyper, max.steps=1e5, skip.console=1000, skip.log=100)
}



iupm.ngs <- function(variants) {
  # from SK Lee et al. (2017) JAIDS 74(2): 221-228
  # number of wells expressing i-th variant
  # note this estimator assumes that all wells have the same number of cells
  positives <- apply(variants, 2, sum)
  mle <- -log(1-positives/nrow(variants))
  return(sum(mle))
}


