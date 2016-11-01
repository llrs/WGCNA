# Function to control the number of threads to use in threaded calculations.

.useNThreads = function(nThreads = 0) {
  if (nThreads == 0) {
    nt.env = Sys.getenv(.threadAllowVar, unset = NA)
    if (is.na(nt.env)) {
        return(1)
    }
    if (nt.env == "") {
        return(1)
    }

    if (nt.env == "ALL_PROCESSORS") {
        return(.nProcessorsOnline())
    }

    nt = suppressWarnings(as.numeric(nt.env))

    if (!is.finite(nt)) {
        return(2)
    }

    return(nt)
  } else
    return(nThreads)
}

.nProcessorsOnline = function() {
  n = detectCores()
  if (!is.numeric(n)) n = 2
  if (!is.finite(n)) n = 2
  if (n<1) n = 2
  n
}


#' Allow and disable multi-threading for certain WGCNA calculations
#'
#' These functions allow and disable multi-threading for WGCNA calculations
#' that can optionally be multi-threaded, which includes all functions using
#' \code{\link{cor}} or \code{\link{bicor}} functions.
#'
#' \code{allowWGCNAThreads} enables parallel calculation within the compiled
#' code in WGCNA, principally for calculation of correlations in the presence
#' of missing data. This function is now deprecated; use
#' \code{enableWGCNAThreads} instead.
#'
#' \code{enableWGCNAThreads} enables parallel calculations within user-level R
#' functions as well as within the compiled code, and registers an appropriate
#' parallel calculation back-end for the operating system/platform.
#'
#' \code{disableWGCNAThreads} disables parallel processing.
#'
#' \code{WGCNAnThreads} returns the number of threads (parallel processes) that
#' WGCNA is currently configured to run with.
#'
#' @aliases allowWGCNAThreads enableWGCNAThreads disableWGCNAThreads
#' WGCNAnThreads
#' @param nThreads Number of threads to allow. If not given, the number of
#' processors online (as reported by system configuration) will be used. There
#' appear to be some cases where the automatically-determined number is wrong
#' please check the output to see that the number of threads makes sense.
#' Except for testing and/or torturing your system, the number of threads
#' should be no more than the number of actual processors/cores.
#' @return \code{allowWGCNAThreads}, \code{enableWGCNAThreads}, and
#' \code{disableWGCNAThreads} return the maximum number of threads WGCNA
#' calculations will be allowed to use.
#' @note Multi-threading within compiled code is not available on Windows; R
#' code parallelization works on all platforms.
#' @author Peter Langfelder
#' @keywords misc
#' @export allowWGCNAThreads
allowWGCNAThreads = function(nThreads = NULL) {
  # Stop any clusters that may be still running
  disableWGCNAThreads()
  # Enable WGCNA threads
  if (is.null(nThreads)) {
      nThreads = .nProcessorsOnline()
  }
  if (!is.numeric(nThreads) || nThreads < 2) {
    stop("nThreads must be numeric and at least 2.")
  }

  if (nThreads > .nProcessorsOnline()) {
    printFlush(paste("Warning in allowWGCNAThreads: Requested number of threads is higher than number\n",
                     "of available processors (or cores). Using too many threads may degrade code",
                     "performance. It is recommended that the number of threads is no more than number\n",
                     "of available processors.\n"))
  }

  printFlush(paste("Allowing multi-threading with up to", nThreads, "threads."))

  pars = list(nThreads)
  names(pars) = .threadAllowVar
  do.call(Sys.setenv, pars)
  invisible(nThreads)
}

#' @rdname allowWGCNAThreads
#' @aliases disableWGCNAThreads
#' @export
disableWGCNAThreads = function() {
  Sys.unsetenv(.threadAllowVar)
  pars = list(1)
  names(pars) = .threadAllowVar
  do.call(Sys.setenv, pars)
  if (exists(".revoDoParCluster", where = ".GlobalEnv")) {
    parallel::stopCluster(get(".revoDoParCluster", pos = ".GlobalEnv"))
  }
  foreach::registerDoSEQ()
}

.checkAvailableMemory = function() {
  size = 0
  res = .C("checkAvailableMemoryForR", size = as.double(size), PACKAGE = "WGCNA")
  res$size
}

#' @rdname allowWGCNAThreads
#' @aliases enableWGCNAThreads
#' @export
enableWGCNAThreads = function(nThreads = NULL) {
  nCores = detectCores()
  if (is.null(nThreads)) {
      if (nCores < 4) {
          nThreads = nCores
      } else {
          nThreads = nCores - 1
      }
  }
  if (!is.numeric(nThreads) || nThreads < 2) {
        stop("nThreads must be numeric and at least 2.")
  }
  if (nThreads > nCores) {
     printFlush(paste("Warning in allowWGCNAThreads: Requested number of threads is higher than number\n",
          "of available processors (or cores). Using too many threads may degrade code",
          "performance. It is recommended that the number of threads is no more than number\n",
          "of available processors.\n"))
  }
  printFlush(paste("Allowing parallel execution with up to", nThreads,
                   "working processes."))
  pars = list(nThreads)
  names(pars) = .threadAllowVar
  do.call(Sys.setenv, pars)

  # Register a parallel backend for foreach
  doParallel::registerDoParallel(nThreads)

  # Return the number of threads invisibly
  invisible(nThreads)
}
#' @rdname allowWGCNAThreads
#' @aliases WGCNAnThreads
#' @export
WGCNAnThreads = function() {
  n = suppressWarnings(as.numeric(as.character(
      Sys.getenv(.threadAllowVar, unset = 1))))
  if (is.na(n)) {
      n = 1
  }
  if (length(n)  ==  0) {
      n = 1
  }
  n
}

# Facilitates multi-threading by producing an even allocation of jobs
# Works even when number of jobs is less than number of threads in which case some components of the
# returned allocation will have length 0.


#' Divide tasks among workers
#'
#' This function calculates an even splitting of a given number of tasks among
#' a given number of workers (threads).
#'
#' Tasks are labeled consecutively 1,2,..., \code{nTasks}. The tasks are split
#' in contiguous blocks as evenly as possible.
#'
#' @param nTasks number of tasks to be divided
#' @param nWorkers number of workers
#' @return A list with one component per worker giving the task indices to be
#' worked on by each worker. If there are more workers than tasks, the tasks
#' for the extra workers are 0-length numeric vectors.
#' @author Peter Langfelder
#' @keywords misc
#' @examples
#'
#' allocateJobs(10, 3)
#' allocateJobs(2,4)
#'
#' @export allocateJobs
allocateJobs = function(nTasks, nWorkers) {
  if (is.na(nWorkers)) {
    warning("In function allocateJobs: 'nWorkers' is NA. Will use 1 worker.")
    nWorkers = 1
  }
  n1 = floor(nTasks/nWorkers)
  n2 = nTasks - nWorkers*n1
  allocation = list()
  start = 1
  for (t in 1:nWorkers) {
    end = start + n1 - 1 + as.numeric(t<=n2)
    if (start > end) {
      allocation[[t]] = numeric(0)
    } else {
        allocation[[t]] = c(start:end)
    }
    start = end+1
  }
  allocation
}
