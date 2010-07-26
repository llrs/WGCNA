# slight re-definition of the bicor function

.pearsonFallbacks = c("none", "individual", "all");

bicor = function(x, y = NULL, robustX = TRUE, robustY = TRUE, use = 'all.obs', maxPOutliers = 1, quick = 0,
                      pearsonFallback = "individual",
                      nThreads = 0, verbose = 0, indent = 0)
{
  Cerrors = c("Memory allocation error")
  nKnownErrors = length(Cerrors);
  na.method = pmatch(use, c("all.obs", "pairwise.complete.obs"))
  if (is.na(na.method))
      stop(paste("Unrecognized parameter 'use'. Recognized values are \n",
          "'all.obs', 'pairwise.complete.obs'"))
  if (na.method==1)
  {
    if (is.null(y))
    {
        sumNAy = 0
    } else sumNAy = sum(is.na(y));
    if (sum(is.na(x)) + sumNAy > 0)
      stop("Missing values present in input data. Consider using use = 'pairwise.complete.obs'.");
  }

  fallback = pmatch(pearsonFallback, .pearsonFallbacks)
  if (is.na(na.method))
      stop(paste("Unrecognized 'pearsonFallback'. Recognized values are (unique abbreviations of)\n",
           paste(.pearsonFallbacks, collapse = ", ")))

  if (quick < 0) stop("quick must be non-negative.");
  if (nThreads < 0) stop("nThreads must be non-negative.");

  x = as.matrix(x);
  nNA = 0;
  err = 0;
  if (is.null(y))
  {
    if (!robustX)
    {
      bi = cor(x, use = use)
    } else {
      bi = matrix(0, ncol(x), ncol(x));
      res = .C("bicor1Fast", x = as.double(x), nrow = as.integer(nrow(x)), ncol = as.integer(ncol(x)),
               maxPOutliers = as.double(maxPOutliers), 
               quick = as.double(quick), 
               fallback = as.integer(fallback),
               res = as.double(bi), nNA = as.integer(nNA),
               err = as.integer(err), nThreads = as.integer(nThreads),
               verbose = as.integer(verbose), indent = as.integer(indent),
               DUP = FALSE, NAOK = TRUE);
    }
    dim(res$res) = dim(bi);
    if (!is.null(dimnames(x)[[2]])) dimnames(res$res) = list(dimnames(x)[[2]],  dimnames(x)[[2]] );
  } else {
    y = as.matrix(y);
    if (nrow(x)!=nrow(y))
      stop("'x' and 'y' have incompatible dimensions (unequal numbers of rows).");
    bi = matrix(0, ncol(x), ncol(y));
    res = .C("bicorFast", x = as.double(x), nrow = as.integer(nrow(x)), ncolx = as.integer(ncol(x)),
             y = as.double(y), ncoly = as.integer(ncol(y)),
             robustX = as.integer(robustX), robustY = as.integer(robustY),
             maxPOutliers = as.double(maxPOutliers), 
             quick = as.double(quick), 
             fallback = as.integer(fallback),
             res = as.double(bi), nNA = as.integer(nNA), err = as.integer(err),
             nThreads = as.integer(nThreads),
             verbose = as.integer(verbose), indent = as.integer(indent), DUP = FALSE, NAOK = TRUE);
    dim(res$res) = dim(bi);
    if (!is.null(dimnames(x)[[2]]) || !is.null(dimnames(y)[[2]]))
        dimnames(res$res) = list(dimnames(x)[[2]], dimnames(y)[[2]]);
  }
  if (res$err > 0)
  {
    if (err > nKnownErrors)
    {
      stop(paste("An error occurred in compiled code. Error code is", err));
    } else {
      stop(paste(Cerrors[err], "occurred in compiled code. "));
    }
  }
  if (res$nNA > 0)
  {
    warning(paste("Missing values generated in calculation of bicor.",
                  "Likely cause: too many missing entries, zero median absolute deviation, or zero variance."));
  }
  res$res;
}

# Code to call my implementation of correlation

cor = function(x, y = NULL, use = "all.obs", method = c("pearson", "kendall", "spearman"),
                quick = 0, nThreads = 0, verbose = 0, indent = 0)
{
    na.method <- pmatch(use, c("all.obs", "complete.obs", "pairwise.complete.obs", 
        "everything", "na.or.complete"), nomatch = 0)
    method <- match.arg(method)

    if ((method=="pearson") && ( (na.method==1) || (na.method==3) ) )
    {
      Cerrors = c("Memory allocation error")
      nKnownErrors = length(Cerrors);
      na.method = pmatch(use, c("all.obs", "pairwise.complete.obs"))
      if (is.na(na.method))
          stop(paste("Unrecognized parameter 'use'. Recognized values are \n",
              "'all.obs', 'pairwise.complete.obs'"))
      if (na.method==1)
      {
        if (sum(is.na(x)) > 0)
          stop("Missing values present in input data. Consider using use = 'pairwise.complete.obs'.");
      }
      
      if (quick < 0) stop("quick must be non-negative.");
      if (nThreads < 0) stop("nThreads must be non-negative.");
    
      x = as.matrix(x);
      nNA = 0;
      err = 0;
      if (is.null(y))
      {
         bi = matrix(0, ncol(x), ncol(x));
         res = .C("cor1Fast", x = as.double(x), nrow = as.integer(nrow(x)), ncol = as.integer(ncol(x)),
                   quick = as.double(quick), res = as.double(bi), nNA = as.integer(nNA),
                   err = as.integer(err), nThreads = as.integer(nThreads), 
                   verbose = as.integer(verbose), indent = as.integer(indent),
                   DUP = FALSE, NAOK = TRUE);
         dim(res$res) = dim(bi);
         if (!is.null(dimnames(x)[[2]])) dimnames(res$res) = list(dimnames(x)[[2]],  dimnames(x)[[2]] );
      } else {
         y = as.matrix(y);
         if (nrow(x)!=nrow(y))
            stop("'x' and 'y' have incompatible dimensions (unequal numbers of rows).");
         bi = matrix(0, ncol(x), ncol(y));
         res = .C("corFast", x = as.double(x), nrow = as.integer(nrow(x)), ncolx = as.integer(ncol(x)),
                 y = as.double(y), ncoly = as.integer(ncol(y)),
                 quick = as.double(quick), res = as.double(bi), nNA = as.integer(nNA), err = as.integer(err),
                 nThreads = as.integer(nThreads),
                 verbose = as.integer(verbose), indent = as.integer(indent), DUP = FALSE, NAOK = TRUE);
         dim(res$res) = dim(bi);
         if (!is.null(dimnames(x)[[2]]) || !is.null(dimnames(y)[[2]]))
            dimnames(res$res) = list(dimnames(x)[[2]], dimnames(y)[[2]]);
      }
      if (res$err > 0)
      {
        if (err > nKnownErrors)
        {
          stop(paste("An error occurred in compiled code. Error code is", err));
        } else {
          stop(paste(Cerrors[err], "occurred in compiled code. "));
        }
      }
      if (res$nNA > 0)
      {
        warning(paste("Missing values generated in calculation of bicor.",
                      "Likely cause: too many missing entries or zero variance."));
      }
      res$res;
    } else {
      stats::cor(x,y, use, method);
    }
}

cor1 = function(x, use = "all.obs", verbose = 0, indent = 0) 
{
   cor(x, use = use, verbose = verbose, indent = indent);
}

# Wrapper for compatibility with older scripts

corFast = function(x, y = NULL, use = "all.obs",
                quick = 0, nThreads = 0, verbose = 0, indent = 0)
{
  cor(x,y, use, method = "pearson", quick, nThreads, verbose, indent)
}


      