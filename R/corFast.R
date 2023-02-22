#' Title     : corFast
#' Objective : Import the implementation of matrix-correlation calculation from the WGCNA package
#' Created by: Nicholas Schmitt
#' Created on: 04.11.22

#' The following code is taken wholesale from the WGCNA package.
#' The original paper can be found here:
#' Langfelder P, Horvath S (2008). “WGCNA: an R package for weighted correlation network analysis.”
#' BMC Bioinformatics, 559. https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-9-559.

#' The reason why we aren't including the WGCNA package as a dependency is because the package relies on a few libraries
#' which aren't available for R >= 4.0.0 as of writing this code.

#'  Copyright (C) 2008 Peter Langfelder; parts based on R by R Development team
#'  This program is free software; you can redistribute it and/or
#'  modify it under the terms of the GNU General Public License
#'  as published by the Free Software Foundation; either version 2
#'  of the License, or (at your option) any later version.
#'  This program is distributed in the hope that it will be useful,
#'  but WITHOUT ANY WARRANTY; without even the implied warranty of
#'  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#'  GNU General Public License for more details.
#'  You should have received a copy of the GNU General Public License
#'  along with this program; if not, write to the Free Software
#'  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

dyn.load("src/corFunctions.so")

print("load successful")

corFast <- function(x, y = NULL, use = "all.obs", method = "pearson",
               weights.x = NULL, weights.y = NULL,
               quick = 0,
               cosine = FALSE,
               cosineX = cosine,
               cosineY = cosine,
               drop = FALSE,
               nThreads = 0, verbose = 0, indent = 0)
{
    na.method <- pmatch(use, c("all.obs", "complete.obs", "pairwise.complete.obs",
        "everything", "na.or.complete"), nomatch = 0)
    method <- match.arg(method)

    if (length(weights.x)==0) weights.x <- NULL
    if (length(weights.y)==0) weights.y <- NULL

    x <- as.matrix(x)
    if (!is.null(y))
    {
      y <- as.matrix(y)
    }

    if ((method=="pearson") && ( (na.method==1) || (na.method==3) ))
    {
      Cerrors <- "Memory allocation error"
      nKnownErrors <- length(Cerrors)
      na.method <- pmatch(use, c("all.obs", "pairwise.complete.obs"))
      if (is.na(na.method))
          stop(paste("Unrecognized parameter 'use'. Recognized values are \n",
              "'all.obs', 'pairwise.complete.obs'"))
      if (na.method==1)
      {
         if (sum(is.na(x))> 0)
           stop("Missing values present in input variable 'x'. Consider using use = 'pairwise.complete.obs'.")
         if (!is.null(y))
         {
           if (sum(is.na(y)) > 0)
             stop("Missing values present in input variable 'y'. Consider using use = 'pairwise.complete.obs'.")
         }
      }

      if (quick < 0) stop("quick must be non-negative.")
      if (nThreads < 0) stop("nThreads must be non-negative.")
      if (is.null(nThreads) || (nThreads==0)) nThreads <- .useNThreads()

      if (prod(dim(x))==0) stop("'x' has a zero dimension.")

      if (!is.null(weights.x))
      {
        if (is.null(dim(weights.x)))
        {
          if (length(weights.x)!=nrow(x))
            stop("When 'weights.x' are given, they must be a vector of length 'nrow(x)' or a matrix\n",
                 "of the same dimensions as 'x'.")
          weights.x <- matrix(weights.x, nrow(x), ncol(x))
        } else
          if (!isTRUE(all.equal(dim(weights.x), dim(x))))
             stop("When 'weights.x' are given, they must be a vector of length 'nrow(x)' or a matrix\n",
                 "of the same dimensions as 'x'.")
        if (any(!is.finite(weights.x)))
        {
          if (verbose > 0)
            warning("cor: found non-finite weights. These will be removed (set to missing), ",
                    "and the corresponding entries in 'x' will be treated as missing.")
          weights.x[!is.finite(weights.x)] <- NA
        }
        if (any(weights.x < 0, na.rm = TRUE))
          stop("All weights must be non-negative.")

        if (!is.null(y) && is.null(weights.y)) weights.y <- matrix(1, nrow(y), ncol(y))
      }

      if (!is.null(weights.y))
      {
        if (is.null(y)) stop("'weights.y' can only be used if 'y' is non-NULL.")
        if (is.null(dim(weights.y)))
        {
          if (length(weights.y)!=nrow(y))
            stop("When 'weights.y' are given, they must be a vector of length 'nrow(y)' or a matrix\n",
                 "of the same dimensions as 'y'.")
          weights.y <- matrix(weights.y, nrow(y), ncol(y))
        } else
          if (!isTRUE(all.equal(dim(weights.y), dim(y))))
             stop("When 'weights.y' are given, they must be a vector of length 'nrow(y)' or a matrix\n",
                 "of the same dimensions as 'y'.")
        if (any(!is.finite(weights.y)))
        {
          if (verbose > 0)
            warning("cor: found non-finite weights. These will be removed (set to missing), ",
                    "and the corresponding entries in 'x' will be treated as missing.")
          weights.y[!is.finite(weights.y)] <- NA
        }
        if (any(weights.y < 0, na.rm = TRUE))
          stop("All weights must be non-negative.")

        if (is.null(weights.x)) weights.x <- matrix(1, nrow(x), ncol(x))
      }

      storage.mode(x) <- "double"
      if (!is.null(weights.x)) storage.mode(weights.x) <- "double"
      if (!is.null(weights.y)) storage.mode(weights.y) <- "double"
      nNA <- 0L
      err <- as.integer(nNA-1 + 1/1)
      cosineX <- as.integer(cosineX)
      nThreads <- as.integer(nThreads)
      verbose <- as.integer(verbose)
      indent <- as.integer(indent)
      if (is.null(y))
      {
         res <- .Call("cor1Fast_call", x, weights.x,
                   quick, cosine,
                   nNA, err, nThreads,
                   verbose, indent)
         if (!is.null(dimnames(x)[[2]])) dimnames(res) <- list(dimnames(x)[[2]],  dimnames(x)[[2]] )
      } else {
         y <- as.matrix(y)
         storage.mode(y)<- "double"
         cosineY <- as.integer(cosineY)
         if (prod(dim(y))==0) stop("'y' has a zero dimension.")
         if (nrow(x)!=nrow(y))
            stop("'x' and 'y' have incompatible dimensions (unequal numbers of rows).")
         res <- .C("cor1Fast_call", x, y,
                 weights.x, weights.y,
                 quick,
                 cosineX,
                 cosineY,
                 nNA, err,
                 nThreads,
                 verbose, indent,
                 PACKAGE = "WGCNA")
         if (!is.null(dimnames(x)[[2]]) || !is.null(dimnames(y)[[2]]))
            dimnames(res) <- list(dimnames(x)[[2]], dimnames(y)[[2]])
      }
      if (err > 0)
      {
        if (err > nKnownErrors)
        {
          stop(paste("An error occurred in compiled code. Error code is", err))
        } else {
          stop(paste(Cerrors[err], "occurred in compiled code. "))
        }
      }
      if (nNA > 0)
      {
        warning(paste("Missing values generated in calculation of cor.",
                      "Likely cause: too many missing entries or zero variance."))
      }
      if (drop) res[, , drop = TRUE] else res
    } else {
      stats::cor(x,y, use, method)
    }
}

.useNThreads <- function(nThreads = 0)
{
  if (nThreads==0)
  {
    nt.env <- Sys.getenv("ALLOW_WGCNA_THREADS", unset = NA)
    if (is.na(nt.env)) return(1)
    if (nt.env=="") return(1)

    if (nt.env=="ALL_PROCESSORS") return (.nProcessorsOnline())

    nt <- suppressWarnings(as.numeric(nt.env))

    if (!is.finite(nt)) return(2)

    return(nt)
  } else
    return (nThreads)
}

.nProcessorsOnline <- function()
{
  n <- detectCores()
  if (!is.numeric(n)) n <- 2
  if (!is.finite(n)) n <- 2
  if (n<1) n <- 2
  n
}
