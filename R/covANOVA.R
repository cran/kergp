setClass("covANOVA",   	
         representation(
             k1Fun1 = "function",
             k1Fun1Char = "character",
             hasGrad = "logical",
             cov = "integer",               ## 0 : corr, 1 : homo
             iso = "integer",               ## iso = 0
             iso1 = "integer",
             label = "character",
             d = "integer",                 ## (spatial) dimension
             inputNames = "optCharacter",   ## spatial var names length d
             parLower = "numeric",          ## lower bound on pars
             parUpper = "numeric",          ## upper bound on pars
             par  = "numeric",              ## params values
             kern1ParN1 = "integer",        ## nbr of 'shape' pars in
                                            ## 'k1Fun1', usually 0
             parN1 = "integer",             ## nbr of 'shape' pars
             parN = "integer",              ## number of pars
             kern1ParNames = "character", 
             kernParNames  = "character"    ## depending on kernel
         ),
         validity = function(object){
             if (length(object@kernParNames) != object@parN) {
                 stop("Incorrect number of parameter names")
             }
           
         }
         )

## *****************************************************************************
## DO NOT Roxygenize THIS!!!
##
##' Creator for the class \code{"covANOVA"}.
##'
##' An ANOVA-product kernel on the \eqn{d}-dimensional Euclidean space
##' takes the form \deqn{K(\mathbf{x},\,\mathbf{x}') = \delta^2
##' \prod_{\elll = 1}^d [ 1 + \tau^2_\ell \kappa(r_\ell)]}{K(x1, x2) = \delta^2 * prod_l
##' [1 + \tau^2[l] k1(r[l])]} where \eqn{\kappa(r)}{k1(r)} is a suitable correlation
##' kernel for a one-dimensional input, and \eqn{r_\ell}{r[l]} is
##' given by \eqn{r_\ell := [x_\ell - x'_\ell] / \theta_\ell}{r[l] = (x2[l] -
##' x2[l]) / \theta[l]} for \eqn{\ell = 1}{l =1} to \eqn{d}.
##' 
##' In this default form, the ANOVA-product kernel depends on \eqn{d +
##' d + 1} parameters: the \emph{ranges} \eqn{\theta_\ell
##' >0}{\theta[l] > 0}, the \emph{weights}
##' \eqn{\tau^2_\ell}{\tau^2[l]}, and the \emph{variance}
##' \eqn{\delta^2}{\delta^2}.
##'
##' An \emph{isotropic} form uses the same range \eqn{\theta}{\theta}
##' for all inputs, i.e. sets \eqn{\theta_\ell =
##' \theta}{\theta[l]=\theta} for all \eqn{\ell}{l}. This is obtained by
##' using \code{iso = TRUE}.
##'
##' A \emph{correlation} version uses \eqn{\delta^2 = 1}{\delta^2 =
##' 1}. This is obtained by using \code{cov = "corr"}.
##'
##' Finally, the correlation kernel \eqn{\kappa(r)}{k1(r)} can depend on
##' a "shape" parameter, e.g. have the form
##' \eqn{\kappa(r;\,\alpha)}{k1(r, \alpha)}. The extra shape parameter
##' \eqn{\alpha}{\alpha} will be considered then as a parameter of the
##' resulting tensor-product kernel, making it possible to estimate it
##' by ML along with the range(s) and the variance.
##' 
##' @title Creator for the Class \code{"covANOVA"}
##' 
##' @param k1Fun1 A kernel function of a \emph{scalar} numeric
##' variable, and possibly of an extra "shape" parameter. This
##' function should return the first-order derivative or the two-first
##' order derivatives as an attibute with name \code{"der"} and with a
##' matrix content. When an extra shape parameter exists, the gradient
##' should also be returned as an attribute with name
##' \code{"gradient"}, see \bold{Examples} later. The name of the
##' function can be given as a character string.
##'
##' @param cov A character string specifying the kind of covariance
##' kernel: correlation kernel ("corr") or kernel of a homoscedastic
##' GP ("homo"). Partial matching is allowed.
##' 
##' @param iso Integer. The value \code{1L} coresponds to an isotropic
##' covariance, with all the inputs sharing the same range value.
##'
##' @param iso1 Integer. This applies only when \code{k1Fun1} contains
##' one or more parameters that can be called 'shape' parameters. For
##' now only one such parameter can be found in \code{k1Fun1} and
##' consequently \code{iso1} must be of length one. With \code{iso1 = 0}
##' the shape parameter in \code{k1Fun1} will generate \code{d}
##' parameters in the \code{covANOVA} object with their name suffixed by
##' the dimension. When \code{iso1} is \code{1} only one shape
##' parameter will be created in the \code{covANOVA} object.
##' 
##' @param hasGrad Integer or logical. Tells if the value returned by
##' the function \code{k1Fun1} has an attribute named \code{"der"}
##' giving the derivative(s).
##' 
##' @param inputs Character. Names of the inputs.
##'
##' @param d Integer. Number of inputs. 
##' 
##' @param parNames Names of the parameters. By default, ranges are
##' prefixed "theta_" in the non-iso case and the range is named
##' "theta" in the iso case.
##' 
##' @param par Numeric values for the parameters. Can be \code{NA}.
##'
##' @param parLower Numeric values for the lower bounds on the
##' parameters. Can be \code{-Inf}.
##'
##' @param parUpper Numeric values for the uppper bounds on the
##' parameters. Can be \code{Inf}.
##'
##' @param label A short description of the kernel object.
##' 
##' @param ... Other arguments passed to the method \code{new}.
##'
##' @param parNames. Names of the parameters. By default, ranges
##' are prefixed \code{"theta_"} in the non-iso case and the range is
##' names \code{"theta"} in the iso case.
##' 
##' @return An object with class \code{"covANOVA"}.
##'
##' @examples
##'
##' 
covANOVA <- function(k1Fun1 = k1Fun1Gauss,
                  cov = c("corr", "homo"),
                  iso = 0,
                  iso1 = 1L,
                  hasGrad = TRUE,
                  inputs = NULL,
                  d = NULL,
                  parNames,
                  par = NULL,
                  parLower = NULL,
                  parUpper = NULL,
                  label = "ANOVA kernel",
                  ...) {

    L <- list(...)
    
    cov <- match.arg(cov)
    cov <- switch(cov, corr = 0, homo = 1)
    
    if (is.null(d)) {
        if (is.null(inputs)) {
            stop ("at least one of the formals 'd' and 'inputs' must be given")
        } else {
            d <- length(inputs)
        }
    } else {
        if (is.null(inputs)) {
            inputs <-  paste("x", 1:d, sep = "")
        } else {
            if (length(inputs) != d) {
                stop("'inputs' must be of length ", d)
            }
        }
    }

    k1Fun1Char <- deparse(substitute(k1Fun1))
    k1Fun1 <- match.fun(k1Fun1)

    ## =========================================================================
    ## If 'k1Fun1' has more than one formal argument, the arguments
    ## having position > 1 become parameters for the kernel, with
    ## their name unchanged ! Caution: possible conflicts here for
    ## some names like "theta" or "delta2".
    ## =========================================================================
    f <- as.list(formals(k1Fun1))
    f[[1]] <- NULL
    kern1ParN1 <- length(f)
    
    if (kern1ParN1) {
        
        if (kern1ParN1 > 1) {
            stop("For now, only one \"shape\" parameter is allowed!")
        }
        
        kern1ParNames <- names(f)
        
        if (length(grep("theta", kern1ParNames)) ||
            length(grep("delta", kern1ParNames))) {
            stop("The names of the formals arguments of 'k1Fun1' ",
                 "are in conflict with default names  \"theta\" and ",
                 "\"delta\". Please rename.")
        }
        
        ## Set the shape parameters to their default value
        par1 <- sapply(f, function(x) if (x == "") return(NA) else return(x))
                       
        if (!iso1) {
            parN1 <- d
            par1 <- rep(par1, length.out = d)
            kern1ParNamesANOVA <- paste(kern1ParNames, 1:d, sep = "_")
        } else {
            parN1 <- kern1ParN1
            kern1ParNamesANOVA <- kern1ParNames
        }
        
    } else {
        kern1ParNames <- kern1ParNamesANOVA <- character(0)
        parN1 <- 0L
    }
    
    parN <- parN1
    
    if (iso) {
        kernParNames <- c(kern1ParNamesANOVA, "theta")
        parN <- parN + 1L
    } else {
        kernParNames <- c(kern1ParNamesANOVA, paste0("theta_", 1:d))
        parN <- parN + d
    }
    
    kernParNames <- c(kernParNames, paste0("tau2_", 1:d))
    parN <- parN + d
    
    if (cov) {
        kernParNames <- c(kernParNames, "delta2")
        parN <- parN + 1L
    }
    
    if (is.null(par)) {
        par <- as.numeric(rep(1.0, parN))
        if (parN1) par[1:parN1] <- par1 
    }
    
    if (is.null(parLower)) {
        parLower <- as.numeric(rep(0, parN))
        if (parN1) parLower[1:parN1] <- -Inf 
    }
    if (is.null(parUpper)) {
        parUpper <- as.numeric(rep(Inf, parN))
    }
    
    if (missing(d) & missing(inputs)) {
        stop("at least one of 'd' or 'inputs' must be provided")
    }
    if (length(inputs) != d) {
        stop("'d' must be equal to 'length(inputs)'")
    }

    ## XXX check that the function is zero. The extra
    ## arguments should have a default value
    if (k1Fun1(0) != 1.0) {
        stop("the function given in 'k1Fun1' must be such that ",
             "k1Fun1(0) = 1.0")
    }
    
    ## XXX check that the derivative is given 
    if (hasGrad) {
        
    }
    
    new("covANOVA", 
        k1Fun1 = k1Fun1,
        k1Fun1Char = k1Fun1Char,
        hasGrad = as.logical(hasGrad),
        cov = as.integer(cov),
        iso = as.integer(iso),
        iso1 = as.integer(iso1),
        label = as.character(label),
        d = as.integer(d),
        inputNames = as.character(inputs),
        par = as.numeric(par),
        parLower = as.numeric(parLower),
        parUpper = as.numeric(parUpper),
        kern1ParN1 = as.integer(kern1ParN1),
        parN1 = as.integer(parN1),
        parN = as.integer(parN),
        kern1ParNames = as.character(kern1ParNames),
        kernParNames = as.character(kernParNames),
        ...)
    
    } 


##*****************************************************************************
## 'covMat' method
##
## When 'Xnew' is not given, i.e. in the symmetric case, efforts are
## made to exploit the symmetry. In particular, the univariate kernel
## function (which can be very costly in some cases) is invoked only
## (n - 1) * n / 2 times.
##
## See the document 'kergp Computing Details' to understand the
## computation of the gradient or of the derivatives.
## 
##******************************************************************************
setMethod("covMat",
          signature = "covANOVA", 
          definition = function(object, X, Xnew,
              compGrad = hasGrad(object),
              deriv = 0,
              checkNames = NULL, ...) {
              
              isXNew <- !is.null(Xnew)
              if (isXNew) compGrad <- FALSE
              
              ## ==============================================================
              ## Specific checks related to 'deriv'.
              ## ==============================================================
              
              if (deriv > 0) {
                  if (compGrad) {
                      stop("For now, 'deriv' can be non-zero only when ",
                           "'compGrad' is FALSE")
                  }
                  if (!isXNew) {
                      stop("When 'Xnew' is not given, 'deriv' can only be 0") 
                  }
                  if (deriv > 1) {
                      stop("Second-order derivatives are not yet available")
                  }
              }
              
              ## ==============================================================
              ## The kernel has only continuous (numeric) inputs...
              ## ==============================================================
              
              if (is.null(checkNames)) {
                  checkNames <- TRUE 
                  if (object@d == 1L) {
                      if (ncol(X) == 1L) checkNames <- FALSE
                  }      
              }
              
              ## ==============================================================
              ## check names and content of 'X', and of 'Xnew' if
              ## provided.
              ## ==============================================================
      
              if (checkNames) X <- checkX(object, X = X)
              if (any(is.na(X))) stop("'X' must not contain NA elements")
              X <- as.matrix(X)
              
              if (isXNew){
                  if (checkNames) Xnew <- checkX(object, X = Xnew)
                  Xnew <- as.matrix(Xnew)
                  if (ncol(X) != ncol(Xnew)) {
                      stop("'X' and 'Xnew' must have the same number of columns")
                  }
                  if (any(is.na(Xnew))) {
                      stop("'Xnew' must not contain NA elements")
                  }
                  nNew <- nrow(Xnew)
              } 

              ## ===============================================================
              ## Some parameters. The vector 'theta' of ranges is
              ## recycled to have length 'd' in the 'iso' case
              ## ===============================================================
              
              d <- object@d
              n <- nrow(X)
              compGrad <- as.integer(compGrad)
              kern1ParN1 <- object@kern1ParN1
              parN1 <- object@parN1
              parN <- object@parN
              par <- coef(object)
              nTheta <- ifelse(object@iso, 1L, d)
              theta <- par[parN1 + (1L:nTheta)]
              theta <- rep(theta, length.out = d)
              tau2 <-  par[parN1 + nTheta + (1:d)]
              if (parN1) {
                  psi <- par[1:parN1]
                  psi <- rep(psi, length.out = d)
              }
              
              if (object@cov) delta2 <- par[parN1 + nTheta + d + 1L]
              else delta2 <- 1.0
              
              if (compGrad) {
                  
                  if (!object@hasGrad) {
                      stop("'object' does not allow the computation of the",
                           " gradient")
                  }
                  if (isXNew) {
                      stop("the computation of the gradient is not allowed",
                           " for now when 'Xnew' is given (only the symmetric",
                           " case is allowed)")
                  }
                  
                  sI <- symIndices(n)
                  
                  ## ==========================================================
                  ## We need to store the quantities 'R' used in the
                  ## gradient. We begin by working with (helf-)
                  ## vectorized versions of the symmetric matrices
                  ## ==========================================================
                  
                  m <- (n - 1L) * n / 2
                  Rs <- Cs <- Kappas <- dKappas <- array(NA, dim = c(m, d))
                  
                  ## gradients (w.r.t) parameters. We must add a dimension
                  ## if kern1ParN1 > 1
                  if (kern1ParN1 > 0L) gKappas <- array(NA, dim = c(m, d))
                  
                  for (ell in 1:d) {                      
                      
                      Rs[ , ell] <- abs(X[sI$i, ell] - X[sI$j, ell]) / theta[ell]
                      
                      if (parN1) {
                          Kappaell <- do.call(object@k1Fun1,
                                              list(Rs[ , ell], psi[ell]))
                      } else {
                          Kappaell <- do.call(object@k1Fun1,
                                              list(Rs[ , ell]))
                      }
                      
                      Kappas[ , ell] <- Kappaell
                      dKappas[ , ell] <- attr(Kappaell, "der")[ , 1L]
                      Cs[ , ell] <- 1.0 + tau2[ell] * Kappaell
                      
                      if (kern1ParN1 > 0L) {
                          gKappas[ , ell]  <- attr(Kappaell, "gradient")[ , 1L]
                      }
                      
                      if (ell == 1) {
                          C <- Cs[ , ell]
                      } else {
                          C <- C * Cs[ , ell]
                      }
                  }
                  
                  ## ==========================================================
                  ## We now make
                  ##
                  ## o 'KappaSym' to be a n * n symmetric matrix and
                  ## 
                  ## o 'GradSym' an array with dim c(n, n, parN) and
                  ## with symmetric n * n slices GradSym[ , , k].
                  ##
                  ## Note that the temporary array 'Grad' is needed
                  ## because it is not possible to use a two-indice
                  ## style of indexation for a 3-dimensional
                  ## array. 'Grad' keeps a zero diagonal all along the
                  ## loop.
                  ## ==========================================================
                  
                  prodTau <- prod(1 + tau2)
                  sigma2 <- delta2 * prodTau

                  CSym <- array(sigma2, dim = c(n, n))
                  CSym[sI$kL] <- CSym[sI$kU] <- C
                  diag(CSym) <- prodTau

                  
                  GradSym <- array(0.0, dim = c(n, n, parN))
                  dimnames(GradSym) <- list(NULL, NULL, object@kernParNames)
                  Grad <- array(0.0, dim = c(n, n))
                  
                  ## ==========================================================
                  ## 'shape' parameters of the kernel For now, 'kernParN1'
                  ## can only be 0 or 1...
                  ## ==========================================================
                  
                  if (kern1ParN1 > 0L) {
                      
                    if (kern1ParN1 > 1L) {
                      stop("Using a kernel with more than one shape parameter is ",
                          " not yet possible when 'compGrad' is TRUE")
                    }
                      
                      if (!object@iso1) {                          
                          
                          for (ell in 1L:d) {
                              
                              gKappas[ , ell] <- delta2 * tau2[ell] * gKappas[ , ell] /
                                  Cs[ , ell] * C
                              
                              Grad[sI$kL] <- Grad[sI$kU] <- gKappas[ , ell]
                              GradSym[ , , ell] <- array(Grad, dim = c(n, n))
                              
                          }
                          
                      } else {                          
                          
                          for (ell in 1L:d) {
                              
                              gKappas[ , ell] <- delta2 * tau2[ell] * gKappas[ , ell] /
                                  Cs[ , ell] * C
                              
                              if (ell == 1) gKappa <- gKappas[ , ell]
                              else gKappa <- gKappa + gKappas[ , ell]
                              
                          }
                          
                          Grad[sI$kL] <- Grad[sI$kU] <- gKappa
                          GradSym[ , , 1L] <- array(Grad, dim = c(n, n))
                          
                      }
                      
                  }
                  
                  ## ==========================================================
                  ## range parameter(s)
                  ## ==========================================================
                  
                  if (!object@iso) {
                      
                      for (ell in 1L:d) {
                          
                          dKappas[ , ell] <- - delta2 * tau2[ell] * Rs[ , ell] *
                              dKappas[ , ell] / Cs[ , ell] *
                                  C / theta[ell]
                              
                          Grad[sI$kL] <- Grad[sI$kU] <- dKappas[ , ell]
                          GradSym[ , , parN1 + ell] <- array(Grad, dim = c(n, n))
                          
                      }
                      
                  } else {
                      
                      dKappas  <- - delta2 * Rs * dKappas  / Cs * C /
                          theta[1L]
                      
                      dKappas <- sweep(dKappas, MARGIN = 2, FUN = "*",
                                       STATS = tau2)
                      
                      Grad[sI$kL] <- Grad[sI$kU] <- apply(dKappas, 1L, sum)
                      GradSym[ , , parN1 + 1L] <- array(Grad, dim = c(n, n))
                  }
                  
                  ## ==========================================================
                  ## Weights parameters tau2
                  ## ==========================================================
                  
                  for (ell in 1L:d) {

                      ## cat(sprintf("tau2[%d] = %6.4f\n", ell, tau2[ell]))
                      
                      dKappas[ , ell] <- delta2 * Kappas[ , ell] / Cs[ , ell] * C 
                      
                      Grad[sI$kL] <- Grad[sI$kU] <- dKappas[ , ell]
                      diag(Grad) <- sigma2 / (1.0 + tau2[ell])
                      GradSym[ , , parN1 + nTheta + ell] <- array(Grad, dim = c(n, n))
                      
                  }
                  
                  ## delta2 parameter Now the diagonal of the slice is 1.0
                  if (object@cov) {
                      GradSym[ , , parN] <- CSym
                      CSym <- delta2 * CSym
                  }
                  
                  attr(CSym, "gradient") <- GradSym
                  return(CSym)
                  
              } else {                  
                  
                  ## ==========================================================
                  ## General (unsymmetric) case.
                  ## ==========================================================
                  
                  if (isXNew) {

                      if (deriv == 0) {

                          ## ==================================================
                          ## No need to store anything
                          ## ==================================================
                      
                          for (ell in 1:d) {
                              
                              Rell <- outer(X = X[ , ell], Y = Xnew[ , ell], "-") /
                                  theta[ell]
                                       
                              if (parN1) {
                                  Kappaell <- do.call(object@k1Fun1,
                                                   list(Rell, psi[ell]))
                              } else {
                                  Kappaell <- do.call(object@k1Fun1, list(Rell))
                              }
                              
                              attr(Kappaell, "der") <- NULL
                              
                              if (ell == 1) C <- 1 + tau2[1L] * Kappaell
                              else C <- C * (1.0 + tau2[ell] * Kappaell)
                          }    
                          
                          if (object@cov) C <- delta2 * C

                          dim(C) <- c(n, nNew)
                          return(C)

                          
                      } else {

                          ## ==================================================
                          ## We use a vectorized form of the matrices
                          ## 'K_ell' and their derivatives 'dK_ell /
                          ## dx_ell'
                          ## ==================================================
                          Rs <- Kappas <- dKappas <- Cs <- 
                              array(NA, dim = c(n * nNew, d))

                          ## XXX why not try an 'apply' here???
                          
                          for (ell in 1:d) {
                              
                              Rs[ , ell] <-
                                  outer(X[ , ell], Xnew[ , ell], "-") /
                                      theta[ell]
                              
                              if (parN1) {
                                  Kappaell <-
                                      do.call(object@k1Fun1,
                                              list(Rs[ , ell], psi[ell]))
                              } else {
                                  Kappaell <- do.call(object@k1Fun1,
                                                      list(Rs[ , ell]))
                              }
                              
                              Kappas[ , ell] <- Kappaell
                              dKappas[ , ell] <-
                                  array(attr(Kappaell, "der")[ , 1L],
                                        dim = c(n, nNew))
                              Cs[ , ell] <- 1.0 + tau2[ell] * Kappaell
                              if (ell == 1) C <- Cs[ , 1L]
                              else C <- C * Cs[ , ell]
                          }
                  
                          C <- array(C, dim = c(n , nNew))
                          
                          ## ==================================================
                          ## Derivatives
                          ## ==================================================

                          der <- array(0.0, dim = c(n, nNew, d))
                          dimnames(der) <- list(NULL, NULL, colnames(X))
                          
                          for (ell in 1L:d) {
                              
                              dKappas[ , ell] <- - delta2 * tau2[ell] *
                                  dKappas[ , ell] / Cs[ , ell] * C  / theta[ell]
                              
                              der[ , , ell] <- dKappas[ , ell] 
                              
                          }
                          
                          if (object@cov) {
                              C <- delta2 * C
                          }
                          
                          ## note that 'kappa' already has a 'der' attribute! 
                          attr(C, "der") <- der
                          return(C)
                          
                      }
                          
                  } else {

                      sI <- symIndices(n)
                      
                      ## ======================================================
                      ## We need to store the quantities 'R' used in the
                      ## gradient. We begin by working with (helf-)
                      ## vectorized versions of the symmetric matrices
                      ## ======================================================
                  
                      m <- (n - 1L) * n / 2
                      
                      for (ell in 1:d) {                      
                          
                          Rell <- (X[sI$i, ell] - X[sI$j, ell]) / theta[ell]
                          
                          if (parN1) {
                              Kappaell <- do.call(object@k1Fun1,
                                                  list(Rell, psi[ell]))
                          } else {
                              Kappaell <- do.call(object@k1Fun1, list(Rell))
                          }
                          
                          if (ell == 1) {
                              C  <- 1.0 + tau2[ell] * Kappaell
                          } else {
                              C <- C * (1.0 + tau2[ell] * Kappaell)
                          }
                      }

                      ## =======================================================
                      ## Reshape to a symmetric matrix. We assume that
                      ## the kernel function is a correlation function
                      ## so the diagonal of the returned matrix
                      ## contains ones.
                      ## =======================================================
                      
                      if (object@cov) C <- delta2 * C
                      sigma2 <- delta2 * prod(1 + tau2)
                      CSym <- array(sigma2, dim = c(n, n))
                      CSym[sI$kL] <- CSym[sI$kU] <- C                   

                      return(CSym)
                      
                  }

              }
              
          })



## *************************************************************************
## 'varVec' method: compute the variance vector.
## *************************************************************************
setMethod("varVec",
          signature = "covANOVA", 
          definition = function(object, X, compGrad = FALSE, 
                                checkNames = NULL, index = -1L, ...) {
              
              if (is.null(checkNames)) {
                  checkNames <- TRUE 
                  if (object@d == 1L) {
                      if (ncol(X) == 1L) checkNames <- FALSE
                  }      
              }
              
              if (checkNames) X <- checkX(object, X = X)
              if (any(is.na(X))) stop("'X' must not contain NA elements")
              
              if (compGrad) {
                  warning("'compGrad' is not yet implemented. Ignored.")
              }
              
              if (object@cov) {
                  Var <- rep(unname(coef(object)[object@parN]), nrow(X))
              } else {
                  Var <- rep(1.0, nrow(X))
              }
             
              
              return(Var) 
              
          })

## *************************************************************************
##' npar method for class "covANOVA".
##'
##' npar method for the "covANOVA" class
##'
##' @param object An object with class "covANOVA"
##' 
##' @return The number of free parmaeters in a `covANOVA` covariance.
##'
##' @docType methods
##'
##' @rdname covANOVA-methods
##'
setMethod("npar",
          signature = signature(object = "covANOVA"),
          definition = function(object,  ...){
            object@parN
          })


setMethod("coef", 
          signature = signature(object = "covANOVA"), 
          definition = function(object){         
            res <- object@par
            names(res) <- object@kernParNames
            res
          })

setMethod("coef<-", 
          signature = signature(object = "covANOVA", value = "numeric"),
          definition = function(object, ..., value){
            if (length(value) != object@parN) {
              stop(sprintf("'value' has length %d but must have length %d",
                           length(value), object@parN))
            }
            object@par[] <- value
            object
          })

##***********************************************************************
##***********************************************************************
setMethod("coefLower", 
          signature = signature(object = "covANOVA"),
          definition = function(object){
              res <- object@parLower
              names(res) <- object@kernParNames
              res
          })

setMethod("coefLower<-",
          signature = signature(object = "covANOVA"),
          definition = function(object, ..., value){
            if (length(value) != object@parN) {
              stop(sprintf("'value' must have length %d", object@parN))
            }
            object@parLower[] <- value
            object   
          })

setMethod("coefUpper", 
          signature = signature(object = "covANOVA"),
          definition = function(object){
              res <- object@parUpper
              names(res) <- object@kernParNames
              res            
          })

setMethod("coefUpper<-",
          signature = signature(object = "covANOVA"),
          definition = function(object, ..., value){
            if (length(value) != object@parN) {
              stop(sprintf("'value' must have length %d", object@parN))
            }
            object@parUpper[] <- value
            object   
          })

##***********************************************************************
## coercion method to cleanly extract the kernel slot
## XXX to be written
##***********************************************************************

## setAs("covANOVA", "function", function(from) from@kernel)

##***********************************************************************
## scores method 
## 
## Note that the scores method can only be used when the weight matrix
## is known.
##***********************************************************************
setMethod("scores",
          signature = "covANOVA", 
          definition = function(object, X, weights, ...) {
              
              X <- as.matrix(X)
              n <- nrow(X)
              d <- ncol(X)
              if (any(is.na(X))) stop("'X' must not contain NA elements") 
              Cov <- covMat(object = object, X = X, compGrad = TRUE)
              
              dCov <- attr(Cov, "gradient")
              if (any(dim(dCov) != c(n, n, object@parN))) {
                  stop("\"gradient\" attribute with wrong dimension")  
              }
              lt <- lower.tri(matrix(NA, nrow = n , ncol = n), diag = TRUE)
              agFun <- function(mat) sum(weights * mat[lt])
              scores <- apply(dCov, MARGIN = 3L, FUN = agFun)
              return(scores)
                            
          })
          


##***********************************************************************
## The 'show' method must show the kernel name and parameter structure.
## It should also provide information of the parameterisation of the
## structure itself (sharing of the parameters across inputs).
##
##' show method for class "covANOVA"
##' @aliases show,covANOVA-method
##'
##' @param object XXX
##' 
##' @docType methods
##'
##' @rdname covANOVA-methods
##'
##***********************************************************************
setMethod("show", 
          signature = signature(object = "covANOVA"), 
          definition = function(object){
            cat("ANOVA product covariance kernel\n")
            cat(paste("o Description:"), object@label, "\n")
            cat(sprintf("o Dimension 'd' (nb of inputs): %d\n", object@d))
            cat(sprintf("o 'k1' Function used: %s\n", object@k1Fun1Char))
#             cat(paste("o Kernel depending on: \"", 
#                       argNames[1], "\", \"", argNames[2], "\"\n", sep=""))
            cat(paste("o Parameters: ",
                      paste(sprintf("\"%s\"", object@kernParNames),
                              collapse = ", "),
                      "\n",sep = ""))
            cat(sprintf("o Number of parameters: %d\n", object@parN))
            if (object@hasGrad) {
                cat("o Analytical gradient is provided.\n")
            }
            cat("o Param. values: \n")
            co <- cbind(coef(object), coefLower(object), coefUpper(object))
            colnames(co) <- c("Value", "Lower", "Upper")
            print(co)
          })


