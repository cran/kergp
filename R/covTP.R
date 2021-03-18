setClass("covTP",   	
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
##' Creator for the class \code{"covTP"}.
##'
##' A tensor-product kernel on the \eqn{d}-dimensional euclidean space
##' takes the form \deqn{K(\mathbf{x},\,\mathbf{x}') = \sigma^2
##' \prod_{\elll = 1}^d \kappa(r_\ell)}{K(x1, x2) = \sigma^2 * prod_l
##' k1(r[l])} where \eqn{\kappa(r)}{k1(r)} is a suitable correlation
##' kernel for a one-dimensional input, and \eqn{r_\ell}{r[l]} is
##' given by \eqn{r_\ell := [x_\ell - x'_\ell] / \theta_\ell}{r[l] = (x2[l] -
##' x2[l]) / \theta[l]} for \eqn{\ell = 1}{l =1} to \eqn{d}.
##' 
##' In this default form, the tensor-product kernel depends on \eqn{d
##' + 1} parameters: the \emph{ranges} \eqn{\theta_\ell >0}{\theta[l] >
##' 0} and the variance \eqn{\sigma^2}{\sigma^2}.
##'
##' An \emph{isotropic} form uses the same range \eqn{\theta}{\theta}
##' for all inputs, i.e. sets \eqn{\theta_\ell =
##' \theta}{\theta[l]=\theta} for all \eqn{\ell}{l}. This is obtained by
##' using \code{iso = TRUE}.
##'
##' A \emph{correlation} version uses \eqn{\sigma^2 = 1}{\sigma^2 =
##' 1}. This is obtained by using \code{cov = "corr"}.
##'
##' Finally, the correlation kernel \eqn{\kappa(r)}{k1(r)} can depend on
##' a "shape" parameter, e.g. have the form
##' \eqn{\kappa(r;\,\alpha)}{k1(r, \alpha)}. The extra shape parameter
##' \eqn{\alpha}{\alpha} will be considered then as a parameter of the
##' resulting tensor-product kernel, making it possible to estimate it
##' by ML along with the range(s) and the variance.
##' 
##' @title Creator for the Class \code{"covTP"}
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
##' parameters in the \code{covTP} object with their name suffixed by
##' the dimension. When \code{iso1} is \code{1} only one shape
##' parameter will be created in the \code{covTP} object.
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
##' @return An object with class \code{"covTP"}.
##'
##' %% @note When \code{k1Fun1} has more than one formal argument, the
##' %% arguments with position \code{> 1} are assumed to be parameters of
##' %% the model. Examples are functions with formals \code{function(x,
##' %% shape = 1.0)} or \code{function(x, alpha = 2.0, beta = 3.0)},
##' %% corresponding to vector of parameter names \code{c("shape")} and
##' %% \code{c("alpha", "beta")}.
##'
##' @examples
##'
##' if (require(DiceKring)) {
##'     ## a 16-points factorial design and the corresponding response
##'     d <- 2; n <- 16; x <- seq(from = 0.0, to = 1.0, length.out = 4)
##'     X <- expand.grid(x1 = x, x2 = x)
##'     y <- apply(X, 1, DiceKriging::branin)
##'
##'     ## kriging model with matern5_2 covariance structure, constant
##'     ## trend. A crucial point is to set the upper bounds!
##'     mycov <- covTP(k1Fun1 = k1Fun1Matern5_2, d = 2, cov = "homo")
##'     coefUpper(mycov) <- c(2.0, 2.0, 1e10)
##'     mygp <- gp(y ~ 1, data = data.frame(X, y),
##'                cov = mycov, multistart = 100, noise = FALSE)
##' 
##'     nGrid <- 50; xGrid <- seq(from = 0, to = 1, length.out = nGrid)
##'     XGrid <- expand.grid(x1 = xGrid, x2 = xGrid)
##'     yGrid <- apply(XGrid, 1, Dicebranin)
##'     pgp <- predict(mygp, XGrid)$mean
##'
##'     mykm <- km(design = X, response = y)
##'     pkm <- predict(mykm, XGrid, "UK")$mean
##'     c("km" = sqrt(mean((yGrid - pkm)^2)),
##'       "gp" = sqrt(mean((yGrid - pgp)^2)))
##' }
##' 
covTP <- function(k1Fun1 = k1Fun1Gauss,
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
                  label = "Tensor product kernel",
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
    ## some names like "theta" or "sigma2".
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
            length(grep("sigma", kern1ParNames))) {
            stop("The names of the formals arguments of 'k1Fun1' ",
                 "are in conflict with default names  \"theta\" and ",
                 "\"sigma\". Please rename.")
        }
        
        ## Set the shape parameters to their default value
        par1 <- sapply(f, function(x) if (x == "") return(NA) else return(x))
                       
        if (!iso1) {
            parN1 <- d
            par1 <- rep(par1, length.out = d)
            kern1ParNamesTP <- paste(kern1ParNames, 1:d, sep = "_")
        } else {
            parN1 <- kern1ParN1
            kern1ParNamesTP <- kern1ParNames
        }
        
    } else {
        kern1ParNames <- kern1ParNamesTP <- character(0)
        parN1 <- 0L
    }
    
    parN <- parN1
    
    if (iso) {
        kernParNames <- c(kern1ParNamesTP, "theta")
        parN <- parN + 1L
    } else {
        kernParNames <- c(kern1ParNamesTP, paste0("theta_", 1:d))
        parN <- parN + d
    }

    if (cov) {
        kernParNames <- c(kernParNames, "sigma2")
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
    
    new("covTP", 
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
          signature = "covTP", 
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

              ## Moved after checkNames, or the matrix may be of mode
              ## "character!!!
              ## X <- as.matrix(X)
              
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
                  ## XXX unclear for now if 'Xnew' being a vector is
                  ## accepted or not
                  ## Xnew <- as.matrix(Xnew)
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
              if (parN1) {
                  psi <- par[1:parN1]
                  psi <- rep(psi, length.out = d)
              }
                  
              if (object@cov) sigma2 <- par[parN1 + nTheta + 1L]
              else sigma2 <- 1.0
      
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
                  Rs <- Kappas <- dKappas <- array(NA, dim = c(m, d))

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

                      if (kern1ParN1 > 0L) {
                          gKappas[ , ell]  <- attr(Kappaell, "gradient")[ , 1L]
                      }
                      
                      if (ell == 1) Kappa <- Kappaell
                      else Kappa <- Kappa * Kappaell
                      
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
                  
                  KappaSym <- array(1.0, dim = c(n, n))
                  KappaSym[sI$kL] <- KappaSym[sI$kU] <- Kappa
                  
                  GradSym <- array(0.0, dim = c(n, n, parN))
                  dimnames(GradSym) <- list(NULL, NULL, object@kernParNames)
                  Grad <- array(0.0, dim = c(n, n))

                  ## ==========================================================
                  ## 'shape' parameters of the kernel For now, 'kernParN1'
                  ## can only be 0 or 1...
                  ## ==========================================================
                  
                  if (kern1ParN1 > 0L) {
                      
                      ## stop("Using a kernel with a shape parameter is ",
                      ##    " not yet possible when 'compGrad' is TRUE")
                      
                      if (!object@iso1) {                          
                                                    
                          for (ell in 1L:d) {
                              
                              ind <- abs(Kappas[ , ell]) > 1e-8
                              
                              gKappas[ind, ell] <- sigma2 * gKappas[ind, ell] /
                                  Kappas[ind, ell] * Kappa[ind]
                              
                              if (any(!ind)) {
                                  gKappas[!ind, ell] <- sigma2 *
                                      gKappas[!ind, ell] *
                                          apply(Kappas[!ind, -ell, drop = FALSE], 1, prod)
                              }
                          
                              Grad[sI$kL] <- Grad[sI$kU] <- gKappas[ , ell]
                              GradSym[ , , ell] <- array(Grad, dim = c(n, n))
                              
                          }
                          
                      } else {                          
                          
                          for (ell in 1L:d) {
                              
                              ind <- abs(Kappas[ , ell]) > 1e-8
                              
                              gKappas[ind, ell] <- sigma2 * gKappas[ind, ell] /
                                  Kappas[ind, ell] * Kappa[ind]
                              
                              if (any(!ind)) {
                                  gKappas[!ind, ell] <- sigma2 *
                                      gKappas[!ind, ell] *
                                          apply(Kappas[!ind, -ell, drop = FALSE], 1, prod)
                              }
                              
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

                          ind <- abs(Kappas[ , ell]) > 1e-8

                          dKappas[ind, ell] <- - sigma2 * Rs[ind, ell] *
                              dKappas[ind, ell] / Kappas[ind, ell] *
                                  Kappa[ind] / theta[ell]
                          
                          if (any(!ind)) {
                              dKappas[!ind, ell] <- - sigma2 * Rs[!ind, ell] *
                                  dKappas[!ind, ell] *
                                      apply(Kappas[!ind, -ell, drop = FALSE], 1, prod) /
                                          theta[ell]
                          }
                              
                          Grad[sI$kL] <- Grad[sI$kU] <- dKappas[ , ell]
                          GradSym[ , , parN1 + ell] <- array(Grad, dim = c(n, n))

                      }
                      
                  } else {
                      
                      ind <- abs(Kappa) > 1e-8

                      dKappas[ind, ] <- - sigma2 * Rs[ind,  ] *
                          dKappas[ind, ] / Kappas[ind, ] *
                              kappa[ind] / theta[1L]
                      
                      if (any(!ind)) {
                          
                          for (ell in 1L:d) {
                              dKappas[!ind, ell] <- - sigma2 * Rs[!ind,  ell] *
                                  dKappas[!ind, ell] *
                                      apply(Kappas[!ind, -ell], 1L, prod) /
                                          theta[1L]   
                          }
                          
                      }
                         
                      Grad[sI$kL] <- Grad[sI$kU] <- apply(dKappas, 1L, sum)
                      GradSym[ , , parN1 + 1L] <- array(Grad, dim = c(n, n))
                  }
                 
                  ## variance parameter Now the diagonal of the slice is 1.0
                  if (object@cov) {
                      GradSym[ , , parN]  <- KappaSym
                      KappaSym <- sigma2 * KappaSym
                  }
                  
                  attr(KappaSym, "gradient") <- GradSym
                  return(KappaSym)
                  
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
                              
                              if (ell == 1) Kappa <- Kappaell
                              else Kappa <- Kappa * Kappaell
                          }    
                         
                          if (object@cov) Kappa <- sigma2 * Kappa

                          dim(Kappa) <- c(n, nNew)
                          return(Kappa)

                          
                      } else {

                          ## ==================================================
                          ## We use a vectorized form of the matrices
                          ## 'K_ell' and their derivatives 'dK_ell /
                          ## dx_ell'
                          ## ==================================================
                          Rs <- Kappas <- dKappas <-
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
                                  Kappaell <-
                                      do.call(object@k1Fun1,
                                              list(Rs[ , ell]))
                              }
                              
                              Kappas[ , ell] <- Kappaell
                              dKappas[ , ell] <-
                                  array(attr(Kappaell, "der")[ , 1L],
                                        dim = c(n, nNew))
                              if (ell == 1) Kappa <- Kappaell
                              else Kappa <- Kappa * Kappaell
                          }
                  
                          Kappa <- array(Kappa, dim = c(n , nNew))
                          
                          ## ==================================================
                          ## ! It can be the case that 'R' contains
                          ## zeros, because a row of 'Xnew' can be
                          ## identical to a row of 'X'.
                          ## ==================================================

                          
                          der <- array(0.0, dim = c(n, nNew, d))
                          dimnames(der) <- list(NULL, NULL, colnames(X))
                          
                          for (ell in 1L:d) {
                              
                              ind <- abs(Kappas[ , ell]) > 1e-8
                              
                              dKappas[ind, ell] <- - sigma2 * dKappas[ind, ell] /
                                  Kappas[ind, ell] * Kappa[ind] / theta[ell]
                              
                              if (any(!ind)) {
                                  ## cat("XXX AAA\n")
                                  ## print(dim(Kappas[!ind, -ell]))
                                  dKappas[!ind, ell] <- - sigma2 * dKappas[!ind, ell] *
                                      apply(Kappas[!ind, -ell, drop = FALSE], 1, prod) /
                                          theta[ell]
                              }

                              ## Mind the derivative of the absolute value!!!
                              ## ind <- Rs[ , ell] < 0
                              ## dKappas[ind , ell] <- -dKappas[ind , ell]

                              der[ , , ell] <- dKappas[ , ell] 
                              
                          }
                          
                          if (object@cov) {
                              Kappa <- sigma2 * Kappa
                          }
                          
                          ## note that 'kappa' already has a 'der' attribute! 
                          attr(Kappa, "der") <- der
                          return(Kappa)
                          
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
                          
                          if (ell == 1) Kappa <- Kappaell
                          else Kappa <- Kappa * Kappaell
                      
                      }

                      ## =======================================================
                      ## Reshape to a symmetric matrix. We assume that
                      ## the kernel function is a correlation function
                      ## so the diagonal of the returned matrix
                      ## contains ones.
                      ## =======================================================
                      
                      if (object@cov) Kappa <- sigma2 * Kappa
                      KappaSym <- array(sigma2, dim = c(n, n))
                      KappaSym[sI$kL] <- KappaSym[sI$kU] <- Kappa                    

                      return(KappaSym)
                      
                  }

              }
              
          })



## *************************************************************************
## 'varVec' method: compute the variance vector.
## *************************************************************************
setMethod("varVec",
          signature = "covTP", 
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
##' npar method for class "covTP".
##'
##' npar method for the "covTP" class
##'
##' @param object An object with class "covTP"
##' 
##' @return The number of free parmaeters in a `covTP` covariance.
##'
##' @docType methods
##'
##' @rdname covTP-methods
##'
setMethod("npar",
          signature = signature(object = "covTP"),
          definition = function(object,  ...){
            object@parN
          })


setMethod("coef", 
          signature = signature(object = "covTP"), 
          definition = function(object){         
            res <- object@par
            names(res) <- object@kernParNames
            res
          })

setMethod("coef<-", 
          signature = signature(object = "covTP", value = "numeric"),
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
          signature = signature(object = "covTP"),
          definition = function(object){
              res <- object@parLower
              names(res) <- object@kernParNames
              res  
          })

setMethod("coefLower<-",
          signature = signature(object = "covTP"),
          definition = function(object, ..., value){
            if (length(value) != object@parN) {
              stop(sprintf("'value' must have length %d", object@parN))
            }
            object@parLower[] <- value
            object   
          })

setMethod("coefUpper", 
          signature = signature(object = "covTP"),
          definition = function(object){
              res <- object@parUpper
              names(res) <- object@kernParNames
              res          
          })

setMethod("coefUpper<-",
          signature = signature(object = "covTP"),
          definition = function(object, ..., value){
            if (length(value) != object@parN) {
              stop(sprintf("'value' must have length %d", object@parN))
            }
            object@parUpper[] <- value
            object   
          })

## Commented out on 2021-03-17. Not really needed inasmuch 'covTP' inherits
## from 'covAll'
## setMethod("checkX",
##           signature = signature(object = "covTP", X = "matrix"),
##           definition = function(object, X, ...){
##               callNextMethod(object, X, ...)
##           })

## setMethod("checkX",
##           signature = signature(object = "covTP", X = "data.frame"),
##           definition = function(object, X, ...){
##               callNextMethod(object, X, ...)
##           })

##***********************************************************************
## coercion method to cleanly extract the kernel slot
## XXX to be written
##***********************************************************************

## setAs("covTP", "function", function(from) from@kernel)

##***********************************************************************
## scores method 
## 
## Note that the scores method can only be used when the weight matrix
## is known.
##***********************************************************************
setMethod("scores",
          signature = "covTP", 
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
##' show method for class "covTP"
##' @aliases show,covTP-method
##'
##' @param object XXX
##' 
##' @docType methods
##'
##' @rdname covTP-methods
##'
##***********************************************************************
setMethod("show", 
          signature = signature(object = "covTP"), 
          definition = function(object){
            cat("Tensor-Product covariance kernel\n")
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


## ****************************************************************************
## DO NOT Roxigenize this!
##
##' Coerce a \code{covTP} object representing a Tensor-Product
##' covariance kernel on the \eqn{d}-dimensional Euclidean space
##' into a list containing \eqn{d} one-dimensional kernels.
##' 
##' @title Coerce a \code{covTP} Object into a List
##'
##' @param x A \code{covTP} object representing a Tensor-Product
##' covariance kernel.
##'
##' @return A list with length \code{d} or \code{d + 1} where \code{d}
##' is the "dimension" slot \code{x@d} of the object \code{x}. The
##' first \code{d} elements of the list are one-dimensional covariance
##' kernel objects with class \code{"covTP"}. When \code{x} is a
##' \emph{covariance} kernel (as opposed to a \emph{correlation}
##' kernel), the list contains one more element which gives the
##' variance.
##'
##' @section Caution: When \code{x} is not a correlation kernel the
##' \code{(d + 1)}-th element of the returned list may be different in
##' future versions: it may be a constant covariance kernel.
##'
##' @seealso \code{\link{covTP}} and \code{\link{covTP-class}}.
##'
##' @examples
##' set.seed(123)
##' d <- 6
##' myCov1 <- covTP(d = d, cov = "corr")
##' coef(myCov1) <- as.vector(simulPar(myCov1, nsim = 1))
##' as.list(myCov1)
##'
##' ## more examples and check the value of a 'covMat'
##' L <- list()
##' myCov <- list()
##'
##' myCov[[1]] <- covTP(d = d, cov = "corr")
##' coef(myCov[[1]]) <- as.vector(simulPar(myCov[[1]], nsim = 1))
##' L[[1]] <- as.list(myCov[[1]])
##'
##' myCov[[2]] <- covTP(k1Fun1 = k1Fun1PowExp, d = d, cov = "corr")
##' coef(myCov[[2]]) <- as.vector(simulPar(myCov[[2]], nsim = 1))
##' L[[2]] <- as.list(myCov[[2]])
##'
##' myCov[[3]] <- covTP(k1Fun1 = k1Fun1PowExp, d = d, iso1 = 0L, cov = "corr")
##' coef(myCov[[3]]) <- as.vector(simulPar(myCov[[3]], nsim = 1))
##' L[[3]] <- as.list(myCov[[3]])
##'
##' n <- 10
##' X <- matrix(runif(n * d), nrow = n,
##'             dimnames = list(NULL, paste("x", 1:d, sep = "")))
##' for (iTest in 1:3) {
##'    C <- covMat(L[[iTest]][[1]], X[ , 1, drop = FALSE])
##'    for (j in 2:d) {
##'       C <- C * covMat(L[[iTest]][[j]], X[ , j, drop = FALSE])
##'    }
##'    CTest <- covMat(myCov[[iTest]], X)
##'    print(max(abs(abs(C - CTest))))
##' }
##' 
as.list.covTP <- function(x) {

    d <- x@d

    if (x@kern1ParN1 > 1) {
        stop("A value > 1 for the '@parN1' slot is not accepted yet") 
    }
    
    parN1 <- x@parN1    
    parN <- x@parN
    par <- coef(x)
    
    nTheta <- ifelse(x@iso, 1L, d)
    theta <- par[parN1 + (1L:nTheta)]
    theta <- rep(theta, length.out = d)
    
    if (x@cov > 0L) iVar <- parN
    else iVar <- integer(0)
    
    L <- list()
 
    for (ell in 1:d) {
                                    
        if (parN1 > 0) {
            if (x@iso1) indell <- c(1L, 1L + ell)
            else indell <- c(ell, ell + d)
            
        } else {
            indell <- ell
        }
        inedll <- c(indell, iVar)
        parN1ell <- length(indell)
        parell <- x@par[indell]
        parLowerell <- x@parLower[indell]
        parUpperell <- x@parUpper[indell]
        kernParNamesell <- c(x@kern1ParNames, "theta")
        L[[ell]] <-  new("covTP", 
                         k1Fun1 = x@k1Fun1,
                         k1Fun1Char = x@k1Fun1Char,
                             hasGrad = x@hasGrad,
                             cov = 0L,
                             iso = 0L,
                             iso1 = 0L,
                             label = as.character(x@label),
                             d = 1L,
                             inputNames = x@inputNames[[ell]],
                             par = as.numeric(parell),
                             parLower = as.numeric(parLowerell),
                             parUpper = as.numeric(parUpperell),
                             kern1ParN1 = as.integer(x@kern1ParN1),
                             parN1 = as.integer(parN1 > 0),
                             parN = as.integer(length(parell)),
                             kern1ParNames = as.character(x@kern1ParNames),
                             kernParNames = as.character(kernParNamesell))
    }
        
    if (x@cov) {
        L[[ell + 1]] <- x@par[parN]
    }

    L
    
}

setMethod(f = "as.list",
          signature = signature(x = "covTP"),
          definition = function(x) as.list.covTP(x))

