#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>  
#define NODEBUG

/*============================================================================*
 * Author: Yves                                                               *
 *                                                                            *
 * In this version, the argument 'par' can be used and passed to the kernel   *
 * function from the C-level. Moreover, the gradient of the matrix with       *
 * respect to ONE of the parameters can be computed if wanted, and then       *
 * returned as the "gradient" attribute of the covariance. Note that          *
 * "gradient" is  NOT the full gradient (which is an huge array) but only one *
 * element of it.                                                             *
 *                                                                            *
 *                                                                            *
 *============================================================================*/

SEXP covMat_covMan(SEXP fun,      // A cov. kernel with TWO vector args + par
		    SEXP Xt,       // Transpose of the spatial design matrix
		    SEXP par,      // Vector of (real) parameters' values
		    SEXP compGrad, // Integer 0 or 1: compute gradient matrix?
		    SEXP index,    // Index for grad. INTEGER 0 <= < length(par) 
		    SEXP rho) {    // An R environment
  
  int  i, j, k, n, d;   //*icompGrad;
  double *rxt = REAL(Xt), *rx1, *rx2, *rCov, *rpar;
  SEXP R_fcall, Cov, x1, x2;
  
  if (!isFunction(fun)) error("'fun' must be a function");
  if (!isMatrix(Xt)) error("'Xt' must be a matrix");
  if(!isEnvironment(rho)) error("'rho' should be an environment");

  /* find the number of rows and cols in 'X'  */
  SEXP dimXt;
  Xt = coerceVector(Xt, REALSXP);
  dimXt = getAttrib(Xt, R_DimSymbol);
  d = INTEGER(dimXt)[0];
  n = INTEGER(dimXt)[1];
  //icompGrad[0] = INTEGER(compGrad)[0];

#ifdef DEBUG 
  Rprintf("d = %d n = %d\n", d, n);
#endif

  PROTECT(Cov = allocMatrix(REALSXP, n, n)); //allocate the n x n matrix  
  PROTECT(x1 = allocVector(REALSXP, d));
  PROTECT(x2 = allocVector(REALSXP, d));

  rCov = REAL(Cov);
  rx1 = REAL(x1);
  rx2 = REAL(x2);
  rpar = REAL(par);

#ifdef DEBUG 
  Rprintf("REAL fait\n");
#endif

  PROTECT(R_fcall = lang4(fun, x1, x2, par));
  SETCADDDR(R_fcall, par);
  // SETCAD4R(R_fcall, compGrad);
  
  if ( INTEGER(compGrad)[0] ) {

    SEXP dCov, kernValue, dkernValue, attrNm;
    double *rdCov;
    int iindex, npar;

    npar = LENGTH(par);

#ifdef DEBUG 
    Rprintf("npar = %d\n", npar);
#endif

    iindex = INTEGER(index)[0];
  
    PROTECT(dCov = allocMatrix(REALSXP, n, n)); //allocate the n x n matrix  
    PROTECT(kernValue = allocVector(REALSXP, 1));
    PROTECT(dkernValue = allocVector(REALSXP, npar));
    
    PROTECT(attrNm = NEW_CHARACTER(1));
    SET_STRING_ELT(attrNm, 0, mkChar("gradient"));
    rdCov = REAL(dCov);

    //Rprintf("Now loop iindex = %d\n", iindex);

#ifdef DEBUG 
    Rprintf("CHECK\n");
#endif

    for (i = 0; i < n; i++) {
      
      for (k = 0; k < d; k++) {
	rx1[k] = rxt[i*d + k];
      }
      
      SETCADR(R_fcall, x1);
      
      for (j = i; j < n; j++) {
	
	//Rprintf("i = %d; j = %d\n", i, j);
	
	for (k = 0; k < d; k++) {
	  rx2[k] = rxt[j*d + k]; 
	}
	
	SETCADDR(R_fcall, x2);
	kernValue = eval(R_fcall, rho);
	// Rprintf("kernValue = %7.4f\n", REAL(kernValue)[0]);

	rCov[i + j*n] = REAL(kernValue)[0]; 
	rCov[j + i*n] = rCov[i + j*n];
	dkernValue = GET_ATTR(kernValue, attrNm);

	//Rprintf("Long attr. %d\n", LENGTH(dkernValue));
	//Rprintf("Elt 1 %7.3f\n", REAL(dkernValue)[0]);

	rdCov[i + j*n] = REAL(dkernValue)[iindex]; 
	rdCov[j + i*n] = rdCov[i + j*n];
      }
      
    }
    
    // attribut "gradient" de 'Cov'.
    SET_ATTR(Cov, attrNm, dCov);

    UNPROTECT(8);
    return(Cov);

  } else {
   
    for (i = 0; i < n; i++) {
      
      for (k = 0; k < d; k++) {
	rx1[k] = rxt[i*d + k];
      }
      
      SETCADR(R_fcall, x1);
      
      for (j = i; j < n; j++) {
	
	//Rprintf("i = %d; j = %d\n", i, j);
	
	for (k = 0; k < d; k++) {
	  rx2[k] = rxt[j*d + k]; 
	}
	
	SETCADDR(R_fcall, x2);
	rCov[i + j*n] = REAL(eval(R_fcall, rho))[0];  
	rCov[j + i*n] = rCov[i + j*n];
      }
      
    }
    
    UNPROTECT(4);
    return(Cov);
    
  }

}

/*============================================================================*
 * covMatMat                                                                *
 *                                                                            *
 * At the time, no derivation is possible.                                    * 
 *                                                                            *
 *============================================================================*/

SEXP covMatMat_covMan(SEXP fun,      // kernel depends on 2 scalar sites + 1 par
		       SEXP X1t,      // Transpose of the spatial design : d*n1
		       SEXP X2t,      // Transpose of the spatial design : d*n2
		       SEXP par,      // vector of 'npar' param for the Kd kernel
		       SEXP compGrad, // NOT USED
		       SEXP index,    // NOT USED
		       SEXP rho) {    // An R environment
  
  int  i, j, k, n1, n2, d, p, npar;
  
  double *rx1t = REAL(X1t),  *rx2t = REAL(X2t),  
    *rx1, *rx2,
    *rCov, *rpar = REAL(par);
  
  SEXP dimX1t, dimX2t, x1, x2, R_fcall, Cov;
  
  if (!isFunction(fun)) error("'fun' must be a function");
  if (!isMatrix(X1t)) error("'X1t' must be a matrix"); 
  if (!isMatrix(X2t)) error("'X2t' must be a matrix");
  if(!isEnvironment(rho)) error("'rho' should be an environment");
  
  /* find the number of rows and cols in 'X'  */
  X1t = coerceVector(X1t, REALSXP);
  PROTECT(dimX1t = getAttrib(X1t, R_DimSymbol));
  d = INTEGER(dimX1t)[0];
  n1 = INTEGER(dimX1t)[1];
 
  X2t = coerceVector(X2t, REALSXP);
  PROTECT(dimX2t = getAttrib(X2t, R_DimSymbol));
  if (INTEGER(dimX2t)[0] != d) {
    error("'X1t' and 'X2t must have the same number of rows (number of inputs)");
  }
  n2 = INTEGER(dimX2t)[1]; 

#ifdef DEBUG 
  Rprintf("d = %d n1 = %d n2 =%d\n", d, n1);
#endif
  
  PROTECT(Cov = allocMatrix(REALSXP, n1, n2)); //allocate the n x n matrix  
  PROTECT(x1 = allocVector(REALSXP, d));
  PROTECT(x2 = allocVector(REALSXP, d));

  rCov = REAL(Cov);
  rx1 = REAL(x1);
  rx2 = REAL(x2);
  rpar = REAL(par);

  PROTECT(R_fcall = lang4(fun, x1, x2, par));
  
  /* check the parameters arrays  */
  npar = LENGTH(par);
  par = coerceVector(par, REALSXP);
  
  if ( INTEGER(compGrad)[0] ) {
    
    UNPROTECT(6);
    error("Gradient computation not implemented for covMatMat");
    
  } else {
    
    SEXP kernValue;
   
    PROTECT(kernValue = allocVector(REALSXP, 1));    
    SETCADDDR(R_fcall, par);

    for (i = 0; i < n1; i++) {
      
      for (k = 0; k < d; k++) {
	rx1[k] = rx1t[i*d + k];
      }
     
      SETCADR(R_fcall, x1);
      
      for (j = 0; j < n2; j++) {

	for (k = 0; k < d; k++) {
	  rx2[k] = rx2t[j*d + k];
	}

	SETCADDR(R_fcall, x2);
	kernValue = eval(R_fcall, rho);
	rCov[i + j * n1] = REAL(kernValue)[0];
      }
    }

    UNPROTECT(7);
    return(Cov);

  }

}
