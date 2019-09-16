##
##  MINPACK-1 Least Squares Fitting Library
##
##  Original public domain version by B. Garbow, K. Hillstrom, J. More'
##    (Argonne National Laboratory, MINPACK project, March 1980)
##
##  Tranlation to C Language by S. Moshier (moshier.net)
##
##  Enhancements and packaging by C. Markwardt
##    (comparable to IDL fitting routine MPFIT
##     see http://cow.physics.wisc.edu/~craigm/idl/idl.html)
##
##  Header file defining constants, data structures and functions of
##    mpfit library
##    $Id: mpfit.h,v 1.16 2016/06/02 19:14:16 craigm Exp $
##

##  This is a C library.  Allow compilation with a C++ compiler

##  MPFIT version string
import math

const
  MPFIT_VERSION* = "1.3"

##  Definition of a parameter constraint structure

type
  mp_par_struct* {.bycopy.} = object
    fixed*: cint               ##  1 = fixed; 0 = free
    limited*: array[2, cint]    ##  1 = low/upper limit; 0 = no limit
    limits*: array[2, cdouble]  ##  lower/upper limit boundary value
    parname*: cstring          ##  Name of parameter, or 0 for none
    step*: cdouble             ##  Step size for finite difference
    relstep*: cdouble          ##  Relative step size for finite difference
    side*: cint ##  Sidedness of finite difference derivative
              ## 		        0 - one-sided derivative computed automatically
              ## 		        1 - one-sided derivative (f(x+h) - f(x)  )/h
              ## 		       -1 - one-sided derivative (f(x)   - f(x-h))/h
              ## 		        2 - two-sided derivative (f(x+h) - f(x-h))/(2*h)
              ## 			3 - user-computed analytical derivatives
              ##
    deriv_debug*: cint ##  Derivative debug mode: 1 = Yes; 0 = No;
                     ##
                     ##                        If yes, compute both analytical and numerical
                     ##                        derivatives and print them to the console for
                     ##                        comparison.
                     ##
                     ## 		       NOTE: when debugging, do *not* set side = 3,
                     ## 		       but rather to the kind of numerical derivative
                     ## 		       you want to compare the user-analytical one to
                     ## 		       (0, 1, -1, or 2).
                     ##
    deriv_reltol*: cdouble     ##  Relative tolerance for derivative debug
                         ## 			  printout
    deriv_abstol*: cdouble     ##  Absolute tolerance for derivative debug
                         ## 			  printout


##  Just a placeholder - do not use!!

type
  mp_iterproc* = proc () {.cdecl.}

##  Definition of MPFIT configuration structure

const
  MP_NO_ITER* = (-1)            ##  No iterations, just checking

type
  mp_config_struct* {.bycopy.} = object
    ftol*: cdouble ##  NOTE: the user may set the value explicitly; OR, if the passed
                 ##      value is zero, then the "Default" value will be substituted by
                 ##      mpfit().
    ##  Relative chi-square convergence criterium Default: 1e-10
    xtol*: cdouble             ##  Relative parameter convergence criterium  Default: 1e-10
    gtol*: cdouble             ##  Orthogonality convergence criterium       Default: 1e-10
    epsfcn*: cdouble           ##  Finite derivative step size               Default: MP_MACHEP0
    stepfactor*: cdouble       ##  Initial step bound                     Default: 100.0
    covtol*: cdouble           ##  Range tolerance for covariance calculation Default: 1e-14
    maxiter*: cint ##  Maximum number of iterations.  If maxiter == MP_NO_ITER,
                 ##                      then basic error checking is done, and parameter
                 ##                      errors/covariances are estimated based on input
                 ##                      parameter values, but no fitting iterations are done.
                 ## 		     Default: 200
                 ##
    maxfev*: cint ##  Maximum number of function evaluations, or 0 for no limit
                ## 		     Default: 0 (no limit)
    nprint*: cint              ##  Default: 1
    douserscale*: cint         ##  Scale variables by user values?
                     ## 		     1 = yes, user scale values in diag;
                     ## 		     0 = no, variables scaled internally (Default)
    nofinitecheck*: cint       ##  Disable check for infinite quantities from user?
                       ## 			0 = do not perform check (Default)
                       ## 			1 = perform check
                       ##
    iterproc*: mp_iterproc     ##  Placeholder pointer - must set to 0


##  Definition of results structure, for when fit completes

type
  mp_result_struct* {.bycopy.} = object
    bestnorm*: cdouble         ##  Final chi^2
    orignorm*: cdouble         ##  Starting value of chi^2
    niter*: cint               ##  Number of iterations
    nfev*: cint                ##  Number of function evaluations
    status*: cint              ##  Fitting status code
    npar*: cint                ##  Total number of parameters
    nfree*: cint               ##  Number of free parameters
    npegged*: cint             ##  Number of pegged parameters
    nfunc*: cint               ##  Number of residuals (= num. of data points)
    resid*: ptr cdouble         ##  Final residuals
                     ## 			  nfunc-vector, or 0 if not desired
    xerror*: ptr cdouble        ##  Final parameter uncertainties (1-sigma)
                      ## 			  npar-vector, or 0 if not desired
    covar*: ptr cdouble         ##  Final parameter covariance matrix
                     ## 			  npar x npar array, or 0 if not desired
    version*: array[20, char]   ##  MPFIT version string


##  Convenience typedefs

type
  mp_par* = mp_par_struct
  mp_config* = mp_config_struct
  mp_result* = mp_result_struct

##  Enforce type of fitting function

type
  mp_func* = proc (m: cint; n: cint; x: ptr cdouble; fvec: ptr cdouble;
                   dvec: ptr ptr cdouble; private_data: var pointer): cint {.cdecl.} ##  Number of functions (elts of fvec)
                                                                      ##  Number of variables (elts of x)
                                                                      ##  I - Parameters
                                                                      ##  O - function values
                                                                      ##  O - function derivatives (optional)

##  I/O - function private data
##  Error codes

const
  MP_ERR_INPUT* = (0)           ##  General input parameter error
  MP_ERR_NAN* = (-16)           ##  User function produced non-finite values
  MP_ERR_FUNC* = (-17)          ##  No user function was supplied
  MP_ERR_NPOINTS* = (-18)       ##  No user data points were supplied
  MP_ERR_NFREE* = (-19)         ##  No free parameters
  MP_ERR_MEMORY* = (-20)        ##  Memory allocation error
  MP_ERR_INITBOUNDS* = (-21)    ##  Initial values inconsistent w constraints
  MP_ERR_BOUNDS* = (-22)        ##  Initial constraints inconsistent
  MP_ERR_PARAM* = (-23)         ##  General input parameter error
  MP_ERR_DOF* = (-24)           ##  Not enough degrees of freedom

##  Potential success status codes

const
  MP_OK_CHI* = (1)              ##  Convergence in chi-square value
  MP_OK_PAR* = (2)              ##  Convergence in parameter value
  MP_OK_BOTH* = (3)             ##  Both MP_OK_PAR and MP_OK_CHI hold
  MP_OK_DIR* = (4)              ##  Convergence in orthogonality
  MP_MAXITER* = (5)             ##  Maximum number of iterations reached
  MP_FTOL* = (6)                ##  ftol is too small; no further improvement
  MP_XTOL* = (7)                ##  xtol is too small; no further improvement
  MP_GTOL* = (8)                ##  gtol is too small; no further improvement

##  Double precision numeric constants

const
  MP_MACHEP0* = 2.220446e-16
  MP_DWARF* = 2.2250739e-308
  MP_GIANT* = 1.7976931e+308

when sizeof(float) == 4:
  ##  Float precision
  const
    MP_MACHEP0* = 1.19209e-07
    MP_DWARF* = 1.17549e-38
    MP_GIANT* = 3.40282e+38
const
  MP_RDWARF* = (sqrt(MP_DWARF * 1.5) * 10)
  MP_RGIANT* = (sqrt(MP_GIANT) * 0.1)

##  External function prototype declarations

proc mpfit*(funct: mp_func; m: cint; npar: cint; xall: ptr cdouble; pars: ptr mp_par;
            config: ptr mp_config; private_data: pointer; result: ptr mp_result): cint {.
    cdecl, importc: "mpfit", dynlib: "libmpfit.so".}
##  C99 uses isfinite() instead of finite()

# when defined(__STDC_VERSION__) and __STDC_VERSION__ >= 199901:
#   template mpfinite*(x: untyped): untyped =
#     isfinite(x)

#   ##  Microsoft C uses _finite(x) instead of finite(x)
# elif defined(_MSC_VER) and _MSC_VER:
#   template mpfinite*(x: untyped): untyped =
#     _finite(x)

#   ##  Default is to assume that compiler/library has finite() function
# else:
template mpfinite*(x: untyped): untyped =
  finite(x)
