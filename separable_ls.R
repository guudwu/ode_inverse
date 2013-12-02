# Separable Least-Squares algorithm to reconstruct a linear ODE system

separable_ls <- function (
  observation
  , eigen_complex = numeric(floor((ncol(observation)-1)/2))
  , eigen_real = numeric((ncol(observation)-1)%%2)
  , lower_bound_complex = -Inf*(1+1i)
  , upper_bound_complex = Inf*(1+1i)
  , lower_bound_real = -Inf
  , upper_bound_real = Inf
  , intercept = TRUE
  , method = c('dfo','lm')
  , dfo.parameter = list (
      npt = 4*length(eigen_complex) + 2*length(eigen_real) + 1
      , rhobeg = 1e-2
      , rhoend = 1e-10
      , maxfun = 1e8
    )
  , lm.parameter = minpack.lm::nls.lm.control (
      maxiter = 1024
      , maxfev = 1e8
    )
  , verbose = 0
)

# INPUT:
# observation: Each row is for one time point.
#   First column is for time points,
#   hence must be of ascending order without identical elements.
#   Each other column is a curve to fit.
# eigen_complex: Each item is for a pair of complex eigen-values
#   of the coefficient matrix.
#   The number of complex eigen-values is not changed in the iteration.
#   The values are used as starting value to find the optimal spectrum.
# eigen_real: Each item is for a real eigen-value of the coefficient
#   matrix.
#   The number of real eigen-values is not changed in the iteration.
#   The values are used as starting value to find the optimal spectrum.
#   Also notice that the total number of complex and real eigen-values
#   must be equal to the number of columns of "observation" minus one
#   (to exclude the time-point column).
# lower_bound_complex,upper_bound_complex: Lower,upper bounds to the
#   real,imaginary part respectively of the complex eigen-values.
#   Notice that complex eigen-values are in conjugate pairs,
#   hence one may set the lower bound of imaginary bound to 0,
#   without affecting the final result.
#   The length of the two vectors should be equal to length of
#   "eigen_complex", but a scalar will be automatically expanded.
# lower_bound_real,upper_bound_real: Lower,upper bounds to the real
#   eigen-values.
#   The length of the two vectors should be equal to length of
#   "eigen_real", but a scalar will be automatically expanded.
# intercept: Whether an intercept term is contained in the system.
# method: Methods used to find the optimal spectrum.
#   "dfo": Derivative-free optimization solver "NEWUOA","BOBYQA".
#   "lm": Levenberg-Marquardt solver with finite difference.
# dfo.parameter: Parameters for "dfo" solver.
#   Refer to manual page of "newuoa" in package "minqa".
# lm.parameter: Parameters for "lm" solver.
#   Refer to manual page of "nls.lm.control" in package "minpack.lm".
# verbose: Control output.
#   0: Disable outputs.
#   >=1: Summary of results for each solver.
#   2-4: verbose-1 will be passed as "iprint" value to "newuoa"
#        or "bobyqa" solver.
#        Refer to manual page of "newuoa" in package "minqa".

# OUTPUT
# rss: Optimal residual sum-of-squares between curve of reconstructed
#   ODE system and "observation".
# curve: Curve of reconstructed ODE system.
# coefficient: Coefficient matrix.
# intercept: Intercept term(NULL if not existed).

{

# Sanity check#{{{

temp <- sort(unique(observation[,1]))
if ( length(temp)!=nrow(observation) || all(temp!=observation[,1]) )
{
  stop('First column of "observation" is for time points, ' ,
    'which must be of ascending order and without identical elements.')
}

eigen_complex <- as.complex(eigen_complex)
num_complex_eigen = length(eigen_complex)

eigen_real <- as.numeric(eigen_real)
num_real_eigen = length(eigen_real)

if ( num_complex_eigen*2+num_real_eigen != ncol(observation)-1 )
{
  stop('Number of eigen-values dismatch "observation".')
}

lower_bound_complex <- as.complex(lower_bound_complex)
upper_bound_complex <- as.complex(upper_bound_complex)
lower_bound_real <- as.numeric(lower_bound_real)
upper_bound_real <- as.numeric(upper_bound_real)

if ( length(lower_bound_complex)==1 )
{
  lower_bound_complex <- rep ( lower_bound_complex , num_complex_eigen )
}
if ( length(upper_bound_complex)==1 )
{
  upper_bound_complex <- rep ( upper_bound_complex , num_complex_eigen )
}
if ( length(lower_bound_real)==1 )
{
  lower_bound_real <- rep ( lower_bound_real , num_real_eigen )
}
if ( length(upper_bound_real)==1 )
{
  upper_bound_real <- rep ( upper_bound_real , num_real_eigen )
}

intercept <- as.logical(intercept)

verbose <- as.integer(verbose)
if ( verbose<0 )
{
  verbose <- 0
}
if ( verbose>4 )
{
  verbose <- 4
}
#}}}

# Construct iteration variable and upper,lower bounds#{{{

parameter <- c (
  Im(eigen_complex)
  , Re(eigen_complex)
  , eigen_real
)

lower_bound <- c (
  Im(lower_bound_complex)
  , Re(lower_bound_complex)
  , lower_bound_real
)
upper_bound <- c (
  Im(upper_bound_complex)
  , Re(upper_bound_complex)
  , upper_bound_real
)
#}}}

# Solver#{{{

dfo_res <- NULL
lm_res <- NULL

source('objective.R')
for ( item in method )
{
# DFO solver#{{{
  if ( item=='dfo' )
  {
    if ( verbose >= 2 )
    {
      dfo.parameter = c ( dfo.parameter , iprint = verbose-1 )
    }
    require('minqa')
    if ( all(lower_bound==(-Inf)) && all(upper_bound==Inf) )
    {
      time_dfo <- system.time (
        dfo_res <- minqa::newuoa (
          parameter
          , fn = .objective
          , control = dfo.parameter
          , num_complex_eigen = num_complex_eigen
          , num_real_eigen = num_real_eigen
          , observation = observation
          , intercept = intercept
          , type = 'dfo'
        )
      )
    }
    else
    {
      time_dfo <- system.time (
        dfo_res <- minqa::bobyqa (
          parameter
          , fn = .objective
          , lower = lower_bound
          , upper = upper_bound
          , control = dfo.parameter
          , num_complex_eigen = num_complex_eigen
          , num_real_eigen = num_real_eigen
          , observation = observation
          , intercept = intercept
          , type = 'dfo'
        )
      )
    }
    if ( verbose >= 1 )
    {
      cat('DFO solver:\nStatus: ')
      if ( dfo_res$ierr==0 )
      {
        cat('Success.\n')
      }
      else
      {
        cat('Early termination.\n')
      }
      cat('Time:\n')
      print(time_dfo)
      cat('\nResidual sum-of-squares: ')
      cat(dfo_res$fval)
      cat('\nNumber of function evaluations: ')
      cat(dfo_res$feval)
      cat('\n')
    }
  }
#}}}
# LM solver#{{{
  if ( item=='lm' )
  {
    require('minpack.lm')
    time_lm <- system.time (
      lm_res <- minpack.lm::nls.lm (
        parameter
        , lower = lower_bound
        , upper = upper_bound
        , fn = .objective
        , jac = NULL
        , control = lm.parameter
        , num_complex_eigen = num_complex_eigen
        , num_real_eigen = num_real_eigen
        , observation = observation
        , intercept = intercept
        , type = 'lm'
      )
    )
    if ( verbose >= 1 )
    {
      cat('LM solver:\nStatus: ')
      if ( lm_res$info==0 )
      {
        cat('Improper input arguments.\n')
      }
      else if ( lm_res$info==5 || lm_res$info==9 )
      {
        cat('Early termination.\n')
      }
      else
      {
        cat('Success.\n')
      }
      cat('Time:\n')
      print(time_lm)
      cat('\nResidual sum-of-squares: ')
      cat(lm_res$deviance)
      cat('\nNumber of function evaluations: ')
      cat((lm_res$niter+1)*(2*num_complex_eigen+num_real_eigen+1))
      cat('\n')
    }
  }
#}}}
}
#}}}

# Post process#{{{

if ( !is.null(dfo_res) )
{
  if ( is.null(lm_res) || dfo_res$fval<lm_res$deviance )
  {
    par_optim <- dfo_res$par
    rss_optim <- dfo_res$fval
  }
  else
  {
    par_optim <- lm_res$par
    rss_optim <- lm_res$deviance
  }
}
else
{
  if ( is.null(lm_res) )
  {
    par_optim <- NULL
    rss_optim <- NULL
  }
  else
  {
    par_optim <- lm_res$par
    rss_optim <- lm_res$deviance
  }
}

if ( !is.null(par_optim) )
{
  res_post <- .objective (
    par_optim
    , num_complex_eigen = num_complex_eigen
    , num_real_eigen = num_real_eigen
    , observation = observation
    , intercept = intercept
    , type = 'post'
  )
}
#}}}

# Return#{{{

ret <- list (
  rss = rss_optim
  , curve = res_post$curve
  , coefficient = res_post$coefficient
)

if (!is.null(res_post$intercept))
{
  ret$intercept <- res_post$intercept
}

return(ret)
#}}}

}
