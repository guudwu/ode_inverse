# From a non-sparse estimation,
# to compute a sparse solution,
# by solving LASSO problem based on the linearization at the estimation.

lasso_linearize <- function (
  estimation
  , observation
)

# INPUT:
# estimation: A list containing estimated parameters
#   of a linear ODE system, which consists of following items:
#   coefficient: Coefficient matrix.
#   initial: Initial condition.
#   intercept: Intercept term.
#     Might be NULL or non-existent for a system without intercpet term.
#   curve: Numerical solution of the estimated system.
#     Of same structure and same time points as argument "observation".
# observation: Data to be fit.
#   Each row is for a time point.
#   Each column is a curve.
#   First column is time points.

# OUTPUT:
# A matrix which returns the optimal primal variables of the LASSO
#   problem.
# Each column is for a penalty parameter \lambda.
# Variables in each column is of order:
#   coefficient matrix(column-oriented)
#   initial condition
#   intercept term(if existed)
# Notice that only the non-zero structure matters here since
# the LASSO problem is formulated using a linear approximation.
# A further model refinement is needed after determining
# the sparsity structure.

{

  intercept <- !is.null(estimation$intercept)
  dimension <- ncol(observation)-1
  time_point <- observation[,1]

  source('gradient.R')

  temp <-
    list (
      coefficient = estimation$coefficient
      , initial = estimation$curve[1,-1]
      , time_point = time_point
    )
  if ( intercept )
  {
    temp$intercept <- estimation$intercept
  }

  gradient_res <-
    .gradient (
      temp
    )

  lasso_basis <-
    matrix (
      0
      , dimension * length(time_point)
      , dimension**2 + dimension + dimension * as.integer(intercept)
    )

  lapply ( 1:(dimension**2) , function(index)
  {
    lasso_basis[,index] <<-
      as.numeric(t(gradient_res$coefficient[[index]]))
    return()
  } )
  lapply ( 1:dimension , function(index)
  {
    lasso_basis[,(dimension**2+index)] <<-
      as.numeric(t(gradient_res$initial[[index]]))
    return()
  } )
  if ( intercept )
  {
    lapply ( 1:dimension , function(index)
    {
      lasso_basis[,(dimension**2+dimension+index)] <<-
        as.numeric(t(gradient_res$intercept[[index]]))
      return()
    } )
  }

  lasso_residue <-
    as.numeric (
      ( observation - estimation$curve ) [,-1]
    )

  parameter <-
    c (
      as.numeric(estimation$coefficient)
      , as.numeric(estimation$curve[1,-1])
    )
  if ( intercept )
  {
    parameter <- c ( parameter , as.numeric(estimation$intercept) )
  }

  lasso_residue <- lasso_residue + lasso_basis %*% parameter
  print(dim(lasso_basis))
  print(length(lasso_residue))

  require('genlasso')
  genlasso_res <-
    genlasso (
      as.numeric(lasso_residue)
      , lasso_basis
      , diag(length(parameter))[1:(dimension**2),]
    )

  return(genlasso_res$beta)

}
