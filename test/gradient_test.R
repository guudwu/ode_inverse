set.seed(0)
dimension <- 4
time_point <- 0:10
intercept <- c(0,1e-2)

source('linear_ode_generation/linear_ode_generation.R')

linear_ode <-
  linear_ode_generation (
    dimension
    , time_point
    , orthogonal_transformation = list(c(1,dimension))
    , row_column_permutation = FALSE
    , intercept = intercept
  )

linODE <- function ( time , state , pars )
{
  res <- pars[[1]] %*% state
  if ( ! is.null ( pars[[2]] ) )
    res <- res + pars[[2]]
  return ( list(res) )
}

source('gradient.R')

linear_ode$initial <- linear_ode$observation[1,-1]

gradient <-
  .gradient ( linear_ode )

require('deSolve')

step_length <- 1e-3

rel_coefficiet <- lapply ( 1:(dimension**2) , function(index)
{
  perturb <- matrix ( 0 , dimension , dimension )
  perturb[index] <- 1

  perturb_res <-
    ode (
      linear_ode$initial
      , linear_ode$observation[,1]
      , linODE
      , list (
          linear_ode$coefficient + step_length * perturb
          , linear_ode$intercept
        )
    ) [,-1]

  ret <-
    sum (
      ( perturb_res - linear_ode$observation[,-1]
        - step_length * t(gradient$coefficient[[index]])
      ) **2
    ) / sum (
      ( perturb_res - linear_ode$observation[,-1] ) ** 2
    )

  return(ret)
} )

print(rel_coefficiet)

rel_initial <- lapply ( 1:dimension , function(index)
{
  perturb <- numeric(dimension)
  perturb[index] <- step_length

  perturb_res <-
    ode (
      linear_ode$initial + perturb
      , linear_ode$observation[,1]
      , linODE
      , list ( linear_ode$coefficient , linear_ode$intercept )
    ) [,-1]

  ret <-
    sum (
      ( perturb_res - linear_ode$observation[,-1]
        - step_length * t(gradient$initial[[index]])
      ) **2
    ) / sum (
      ( perturb_res - linear_ode$observation[,-1] ) ** 2
    )

  return(ret)
} )

print(rel_initial)

if ( !is.null(linear_ode$intercept) )
{
  rel_intercept <- lapply ( 1:dimension , function(index)
  {
    perturb <- numeric(dimension)
    perturb[index] <- step_length

    perturb_res <-
      ode (
        linear_ode$initial
        , linear_ode$observation[,1]
        , linODE
        , list ( linear_ode$coefficient , linear_ode$intercept + perturb )
      ) [,-1]

    ret <-
      sum (
        ( perturb_res - linear_ode$observation[,-1]
          - step_length * t(gradient$intercept[[index]])
        ) **2
      ) / sum (
        ( perturb_res - linear_ode$observation[,-1] ) ** 2
      )

    return(ret)
  } )

  print(rel_intercept)
}

