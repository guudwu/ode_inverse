set.seed(0)
dimension <- 9
time_point <- 0:100
intercept <- c(0,1e-2)

source('linear_ode_generation/linear_ode_generation.R')

linear_ode <-
  linear_ode_generation (
    dimension
    , time_point
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

rel_rss <- lapply ( 1:dimension , function(index)
{
  step_length <- 1e-5

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
        - step_length * t(gradient$gradient_initial[[index]])
      ) **2
    ) / sum (
      perturb_res ** 2
    )

  return(ret)
} )

print(rel_rss)
