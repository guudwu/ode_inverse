# Function to compute the variance of a linear ODE system parameter
# estimation.

variance_linear_ode <- function (
  estimation
  , observation
  , method = 'fim'
)

{

dimension <- nrow(estimation$coefficient_est)
intercept <- !is.null(estimation$intercept_est)

# Method "fim"#{{{

source('gradient.R')
if ( method=='fim' )
{
  ret <- list()

  dof <- (
    length(observation) - nrow(observation)
    - dimension * ( dimension + 1 )
  )
  if ( intercept )
  {
    dof <- dof - dimension
  }

  if ( dof<=0 )
  {
    stop('Degree of Freedom is not positive.')
  }

  fim <-
    .gradient (
      list (
        coefficient = estimation$coefficient_est
        , initial = estimation$curve_est[1,-1]
        , intercept = estimation$intercept_est
        , time_point = observation[,1]
      )
      , type = 'fim'
    )

  variance <- diag(solve(fim))
  variance <- variance * estimation$rss / dof

  ret$coefficient <-
    matrix (
      variance[1:(dimension**2)]
      , dimension
      , dimension
    )
  ret$initial <- variance[(dimension**2+1):(dimension**2+dimension)]
  if ( intercept )
  {
    ret$intercept <- tail(variance,dimension)
  }
}
#}}}

# Return#{{{

return(ret)
#}}}

}
