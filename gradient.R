# Private function.
# Compute gradient of ODE system solution with respect to
# coefficient matrix, intercept and initial condition, respectively.

# TODO: output type

.gradient <- function (
  linear_ode
  , type = 'gradient'
)

{

#ainv <- solve(linear_ode$coefficient)
#ainv_b <- ainv %*% linear_ode$intercept

ret <- list()

dimension <- nrow(linear_ode$coefficient)
time_point <- linear_ode$observation[,1]

if ( type=='gradient' )
{
  gradient_initial <-
    lapply ( 1:dimension , function(index)
    {
      ret <-
        matrix (
          0
          , dimension
          , length(time_point)
        )
      return(ret)
    } )
}

require('expm')

lapply ( 1:length(time_point) , function(index)
{
  exp_at <-
    expm::expm (
      linear_ode$coefficient * time_point[index]
    )
  lapply ( 1:dimension , function(index2)
  {
    gradient_initial[[index2]][,index] <<- exp_at[,index2]
    return()
  } )
  return()
} )

ret$gradient_initial <- gradient_initial

# Return#{{{

  return(ret)
#}}}

}
