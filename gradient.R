# Private function.
# Compute gradient of ODE system solution with respect to
# coefficient matrix, intercept and initial condition, respectively.

.gradient <- function (
  linear_ode
)

{

# Initialization#{{{

ret <- list()

dimension <- nrow(linear_ode$coefficient)
time_point <- linear_ode$observation[,1]
intercept <- !is.null(linear_ode$intercept)

gradient_coefficient <-
  lapply ( 1:(dimension**2) , function(index)
  {
    ret <-
      matrix (
        0
        , dimension
        , length(time_point)
      )
    return(ret)
  } )

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
#}}}

# Gradient#{{{

if ( intercept )
{
  ainv_b <-
    solve (
      linear_ode$coefficient
      , linear_ode$intercept
    )
  ainv <-
    solve ( linear_ode$coefficient )
}

require('expm')

lapply ( 1:length(time_point) , function(index)
{
  exp_at <-
    expm::expm (
      linear_ode$coefficient * time_point[index]
    )

# Coefficient gradient#{{{

  lapply ( 0:(dimension**2-1) , function(index2)
  {
    num_row <- floor(index2/dimension) + 1
    num_col <- index%%dimension + 1

    direction <- matrix ( 0 , dimension , dimension )
    direction[num_row,num_col] <- 1

    frechet <-
      expm::expmFrechet (
        linear_ode$coefficient * time_point[index]
        , direction * time_point[index]
        , expm = FALSE
      ) $ Lexpm

    if ( intercept )
    {
      partial_ainv_b <-
        solve (
          linear_ode$coefficient
          , direction %*% ainv_b
        )
    }

    if ( intercept )
    {
      gradient_coefficient[[index2+1]][,index] <<- (
        exp_at %*% partial_ainv_b
        + frechet %*% ( linear_ode$initial + ainv_b )
        - partial_ainv_b
      )
    }
    else
    {
      gradient_coefficient[[index2+1]][,index] <<-
        frechet %*% linear_ode$initial
    }

    return()
  } )
#}}}

# Initial condition gradient#{{{

  lapply ( 1:dimension , function(index2)
  {
    gradient_initial[[index2]][,index] <<- exp_at[,index2]
    return()
  } )
#}}}

  return()
} )
#}}}

# Return#{{{

ret$initial <- gradient_initial
ret$coefficient <- gradient_coefficient

return(ret)
#}}}

}
