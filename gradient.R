# Private function.
# Compute gradient of ODE system solution with respect to
# coefficient matrix, intercept and initial condition, respectively.

.gradient <- function (
  linear_ode
  , type = 'gradient'
)

# INPUT:
# linear_ode: A list containing following parameters:
#   coefficient: Estimated coefficient matrix.
#   initial: Estimated initial condition.
#   intercept: Estimated intercept term(if existed).
#   observation: Data to be fit.
#     Each row is for one time point.
#     Each column is a curve.
#     First column is time points.
# type: Determine return value.
#   Can be following value:
#   "gradient": partial derivative matrices.
#   "fim": Fisher Information Matrix.

# OUTPUT for type "gradient":
# coefficient: Derivative matrix for each entry of the coefficient
#   matrix.
#   The order is column-oriented.
# initial: Derivative matrix for each entry of the initial condition
#   vector.
# intercept: Only exists when input "linear_ode" contains an
#   "intercept" item.
#   Derivative matrix for each entry of the intercept term.

# OUTPUT for type "fim":
# Fisher information matrix.
# The order is: the first "dimension**2" elements is for coefficient
# matrix entries, and the following "dimension" elements for initial
# condition vector, and the rest(if existed) for intercept term.

{

# Initialization#{{{

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

if ( intercept )
{
  gradient_intercept <-
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

  lapply ( 1:(dimension**2) , function(index2)
  {
    direction <- matrix ( 0 , dimension , dimension )
    direction[index2] <- 1

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
      gradient_coefficient[[index2]][,index] <<- (
        exp_at %*% partial_ainv_b
        + frechet %*% ( linear_ode$initial + ainv_b )
        - partial_ainv_b
      )
    }
    else
    {
      gradient_coefficient[[index2]][,index] <<-
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

# intercept gradient#{{{

  if ( intercept )
  {
    temp <-
      ( exp_at - diag(dimension) ) %*% ainv
    lapply ( 1:dimension , function(index2)
    {
      gradient_intercept[[index2]][,index] <<- temp[,index2]
      return()
    } )
  }
#}}}

  return()
} )
#}}}

# Return for type "gradient"#{{{

if ( type=='gradient' )
{
  ret <- list()

  ret$initial <- gradient_initial
  ret$coefficient <- gradient_coefficient
  if ( intercept )
  {
    ret$intercept <- gradient_intercept
  }

  return(ret)
}
#}}}

# Return for type "fim"#{{{

if ( type=='fim' )
{
  fim_dim <- dimension * (dimension+1)
  if ( intercept )
  {
    fim_dim <- fim_dim + dimension
  }

  ret <- matrix ( 0 , fim_dim , fim_dim )

  lapply ( 1:fim_dim , function(index)
  {
    if ( index<=dimension**2 )
    {
      temp1 <- gradient_coefficient[[index]]
    }
    else if ( index<=dimension*(dimension+1) )
    {
      temp1 <- gradient_initial[[index-dimension**2]]
    }
    else
    {
      temp1 <- gradient_intercept[[index-dimension*(dimension+1)]]
    }
    lapply ( 1:index , function(index2)
    {
      if ( index2<=dimension**2 )
      {
        temp2 <- gradient_coefficient[[index2]]
      }
      else if ( index2<=dimension*(dimension+1) )
      {
        temp2 <- gradient_initial[[index2-dimension**2]]
      }
      else
      {
        temp2 <- gradient_intercept[[index2-dimension*(dimension+1)]]
      }
      ret[index,index2] <<-
        sum (
          temp1 * temp2
        )
      ret[index2,index] <<- ret[index,index2]
      return()
    } )
    return()
  } )

  return(ret)
}
#}}}

}
