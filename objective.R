# Objective function which gives the residual sum-of-squares
# for a given spectrum.
# Used as objective function for LM,NEWUOA,BOBYQA algorithms.
# Also used to reconstruct an ODE system after the optimal spectrum
# is found.

.objective <- function (
  par
  , num_complex_eigen
  , num_real_eigen
  , observation
  , intercept = TRUE
  , type = 'post'
  , ridge_coefficient = 1e-3
)

# INPUT:
# par: Parameters.
#   Check following item for details.
# num_complex_eigen, num_real_eigen:
#   First "num_complex_eigen" values of "par" are
#   imaginary parts of complex eigen-values.
#   The following "num_complex_eigen" values are corresponding
#   real parts of those complex eigen-values.
#   Notice that it gives "num_complex_eigen" pair of conjugate
#   complex eigen-values.
#   The last "num_real_eigen" values of "par" are real eigen-values.
# observation: Observation matrix.
#   Each row stands for time-point.
#   Each column is a curve.
#   First column is time points.
# intercept: Whether intercept term is contained in the system.
# type: Decide return value.
#   "dfo": For NEWUOA,BOBYQA solver.
#   "lm": For Levenberg-Marquardt solver.
#   "post": Post process.
#     The input "par" is treated as optimal eigen-values.
#     Return corresponding curves, coefficient matrix,
#     and intercept term(if exists).
# ridge_coefficient: Coefficient for Ridge regression.

{

# Generate basis#{{{

dimension <- length(par)
time_point <- observation[,1]

basis <- matrix ( 0 , length(time_point) , dimension )

lapply ( 1 : num_complex_eigen , function(index)
{
  temp <- exp ( par[index+num_complex_eigen] * time_point )
  basis [,2*index-1] <<- temp * sin ( par[index] * time_point )
  basis [,2*index ] <<- temp * cos ( par[index] * time_point )
} )

if ( num_real_eigen > 0 )
{
  lapply ( (dimension-num_real_eigen+1) : dimension , function(index)
  {
    basis [,index] <<- exp ( par[index] * time_point )
  } )
}

if ( intercept )
{
  basis <- cbind ( basis , rep(1,nrow(basis)) )
}
#}}}

# Return value for type "dfo","lm" #{{{

ridge_basis <- rbind (
  basis
  , ridge_coefficient * diag ( ncol(basis) )
)

q_basis <- qr.Q(qr(ridge_basis))

if ( type == 'lm' || type == 'dfo' )
{
  residue <- sapply ( 1:dimension , function(index)
  {
    ret <-
      q_basis %*%
      t(q_basis) %*%
      c(observation[,index+1],numeric(ncol(basis)))
    ret <- observation[,index+1] - ret[1:length(time_point)]
    ret <- sum(ret**2)
    return(ret)
  } )

  if ( type == 'dfo' )
  {
    return(sum(residue))
  }
  else
  {
    return(sqrt(residue))
  }
}
#}}}

# Return value for type "post"#{{{

if ( type == 'post' )
{
  curve_est <- sapply ( 1:dimension , function(index)
  {
    ret <-
      q_basis %*%
      t(q_basis) %*%
      c(observation[,index+1],numeric(ncol(basis)))
    ret <- ret[1:length(time_point)]
    return(ret)
  } )
  curve_est <- cbind ( time_point , curve_est )

  if ( intercept )
  {
    intercept_est <- numeric(dimension)
  }
  else
  {
    intercept_est <- NULL
  }

  q_optim <- sapply ( 1:dimension , function(index)
  {
    bt_o <- t(ridge_basis) %*% c(observation[,index+1],numeric(ncol(basis)))
    ret <- solve (
      t(ridge_basis)%*%ridge_basis + ridge_coefficient * diag(ncol(basis))
      , bt_o
    )
    if ( intercept )
    {
      intercept_est[index] <<- tail(ret,1)
      ret <- ret[1:dimension]
    }
    return(ret)
  } )

  coefficient_est <- matrix ( 0 , dimension , dimension )

  index <- (1:num_complex_eigen) * 2
  coefficient_est[cbind(index,index)] <-
    par[(num_complex_eigen+1):(2*num_complex_eigen)]
  coefficient_est[cbind(index-1,index-1)] <-
    par[(num_complex_eigen+1):(2*num_complex_eigen)]
  coefficient_est[cbind(index,index-1)] <-
    par[1:num_complex_eigen]
  coefficient_est[cbind(index-1,index)] <-
    -par[1:num_complex_eigen]

  if ( num_real_eigen > 0 )
  {
    index <- (dimension-num_real_eigen+1):dimension
    coefficient_est[cbind(index,index)] <- tail(par,num_real_eigen)
  }

  coefficient_est <- solve ( q_optim , coefficient_est )
  coefficient_est <- coefficient_est %*% q_optim
  coefficient_est <- t(coefficient_est)

  if ( intercept )
  {
    intercept_est <- - coefficient_est %*% intercept_est
    intercept_est <- as.numeric(intercept_est)
  }

  return ( list (
    curve_est = curve_est
    , coefficient_est = coefficient_est
    , intercept_est = intercept_est
  ) )
}
#}}}

}
