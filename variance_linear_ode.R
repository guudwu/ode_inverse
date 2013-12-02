# Function to compute the variance of a linear ODE system parameter
# estimation.

variance_linear_ode <- function (
  estimation
  , observation
  , method = 'fim'
  , latex = NULL
)

{

dimension <- nrow(estimation$coefficient)
intercept <- !is.null(estimation$intercept)

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
        coefficient = estimation$coefficient
        , initial = estimation$curve[1,-1]
        , intercept = estimation$intercept
        , time_point = observation[,1]
      )
      , type = 'fim'
    )

# Numerical adjustment for singular FIM
  fim_svd <- svd(fim)$d
  if ( fim_svd[1]/tail(fim_svd,1)>1e15 )
  {
    warning('FIM computationally singular.')
    fim <- (
      fim
      + 1e-15 * fim_svd[1] * diag(nrow(fim))
    )
  }

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

# Latex output#{{{

if ( !is.null(latex) )
{
  sink (
    latex
    , append = TRUE
  )

  cat('\\begin{equation}\n')
  cat('\\begin{cases}\n')
  cat('x\'(t)=\\begin{bmatrix}\n')
  for ( index_row in 1:dimension )
  {
    for ( index_col in 1:dimension )
    {
      cat(
        sprintf("%.1e",estimation$coefficient[index_row,index_col])
      )
      cat('(')
      cat(
        sprintf("%.1e",sqrt(ret$coefficient[index_row,index_col]))
      )
      cat(')')
      if ( index_col<dimension )
      {
        cat('&')
      }
    }
    if ( index_row<dimension )
    {
      cat('\\\\')
    }
  }
  cat('\n\\end{bmatrix}x(t)\n')
  if ( intercept )
  {
    cat('+\\begin{bmatrix}\n')
    for ( index in 1:dimension )
    {
      cat(
        sprintf("%.1e",estimation$intercept[index])
      )
      cat('(')
      cat(
        sprintf("%.1e",sqrt(ret$intercept[index]))
      )
      cat(')')
      if ( index<dimension )
      {
        cat('\\\\')
      }
    }
    cat('\n\\end{bmatrix}\n')
  }
  cat('&\\\\\nx(0)=')
  cat('\\begin{bmatrix}\n')
  for ( index in 1:dimension )
  {
    cat(
      sprintf("%.1e",estimation$curve[1,index+1])
    )
    cat('(')
    cat(
      sprintf("%.1e",sqrt(ret$initial[index]))
    )
    cat(')')
    if ( index<dimension )
    {
      cat('\\\\')
    }
  }
  cat('\n\\end{bmatrix}&\n')
  cat('\\end{cases}\n')
  cat('\\end{equation}\n')

  sink(NULL)
}
#}}}

# Return#{{{

return(ret)
#}}}

}
