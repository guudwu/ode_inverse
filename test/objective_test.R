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
#    , intercept = intercept
  )

source('objective.R')
rss <- .objective (
  c ( linear_ode$eigen_imaginary , linear_ode$eigen_real )
  , num_complex_eigen = floor(dimension/2)
  , num_real_eigen = dimension%%2
  , observation = linear_ode$observation
  , intercept = !is.null(linear_ode$intercept)
  , type = 'dfo'
)

noise <- matrix (
  rnorm (
    dimension * length(time_point)
    , 0
    , 0.1 * sd(linear_ode$observation[,-1])
  )
  , length(time_point)
  , dimension
)
observation_noise <- linear_ode$observation
observation_noise[,-1] <- observation_noise[,-1] + noise

source('objective.R')
rss_noise <- .objective (
  c ( linear_ode$eigen_imaginary , linear_ode$eigen_real )
  , num_complex_eigen = floor(dimension/2)
  , num_real_eigen = dimension%%2
  , observation = observation_noise
  , intercept = !is.null(linear_ode$intercept)
  , type = 'dfo'
)

print(rss)
print(rss_noise)
print(sum(noise**2))
