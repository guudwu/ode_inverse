dimension <- 9
time_point <- 0:100
ITER <- 1
seed <- 0

print('Generate simulation data:')
repeat
{
  print(seed)
  set.seed(seed)
  source('linear_ode_generation/linear_ode_generation.R')
  linear_ode <- linear_ode_generation (
    dimension
    , time_point
    , row_column_permutation = FALSE
    , intercept = c(0,1e-2)
  )
  source('linear_ode_generation/check.linear_ode.R')
  ode_property <- check.linear_ode(linear_ode)
  if ( ode_property$max_correlation<0.3 && ode_property$vif<5 )
  {
    break
  }
  seed <- seed + ITER
}
print('Start of iteration:')
print(seed)

noise_level <- 0
#noise_level <- 0.3
noise <- matrix (
  rnorm (
    length(time_point) * dimension
    , 0
    , noise_level * sd(linear_ode$observation[,-1])
  )
  , length(time_point)
  , dimension
)
observation_noise <- linear_ode$observation
observation_noise[,-1] <- observation_noise[,-1] + noise

initial <- numeric(dimension)

source('separable_ls.R')
res <- separable_ls (
  observation_noise
  , upper_bound_complex = 1+1i
  , lower_bound_complex = -1-1i
#  , method = 'dfo'
  , lm.parameter = nls.lm.control(maxiter=1024,maxfev=1e8)
  , verbose = 1
)
