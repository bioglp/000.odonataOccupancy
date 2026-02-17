# Load libraries
library(rjags)

# Example data
set.seed(123)
x <- rnorm(100)
y <- 3 + 2 * x + rnorm(100, 0, 1)

# JAGS model
model_string <- "
  model {
    for (i in 1:N) {
      y[i] ~ dnorm(mu[i], tau)
      mu[i] <- alpha + beta * x[i]
    }
    alpha ~ dnorm(0, 0.001)
    beta ~ dnorm(0, 0.001)
    tau ~ dgamma(0.01, 0.01)
    sigma <- 1 / sqrt(tau)
  }
"

# Data for JAGS
data_jags <- list(
  N = length(y),
  y = y,
  x = x
)

# Initial values
inits <- list(
  list(alpha = 0, beta = 0, tau = 1),
  list(alpha = 1, beta = 1, tau = 2)
)

# Compile model
model <- jags.model(textConnection(model_string), data = data_jags, inits = inits, n.chains = 2)

# Burn-in phase
update(model, 1000)

# Sample from posterior
samples <- coda.samples(model, variable.names = c("alpha", "beta", "tau", "sigma"), n.iter = 5000)

# Summary and plot
summary(samples)
plot(samples)

# Convergence diagnostics
gelman.diag(samples)
