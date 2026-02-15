library(nimble)
library(coda)

# 1. PREPARA I DATI
nomeSp='Trithemis annulata'

# celle delle specie
df_occu <- dfo |>
  filter(species == nomeSp) |>
  select(date_year, species, cellcode20) |>
  as_tibble() |>
  count(cellcode20, species, date_year, name = 'n_occurrences')

y_matrix <- df_occu %>% 
  complete(cellcode20, date_year, 
           fill = list(n_occurrences = 0)) %>% 
  mutate(detection = ifelse(n_occurrences > 0, 1, 0)) %>%
  pivot_wider(id_cols = cellcode20, 
              names_from = date_year, 
              values_from = detection,
              values_fill = 0) %>%
  column_to_rownames("cellcode20") %>%
  as.matrix()

effort_matrix <- dfo %>%
  filter(cellcode20 %in% df_occu$cellcode20) %>% 
  select(date_year, cellcode20) |>
  count(date_year, cellcode20, name = 'tot_occurrences') %>% 
  pivot_wider(id_cols = cellcode20,
              names_from = date_year,
              values_from = tot_occurrences,
              values_fill = 0) %>%
  column_to_rownames("cellcode20") %>%
  as.matrix()

# Standardizza lo sforzo
effort_std <- scale(effort_matrix)

# 3. Crea covariata temporale (anno)
years <- as.numeric(colnames(y_matrix))
year_std <- scale(years)

# Standardizza anno
year_std <- scale(as.numeric(colnames(y_matrix)))[,1]

# Dati per NIMBLE
nim_data <- list(
  y = y_matrix,
  effort = effort_std,
  year = year_std
)

# Costanti
nim_constants <- list(
  n_sites = nrow(y_matrix),
  n_years = ncol(y_matrix)
)

# 2. DEFINISCI IL MODELLO
code <- nimbleCode({
  
  # PRIORS
  # Occupancy
  alpha_psi ~ dnorm(0, sd = 2)        # intercetta occupancy
  beta_year_psi ~ dnorm(0, sd = 2)    # effetto trend su occupancy
  
  # Detection
  alpha_p ~ dnorm(0, sd = 2)          # intercetta detection
  beta_effort ~ dnorm(0, sd = 2)      # effetto sforzo
  beta_year_p ~ dnorm(0, sd = 2)      # effetto anno su detection
  
  # LIKELIHOOD
  for(i in 1:n_sites) {
    
    # Occupancy per ogni sito (può variare nel tempo)
    for(t in 1:n_years) {
      logit(psi[i,t]) <- alpha_psi + beta_year_psi * year[t]
      z[i,t] ~ dbern(psi[i,t])  # true occupancy state
    }
    
    # Detection dato occupancy
    for(t in 1:n_years) {
      logit(p[i,t]) <- alpha_p + beta_effort * effort[i,t] + 
        beta_year_p * year[t]
      mu[i,t] <- z[i,t] * p[i,t]
      y[i,t] ~ dbern(mu[i,t])
    }
  }
  
  # QUANTITÀ DERIVATE
  # Occupancy media per anno
  for(t in 1:n_years) {
    mean_psi[t] <- ilogit(alpha_psi + beta_year_psi * year[t])
  }
  
  # Detection probability media
  mean_p <- ilogit(alpha_p)
  
})

# 3. INIZIALIZZA
# Valori iniziali per z (true state)
z_init <- ifelse(y_matrix == 1, 1, 0)
# Assicurati che z=1 dove c'è detection
for(i in 1:nrow(y_matrix)) {
  for(j in 1:ncol(y_matrix)) {
    if(is.na(z_init[i,j])) z_init[i,j] <- 0
  }
}

inits <- list(
  alpha_psi = 0,
  beta_year_psi = 0,
  alpha_p = 0,
  beta_effort = 0,
  beta_year_p = 0,
  z = z_init
)

# 4. COMPILA ED ESEGUI
model <- nimbleModel(
  code = code,
  constants = nim_constants,
  data = nim_data,
  inits = inits
)

# Compila
c_model <- compileNimble(model,
                         showCompilerOutput = T)

# Configura MCMC
mcmc_conf <- configureMCMC(model, monitors = c(
  'alpha_psi', 'beta_year_psi', 
  'alpha_p', 'beta_effort', 'beta_year_p',
  'mean_psi', 'mean_p'
))

mcmc <- buildMCMC(mcmc_conf)
c_mcmc <- compileNimble(mcmc, project = model)

# 5. RUN MCMC
set.seed(123)

samples <- runMCMC(
  c_mcmc,
  niter = 10000,
  nburnin = 2000,
  nchains = 3,
  thin = 5,
  samplesAsCodaMCMC = TRUE
)

# 6. DIAGNOSTICA
library(coda)
effectiveSize(samples)
plot(samples[, c('alpha_psi', 'beta_year_psi')])
main_params <- c('alpha_psi', 'beta_year_psi')
autocorr.plot(samples[[1]][, main_params])

# 7. RISULTATI
summary(samples)

# Estrai occupancy per anno
mean_psi_samples <- as.matrix(samples)[, grep("mean_psi", colnames(as.matrix(samples)))]
mean_psi_est <- colMeans(mean_psi_samples)
mean_psi_ci <- apply(mean_psi_samples, 2, quantile, probs = c(0.025, 0.975))

# 8. PLOT TREND
anni <- as.numeric(colnames(y_matrix))

plot(anni, mean_psi_est, type = "l", lwd = 2,
     ylim = c(0, 1), xlab = "Anno", ylab = "Occupancy",
     main = paste("Trend",nomeSp,"(Bayesiano)"))
lines(anni, mean_psi_ci[1,], lty = 2)
lines(anni, mean_psi_ci[2,], lty = 2)
polygon(c(anni, rev(anni)), 
        c(mean_psi_ci[1,], rev(mean_psi_ci[2,])),
        col = rgb(0,0,1,0.2), border = NA)
