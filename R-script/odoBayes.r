# funzione odoBayes
# per il calcolo dell'occupancy
# richiede la definizione di un campo cellcodeX
odoBayes <- function(nomeSp, nk = 10) {
  cat(nomeSp, '\n')

  # 1. PREPARAZIONE MATRICE Y
  y_matrix <- dfo %>%
    filter(species == nomeSp) %>%
    count(cellcodeX, date_year, name = 'n_occurrences') %>%
    right_join(
      dfo %>% count(cellcodeX, date_year, name = 'tot_occurrences'),
      by = c("cellcodeX", "date_year")
    ) %>%
    complete(
      cellcodeX = all_cells,
      date_year = all_years,
      fill = list(n_occurrences = NA, tot_occurrences = NA)
    ) %>%
    mutate(
      detection = case_when(
        is.na(tot_occurrences) ~ NA_real_,
        n_occurrences > 0 ~ 1,
        tot_occurrences > 0 ~ 0
      )
    ) %>%
    pivot_wider(
      id_cols = cellcodeX,
      names_from = date_year,
      values_from = detection
    ) %>%
    column_to_rownames("cellcodeX") %>%
    as.matrix()

  # Standardizzazione effort
  effort_vec <- as.numeric(effort_matrix)
  effort_std_mat <- matrix(
    (effort_vec - mean(effort_vec, na.rm = TRUE)) /
      sd(effort_vec, na.rm = TRUE),
    nrow = nrow(y_matrix),
    ncol = ncol(y_matrix)
  )

  anni_num <- as.numeric(colnames(y_matrix))
  year_std <- as.vector(scale(anni_num))

  # Covariata sito (tin_1)
  cells_in_y <- rownames(y_matrix)
  tin_data <- dfo %>%
    select(cellcodeX, tin_1) %>%
    distinct() %>%
    mutate(cellcodeX = as.character(cellcodeX)) %>%
    filter(cellcodeX %in% cells_in_y) %>%
    group_by(cellcodeX) %>%
    summarise(tin_1 = mean(round(tin_1), na.rm = TRUE))

  tin_df <- data.frame(cellcodeX = cells_in_y) %>%
    left_join(tin_data, by = "cellcodeX")

  tin_std <- as.vector(scale(tin_df$tin_1))

  # Preparazione dati osservati (long format per NIMBLE)
  obs_data <- data.frame(
    site_idx = integer(),
    year_idx = integer(),
    y_value = numeric(),
    effort_value = numeric(),
    year_value = numeric()
  )

  for (i in 1:nrow(y_matrix)) {
    for (t in 1:ncol(y_matrix)) {
      if (!is.na(y_matrix[i, t])) {
        obs_data <- rbind(
          obs_data,
          data.frame(
            site_idx = i,
            year_idx = t,
            y_value = y_matrix[i, t],
            effort_value = effort_std_mat[i, t],
            year_value = year_std[t]
          )
        )
      }
    }
  }

  # 2. DEFINIZIONE MODELLO DINAMICO
  code <- nimbleCode({
    # --- PRIORS ---
    alpha_psi1 ~ dnorm(0, sd = 2)
    beta_tin_psi1 ~ dnorm(0, sd = 2)
    alpha_gamma ~ dnorm(0, sd = 2)
    beta_year_gamma ~ dnorm(0, sd = 2)
    alpha_eps ~ dnorm(0, sd = 2)
    beta_year_eps ~ dnorm(0, sd = 2)
    alpha_p ~ dnorm(0, sd = 2)
    beta_effort ~ dnorm(0, sd = 2)
    beta_year_p ~ dnorm(0, sd = 2)

    # --- LIKELIHOOD ---
    for (i in 1:n_sites) {
      # Anno 1
      logit(psi1[i]) <- alpha_psi1 + beta_tin_psi1 * tin_1[i]
      z[i, 1] ~ dbern(psi1[i])
      # Anni successivi
      for (t in 2:n_years) {
        logit(gamma[i, t - 1]) <- alpha_gamma + beta_year_gamma * year[t]
        logit(eps[i, t - 1]) <- alpha_eps + beta_year_eps * year[t]
        prob_occ[i, t] <- z[i, t - 1] *
          (1 - eps[i, t - 1]) +
          (1 - z[i, t - 1]) * gamma[i, t - 1]
        z[i, t] ~ dbern(prob_occ[i, t])
      }
    }
    # Rilevamento
    for (k in 1:n_obs) {
      logit(p[k]) <- alpha_p +
        beta_effort * effort_obs[k] +
        beta_year_p * year_value_obs[k]
      mu[k] <- z[site_obs[k], year_idx_obs[k]] * p[k]
      y_obs[k] ~ dbern(mu[k])
    }
    # Quantità derivate
    for (t in 1:n_years) {
      mean_psi[t] <- sum(z[1:n_sites, t]) / n_sites
    }
  })

  # 3. INIZIALIZZAZIONE
  # Cruciale: z non deve avere NA
  prob_osservata <- mean(y_matrix, na.rm = TRUE)

  z_init <- y_matrix
  z_init[is.na(z_init)] <- rbinom(sum(is.na(y_matrix)), 1, prob_osservata)
  z_init[y_matrix == 1] <- 1
  z_init[y_matrix == 0] <- 0

  nim_constants <- list(
    n_sites = nrow(y_matrix),
    n_years = ncol(y_matrix),
    n_obs = nrow(obs_data),
    site_obs = obs_data$site_idx,
    year_idx_obs = obs_data$year_idx,
    tin_1 = tin_std,
    year = year_std
  )

  nim_data <- list(
    y_obs = obs_data$y_value,
    effort_obs = obs_data$effort_value,
    year_value_obs = obs_data$year_value
  )

  inits <- list(
    alpha_psi1 = 0,
    beta_tin_psi1 = 0,
    alpha_gamma = 0,
    beta_year_gamma = 0,
    alpha_eps = 0,
    beta_year_eps = 0,
    alpha_p = 0,
    beta_effort = 0,
    beta_year_p = 0,
    z = z_init
  )

  # 4. ESECUZIONE
  model <- nimbleModel(
    code,
    constants = nim_constants,
    data = nim_data,
    inits = inits
  )
  c_model <- compileNimble(model)

  mcmc_conf <- configureMCMC(
    model,
    monitors = c(
      'alpha_psi1',
      'alpha_gamma',
      'alpha_eps',
      'alpha_p',
      'mean_psi'
    )
  )
  mcmc <- buildMCMC(mcmc_conf)
  c_mcmc <- compileNimble(mcmc, project = model)

  samples <- runMCMC(
    c_mcmc,
    niter = 20000,
    nburnin = 5000,
    nchains = 3,
    thin = 10,
    samplesAsCodaMCMC = T,
    progressBar = F
  )

  # 5. DIAGNOSTICHE E RISULTATI
  rhat_df <- as.data.frame(gelman.diag(samples, multivariate = FALSE)$psrf)

  mean_psi_cols <- grep("mean_psi", colnames(as.matrix(samples)))
  mean_psi_samples <- as.matrix(samples)[, mean_psi_cols]

  df_plot <- data.frame(
    Anno = anni_num,
    Occupancy = colMeans(mean_psi_samples),
    Lower = apply(mean_psi_samples, 2, quantile, probs = 0.025),
    Upper = apply(mean_psi_samples, 2, quantile, probs = 0.975)
  )

  # 6. PLOT E SALVATAGGIO
  pOccu <- ggplot(df_plot, aes(x = Anno, y = Occupancy)) +
    geom_ribbon(
      aes(ymin = Lower, ymax = Upper),
      fill = "forestGreen",
      alpha = 0.2
    ) +
    geom_line(color = "forestGreen", linewidth = 1) +
    geom_point(color = "forestGreen", size = 3) +
    scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
    labs(
      title = nomeSp,
      subtitle = paste("R-hat:", round(mean(rhat_df[, 1]), 2)),
      x = '',
      y = "Occupancy",
      caption = "nimble (dynamic)"
    ) +
    theme_bw(base_size = 14) +
    theme(plot.title = element_text(face = "italic"))

  print(pOccu)

  ggsave(
    paste0('output_bayes/', nomeSp, '_', today(), '.pdf'),
    width = 9,
    height = 6
  )

  saveRDS(
    list(samples = samples, plot_data = df_plot, rhat = rhat_df),
    paste0("output_bayes/dyn_", nomeSp, "_", today(), ".rds")
  )

  return(range(df_plot$Occupancy))
}
