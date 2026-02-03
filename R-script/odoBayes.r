# funzione odoBayes
# per il calcolo dell'occupancy
# richiede la definizione di un campo cellcodeX

odoBayes <- function(nomeSp, nk = 30) {
  # matrice di rilevamento (y)
  y_matrix <- dfo %>%
    filter(species == nomeSp) %>%
    count(cellcodeX, date_year, name = 'n_occurrences') %>%
    # Join con effort
    right_join(
      dfo %>%
        count(cellcodeX, date_year, name = 'tot_occurrences'),
      by = c("cellcodeX", "date_year")
    ) %>%
    complete(
      cellcodeX = all_cells,
      date_year = all_years,
      fill = list(n_occurrences = NA, tot_occurrences = NA)
    ) %>%
    mutate(
      detection = case_when(
        is.na(tot_occurrences) ~ NA_real_, # Cella non indagata
        n_occurrences > 0 ~ 1, # Specie trovata
        tot_occurrences > 0 ~ 0, # Indagata ma specie non trovata
        TRUE ~ NA_real_ # Cella non indagata
      )
    ) %>%
    pivot_wider(
      id_cols = cellcodeX,
      names_from = date_year,
      values_from = detection
    ) %>%
    column_to_rownames("cellcodeX") %>%
    as.matrix()

  # Covariata Anno (vettore atomico per NIMBLE)
  anni_num <- as.numeric(colnames(y_matrix))
  year_std <- as.vector(scale(anni_num))

  # Standardizzazione (vettoriale per evitare attributi scale)
  effort_vec <- as.numeric(effort_matrix)
  effort_std_mat <- matrix(
    (effort_vec - mean(effort_vec, na.rm = TRUE)) /
      sd(effort_vec, na.rm = TRUE),
    nrow = nrow(y_matrix),
    ncol = ncol(y_matrix)
  )

  cells_in_y <- rownames(y_matrix)
  tin_data <- dfo %>%
    select(cellcodeX, tin_1) %>%
    distinct() %>%
    mutate(cellcodeX = as.character(cellcodeX)) %>%
    filter(cellcodeX %in% cells_in_y) %>%
    group_by(cellcodeX) %>%
    summarise(tin_1 = mean(round(tin_1), na.rm = TRUE))

  # Crea un dataframe con tutte le celle nell'ordine corretto
  tin_df <- data.frame(cellcodeX = cells_in_y) %>%
    left_join(tin_data, by = "cellcodeX")

  # Standardizzazione di tin_1
  tin_mean <- mean(tin_df$tin_1, na.rm = TRUE)
  tin_sd <- sd(tin_df$tin_1, na.rm = TRUE)
  tin_std <- as.vector((tin_df$tin_1 - tin_mean) / tin_sd)

  # Dati presenza e effort
  nim_data <- list(
    y = y_matrix,
    effort = effort_std_mat,
    year = year_std,
    tin_1 = tin_std
  )

  nim_constants <- list(
    n_sites = nrow(y_matrix),
    n_years = ncol(y_matrix)
  )

  # 2. DEFINIZIONE MODELLO
  code <- nimbleCode({
    # PRIORS
    alpha_psi ~ dnorm(0, sd = 2)
    beta_year_psi ~ dnorm(0, sd = 2)
    beta_tin_psi ~ dnorm(0, sd = 2)
    alpha_p ~ dnorm(0, sd = 2)
    beta_effort ~ dnorm(0, sd = 2)
    beta_year_p ~ dnorm(0, sd = 2)

    # LIKELIHOOD
    for (i in 1:n_sites) {
      for (t in 1:n_years) {
        # Occupancy
        logit(psi[i, t]) <- alpha_psi +
          beta_year_psi * year[t] +
          beta_tin_psi * tin_1[i]
        z[i, t] ~ dbern(psi[i, t])

        # Detection
        logit(p[i, t]) <- alpha_p +
          beta_effort * effort[i, t] +
          beta_year_p * year[t]
        mu[i, t] <- z[i, t] * p[i, t]
        y[i, t] ~ dbern(mu[i, t])
      }
    }

    # DERIVED QUANTITIES
    # Media occupancy per anno (calcolata al valore medio di tin_1 = 0)
    for (t in 1:n_years) {
      mean_psi[t] <- ilogit(alpha_psi + beta_year_psi * year[t])
    }
    mean_p <- ilogit(alpha_p)
  })

  # 3. INIZIALIZZAZIONE
  # Importante: z deve essere 1 se y è 1
  z_init <- y_matrix
  z_init[is.na(z_init)] <- 0

  inits <- list(
    alpha_psi = 0,
    beta_year_psi = 0,
    beta_tin_psi = 0,
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
      'alpha_psi',
      'beta_year_psi',
      'beta_tin_psi',
      'alpha_p',
      'beta_effort',
      'mean_psi'
    )
  )
  mcmc <- buildMCMC(mcmc_conf)
  c_mcmc <- compileNimble(mcmc, project = model, showCompilerOutput = F)

  # 5. RUN MCMC
  set.seed(32)
  samples <- runMCMC(
    c_mcmc,
    niter = 20000,
    nburnin = 5000,
    nchains = 3,
    thin = 10,
    samplesAsCodaMCMC = TRUE
  )

  # 6. CALCOLO DIAGNOSTICHE (R-hat)
  rhat_values <- gelman.diag(samples, multivariate = FALSE)
  rhat_df <- as.data.frame(rhat_values$psrf)

  # 7. RISULTATI E SALVATAGGIO
  risultati_finali <- list(
    specie = nomeSp,
    summary = summary(samples),
    rhat = rhat_df,
    samples = samples
  )

  # Salva in formato RDS (il nome del file include il nome della specie)
  file_rds <- paste0("output_bayes/risultati", nk, "k_", nomeSp, ".rds")
  saveRDS(risultati_finali, file = file_rds)

  message(paste("Statistiche salvate in:", file_rds))

  # 8. PLOT TREND
  anni_num <- as.numeric(colnames(y_matrix))
  mean_psi_cols <- grep("mean_psi", colnames(as.matrix(samples)))
  mean_psi_samples <- as.matrix(samples)[, mean_psi_cols]
  mean_rhat <- round(mean(rhat_df$`Point est.`), 2)
  est <- colMeans(mean_psi_samples)
  ci <- apply(mean_psi_samples, 2, quantile, probs = c(0.025, 0.975))

  ## data.frame
  df_plot <- data.frame(
    Anno = anni_num,
    Occupancy = est,
    Lower = ci["2.5%", ],
    Upper = ci["97.5%", ]
  )

  ## plot
  pOccu <- df_plot |>
    ggplot(aes(x = Anno, y = Occupancy)) +
    geom_ribbon(
      aes(ymin = Lower, ymax = Upper),
      fill = "forestGreen",
      alpha = 0.2
    ) +
    geom_line(color = "forestGreen", linewidth = 1) +
    geom_point(color = "forestGreen", size = 3) +
    geom_smooth(
      method = "loess",
      color = "red",
      linetype = "dashed",
      se = F,
      alpha = 0.5
    ) +
    scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
    labs(
      title = paste(nomeSp),
      x = "",
      y = "Occupancy",
      tag = "nimble"
    ) +
    theme_bw(base_size = 14) +
    theme(plot.title = element_text(face = "italic"))

  print(pOccu)

  beta_tin_samples <- as.matrix(samples)[, "beta_tin_psi"]
  cat("Media posteriore:", mean(beta_tin_samples), "\n")
  cat("IC 95%:", quantile(beta_tin_samples, c(0.025, 0.975)), "\n")

  return(risultati_finali)
}
