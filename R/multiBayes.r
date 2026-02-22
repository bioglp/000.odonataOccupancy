#' funzione per il calcolo multispecie
#'
#'

odoMultiSpecies_AllGroups <- function(species_index_df, nk = 10, ncores = 7) {
  cat("Preparazione dati multispecie per TUTTI i gruppi...\n")

  species_subset <- species_index_df %>%
    select(species, gindex) %>%
    distinct()

  all_species <- unique(species_subset$species)
  n_species <- length(all_species)

  cat("Numero totale di specie:", n_species, "\n")
  cat("Distribuzione per gruppo:\n")
  print(table(species_subset$gindex))

  # Creare matrice Y per ogni specie
  y_list <- list()
  tin_std_list <- list()

  for (sp in all_species) {
    cat("Preparando dati per", sp, "\n")

    y_matrix <- dfo %>%
      filter(species == sp) %>%
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

    y_list[[sp]] <- y_matrix

    # Covariata sito (tin_1)
    cells_in_y <- rownames(y_matrix)
    tin_data <- dfo %>%
      select(cellcodeX, tin_1) %>%
      distinct() %>%
      mutate(cellcodeX = as.character(cellcodeX)) %>%
      filter(cellcodeX %in% cells_in_y) %>%
      group_by(cellcodeX) %>%
      summarise(tin_1 = mean(round(tin_1), na.rm = TRUE), .groups = 'drop')

    tin_df <- data.frame(cellcodeX = cells_in_y) %>%
      left_join(tin_data, by = "cellcodeX")

    tin_std_list[[sp]] <- as.vector(scale(tin_df$tin_1))
  }

  # Standardizzazione effort (uguale per tutte le specie)
  effort_vec <- as.numeric(effort_matrix)
  effort_std_mat <- matrix(
    (effort_vec - mean(effort_vec, na.rm = TRUE)) /
      sd(effort_vec, na.rm = TRUE),
    nrow = nrow(y_list[[1]]),
    ncol = ncol(y_list[[1]])
  )

  anni_num <- as.numeric(colnames(y_list[[1]]))
  year_std <- as.vector(scale(anni_num))

  # Preparazione dati osservati per tutte le specie
  obs_data <- data.frame(
    species_idx = integer(),
    site_idx = integer(),
    year_idx = integer(),
    y_value = numeric(),
    effort_value = numeric(),
    year_value = numeric()
  )

  for (s in 1:n_species) {
    sp <- all_species[s]
    y_matrix <- y_list[[sp]]

    for (i in 1:nrow(y_matrix)) {
      for (t in 1:ncol(y_matrix)) {
        if (!is.na(y_matrix[i, t])) {
          obs_data <- rbind(
            obs_data,
            data.frame(
              species_idx = s,
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
  }

  cat("N. osservazioni totali:", nrow(obs_data), "\n")

  # Mappare gindex per ogni specie
  gindex_vec <- species_subset$gindex[match(
    all_species,
    species_subset$species
  )]

  cat("Verifica gindex:\n")
  print(table(gindex_vec))

  # Numero gruppi
  n_groups <- length(unique(gindex_vec))
  cat("Numero gruppi:", n_groups, "\n")

  if (!all(sort(unique(gindex_vec)) == 1:n_groups)) {
    stop(
      "gindex deve contenere valori contigui da 1 a ",
      n_groups,
      ". Valori trovati: ",
      paste(sort(unique(gindex_vec)), collapse = ", ")
    )
  }

  # 2. DEFINIZIONE MODELLO
  code <- nimbleCode({
    # --- HYPERPRIORS PER OGNI GRUPPO ---
    for (tg in 1:n_groups) {
      mu_alpha_psi1[tg] ~ dnorm(0, sd = 2)
      mu_alpha_gamma[tg] ~ dnorm(0, sd = 2)
      mu_alpha_eps[tg] ~ dnorm(0, sd = 2)
      mu_alpha_p[tg] ~ dnorm(0, sd = 2)

      sigma_alpha_psi1[tg] ~ dexp(1)
      sigma_alpha_gamma[tg] ~ dexp(1)
      sigma_alpha_eps[tg] ~ dexp(1)
      sigma_alpha_p[tg] ~ dexp(1)
    }

    # --- PRIORS COMUNI ---
    beta_tin_psi1 ~ dnorm(0, sd = 2)
    beta_year_gamma ~ dnorm(0, sd = 2)
    beta_tin_gamma ~ dnorm(0, sd = 2)
    beta_year_eps ~ dnorm(0, sd = 2)
    beta_effort ~ dnorm(0, sd = 2)
    beta_year_p ~ dnorm(0, sd = 2)

    # --- PRIORS SPECIFICI PER SPECIE ---
    for (s in 1:n_species) {
      alpha_psi1[s] ~ dnorm(
        mu_alpha_psi1[gindex[s]],
        sd = sigma_alpha_psi1[gindex[s]]
      )
      alpha_gamma[s] ~ dnorm(
        mu_alpha_gamma[gindex[s]],
        sd = sigma_alpha_gamma[gindex[s]]
      )
      alpha_eps[s] ~ dnorm(
        mu_alpha_eps[gindex[s]],
        sd = sigma_alpha_eps[gindex[s]]
      )
      alpha_p[s] ~ dnorm(mu_alpha_p[gindex[s]], sd = sigma_alpha_p[gindex[s]])
    }

    # --- LIKELIHOOD ---
    for (s in 1:n_species) {
      for (i in 1:n_sites) {
        # Anno 1
        logit(psi1[s, i]) <- alpha_psi1[s] + beta_tin_psi1 * tin_1[s, i]
        z[s, i, 1] ~ dbern(psi1[s, i])

        # Anni successivi
        for (t in 2:n_years) {
          logit(gamma[s, i, t - 1]) <- alpha_gamma[s] +
            beta_year_gamma * year[t] +
            beta_tin_gamma * tin_1[s, i]
          logit(eps[s, i, t - 1]) <- alpha_eps[s] + beta_year_eps * year[t]
          prob_occ[s, i, t] <- z[s, i, t - 1] *
            (1 - eps[s, i, t - 1]) +
            (1 - z[s, i, t - 1]) * gamma[s, i, t - 1]
          z[s, i, t] ~ dbern(prob_occ[s, i, t])
        }
      }
    }

    # Rilevamento
    for (k in 1:n_obs) {
      logit(p[k]) <- alpha_p[species_obs[k]] +
        beta_effort * effort_obs[k] +
        beta_year_p * year_value_obs[k]
      mu[k] <- z[species_obs[k], site_obs[k], year_idx_obs[k]] * p[k]
      y_obs[k] ~ dbern(mu[k])
    }

    # Quantità derivate per singola specie
    for (s in 1:n_species) {
      for (t in 1:n_years) {
        mean_psi[s, t] <- sum(z[s, 1:n_sites, t]) / n_sites
      }
    }
  })

  # 3. INIZIALIZZAZIONE
  z_init <- array(NA, dim = c(n_species, nrow(y_list[[1]]), ncol(y_list[[1]])))

  for (s in 1:n_species) {
    sp <- all_species[s]
    y_matrix <- y_list[[sp]]
    prob_osservata <- mean(y_matrix, na.rm = TRUE)
    if (is.na(prob_osservata) || prob_osservata == 0) {
      prob_osservata <- 0.5
    }

    z_temp <- matrix(NA, nrow = nrow(y_matrix), ncol = ncol(y_matrix))

    for (i in 1:nrow(y_matrix)) {
      for (t in 1:ncol(y_matrix)) {
        if (!is.na(y_matrix[i, t])) {
          if (y_matrix[i, t] == 1) {
            z_temp[i, t] <- 1
          } else {
            z_temp[i, t] <- rbinom(1, 1, 0.3)
          }
        } else {
          z_temp[i, t] <- rbinom(1, 1, prob_osservata)
        }
      }
    }

    z_init[s, , ] <- z_temp
  }

  # Preparare tin_1 come matrice
  tin_matrix <- matrix(NA, nrow = n_species, ncol = nrow(y_list[[1]]))
  for (s in 1:n_species) {
    sp <- all_species[s]
    tin_matrix[s, ] <- tin_std_list[[sp]]
  }

  nim_constants <- list(
    n_sites = nrow(y_list[[1]]),
    n_years = ncol(y_list[[1]]),
    n_species = n_species,
    n_groups = n_groups,
    n_obs = nrow(obs_data),
    species_obs = obs_data$species_idx,
    site_obs = obs_data$site_idx,
    year_idx_obs = obs_data$year_idx,
    tin_1 = tin_matrix,
    year = year_std,
    gindex = gindex_vec
  )

  nim_data <- list(
    y_obs = obs_data$y_value,
    effort_obs = obs_data$effort_value,
    year_value_obs = obs_data$year_value
  )

  # 4. ESECUZIONE CON PARALLELIZZAZIONE
  cat("\nEsecuzione MCMC in parallelo con", ncores, "core...\n")

  run_chain <- function(
    chain_id,
    seed_val,
    model_code,
    constants,
    data,
    n_sp,
    n_groups,
    y_matrices
  ) {
    set.seed(seed_val)

    z_init_local <- array(
      NA,
      dim = c(n_sp, constants$n_sites, constants$n_years)
    )

    for (s in 1:n_sp) {
      y_mat <- y_matrices[[s]]
      prob_init <- runif(1, 0.4, 0.6)

      for (i in 1:constants$n_sites) {
        for (t in 1:constants$n_years) {
          if (!is.na(y_mat[i, t])) {
            if (y_mat[i, t] == 1) {
              z_init_local[s, i, t] <- 1
            } else {
              z_init_local[s, i, t] <- rbinom(1, 1, 0.3)
            }
          } else {
            z_init_local[s, i, t] <- rbinom(1, 1, prob_init)
          }
        }
      }
    }

    initial_values <- list(
      mu_alpha_psi1 = rnorm(n_groups, 0, 0.5),
      mu_alpha_gamma = rnorm(n_groups, 0, 0.5),
      mu_alpha_eps = rnorm(n_groups, 0, 0.5),
      mu_alpha_p = rnorm(n_groups, 0, 0.5),
      sigma_alpha_psi1 = rexp(n_groups, 1),
      sigma_alpha_gamma = rexp(n_groups, 1),
      sigma_alpha_eps = rexp(n_groups, 1),
      sigma_alpha_p = rexp(n_groups, 1),
      alpha_psi1 = rnorm(n_sp, 0, 0.5),
      alpha_gamma = rnorm(n_sp, 0, 0.5),
      alpha_eps = rnorm(n_sp, 0, 0.5),
      alpha_p = rnorm(n_sp, 0, 0.5),
      beta_tin_psi1 = rnorm(1, 0, 0.5),
      beta_year_gamma = rnorm(1, 0, 0.5),
      beta_tin_gamma = rnorm(1, 0, 0.5),
      beta_year_eps = rnorm(1, 0, 0.5),
      beta_effort = rnorm(1, 0, 0.5),
      beta_year_p = rnorm(1, 0, 0.5),
      z = z_init_local
    )

    suppressMessages({
      model <- nimbleModel(
        model_code,
        constants = constants,
        data = data,
        inits = initial_values
      )
      c_model <- compileNimble(model, showCompilerOutput = FALSE)

      monitors_vec <- c(
        'mu_alpha_psi1',
        'mu_alpha_gamma',
        'mu_alpha_eps',
        'mu_alpha_p',
        'sigma_alpha_psi1',
        'sigma_alpha_gamma',
        'sigma_alpha_eps',
        'sigma_alpha_p',
        'mean_psi'
      )

      mcmc_conf <- configureMCMC(model, monitors = monitors_vec)
      mcmc <- buildMCMC(mcmc_conf)
      c_mcmc <- compileNimble(mcmc, project = model, showCompilerOutput = FALSE)
    })

    samples <- runMCMC(
      c_mcmc,
      niter = 20000,
      nburnin = 5000,
      thin = 10,
      samplesAsCodaMCMC = TRUE,
      progressBar = TRUE
    )

    return(samples)
  }

  cl <- makeCluster(ncores, outfile = "")
  invisible(clusterEvalQ(cl, {
    library(nimble)
    library(coda)
  }))

  samples_list <- tryCatch(
    {
      parLapply(cl, 1:3, function(i) {
        run_chain(
          i,
          1000 + i * 123,
          code,
          nim_constants,
          nim_data,
          n_species,
          n_groups,
          y_list
        )
      })
    },
    finally = {
      stopCluster(cl)
    }
  )

  samples <- as.mcmc.list(samples_list)

  # Estrarre colonne mean_psi una sola volta
  samples_mat <- as.matrix(samples)
  mean_psi_samples <- samples_mat[, grep("mean_psi\\[", colnames(samples_mat))]

  # 5. CALCOLARE group_mean_psi POST-HOC (per gruppo)
  cat("\nCalcolando statistiche di gruppo...\n")

  group_results <- list()
  for (g in 1:n_groups) {
    species_in_group <- which(gindex_vec == g)

    if (length(species_in_group) == 0) {
      cat("Nessuna specie nel gruppo", g, "\n")
      next
    }

    group_mean_by_year <- matrix(
      NA,
      nrow = nrow(mean_psi_samples),
      ncol = length(anni_num)
    )

    for (t in 1:length(anni_num)) {
      cols_year_t <- grep(paste0(", ", t, "\\]"), colnames(mean_psi_samples))
      cols_group_year <- cols_year_t[sapply(cols_year_t, function(col) {
        sp_idx <- as.numeric(gsub(
          ".*\\[(\\d+),.*",
          "\\1",
          colnames(mean_psi_samples)[col]
        ))
        sp_idx %in% species_in_group
      })]

      if (length(cols_group_year) > 0) {
        group_mean_by_year[, t] <- rowMeans(mean_psi_samples[,
          cols_group_year,
          drop = FALSE
        ])
      }
    }

    df_plot <- data.frame(
      Anno = anni_num,
      Group = g,
      Occupancy = colMeans(group_mean_by_year, na.rm = TRUE),
      Lower = apply(
        group_mean_by_year,
        2,
        quantile,
        probs = 0.025,
        na.rm = TRUE
      ),
      Upper = apply(
        group_mean_by_year,
        2,
        quantile,
        probs = 0.975,
        na.rm = TRUE
      )
    )

    group_results[[g]] <- list(
      group = g,
      n_species = length(species_in_group),
      plot_data = df_plot,
      samples = group_mean_by_year
    )
  }

  # 6. CALCOLO STIMA GLOBALE (tutte le specie, post-hoc)
  cat("\nCalcolando stima globale (tutte le specie)...\n")

  global_mean_by_year <- matrix(
    NA,
    nrow = nrow(mean_psi_samples),
    ncol = length(anni_num)
  )

  for (t in 1:length(anni_num)) {
    # Tutte le colonne per l'anno t, indipendentemente dalla specie
    cols_year_t <- grep(paste0(", ", t, "\\]"), colnames(mean_psi_samples))
    if (length(cols_year_t) > 0) {
      global_mean_by_year[, t] <- rowMeans(
        mean_psi_samples[, cols_year_t, drop = FALSE]
      )
    }
  }

  global_plot_data <- data.frame(
    Anno = anni_num,
    Occupancy = colMeans(global_mean_by_year, na.rm = TRUE),
    Lower = apply(
      global_mean_by_year,
      2,
      quantile,
      probs = 0.025,
      na.rm = TRUE
    ),
    Upper = apply(global_mean_by_year, 2, quantile, probs = 0.975, na.rm = TRUE)
  )

  global_results <- list(
    n_species = n_species,
    plot_data = global_plot_data,
    samples = global_mean_by_year
  )

  # 7. CALCOLO STIMA PER SINGOLA SPECIE (post-hoc)
  cat("\nCalcolando stima per singola specie...\n")

  species_results <- list()
  for (s in 1:n_species) {
    sp <- all_species[s]

    species_mean_by_year <- matrix(
      NA,
      nrow = nrow(mean_psi_samples),
      ncol = length(anni_num)
    )

    for (t in 1:length(anni_num)) {
      col_name <- paste0("mean_psi[", s, ", ", t, "]")
      if (col_name %in% colnames(mean_psi_samples)) {
        species_mean_by_year[, t] <- mean_psi_samples[, col_name]
      }
    }

    species_plot_data <- data.frame(
      Anno = anni_num,
      Species = sp,
      Group = gindex_vec[s],
      Occupancy = colMeans(species_mean_by_year, na.rm = TRUE),
      Lower = apply(
        species_mean_by_year,
        2,
        quantile,
        probs = 0.025,
        na.rm = TRUE
      ),
      Upper = apply(
        species_mean_by_year,
        2,
        quantile,
        probs = 0.975,
        na.rm = TRUE
      )
    )

    species_results[[sp]] <- list(
      species = sp,
      group = gindex_vec[s],
      plot_data = species_plot_data,
      samples = species_mean_by_year
    )
  }

  # 8. PLOT: confronto tra gruppi (come prima)
  all_group_data <- bind_rows(lapply(group_results, function(x) x$plot_data))

  pGroups <- all_group_data %>%
    filter(Anno > 2004) %>%
    ggplot(aes(
      x = Anno,
      y = Occupancy,
      color = factor(Group),
      fill = factor(Group)
    )) +
    geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.2, color = NA) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 2) +
    scale_color_viridis_d(
      name = "Group",
      labels = paste("Group", 1:n_groups)
    ) +
    scale_fill_viridis_d(
      name = "Group",
      labels = paste("Group", 1:n_groups)
    ) +
    scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
    labs(
      title = "Occupancy Comparison across Groups",
      x = "Year",
      y = "Mean Occupancy"
    ) +
    theme_bw(base_size = 14) +
    theme(legend.position = "bottom")

  print(pGroups)

  ggsave(
    paste0("output_bayes/ALL_GROUPS_comparison_", today(), ".pdf"),
    width = 12,
    height = 8
  )

  # 9. PLOT: stima globale (tutte le specie)
  pGlobal <- global_plot_data %>%
    filter(Anno > 2004) %>%
    ggplot(aes(x = Anno, y = Occupancy)) +
    geom_ribbon(
      aes(ymin = Lower, ymax = Upper),
      alpha = 0.25,
      fill = "steelblue"
    ) +
    geom_line(linewidth = 1.4, color = "steelblue4") +
    geom_point(size = 2.5, color = "steelblue4") +
    scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
    labs(
      title = paste0(
        "Global Mean Occupancy — All Species (n = ",
        n_species,
        ")"
      ),
      x = "Year",
      y = "Mean Occupancy"
    ) +
    theme_bw(base_size = 14)

  print(pGlobal)

  ggsave(
    paste0("output_bayes/GLOBAL_occupancy_", today(), ".pdf"),
    plot = pGlobal,
    width = 10,
    height = 6
  )

  # 10. DIAGNOSTICHE
  rhat_df <- tryCatch(
    {
      as.data.frame(gelman.diag(samples, multivariate = FALSE)$psrf)
    },
    error = function(e) {
      data.frame(Point.est = NA, Upper.C.I. = NA)
    }
  )

  # 11. CONFRONTI STATISTICI: tra gruppi
  cat("\n=== CONFRONTI TRA GRUPPI ===\n")

  mu_psi1_samples <- samples_mat[, grep(
    "^mu_alpha_psi1\\[",
    colnames(samples_mat)
  )]
  mu_gamma_samples <- samples_mat[, grep(
    "^mu_alpha_gamma\\[",
    colnames(samples_mat)
  )]

  comparisons <- list()
  for (g1 in 1:(n_groups - 1)) {
    for (g2 in (g1 + 1):n_groups) {
      diff_psi1 <- mu_psi1_samples[, g2] - mu_psi1_samples[, g1]
      diff_gamma <- mu_gamma_samples[, g2] - mu_gamma_samples[, g1]

      cat("\n--- Gruppo", g2, "vs Gruppo", g1, "---\n")
      cat("Occupancy iniziale:\n")
      cat("  Prob(G", g2, "> G", g1, "):", mean(diff_psi1 > 0), "\n")
      cat("Colonizzazione:\n")
      cat("  Prob(G", g2, "> G", g1, "):", mean(diff_gamma > 0), "\n")

      comparisons[[paste0("G", g2, "_vs_G", g1)]] <- list(
        diff_psi1 = diff_psi1,
        diff_gamma = diff_gamma
      )
    }
  }

  # 12. CONFRONTI STATISTICI: ogni gruppo vs stima globale
  cat("\n=== CONFRONTI GRUPPI vs GLOBALE ===\n")

  group_vs_global <- list()
  for (g in 1:n_groups) {
    if (is.null(group_results[[g]])) {
      next
    }

    # Differenza tra media annuale del gruppo e media globale, per ogni campione MCMC
    diff_samples <- group_results[[g]]$samples - global_mean_by_year

    # Prob che il gruppo sia > globale (aggregato su anni)
    prob_above <- mean(diff_samples > 0, na.rm = TRUE)

    # Differenza media e CI per anno
    diff_by_year <- data.frame(
      Anno = anni_num,
      Group = g,
      Diff_mean = colMeans(diff_samples, na.rm = TRUE),
      Diff_lower = apply(
        diff_samples,
        2,
        quantile,
        probs = 0.025,
        na.rm = TRUE
      ),
      Diff_upper = apply(diff_samples, 2, quantile, probs = 0.975, na.rm = TRUE)
    )

    cat(
      "\n--- Gruppo",
      g,
      "(n =",
      group_results[[g]]$n_species,
      "specie) vs Globale ---\n"
    )
    cat("  Prob(Gruppo > Globale):", round(prob_above, 3), "\n")
    cat(
      "  Differenza media (pooled):",
      round(mean(diff_samples, na.rm = TRUE), 3),
      "\n"
    )

    group_vs_global[[paste0("G", g, "_vs_global")]] <- list(
      group = g,
      prob_above_global = prob_above,
      diff_by_year = diff_by_year,
      diff_samples = diff_samples
    )
  }

  # Plot dei confronti gruppo vs globale
  diff_all <- bind_rows(lapply(group_vs_global, function(x) x$diff_by_year))

  pDiff <- diff_all %>%
    filter(Anno > 2004) %>%
    ggplot(aes(
      x = Anno,
      y = Diff_mean,
      color = factor(Group),
      fill = factor(Group)
    )) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
    geom_ribbon(
      aes(ymin = Diff_lower, ymax = Diff_upper),
      alpha = 0.2,
      color = NA
    ) +
    geom_line(linewidth = 1.2) +
    scale_color_viridis_d(
      name = "Group",
      labels = paste("Group", 1:n_groups)
    ) +
    scale_fill_viridis_d(
      name = "Group",
      labels = paste("Group", 1:n_groups)
    ) +
    labs(
      title = "Difference in Occupancy: Each Group vs Global Mean",
      x = "Year",
      y = "Difference (Group – Global)"
    ) +
    theme_bw(base_size = 14) +
    theme(legend.position = "bottom")

  print(pDiff)

  ggsave(
    paste0("output_bayes/GROUP_vs_GLOBAL_diff_", today(), ".pdf"),
    plot = pDiff,
    width = 12,
    height = 7
  )

  # 13. SALVATAGGIO
  saveRDS(
    list(
      n_species = n_species,
      species = all_species,
      gindex = gindex_vec,
      samples = samples,
      group_results = group_results,
      global_results = global_results,
      species_results = species_results,
      comparisons = comparisons,
      group_vs_global = group_vs_global,
      rhat = rhat_df
    ),
    paste0("output_bayes/MULTISPECIES_ALL_GROUPS_", today(), ".rds")
  )

  cat("\nAnalisi completata!\n")

  return(list(
    group_results = group_results,
    global_results = global_results,
    species_results = species_results,
    comparisons = comparisons,
    group_vs_global = group_vs_global,
    rhat = rhat_df
  ))
}
