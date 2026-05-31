#' Modello multispecie con stima della ricchezza (Tingley et al. 2020)
#' family e Habitat
#'
#' === MODIFICHE RISPETTO ALLA VERSIONE ORIGINALE ===
#'
#' Integrazione dell'approccio di Tingley et al. (2020) per la stima della
#' ricchezza di specie tramite data augmentation e parametro omega.
#'
#' LOGICA:
#'   1. Al dataset reale si aggiungono n_zero specie "fantasma" (all-zero).
#'   2. Ogni specie (reale + fantasma) riceve w[s] ~ dbern(omega), che indica
#'      se è un vero membro della comunità regionale.
#'   3. La probabilità di occupancy è condizionata a w[s]: solo le specie con
#'      w[s]=1 possono occupare siti.
#'   4. La ricchezza per anno t è derivata in-model come:
#'         richness[t] = Σ_s  w[s] * I(Σ_i z[s,i,t] > 0)
#'   5. Le specie augmentate attingono da un grand-mean hyperprior comune
#'      anziché dalla struttura per gruppo (gindex), che è specifica delle
#'      specie osservate.
#'
#' PARAMETRI NUOVI:
#'   omega               ~ dunif(0, 1)   — proporzione specie nel pool regionale
#'   w[s]                ~ dbern(omega)  — appartenenza al pool (tutte le specie)
#'   mu_alpha_*_grand    ~ dnorm(0, 2)   — grand-mean hyperprior per specie augmentate
#'   sigma_alpha_*_grand ~ dexp(1)
#'   richness[t]         (nodo deterministico) — ricchezza realizzata per anno
#'
#' QUANTITÀ MONITORATE AGGIUNTIVE:
#'   "omega", "richness"
#'
#' NOTA SUL MODELLO DINAMICO:
#'   Nel processo di stato, prob_occ è moltiplicata per w[s] anche negli anni
#'   successivi, garantendo che le specie fantasma con w[s]=0 rimangano z=0.

odoMultiSpecies_AllGroups <- function(species_index_df,
                                      nk = 10,
                                      ncores = 7,
                                      n_zero = 50) {   # <<< NUOVO parametro

  cat("Preparazione dati multispecie per TUTTI i gruppi...\n")

  species_subset <- species_index_df %>%
    select(species, gindex) %>%
    distinct()

  all_species <- unique(species_subset$species)
  n_species   <- length(all_species)
  n_total     <- n_species + n_zero   # <<< NUOVO: specie osservate + fantasma

  cat("Numero totale di specie osservate:", n_species, "\n")
  cat("Specie augmentate (n_zero):", n_zero, "\n")
  cat("Pool totale (n_total):", n_total, "\n")
  cat("Distribuzione per gruppo:\n")
  print(table(species_subset$gindex))

  # --------------------------------------------------------------------------
  # Preparazione matrice Y e covariata tin per ogni specie osservata
  # --------------------------------------------------------------------------
  y_list         <- list()
  tin_std_list   <- list()

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
          n_occurrences > 0     ~ 1,
          tot_occurrences > 0   ~ 0
        )
      ) %>%
      pivot_wider(
        id_cols    = cellcodeX,
        names_from = date_year,
        values_from = detection
      ) %>%
      column_to_rownames("cellcodeX") %>%
      as.matrix()

    y_list[[sp]] <- y_matrix

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

  # Preparazione obs_data (solo specie osservate, invariato)
  obs_data <- data.frame(
    species_idx = integer(),
    site_idx    = integer(),
    year_idx    = integer(),
    y_value     = numeric(),
    effort_value = numeric(),
    year_value  = numeric()
  )

  for (s in 1:n_species) {
    sp       <- all_species[s]
    y_matrix <- y_list[[sp]]

    for (i in 1:nrow(y_matrix)) {
      for (t in 1:ncol(y_matrix)) {
        if (!is.na(y_matrix[i, t])) {
          obs_data <- rbind(
            obs_data,
            data.frame(
              species_idx  = s,
              site_idx     = i,
              year_idx     = t,
              y_value      = y_matrix[i, t],
              effort_value = effort_std_mat[i, t],
              year_value   = year_std[t]
            )
          )
        }
      }
    }
  }

  cat("N. osservazioni totali:", nrow(obs_data), "\n")

  gindex_vec <- species_subset$gindex[match(all_species, species_subset$species)]
  n_groups   <- length(unique(gindex_vec))

  cat("Numero gruppi:", n_groups, "\n")

  if (!all(sort(unique(gindex_vec)) == 1:n_groups)) {
    stop(
      "gindex deve contenere valori contigui da 1 a ", n_groups,
      ". Valori trovati: ", paste(sort(unique(gindex_vec)), collapse = ", ")
    )
  }

  family_lookup <- dfo %>%
    select(species, family) %>%
    distinct() %>%
    filter(species %in% all_species)

  family_vec  <- family_lookup$family[match(all_species, family_lookup$species)]
  all_families <- sort(unique(family_vec))
  n_families  <- length(all_families)

  cat("Famiglie trovate:", n_families, "\n")

  # --------------------------------------------------------------------------
  # 2. DEFINIZIONE MODELLO NIMBLE (con integrazione Tingley)
  # --------------------------------------------------------------------------
  code <- nimbleCode({

    # --- OMEGA: proporzione specie nel pool regionale --- <<< NUOVO
    omega ~ dunif(0, 1)

    # --- W: indicatore di appartenenza alla comunità --- <<< NUOVO
    # Tutte le specie (osservate + augmentate) hanno w[s] ~ dbern(omega)
    for (s in 1:n_total) {
      w[s] ~ dbern(omega)
    }

    # --- HYPERPRIORS PER OGNI GRUPPO (invariati, per specie osservate) ---
    for (tg in 1:n_groups) {
      mu_alpha_psi1[tg]  ~ dnorm(0, sd = 2)
      mu_alpha_gamma[tg] ~ dnorm(0, sd = 2)
      mu_alpha_eps[tg]   ~ dnorm(0, sd = 2)
      mu_alpha_p[tg]     ~ dnorm(0, sd = 2)

      sigma_alpha_psi1[tg]  ~ dexp(1)
      sigma_alpha_gamma[tg] ~ dexp(1)
      sigma_alpha_eps[tg]   ~ dexp(1)
      sigma_alpha_p[tg]     ~ dexp(1)
    }

    # --- GRAND-MEAN HYPERPRIOR (per specie augmentate) --- <<< NUOVO
    # Le specie fantasma non appartengono a nessun gruppo osservato;
    # attingono da una distribuzione comunitaria "grand-mean" separata.
    mu_alpha_psi1_grand  ~ dnorm(0, sd = 2)
    mu_alpha_gamma_grand ~ dnorm(0, sd = 2)
    mu_alpha_eps_grand   ~ dnorm(0, sd = 2)
    mu_alpha_p_grand     ~ dnorm(0, sd = 2)

    sigma_alpha_psi1_grand  ~ dexp(1)
    sigma_alpha_gamma_grand ~ dexp(1)
    sigma_alpha_eps_grand   ~ dexp(1)
    sigma_alpha_p_grand     ~ dexp(1)

    # --- PRIORS COMUNI (invariati) ---
    beta_tin_psi1  ~ dnorm(0, sd = 2)
    beta_year_gamma ~ dnorm(0, sd = 2)
    beta_tin_gamma  ~ dnorm(0, sd = 2)
    beta_year_eps   ~ dnorm(0, sd = 2)
    beta_effort     ~ dnorm(0, sd = 2)
    beta_year_p     ~ dnorm(0, sd = 2)

    # --- PRIORS SPECIFICI PER SPECIE OSSERVATE (invariati) ---
    for (s in 1:n_species) {
      alpha_psi1[s]  ~ dnorm(mu_alpha_psi1[gindex[s]],  sd = sigma_alpha_psi1[gindex[s]])
      alpha_gamma[s] ~ dnorm(mu_alpha_gamma[gindex[s]], sd = sigma_alpha_gamma[gindex[s]])
      alpha_eps[s]   ~ dnorm(mu_alpha_eps[gindex[s]],   sd = sigma_alpha_eps[gindex[s]])
      alpha_p[s]     ~ dnorm(mu_alpha_p[gindex[s]],     sd = sigma_alpha_p[gindex[s]])
    }

    # --- PRIORS PER SPECIE AUGMENTATE --- <<< NUOVO
    # Attingono dal grand-mean hyperprior, non dalla struttura per gruppo.
    # alpha_p non compare nella likelihood (nessuna osservazione per queste
    # specie) ma viene definito per completezza del superpopulation model.
    for (s in (n_species + 1):n_total) {
      alpha_psi1[s]  ~ dnorm(mu_alpha_psi1_grand,  sd = sigma_alpha_psi1_grand)
      alpha_gamma[s] ~ dnorm(mu_alpha_gamma_grand, sd = sigma_alpha_gamma_grand)
      alpha_eps[s]   ~ dnorm(mu_alpha_eps_grand,   sd = sigma_alpha_eps_grand)
      alpha_p[s]     ~ dnorm(mu_alpha_p_grand,     sd = sigma_alpha_p_grand)
    }

    # --- PROCESSO DI STATO (modificato per condizionare su w[s]) --- <<< MODIFICATO
    # Il loop ora gira su n_total anziché n_species.
    # Per specie con w[s]=0, mu_psi1=0 e prob_occ=0, quindi z rimane 0.
    for (s in 1:n_total) {
      for (i in 1:n_sites) {

        # Anno 1: occupancy iniziale condizionata su w[s]
        logit(psi1[s, i]) <- alpha_psi1[s] + beta_tin_psi1 * tin_1[s, i]
        mu_psi1[s, i]     <- psi1[s, i] * w[s]   # <<< NUOVO: w[s] come cancello
        z[s, i, 1]        ~ dbern(mu_psi1[s, i])

        # Anni successivi: gamma/eps anch'essi condizionati su w[s]
        for (t in 2:n_years) {
          logit(gamma[s, i, t - 1]) <- alpha_gamma[s] +
            beta_year_gamma * year[t] +
            beta_tin_gamma  * tin_1[s, i]
          logit(eps[s, i, t - 1]) <- alpha_eps[s] + beta_year_eps * year[t]

          # <<< MODIFICATO: w[s] moltiplica prob_occ → specie con w=0 restano z=0
          prob_occ[s, i, t] <- w[s] * (
            z[s, i, t - 1] * (1 - eps[s, i, t - 1]) +
            (1 - z[s, i, t - 1]) * gamma[s, i, t - 1]
          )
          z[s, i, t] ~ dbern(prob_occ[s, i, t])
        }
      }
    }

    # --- PROCESSO DI RILEVAMENTO (invariato: solo specie osservate) ---
    # species_obs[k] ∈ {1, ..., n_species}, mai > n_species, quindi
    # le specie augmentate non entrano mai qui.
    for (k in 1:n_obs) {
      logit(p[k]) <- alpha_p[species_obs[k]] +
        beta_effort * effort_obs[k] +
        beta_year_p * year_value_obs[k]
      mu[k]      <- z[species_obs[k], site_obs[k], year_idx_obs[k]] * p[k]
      y_obs[k]   ~ dbern(mu[k])
    }

    # --- QUANTITÀ DERIVATE ---

    # Occupancy media per specie osservata (invariato)
    for (s in 1:n_species) {
      for (t in 1:n_years) {
        mean_psi[s, t] <- sum(z[s, 1:n_sites, t]) / n_sites
      }
    }

    # Ricchezza realizzata per anno <<< NUOVO
    # sp_present[s,t] = 1 se la specie s è nel pool (w=1) E occupa almeno un sito
    # richness[t] = numero di specie presenti in almeno un sito nell'anno t
    for (t in 1:n_years) {
      for (s in 1:n_total) {
        sp_present[s, t] <- w[s] * step(sum(z[s, 1:n_sites, t]) - 0.5)
      }
      richness[t] <- sum(sp_present[1:n_total, t])
    }

    # Pool size: stima del numero totale di specie nel pool regionale <<< NUOVO
    # Equivalente a omega * n_total (calcolabile post-hoc da omega)
    N_pool <- sum(w[1:n_total])

  })

  # --------------------------------------------------------------------------
  # 3. PREPARAZIONE DATI PER NIMBLE
  # --------------------------------------------------------------------------
  n_sites <- nrow(y_list[[1]])
  n_years <- ncol(y_list[[1]])

  # Matrice tin_1: per specie osservate usa valori reali,
  # per specie augmentate usa 0 (= media sulla scala standardizzata)
  tin_matrix_obs <- matrix(NA, nrow = n_species, ncol = n_sites)
  for (s in 1:n_species) {
    sp <- all_species[s]
    tin_matrix_obs[s, ] <- tin_std_list[[sp]]
  }

  # <<< NUOVO: estendi tin_matrix con righe di zeri per le specie augmentate
  tin_matrix_aug <- rbind(
    tin_matrix_obs,
    matrix(0, nrow = n_zero, ncol = n_sites)
  )

  nim_constants <- list(
    n_sites   = n_sites,
    n_years   = n_years,
    n_species = n_species,
    n_zero    = n_zero,       # <<< NUOVO
    n_total   = n_total,      # <<< NUOVO
    n_groups  = n_groups,
    n_obs     = nrow(obs_data),
    species_obs    = obs_data$species_idx,
    site_obs       = obs_data$site_idx,
    year_idx_obs   = obs_data$year_idx,
    tin_1     = tin_matrix_aug,   # <<< MODIFICATO: ora n_total × n_sites
    year      = year_std,
    gindex    = gindex_vec        # solo per s in 1:n_species, invariato
  )

  nim_data <- list(
    y_obs          = obs_data$y_value,
    effort_obs     = obs_data$effort_value,
    year_value_obs = obs_data$year_value
  )

  # --------------------------------------------------------------------------
  # 4. INIZIALIZZAZIONE
  # --------------------------------------------------------------------------

  # z_init per specie osservate (come prima)
  z_init_obs <- array(NA, dim = c(n_species, n_sites, n_years))
  for (s in 1:n_species) {
    sp <- all_species[s]
    y_matrix <- y_list[[sp]]
    prob_osservata <- mean(y_matrix, na.rm = TRUE)
    if (is.na(prob_osservata) || prob_osservata == 0) prob_osservata <- 0.5

    for (i in 1:n_sites) {
      for (t in 1:n_years) {
        if (!is.na(y_matrix[i, t])) {
          z_init_obs[s, i, t] <- if (y_matrix[i, t] == 1) 1L else rbinom(1, 1, 0.3)
        } else {
          z_init_obs[s, i, t] <- rbinom(1, 1, prob_osservata)
        }
      }
    }
  }

  # <<< NUOVO: z_init per specie augmentate: tutto 0 (non osservate)
  z_init_aug <- array(0L, dim = c(n_zero, n_sites, n_years))

  # Combina: array [n_total, n_sites, n_years]
  z_init_full <- array(
    c(z_init_obs, z_init_aug),
    dim = c(n_total, n_sites, n_years)
  )

  # --------------------------------------------------------------------------
  # 5. ESECUZIONE CON PARALLELIZZAZIONE
  # --------------------------------------------------------------------------
  cat("\nEsecuzione MCMC in parallelo con", ncores, "core...\n")

  run_chain <- function(chain_id, seed_val, model_code, constants, data,
                        n_sp, n_zero_in, n_total_in, n_groups,
                        z_init_full_in, y_matrices) {
    set.seed(seed_val)

    n_sites_in <- constants$n_sites
    n_years_in <- constants$n_years

    # z iniziali locali (variazione per catena)
    z_init_local <- z_init_full_in
    for (s in 1:n_sp) {
      y_mat     <- y_matrices[[s]]
      prob_init <- runif(1, 0.4, 0.6)
      for (i in 1:n_sites_in) {
        for (t in 1:n_years_in) {
          if (!is.na(y_mat[i, t])) {
            z_init_local[s, i, t] <- if (y_mat[i, t] == 1) 1L else rbinom(1, 1, 0.3)
          } else {
            z_init_local[s, i, t] <- rbinom(1, 1, prob_init)
          }
        }
      }
    }
    # specie augmentate: z = 0
    z_init_local[(n_sp + 1):n_total_in, , ] <- 0L

    initial_values <- list(
      # <<< NUOVO: omega e w
      omega = runif(1, 0.3, 0.7),
      w     = c(rep(1L, n_sp), rep(0L, n_zero_in)),   # w=1 per osservate, 0 per fantasma

      # Hyperpriors per gruppo (invariati)
      mu_alpha_psi1  = rnorm(n_groups, 0, 0.5),
      mu_alpha_gamma = rnorm(n_groups, 0, 0.5),
      mu_alpha_eps   = rnorm(n_groups, 0, 0.5),
      mu_alpha_p     = rnorm(n_groups, 0, 0.5),

      sigma_alpha_psi1  = rexp(n_groups, 1),
      sigma_alpha_gamma = rexp(n_groups, 1),
      sigma_alpha_eps   = rexp(n_groups, 1),
      sigma_alpha_p     = rexp(n_groups, 1),

      # <<< NUOVO: grand-mean hyperpriors per specie augmentate
      mu_alpha_psi1_grand  = rnorm(1, 0, 0.5),
      mu_alpha_gamma_grand = rnorm(1, 0, 0.5),
      mu_alpha_eps_grand   = rnorm(1, 0, 0.5),
      mu_alpha_p_grand     = rnorm(1, 0, 0.5),

      sigma_alpha_psi1_grand  = rexp(1, 1),
      sigma_alpha_gamma_grand = rexp(1, 1),
      sigma_alpha_eps_grand   = rexp(1, 1),
      sigma_alpha_p_grand     = rexp(1, 1),

      # Parametri specie-specifici: ora di lunghezza n_total
      alpha_psi1  = rnorm(n_total_in, 0, 0.5),
      alpha_gamma = rnorm(n_total_in, 0, 0.5),
      alpha_eps   = rnorm(n_total_in, 0, 0.5),
      alpha_p     = rnorm(n_total_in, 0, 0.5),

      beta_tin_psi1   = rnorm(1, 0, 0.5),
      beta_year_gamma = rnorm(1, 0, 0.5),
      beta_tin_gamma  = rnorm(1, 0, 0.5),
      beta_year_eps   = rnorm(1, 0, 0.5),
      beta_effort     = rnorm(1, 0, 0.5),
      beta_year_p     = rnorm(1, 0, 0.5),

      z = z_init_local
    )

    suppressMessages({
      model   <- nimbleModel(model_code, constants = constants,
                              data = data, inits = initial_values)
      c_model <- compileNimble(model, showCompilerOutput = FALSE)

      monitors_vec <- c(
        # Parametri originali
        'mu_alpha_psi1', 'mu_alpha_gamma', 'mu_alpha_eps', 'mu_alpha_p',
        'sigma_alpha_psi1', 'sigma_alpha_gamma', 'sigma_alpha_eps', 'sigma_alpha_p',
        'mean_psi',
        # <<< NUOVO: quantità per richness
        'omega',          # proporzione pool regionale
        'richness',       # ricchezza realizzata per anno (vettore n_years)
        'N_pool'          # dimensione pool = sum(w)
      )

      mcmc_conf <- configureMCMC(model, monitors = monitors_vec)
      mcmc      <- buildMCMC(mcmc_conf)
      c_mcmc    <- compileNimble(mcmc, project = model, showCompilerOutput = FALSE)
    })

    samples <- runMCMC(
      c_mcmc,
      niter   = 20000,
      nburnin = 5000,
      thin    = 10,
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
          chain_id      = i,
          seed_val      = 1000 + i * 123,
          model_code    = code,
          constants     = nim_constants,
          data          = nim_data,
          n_sp          = n_species,
          n_zero_in     = n_zero,
          n_total_in    = n_total,
          n_groups      = n_groups,
          z_init_full_in = z_init_full,
          y_matrices    = y_list
        )
      })
    },
    finally = { stopCluster(cl) }
  )

  samples     <- as.mcmc.list(samples_list)
  samples_mat <- as.matrix(samples)

  # --------------------------------------------------------------------------
  # 6. ESTRAZIONE mean_psi (invariato)
  # --------------------------------------------------------------------------
  mean_psi_samples <- samples_mat[, grep("mean_psi\\[", colnames(samples_mat))]

  # --------------------------------------------------------------------------
  # 7. ESTRAZIONE RICHNESS (post-hoc) <<< NUOVO
  # --------------------------------------------------------------------------
  cat("\nEstraendo stime di ricchezza...\n")

  richness_samples <- samples_mat[, grep("^richness\\[", colnames(samples_mat))]
  omega_samples    <- samples_mat[, "omega"]
  N_pool_samples   <- samples_mat[, "N_pool"]

  # Riepilogo per anno
  richness_summary <- data.frame(
    Anno     = anni_num,
    richness_mean  = colMeans(richness_samples),
    richness_lower = apply(richness_samples, 2, quantile, probs = 0.025),
    richness_upper = apply(richness_samples, 2, quantile, probs = 0.975),
    richness_sd    = apply(richness_samples, 2, sd)
  )

  # Stima omega (proporzione pool regionale)
  omega_summary <- c(
    mean  = mean(omega_samples),
    lower = quantile(omega_samples, 0.025),
    upper = quantile(omega_samples, 0.975)
  )

  # Pool size (= n_total * omega a posteriori)
  N_pool_summary <- c(
    mean  = mean(N_pool_samples),
    lower = quantile(N_pool_samples, 0.025),
    upper = quantile(N_pool_samples, 0.975)
  )

  cat("\n=== STIMA RICHNESS ===\n")
  cat("omega (proporzione pool):", round(omega_summary["mean"], 3),
      "[", round(omega_summary["lower.2.5%"], 3), "-",
      round(omega_summary["upper.97.5%"], 3), "]\n")
  cat("N_pool (dimensione pool):", round(N_pool_summary["mean"], 1),
      "[", round(N_pool_summary["lower.2.5%"], 1), "-",
      round(N_pool_summary["upper.97.5%"], 1), "]\n")
  cat("Richness per anno:\n")
  print(richness_summary)

  # --------------------------------------------------------------------------
  # 8. CALCOLO RISULTATI PER GRUPPO, SPECIE, FAMIGLIA (invariati)
  # --------------------------------------------------------------------------

  # --- Gruppi ---
  group_results <- list()
  for (g in 1:n_groups) {
    species_in_group <- which(gindex_vec == g)
    if (length(species_in_group) == 0) { cat("Nessuna specie nel gruppo", g, "\n"); next }

    group_mean_by_year <- matrix(NA, nrow = nrow(mean_psi_samples), ncol = length(anni_num))
    for (t in 1:length(anni_num)) {
      cols_year_t    <- grep(paste0(", ", t, "\\]"), colnames(mean_psi_samples))
      cols_group_year <- cols_year_t[sapply(cols_year_t, function(col) {
        sp_idx <- as.numeric(gsub(".*\\[(\\d+),.*", "\\1", colnames(mean_psi_samples)[col]))
        sp_idx %in% species_in_group
      })]
      if (length(cols_group_year) > 0)
        group_mean_by_year[, t] <- rowMeans(mean_psi_samples[, cols_group_year, drop = FALSE])
    }

    df_plot <- data.frame(
      Anno      = anni_num,
      Group     = g,
      Occupancy = colMeans(group_mean_by_year, na.rm = TRUE),
      Lower     = apply(group_mean_by_year, 2, quantile, probs = 0.025, na.rm = TRUE),
      Upper     = apply(group_mean_by_year, 2, quantile, probs = 0.975, na.rm = TRUE)
    )

    group_results[[g]] <- list(group = g, n_species = length(species_in_group),
                                plot_data = df_plot, samples = group_mean_by_year)
  }

  # --- Specie ---
  species_results <- list()
  for (s in 1:n_species) {
    sp <- all_species[s]
    species_mean_by_year <- matrix(NA, nrow = nrow(mean_psi_samples), ncol = length(anni_num))
    for (t in 1:length(anni_num)) {
      col_name <- paste0("mean_psi[", s, ", ", t, "]")
      if (col_name %in% colnames(mean_psi_samples))
        species_mean_by_year[, t] <- mean_psi_samples[, col_name]
    }
    species_plot_data <- data.frame(
      Anno      = anni_num, Species = sp, Group = gindex_vec[s], Family = family_vec[s],
      Occupancy = colMeans(species_mean_by_year, na.rm = TRUE),
      Lower     = apply(species_mean_by_year, 2, quantile, probs = 0.025, na.rm = TRUE),
      Upper     = apply(species_mean_by_year, 2, quantile, probs = 0.975, na.rm = TRUE)
    )
    species_results[[sp]] <- list(species = sp, group = gindex_vec[s], family = family_vec[s],
                                   plot_data = species_plot_data, samples = species_mean_by_year)
  }

  # --- Famiglie ---
  family_results <- list()
  for (fam in all_families) {
    species_in_family <- which(family_vec == fam)
    if (length(species_in_family) == 0) next

    family_mean_by_year <- matrix(NA, nrow = nrow(mean_psi_samples), ncol = length(anni_num))
    for (t in 1:length(anni_num)) {
      cols_year_t     <- grep(paste0(", ", t, "\\]"), colnames(mean_psi_samples))
      cols_family_year <- cols_year_t[sapply(cols_year_t, function(col) {
        sp_idx <- as.numeric(gsub(".*\\[(\\d+),.*", "\\1", colnames(mean_psi_samples)[col]))
        sp_idx %in% species_in_family
      })]
      if (length(cols_family_year) > 0)
        family_mean_by_year[, t] <- rowMeans(mean_psi_samples[, cols_family_year, drop = FALSE])
    }

    family_plot_data <- data.frame(
      Anno      = anni_num, Family = fam,
      Occupancy = colMeans(family_mean_by_year, na.rm = TRUE),
      Lower     = apply(family_mean_by_year, 2, quantile, probs = 0.025, na.rm = TRUE),
      Upper     = apply(family_mean_by_year, 2, quantile, probs = 0.975, na.rm = TRUE)
    )
    family_results[[fam]] <- list(family = fam, n_species = length(species_in_family),
                                   plot_data = family_plot_data, samples = family_mean_by_year)
  }

  # --------------------------------------------------------------------------
  # 9. PLOT OCCUPANCY (invariati)
  # --------------------------------------------------------------------------
  all_group_data <- bind_rows(lapply(group_results, function(x) x$plot_data))

  pGroups <- all_group_data %>%
    filter(Anno > 2004) %>%
    ggplot(aes(x = Anno, y = Occupancy, color = factor(Group), fill = factor(Group))) +
    geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.2, color = NA) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 2) +
    scale_color_viridis_d(name = "Group", labels = paste("Group", 1:n_groups)) +
    scale_fill_viridis_d(name  = "Group", labels = paste("Group", 1:n_groups)) +
    scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
    labs(title = "Occupancy across HV Groups", x = "Year", y = "Mean Occupancy") +
    theme_bw(base_size = 14) + theme(legend.position = "bottom")
  print(pGroups)

  all_family_data <- bind_rows(lapply(family_results, function(x) x$plot_data))

  pFamilies <- all_family_data %>%
    filter(Anno > 2004) %>%
    ggplot(aes(x = Anno, y = Occupancy, color = Family, fill = Family)) +
    geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.2, color = NA) +
    geom_line(linewidth = 1.2) + geom_point(size = 2) +
    scale_color_viridis_d(name = "Family") + scale_fill_viridis_d(name = "Family") +
    scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
    labs(title = "Occupancy across Families", x = "Year", y = "Mean Occupancy") +
    theme_bw(base_size = 14) + theme(legend.position = "bottom")
  print(pFamilies)

  # --------------------------------------------------------------------------
  # 10. PLOT RICHNESS <<< NUOVO
  # --------------------------------------------------------------------------
  cat("\nGenerando plot richness...\n")

  pRichness <- richness_summary %>%
    filter(Anno > 2004) %>%
    ggplot(aes(x = Anno, y = richness_mean)) +
    geom_ribbon(aes(ymin = richness_lower, ymax = richness_upper),
                alpha = 0.25, fill = "steelblue") +
    geom_line(linewidth = 1.4, color = "steelblue4") +
    geom_point(size = 2.5, color = "steelblue4") +
    # linea orizzontale: specie osservate (lower bound)
    geom_hline(yintercept = n_species, linetype = "dashed",
               color = "firebrick", linewidth = 0.8) +
    annotate("text", x = min(richness_summary$Anno[richness_summary$Anno > 2004]),
             y = n_species + 0.5, label = paste0("Osservate: ", n_species),
             hjust = 0, color = "firebrick", size = 3.5) +
    labs(
      title = paste0("Richness stimata (MSOM, n_zero = ", n_zero, ")"),
      subtitle = paste0("omega = ", round(omega_summary["mean"], 3),
                        " [", round(omega_summary["lower.2.5%"], 3),
                        " - ", round(omega_summary["upper.97.5%"], 3), "]"),
      x = "Anno", y = "Numero di specie"
    ) +
    theme_bw(base_size = 14)
  print(pRichness)

  ggsave(
    paste0("output/output_bayes/multispecies/RICHNESS_", today(), ".pdf"),
    plot = pRichness, width = 10, height = 6
  )

  ggsave(
    paste0("output/output_bayes/multispecies/MULTISPECIES_ALL_GROUPS_HV", today(), ".pdf"),
    plot = pGroups, width = 12, height = 8
  )
  ggsave(
    paste0("output/output_bayes/multispecies/Families", today(), ".pdf"),
    plot = pFamilies, width = 12, height = 8
  )

  # --------------------------------------------------------------------------
  # 11. DIAGNOSTICHE
  # --------------------------------------------------------------------------
  rhat_df <- tryCatch(
    as.data.frame(gelman.diag(samples, multivariate = FALSE)$psrf),
    error = function(e) data.frame(Point.est = NA, Upper.C.I. = NA)
  )

  # --------------------------------------------------------------------------
  # 12. CONFRONTI STATISTICI: tra gruppi (invariato)
  # --------------------------------------------------------------------------
  cat("\n=== CONFRONTI TRA GRUPPI ===\n")

  mu_psi1_samples  <- samples_mat[, grep("^mu_alpha_psi1\\[",  colnames(samples_mat))]
  mu_gamma_samples <- samples_mat[, grep("^mu_alpha_gamma\\[", colnames(samples_mat))]

  comparisons <- list()
  for (g1 in 1:(n_groups - 1)) {
    for (g2 in (g1 + 1):n_groups) {
      diff_psi1  <- mu_psi1_samples[, g2]  - mu_psi1_samples[, g1]
      diff_gamma <- mu_gamma_samples[, g2] - mu_gamma_samples[, g1]

      cat("\n--- Gruppo", g2, "vs Gruppo", g1, "---\n")
      cat("Occupancy iniziale:\n")
      cat("  Prob(G", g2, "> G", g1, "):", mean(diff_psi1 > 0), "\n")
      cat("Colonizzazione:\n")
      cat("  Prob(G", g2, "> G", g1, "):", mean(diff_gamma > 0), "\n")

      comparisons[[paste0("G", g2, "_vs_G", g1)]] <- list(
        diff_psi1 = diff_psi1, diff_gamma = diff_gamma
      )
    }
  }

  # --------------------------------------------------------------------------
  # 13. SALVATAGGIO
  # --------------------------------------------------------------------------
  saveRDS(
    list(
      n_species        = n_species,
      n_zero           = n_zero,
      n_total          = n_total,
      species          = all_species,
      gindex           = gindex_vec,
      family           = family_vec,
      samples          = samples,
      group_results    = group_results,
      family_results   = family_results,
      species_results  = species_results,
      comparisons      = comparisons,
      # <<< NUOVO: risultati richness
      richness_summary = richness_summary,
      richness_samples = richness_samples,
      omega_summary    = omega_summary,
      N_pool_summary   = N_pool_summary,
      rhat             = rhat_df
    ),
    paste0("output/output_bayes/multispecies/MULTISPECIES_ALL_GROUPS_HV", today(), ".rds")
  )

  cat("\nAnalisi completata!\n")

  return(list(
    group_results    = group_results,
    family_results   = family_results,
    species_results  = species_results,
    comparisons      = comparisons,
    richness_summary = richness_summary,
    richness_samples = richness_samples,
    omega_summary    = omega_summary,
    N_pool_summary   = N_pool_summary,
    rhat             = rhat_df
  ))
}
