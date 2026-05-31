# odoRichness_robust.r  — versione con selezione siti integrata
#
# [S1] Selezione stratificata dei siti
#        — specie rare: siti forzati indipendentemente dall'effort
#        — siti restanti: campionamento proporzionale per strato
#          con peso sqrt(effort) (privilegia informatività senza
#          escludere siti a medio-basso sforzo)
#        — stratificazione su gradiente termico tin_1
# [S2] Diagnostica selezione: log + plot + metadati per Metodi paper
# [P1] Catene SEQUENZIALI → progressBar visibile, RAM ridotta
# [P2] CHECKPOINT per catena → crash non perde il lavoro fatto
# [P3] RESUME automatico → rileva checkpoint esistenti
# [P4] LOG su file → tail -f richness_run.log
# [P5] TEST-RUN → stima tempo prima della run completa
#
# uso:
#   res <- odoRichness_robust(species_index_df, all_cells, all_years)
#   tail -f output/output_bayes/multispecies/richness_run.log

odoRichness_robust <- function(species_index_df,
                               all_cells,
                               all_years,
                               # selezione siti
                               n_sites_target  = 100,
                               rare_threshold  = 5,
                               n_strata        = 4,
                               site_seed       = 42,
                               # augmentation
                               n_zero          = 15,
                               # MCMC
                               n_chains        = 3,
                               niter           = 20000,
                               nburnin         = 5000,
                               thin            = 10,
                               test_run        = TRUE,
                               checkpoint_dir  = "output/output_bayes/multispecies/checkpoints") {

  # --------------------------------------------------------------------------
  # UTILITY: log su file + console
  # --------------------------------------------------------------------------
  out_dir  <- dirname(checkpoint_dir)
  log_file <- file.path(out_dir, "richness_run.log")
  dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)

  log_msg <- function(...) {
    msg <- paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", ..., "\n")
    cat(msg)
    cat(msg, file = log_file, append = TRUE)
  }

  log_msg("======================================================")
  log_msg("odoRichness_robust — avvio")
  log_msg("Log:            ", log_file)
  log_msg("Checkpoint dir: ", checkpoint_dir)

  all_species <- unique(species_index_df$species)
  n_species   <- length(all_species)
  log_msg("Specie osservate: ", n_species,
          " | n_zero: ", n_zero,
          " | n_total: ", n_species + n_zero)

  # ==========================================================================
  # SEZIONE 1 — SELEZIONE STRATIFICATA DEI SITI
  # ==========================================================================
  log_msg("------------------------------------------------------")
  log_msg("Sezione 1 — Selezione stratificata siti")
  log_msg("  Siti disponibili: ", length(all_cells))
  log_msg("  Target:           ", n_sites_target)
  log_msg("  Soglia rarità:    presenza < ", rare_threshold, " siti")
  log_msg("  Strati (tin_1):   ", n_strata)

  # 1a. effort e gradiente termico per ogni sito
  site_meta <- dfo %>%
    mutate(cellcodeX = as.character(cellcodeX)) %>%
    filter(cellcodeX %in% as.character(all_cells)) %>%
    group_by(cellcodeX) %>%
    summarise(
      effort   = n_distinct(date_year[!is.na(date_year)]),
      mean_tin = mean(tin_1, na.rm = TRUE),
      .groups  = "drop"
    ) %>%
    mutate(
      stratum = as.integer(cut(mean_tin, breaks = n_strata,
                               labels = FALSE, include.lowest = TRUE))
    )

  eff_q <- quantile(site_meta$effort, c(0, 0.25, 0.5, 0.75, 1))
  log_msg("  Effort (anni survey) — min:", eff_q[1], " Q1:", eff_q[2],
          " med:", eff_q[3], " Q3:", eff_q[4], " max:", eff_q[5])

  # 1b. specie rare e siti che le contengono (da forzare)
  species_site_counts <- dfo %>%
    mutate(cellcodeX = as.character(cellcodeX)) %>%
    filter(cellcodeX %in% as.character(all_cells),
           species %in% all_species) %>%
    group_by(species) %>%
    summarise(n_siti = n_distinct(cellcodeX), .groups = "drop")

  rare_species <- species_site_counts %>%
    filter(n_siti < rare_threshold) %>%
    pull(species)

  sites_with_rare <- dfo %>%
    mutate(cellcodeX = as.character(cellcodeX)) %>%
    filter(species %in% rare_species,
           cellcodeX %in% as.character(all_cells)) %>%
    pull(cellcodeX) %>%
    unique()

  sites_forced <- intersect(site_meta$cellcodeX, sites_with_rare)

  log_msg("  Specie rare (<", rare_threshold, " siti): ",
          length(rare_species), " — ",
          paste(rare_species, collapse = ", "))
  log_msg("  Siti forzati (contengono specie rare): ", length(sites_forced))

  # 1c. allocazione proporzionale nei restanti slot per strato
  n_remaining <- max(0L, n_sites_target - length(sites_forced))

  pool <- site_meta %>%
    filter(!cellcodeX %in% sites_forced, !is.na(stratum))

  stratum_alloc <- pool %>%
    count(stratum) %>%
    mutate(n_alloc = round(n_remaining * n / sum(n)))

  # correggi arrotondamento
  diff_r <- n_remaining - sum(stratum_alloc$n_alloc)
  if (diff_r != 0L)
    stratum_alloc$n_alloc[which.max(stratum_alloc$n)] <-
      stratum_alloc$n_alloc[which.max(stratum_alloc$n)] + diff_r

  # non chiedere più siti di quanti ne esistano
  stratum_alloc <- stratum_alloc %>%
    mutate(n_alloc = pmin(n_alloc, n))

  log_msg("  Allocazione per strato:")
  for (i in seq_len(nrow(stratum_alloc)))
    log_msg(sprintf("    Strato %d (%d disponibili) → %d selezionati",
                    stratum_alloc$stratum[i],
                    stratum_alloc$n[i],
                    stratum_alloc$n_alloc[i]))

  # campionamento con peso sqrt(effort) — loop per strato (n deve essere scalare)
  set.seed(site_seed)
  sites_stratified <- do.call(c, lapply(stratum_alloc$stratum, function(st) {
    pool_st <- pool %>% filter(stratum == st)
    n_st    <- stratum_alloc$n_alloc[stratum_alloc$stratum == st]
    if (n_st == 0 || nrow(pool_st) == 0) return(character(0))
    pool_st %>%
      slice_sample(n = n_st, weight_by = sqrt(effort + 0.5)) %>%
      pull(cellcodeX)
  }))

  selected_cells <- c(sites_forced, sites_stratified)
  n_sites_sel    <- length(selected_cells)

  log_msg("  Siti selezionati: ", n_sites_sel,
          " (", length(sites_forced), " forzati + ",
          length(sites_stratified), " stratificati)")

  eff_sel <- site_meta %>% filter(cellcodeX %in% selected_cells) %>% pull(effort)
  eff_exc <- site_meta %>% filter(!cellcodeX %in% selected_cells) %>% pull(effort)
  log_msg(sprintf("  Effort medio — selezionati: %.1f  |  esclusi: %.1f anni",
                  mean(eff_sel), mean(eff_exc)))

  # 1d. salva metadati selezione (riportabili nella sezione Metodi)
  site_selection_meta <- list(
    n_sites_available  = length(all_cells),
    n_sites_selected   = n_sites_sel,
    n_sites_forced     = length(sites_forced),
    n_sites_stratified = length(sites_stratified),
    rare_species       = rare_species,
    rare_threshold     = rare_threshold,
    n_strata           = n_strata,
    stratum_allocation = stratum_alloc,
    selected_cells     = selected_cells,
    forced_cells       = sites_forced,
    site_seed          = site_seed,
    effort_selected    = summary(eff_sel),
    effort_excluded    = summary(eff_exc)
  )
  saveRDS(site_selection_meta,
          file.path(out_dir, paste0("site_selection_", today(), ".rds")))
  log_msg("  Metadati selezione salvati.")

  # aggiorna all_cells
  all_cells <- selected_cells

  # ==========================================================================
  # SEZIONE 2 — PREPARAZIONE DATI
  # ==========================================================================
  log_msg("------------------------------------------------------")
  log_msg("Sezione 2 — Preparazione dati")

  y_list       <- list()
  tin_std_list <- list()

  for (sp in all_species) {
    y_matrix <- dfo %>%
      filter(species == sp) %>%
      count(cellcodeX, date_year, name = "n_occurrences") %>%
      right_join(
        dfo %>% count(cellcodeX, date_year, name = "tot_occurrences"),
        by = c("cellcodeX", "date_year")
      ) %>%
      complete(
        cellcodeX = all_cells, date_year = all_years,
        fill = list(n_occurrences = NA, tot_occurrences = NA)
      ) %>%
      mutate(detection = case_when(
        is.na(tot_occurrences) ~ NA_real_,
        n_occurrences > 0     ~ 1,
        tot_occurrences > 0   ~ 0
      )) %>%
      pivot_wider(id_cols = cellcodeX, names_from = date_year,
                  values_from = detection) %>%
      column_to_rownames("cellcodeX") %>%
      as.matrix()

    y_list[[sp]] <- y_matrix

    cells_in_y <- rownames(y_matrix)
    tin_data   <- dfo %>%
      select(cellcodeX, tin_1) %>% distinct() %>%
      mutate(cellcodeX = as.character(cellcodeX)) %>%
      filter(cellcodeX %in% cells_in_y) %>%
      group_by(cellcodeX) %>%
      summarise(tin_1 = mean(round(tin_1), na.rm = TRUE), .groups = "drop")
    tin_df <- data.frame(cellcodeX = cells_in_y) %>%
      left_join(tin_data, by = "cellcodeX")
    tin_std_list[[sp]] <- as.vector(scale(tin_df$tin_1))
  }

  n_sites  <- nrow(y_list[[1]])
  n_years  <- ncol(y_list[[1]])
  n_total  <- n_species + n_zero
  anni_num <- as.numeric(colnames(y_list[[1]]))
  year_std <- as.vector(scale(anni_num))

  log_msg("Dimensioni: n_sites=", n_sites,
          " | n_years=", n_years, " | n_total=", n_total)
  log_msg("Nodi processo di stato: ",
          n_total, " x ", n_sites, " x ", n_years,
          " = ", n_total * n_sites * n_years)

  # allinea effort_matrix ai siti selezionati nello stesso ordine di y_list
  effort_subset <- effort_matrix[match(as.character(all_cells),
                                       rownames(effort_matrix)), ,
                                  drop = FALSE]
  effort_global_mean <- mean(as.numeric(effort_matrix), na.rm = TRUE)
  effort_global_sd   <- sd(as.numeric(effort_matrix),   na.rm = TRUE)
  effort_std_raw <- (as.numeric(effort_subset) - effort_global_mean) /
                     effort_global_sd
  effort_std_raw[is.na(effort_std_raw)] <- 0   # NA → media standardizzata
  effort_std_mat <- matrix(effort_std_raw, nrow = n_sites, ncol = n_years)

  obs_data <- data.frame(
    species_idx=integer(), site_idx=integer(),
    year_idx=integer(), y_value=numeric(), effort_value=numeric()
  )
  for (s in seq_len(n_species)) {
    y_matrix <- y_list[[all_species[s]]]
    for (i in seq_len(nrow(y_matrix)))
      for (t in seq_len(ncol(y_matrix)))
        if (!is.na(y_matrix[i, t]))
          obs_data <- rbind(obs_data, data.frame(
            species_idx=s, site_idx=i, year_idx=t,
            y_value=y_matrix[i, t], effort_value=effort_std_mat[i, t]
          ))
  }
  log_msg("Osservazioni totali: ", nrow(obs_data))

  tin_matrix <- rbind(
    do.call(rbind, lapply(all_species, function(sp) tin_std_list[[sp]])),
    matrix(0, nrow = n_zero, ncol = n_sites)
  )

  # ==========================================================================
  # SEZIONE 3 — MODELLO NIMBLE
  # ==========================================================================
  code <- nimbleCode({
    omega ~ dunif(0, 1)
    mu_psi1   ~ dnorm(0, sd=2);  sigma_psi1  ~ dexp(1)
    mu_gamma  ~ dnorm(0, sd=2);  sigma_gamma ~ dexp(1)
    mu_eps    ~ dnorm(0, sd=2);  sigma_eps   ~ dexp(1)
    mu_p      ~ dnorm(0, sd=2);  sigma_p     ~ dexp(1)
    beta_tin_psi1   ~ dnorm(0, sd=2)
    beta_year_gamma ~ dnorm(0, sd=2)
    beta_year_eps   ~ dnorm(0, sd=2)
    beta_effort     ~ dnorm(0, sd=2)

    for (s in 1:n_total) {
      w[s]           ~ dbern(omega)
      alpha_psi1[s]  ~ dnorm(mu_psi1,  sd=sigma_psi1)
      alpha_gamma[s] ~ dnorm(mu_gamma, sd=sigma_gamma)
      alpha_eps[s]   ~ dnorm(mu_eps,   sd=sigma_eps)
      alpha_p[s]     ~ dnorm(mu_p,     sd=sigma_p)
    }

    for (s in 1:n_total) {
      for (i in 1:n_sites) {
        logit(psi1[s,i]) <- alpha_psi1[s] + beta_tin_psi1 * tin_1[s,i]
        mu_psi1n[s,i]    <- psi1[s,i] * w[s]
        z[s,i,1]         ~ dbern(mu_psi1n[s,i])
        for (t in 2:n_years) {
          logit(gamma[s,i,t-1]) <- alpha_gamma[s] + beta_year_gamma * year[t]
          logit(eps[s,i,t-1])   <- alpha_eps[s]   + beta_year_eps   * year[t]
          prob_occ[s,i,t] <- w[s] * (
            z[s,i,t-1] * (1 - eps[s,i,t-1]) +
            (1 - z[s,i,t-1]) * gamma[s,i,t-1]
          )
          z[s,i,t] ~ dbern(prob_occ[s,i,t])
        }
      }
    }

    for (k in 1:n_obs) {
      logit(p[k]) <- alpha_p[species_obs[k]] + beta_effort * effort_obs[k]
      mu[k]       <- z[species_obs[k], site_obs[k], year_idx_obs[k]] * p[k]
      y_obs[k]    ~ dbern(mu[k])
    }

    for (t in 1:n_years) {
      for (s in 1:n_total) {
        sp_present[s,t] <- w[s] * step(sum(z[s,1:n_sites,t]) - 0.5)
      }
      richness[t] <- sum(sp_present[1:n_total,t])
    }

    N_pool <- sum(w[1:n_total])
  })

  # ==========================================================================
  # SEZIONE 4 — COSTANTI, DATI, INITIAL VALUES
  # ==========================================================================
  nim_constants <- list(
    n_sites=n_sites, n_years=n_years, n_species=n_species,
    n_total=n_total, n_obs=nrow(obs_data),
    species_obs=obs_data$species_idx, site_obs=obs_data$site_idx,
    year_idx_obs=obs_data$year_idx, tin_1=tin_matrix, year=year_std
  )
  nim_data <- list(y_obs=obs_data$y_value, effort_obs=obs_data$effort_value)

  make_inits <- function(seed_val) {
    set.seed(seed_val)
    z_init <- array(0L, dim = c(n_total, n_sites, n_years))

    for (s in seq_len(n_species)) {
      y_mat <- y_list[[all_species[s]]]

      # Passo 1: z=1 ovunque sia stato rilevato (obbligatorio per consistenza)
      for (i in seq_len(n_sites))
        for (t in seq_len(n_years))
          if (!is.na(y_mat[i, t]) && y_mat[i, t] == 1L)
            z_init[s, i, t] <- 1L

      # Passo 2: back-propagazione — se z[t]=1, z[t-1] deve poter esserlo
      # (evita la situazione z[t-1]=0, gamma≈0, z[t]=1 → prob=0)
      for (i in seq_len(n_sites))
        for (t in n_years:2)
          if (z_init[s, i, t] == 1L && z_init[s, i, t - 1] == 0L)
            z_init[s, i, t - 1] <- 1L

      # Passo 3: anni rimanenti (non osservati, non back-propagati) → bassa prob
      for (i in seq_len(n_sites))
        for (t in seq_len(n_years))
          if (z_init[s, i, t] == 0L)
            z_init[s, i, t] <- rbinom(1, 1, 0.15)
    }
    # specie augmentate: z = 0 (w=0 → z forced to 0)

    list(
      omega         = runif(1, 0.7, 0.95),   # aspettiamo omega alto (≥10 specie non rilevate/68)
      w             = c(rep(1L, n_species), rep(0L, n_zero)),
      mu_psi1       = 0,    sigma_psi1  = 0.5,
      mu_gamma      = -1,   sigma_gamma = 0.5,   # gamma moderata ~0.27
      mu_eps        = -1,   sigma_eps   = 0.5,   # eps moderata ~0.27
      mu_p          = 0,    sigma_p     = 0.5,
      alpha_psi1    = rnorm(n_total,  0,  0.3),  # range stretto: psi ∈ (0.4, 0.6)
      alpha_gamma   = rnorm(n_total, -1,  0.3),
      alpha_eps     = rnorm(n_total, -1,  0.3),
      alpha_p       = rnorm(n_total,  0,  0.3),
      beta_tin_psi1   = 0,   # beta a 0: punto di partenza neutro
      beta_year_gamma = 0,
      beta_year_eps   = 0,
      beta_effort     = 0,
      z = z_init
    )
  }

  monitors_vec <- c("omega", "richness", "N_pool")

  # ==========================================================================
  # SEZIONE 5 — COMPILAZIONE
  # ==========================================================================
  log_msg("------------------------------------------------------")
  log_msg("Sezione 5 — Compilazione NIMBLE")

  log_msg("Fase 1/3 — build modello...")
  t0 <- proc.time()
  model   <- nimbleModel(code, nim_constants, nim_data, make_inits(42))
  lp_init <- model$calculate()
  log_msg("  ", round((proc.time()-t0)["elapsed"]/60,1), " min")
  log_msg("  Log-prob valori iniziali: ", round(lp_init, 1))
  if (!is.finite(lp_init))
    stop("Log-prob iniziale = -Inf o NaN: valori iniziali incompatibili. ",
         "Controlla effort_matrix e y_list.")

  log_msg("Fase 2/3 — compilazione C++...")
  t0 <- proc.time()
  c_model <- compileNimble(model, showCompilerOutput=FALSE)
  log_msg("  ", round((proc.time()-t0)["elapsed"]/60,1), " min")

  log_msg("Fase 3/3 — build e compilazione MCMC...")
  t0 <- proc.time()
  mcmc_conf <- configureMCMC(model, monitors=monitors_vec)
  mcmc      <- buildMCMC(mcmc_conf)
  c_mcmc    <- compileNimble(mcmc, project=model, showCompilerOutput=FALSE)
  log_msg("  ", round((proc.time()-t0)["elapsed"]/60,1), " min")

  # ==========================================================================
  # SEZIONE 6 — TEST RUN
  # ==========================================================================
  if (test_run) {
    log_msg("------------------------------------------------------")
    log_msg("Sezione 6 — Test run (50 iterazioni)")
    t0 <- proc.time()
    tryCatch({
      runMCMC(c_mcmc, niter=50, nburnin=0, thin=1,
              samplesAsCodaMCMC=FALSE, progressBar=TRUE)
      elapsed      <- (proc.time()-t0)["elapsed"]
      iter_per_min <- round(50/elapsed*60, 1)
      stima_ore    <- round((niter*n_chains)/(iter_per_min*60), 1)
      log_msg("  iter/min: ", iter_per_min)
      log_msg("  Stima tempo totale: ~", stima_ore, " ore")
      if (stima_ore > 48)
        log_msg("  ATTENZIONE: run >48h. ",
                "Considera n_sites_target < ", n_sites_target,
                " o marginalizzazione di z.")
    }, error = function(e) {
      log_msg("  ERRORE: ", conditionMessage(e))
      stop("Test run fallito.")
    })
  }

  # ==========================================================================
  # SEZIONE 7 — CATENE SEQUENZIALI CON CHECKPOINT
  # ==========================================================================
  log_msg("------------------------------------------------------")
  log_msg("Sezione 7 — Campionamento MCMC")
  log_msg("  n_chains=", n_chains, " | niter=", niter,
          " | nburnin=", nburnin, " | thin=", thin)
  log_msg("  Campioni per catena: ", (niter-nburnin)/thin)

  samples_list <- vector("list", n_chains)

  for (chain_id in seq_len(n_chains)) {

    ckpt_file <- file.path(checkpoint_dir, paste0("chain_", chain_id, ".rds"))

    if (file.exists(ckpt_file)) {
      log_msg("Catena ", chain_id, "/", n_chains,
              " — checkpoint trovato, carico da disco")
      samples_list[[chain_id]] <- readRDS(ckpt_file)
      next
    }

    log_msg("Catena ", chain_id, "/", n_chains, " — avvio")
    t_chain <- proc.time()
    model$setInits(make_inits(seed_val = 1000 + chain_id*123))

    tryCatch({
      chain_samples <- runMCMC(
        c_mcmc, niter=niter, nburnin=nburnin, thin=thin,
        samplesAsCodaMCMC=TRUE, progressBar=TRUE
      )
      elapsed_chain <- round((proc.time()-t_chain)["elapsed"]/60, 1)
      log_msg("  completata in ", elapsed_chain, " min")
      samples_list[[chain_id]] <- chain_samples
      saveRDS(chain_samples, ckpt_file)
      log_msg("  checkpoint: ", ckpt_file)
    }, error = function(e) {
      log_msg("  ERRORE: ", conditionMessage(e))
      stop("Errore catena ", chain_id, ". Riesegui per riprendere.")
    })
  }

  # ==========================================================================
  # SEZIONE 8 — RISULTATI
  # ==========================================================================
  log_msg("------------------------------------------------------")
  log_msg("Sezione 8 — Risultati")

  samples     <- as.mcmc.list(Filter(Negate(is.null), samples_list))
  samples_mat <- as.matrix(samples)

  richness_samples <- samples_mat[, grep("^richness\\[", colnames(samples_mat)),
                                   drop=FALSE]
  omega_samples  <- samples_mat[, "omega"]
  N_pool_samples <- samples_mat[, "N_pool"]

  richness_summary <- data.frame(
    Anno           = anni_num,
    richness_mean  = colMeans(richness_samples),
    richness_lower = apply(richness_samples, 2, quantile, probs=0.025),
    richness_upper = apply(richness_samples, 2, quantile, probs=0.975),
    richness_sd    = apply(richness_samples, 2, sd)
  )

  omega_summary <- c(mean  = mean(omega_samples),
                     lower = quantile(omega_samples, 0.025),
                     upper = quantile(omega_samples, 0.975))

  N_pool_summary <- c(mean  = mean(N_pool_samples),
                      lower = quantile(N_pool_samples, 0.025),
                      upper = quantile(N_pool_samples, 0.975))

  rhat_df <- tryCatch(
    as.data.frame(gelman.diag(samples, multivariate=FALSE)$psrf),
    error = function(e) data.frame(Point.est=NA, Upper.C.I.=NA)
  )

  if (omega_summary["upper.97.5%"] > 0.97)
    log_msg("ATTENZIONE: omega > 0.97 → n_zero troppo piccolo. ",
            "Ripeti con n_zero = ", n_zero + 10)

  log_msg(sprintf("omega:  %.3f [%.3f - %.3f]",
                  omega_summary["mean"],
                  omega_summary["lower.2.5%"],
                  omega_summary["upper.97.5%"]))
  log_msg(sprintf("N_pool: %.1f [%.1f - %.1f]  (osservate: %d | stimate non rilevate: %.1f)",
                  N_pool_summary["mean"],
                  N_pool_summary["lower.2.5%"],
                  N_pool_summary["upper.97.5%"],
                  n_species,
                  N_pool_summary["mean"] - n_species))

  log_msg("Richness per anno:")
  for (i in seq_len(nrow(richness_summary)))
    log_msg(sprintf("  %d: %.1f [%.1f - %.1f]",
                    richness_summary$Anno[i],
                    richness_summary$richness_mean[i],
                    richness_summary$richness_lower[i],
                    richness_summary$richness_upper[i]))

  # ==========================================================================
  # SEZIONE 9 — PLOT
  # ==========================================================================
  out_prefix <- file.path(out_dir, paste0("RICHNESS_", today()))

  # 9a. richness
  pRichness <- richness_summary %>%
    filter(Anno > 2004) %>%
    ggplot(aes(x=Anno, y=richness_mean)) +
    geom_ribbon(aes(ymin=richness_lower, ymax=richness_upper),
                alpha=0.25, fill="steelblue") +
    geom_line(linewidth=1.4, color="steelblue4") +
    geom_point(size=2.5, color="steelblue4") +
    geom_hline(yintercept=n_species, linetype="dashed",
               color="firebrick", linewidth=0.8) +
    annotate("text",
             x=min(richness_summary$Anno[richness_summary$Anno > 2004]),
             y=n_species+0.8,
             label=paste0("Osservate: ", n_species),
             hjust=0, color="firebrick", size=3.5) +
    labs(
      title    = paste0("Richness stimata — MSOM (n_zero=", n_zero, ")"),
      subtitle = sprintf(
        "omega=%.3f [%.3f-%.3f]  |  N_pool=%.1f [%.1f-%.1f]  |  siti=%d",
        omega_summary["mean"],
        omega_summary["lower.2.5%"], omega_summary["upper.97.5%"],
        N_pool_summary["mean"],
        N_pool_summary["lower.2.5%"], N_pool_summary["upper.97.5%"],
        n_sites
      ),
      x="Anno", y="Numero di specie"
    ) +
    theme_bw(base_size=14)
  print(pRichness)
  ggsave(paste0(out_prefix, ".pdf"), plot=pRichness, width=10, height=6)

  # 9b. diagnostico selezione siti
  site_meta_plot <- site_meta %>%
    mutate(tipo = case_when(
      cellcodeX %in% sites_forced     ~ "Forzato (specie rare)",
      cellcodeX %in% sites_stratified ~ "Stratificato",
      TRUE                            ~ "Escluso"
    ),
    tipo = factor(tipo, levels=c("Forzato (specie rare)",
                                  "Stratificato", "Escluso")))

  pSiti <- ggplot(site_meta_plot,
                  aes(x=mean_tin, y=effort, color=tipo, shape=tipo)) +
    geom_point(alpha=0.7, size=2) +
    scale_color_manual(values=c("firebrick","steelblue4","gray70")) +
    scale_shape_manual(values=c(17, 16, 1)) +
    labs(
      title    = "Selezione siti: stratificata + specie rare forzate",
      subtitle = sprintf(
        "%d selezionati / %d totali  |  %d forzati + %d stratificati",
        n_sites_sel, nrow(site_meta),
        length(sites_forced), length(sites_stratified)
      ),
      x="Gradiente termico (tin_1 medio)",
      y="Effort (anni con survey)",
      color=NULL, shape=NULL
    ) +
    theme_bw(base_size=13) +
    theme(legend.position="bottom")
  print(pSiti)
  ggsave(paste0(out_prefix, "_siti.pdf"), plot=pSiti, width=10, height=6)

  # ==========================================================================
  # SEZIONE 10 — SALVATAGGIO
  # ==========================================================================
  saveRDS(
    list(
      n_species=n_species, n_zero=n_zero, n_total=n_total,
      species=all_species, n_sites=n_sites, n_years=n_years,
      site_selection=site_selection_meta,
      richness_summary=richness_summary, richness_samples=richness_samples,
      omega_summary=omega_summary, N_pool_summary=N_pool_summary,
      samples=samples, rhat=rhat_df
    ),
    paste0(out_prefix, ".rds")
  )
  log_msg("Salvato: ", out_prefix, ".rds / .pdf / _siti.pdf")
  log_msg("Analisi completata.")

  return(list(
    richness_summary = richness_summary,
    richness_samples = richness_samples,
    omega_summary    = omega_summary,
    N_pool_summary   = N_pool_summary,
    site_selection   = site_selection_meta,
    rhat             = rhat_df
  ))
}
