# funzione occuProb
# per il calcolo dell'occupancy
# il campo cellcodeX è definito a monte

occuProb <- function(nomeSp) {
  # nomeSp
  cat(nomeSp, '\n')

  # matrice di rilevamento
  y_matrix <- dfo |>
    filter(species == nomeSp) |>
    count(cellcodeX, date_year, name = 'n_occurrences') |>
    # Join con effort
    right_join(
      dfo |>
        count(cellcodeX, date_year, name = 'tot_occurrences'),
      by = c("cellcodeX", "date_year")
    ) |>
    complete(
      cellcodeX = all_cells,
      date_year = all_years,
      fill = list(n_occurrences = NA, tot_occurrences = NA)
    ) |>
    mutate(
      detection = case_when(
        is.na(tot_occurrences) ~ NA_real_, # Cella non indagata
        n_occurrences > 0 ~ 1, # Specie trovata
        tot_occurrences > 0 ~ 0, # Indagata ma specie non trovata
        TRUE ~ NA_real_ # Cella non indagata
      )
    ) |>
    pivot_wider(
      id_cols = cellcodeX,
      names_from = date_year,
      values_from = detection
    ) |>
    column_to_rownames("cellcodeX") |>
    as.matrix()

  # Standardizzazione effort
  effort_vec <- as.numeric(effort_matrix)
  effort_std_mat <- matrix(
    (effort_vec - mean(effort_vec, na.rm = TRUE)) /
      sd(effort_vec, na.rm = TRUE),
    nrow = nrow(y_matrix),
    ncol = ncol(y_matrix)
  )

  # Covariata Anno
  anni_num <- as.numeric(colnames(y_matrix))
  year_std <- as.vector(scale(anni_num))
  n_siti <- nrow(y_matrix)
  n_anni <- ncol(y_matrix)

  # matrice degli anni (standardizzati) con dimensioni siti x anni
  year_mat <- matrix(rep(year_std, each = n_siti), nrow = n_siti, ncol = n_anni)

  # tin_1 medio per cella 10k
  cells_in_y <- rownames(y_matrix)
  tin_data <- dfo |>
    select(cellcodeX, tin_1) |>
    distinct() |>
    mutate(cellcodeX = as.character(cellcodeX)) |>
    filter(cellcodeX %in% cells_in_y) |>
    group_by(cellcodeX) |>
    summarise(tin_1 = mean(round(tin_1), rm.na = T))

  # Crea un dataframe con tutte le celle nell'ordine corretto
  tin_df <- data.frame(cellcodeX = cells_in_y) |>
    left_join(tin_data, by = "cellcodeX")

  # Standardizzazione di tin_1
  tin_std <- scale(tin_df$tin_1)

  # 2. creazione umf
  umf <- unmarkedMultFrame(
    y = y_matrix,
    siteCovs = data.frame(tin_1 = as.vector(tin_std)),
    yearlySiteCovs = list(
      effort = effort_std_mat,
      year_idx = year_mat
    ),
    numPrimary = n_anni
  )

  # 3. Esecuzione del modello
  m1 <- colext(
    psiformula = ~tin_1,
    gammaformula = ~ year_idx + tin_1,
    epsilonformula = ~year_idx,
    pformula = ~effort,
    data = umf
  )

  sm <- smoothed(m1)
  occ_values <- as.numeric(sm["occupied", ])

  # Calcola SE tramite parametric bootstrap
  set.seed(32)
  pb <- parboot(
    m1,
    statistic = function(fm) smoothed(fm)[2, ],
    nsim = 700,
    parallel = T
  )

  occ_lower <- apply(pb@t.star, 2, quantile, probs = 0.025, na.rm = TRUE)
  occ_upper <- apply(pb@t.star, 2, quantile, probs = 0.975, na.rm = TRUE)

  df_plot <- data.frame(
    Species = nomeSp,
    Anno = anni_num,
    Occupancy = occ_values,
    Lower = occ_lower,
    Upper = occ_upper
  )

  pOccu <- df_plot |>
    filter(Anno > 2004) |>
    ggplot(aes(x = Anno, y = Occupancy)) +
    geom_ribbon(
      aes(ymin = Lower, ymax = Upper),
      fill = "steelblue",
      alpha = 0.2
    ) +
    geom_line(color = "steelblue", linewidth = 1) +
    geom_point(color = "steelblue", size = 3) +
    # geom_smooth(
    #   method = "loess",
    #   color = "red",
    #   linetype = "dashed",
    #   se = F,
    #   alpha = 0.3
    # ) +
    scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
    labs(
      title = nomeSp,
      subtitle = paste("Convergence:", round(m1@opt$convergence, 1)),
      x = "",
      y = "Occupancy",
      caption = "unmarked"
    ) +
    theme_bw(base_size = 14) +
    theme(plot.title = element_text(face = "italic"))

  print(pOccu)

  ggsave(
    paste0('output_unmarked/', nomeSp, '_', today(), '.pdf'),
    width = 9,
    height = 6
  )

  saveRDS(
    list(
      species = nomeSp,
      range = range(df_plot$Occupancy),
      sm = sm,
      df_plot = df_plot
    ),
    paste0('output_unmarked/', nomeSp, "_", today(), ".rds")
  )

  return(range(df_plot$Occupancy))
}
