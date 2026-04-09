# stat per specie
# libs
library(pacman)
p_load(tidyverse, rstatix, scales)
theme_set(theme_bw(base_size = 14))

# Lista delle specie
rds_files <- list.files(
  "output/output_bayes",
  pattern = "\\.rds$",
  full.names = TRUE
)

# combina in una tabella
sp_summary <- map_dfr(rds_files, function(file) {
  tryCatch(
    {
      spRDS <- readRDS(file)

      sp <- spRDS$plot_data
      sp <- sp[sp$Anno != 2004, ]

      sp_stat <- sp |>
        get_summary_stats(Occupancy, type = "common")

      sp_stat$species <- spRDS$species
      sp_stat$y2005 <- round(sp$Occupancy[1], 3)
      sp_stat$y2025 <- round(sp$Occupancy[21], 3)
      sp_stat$delta <- round(sp_stat$y2025 - sp_stat$y2005, 3)
      sp_stat
    },
    error = function(e) {
      message("Errore nel file: ", file, " -> ", e$message)
      NULL # salta il file problematico
    }
  )
})

# Riordina con species come prima colonna
sp_summary <- sp_summary |> relocate(species)

sp_summary
write.csv(
  sp_summary,
  paste0("output/output_bayes/sp_summaryBayes_", today(), ".csv")
)

# grafico con i range
sp_summary <- sp_summary |>
  mutate(
    trend = y2025 - y2005,
    direzione = ifelse(trend > 0, "Aumento", "Diminuzione")
  )

sp_summary |>
  ggplot(aes(x = reorder(species, median), color = direzione)) +
  geom_linerange(
    aes(ymin = min, ymax = max),
    linewidth = 0.8,
    color = "gray70"
  ) +
  geom_linerange(
    aes(ymin = min, ymax = max),
    linewidth = 1
  ) +
  geom_point(aes(y = median), size = 3) +
  coord_flip() +
  scale_color_manual(
    values = c("Aumento" = "#2196F3", "Diminuzione" = "#F44336")
  ) +
  scale_x_discrete(labels = function(x) {
    parse(text = paste0("italic('", x, "')"))
  }) +
  labs(
    x = NULL,
    y = "\nOccupancy"
  ) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )

ggsave(
  paste0("output/output_bayes/sp_summaryBayes_", today(), ".pdf"),
  width = 8,
  height = 12
)

#####################################################################

# Lista delle specie
rds_files <- list.files(
  "output/output_unmarked",
  pattern = "\\.rds$",
  full.names = TRUE
)

# combina in una tabella
sp_summary <- map_dfr(rds_files, function(file) {
  tryCatch(
    {
      spRDS <- readRDS(file)

      sp <- spRDS$df_plot
      sp <- sp[sp$Anno != 2004, ]

      sp_stat <- sp |>
        get_summary_stats(Occupancy, type = "common")

      sp_stat$species <- spRDS$species
      sp_stat$y2005 <- round(sp$Occupancy[1], 3)
      sp_stat$y2025 <- round(sp$Occupancy[21], 3)
      sp_stat$delta <- round(sp_stat$y2025 - sp_stat$y2005, 3)
      sp_stat
    },
    error = function(e) {
      message("Errore nel file: ", file, " -> ", e$message)
      NULL # salta il file problematico
    }
  )
})

# Riordina con species come prima colonna
sp_summary <- sp_summary |> relocate(species)

sp_summary
write.csv(
  sp_summary,
  paste0("output/output_unmarked/sp_summaryUnmarked_", today(), ".csv")
)

# grafico con i range
sp_summary <- sp_summary |>
  mutate(
    trend = y2025 - y2005,
    direzione = ifelse(trend > 0, "Aumento", "Diminuzione")
  )

sp_summary |>
  ggplot(aes(x = reorder(species, median), color = direzione)) +
  geom_linerange(
    aes(ymin = min, ymax = max),
    linewidth = 0.8,
    color = "gray70"
  ) +
  geom_linerange(
    aes(ymin = min, ymax = max),
    linewidth = 1
  ) +
  geom_point(aes(y = median), size = 3) +
  coord_flip() +
  scale_color_manual(
    values = c("Aumento" = "#2196F3", "Diminuzione" = "#F44336")
  ) +
  scale_x_discrete(labels = function(x) {
    parse(text = paste0("italic('", x, "')"))
  }) +
  labs(
    x = NULL,
    y = "\nOccupancy"
  ) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )

ggsave(
  paste0("output/output_unmarked/sp_summaryUnmarked_", today(), ".pdf"),
  width = 8,
  height = 12
)
