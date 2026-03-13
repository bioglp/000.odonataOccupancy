# libs
library(pacman)
p_load(tidyverse)
theme_set(theme_bw(base_size = 14))

# creazione plot con riferimento di partenza al 100%

# Temp Index
ti <- readRDS(
  'output/output_bayes/multispecies/MULTISPECIES_ALL_GROUPS_Tindex_2026-02-17.rds'
)

ti1 <- ti$group_results[[1]]$plot_data
ti2 <- ti$group_results[[2]]$plot_data
ti3 <- ti$group_results[[3]]$plot_data
ti4 <- ti$group_results[[4]]$plot_data

ti_list <- lapply(ti$group_results, function(g) {
  df <- g$plot_data
  df <- df[df$Anno != 2004, ]
  base_value <- df$Occupancy[df$Anno == 2005]
  df$Occupancy <- df$Occupancy / base_value * 100
  df$Lower <- df$Lower / base_value * 100
  df$Upper <- df$Upper / base_value * 100
  df
})

ti_all <- do.call(rbind, ti_list)
write.csv(
  ti_all,
  paste0('output/output_bayes/multispecies/TI100_', today(), ".csv")
)

# Plot
ti_all |>
  #filter(Group!=4) |>
  ggplot(
    aes(x = Anno, y = Occupancy, color = factor(Group), fill = factor(Group))
  ) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.2, color = NA) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  geom_hline(yintercept = 100, linetype = "dashed", color = "gray40") +
  scale_color_viridis_d(labels = paste("Group", 1:4)) +
  scale_fill_viridis_d(labels = paste("Group", 1:4)) +
  scale_y_log10() +
  labs(
    x = "\nYear",
    y = "Occupancy (2005 = 100)",
    color = "TI",
    fill = "TI"
  ) +
  theme(
    legend.position = "top"
  )

ggsave(
  paste0("output/output_bayes/multispecies/TI100_", today(), ".pdf"),
  width = 12,
  height = 8
)

# Habitat Index
ti <- readRDS(
  'output/output_bayes/multispecies/MULTISPECIES_ALL_GROUPS_HV2026-03-02.rds'
)

ti_list <- lapply(ti$group_results, function(g) {
  df <- g$plot_data
  df <- df[df$Anno != 2004, ]
  base_value <- df$Occupancy[df$Anno == 2005]
  df$Occupancy <- df$Occupancy / base_value * 100
  df$Lower <- df$Lower / base_value * 100
  df$Upper <- df$Upper / base_value * 100
  df
})

ti_all <- do.call(rbind, ti_list)
write.csv(
  ti_all,
  paste0('output/output_bayes/multispecies/HV100_', today(), ".csv")
)

# Plot
ti_all |>
  #filter(Group!=4) |>
  ggplot(
    aes(x = Anno, y = Occupancy, color = factor(Group), fill = factor(Group))
  ) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.2, color = NA) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  geom_hline(yintercept = 100, linetype = "dashed", color = "gray40") +
  scale_color_viridis_d(labels = paste("Group", 1:4)) +
  scale_fill_viridis_d(labels = paste("Group", 1:4)) +
  scale_y_sqrt() +
  labs(
    x = "\nYear",
    y = "Occupancy (2005 = 100)",
    color = "HV",
    fill = "HV"
  ) +
  theme(
    legend.position = "top"
  )


ggsave(
  paste0("output/output_bayes/multispecies/HV100_", today(), ".pdf"),
  width = 12,
  height = 8
)

# IW Index
ti <- readRDS(
  'output/output_bayes/multispecies/MULTISPECIES_ALL_GROUPS_IW2026-02-22.rds'
)

ti_list <- lapply(ti$group_results, function(g) {
  df <- g$plot_data
  df <- df[df$Anno != 2004, ]
  base_value <- df$Occupancy[df$Anno == 2005]
  df$Occupancy <- df$Occupancy / base_value * 100
  df$Lower <- df$Lower / base_value * 100
  df$Upper <- df$Upper / base_value * 100
  df
})

ti_all <- do.call(rbind, ti_list)
write.csv(
  ti_all,
  paste0('output/output_bayes/multispecies/IW100_', today(), ".csv")
)

# Plot
ti_all |>
  #filter(Group!=4) |>
  ggplot(
    aes(x = Anno, y = Occupancy, color = factor(Group), fill = factor(Group))
  ) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.2, color = NA) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  geom_hline(yintercept = 100, linetype = "dashed", color = "gray40") +
  scale_color_viridis_d(labels = paste("Group", 1:5)) +
  scale_fill_viridis_d(labels = paste("Group", 1:5)) +
  scale_y_sqrt() +
  labs(
    x = "\nYear",
    y = "Occupancy (2005 = 100)",
    color = "IW",
    fill = "IW"
  ) +
  theme(
    legend.position = "top"
  )

ggsave(
  paste0("output/output_bayes/multispecies/IW100_", today(), ".pdf"),
  width = 12,
  height = 8
)
