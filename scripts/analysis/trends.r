# dumbell plot + bar chart combinati
library(tidyverse)
library(patchwork)
library(rstatix)

dfb <- read.csv('output/tables/risultatiBayes.csv')
dfu <- read.csv('output/tables/risultatiUnmarked.csv')

# correlazione
cor(dfb$tau, dfu$tau)
data.frame(
  species = dfb$species,
  tu = dfu$tau,
  tb = dfb$tau,
  sensu = dfu$sens,
  sensb = dfb$sens
)

spb_stable <- dfb |> filter(trend == 'stable')
spu_stable <- dfu |> filter(trend == 'stable')
spb_stable[spb_stable$species %in% spu_stable$species, ]

# differenze nella valutazione del trend
dfb |>
  count(trend)
dfu |>
  count(trend)

# starting point
mean(dfb$occupancy_2005 - dfu$occupancy_2005)
sd(dfb$occupancy_2005 - dfu$occupancy_2005)

# end point
mean(dfb$occupancy_2025 - dfu$occupancy_2025)
sd(dfb$occupancy_2025 - dfu$occupancy_2025)

# increasing species -- bayes
df <- dfb
df |>
  filter(trend == 'increase') |>
  dplyr::select(species, occupancy_2005, occupancy_2025, change, tau, sens) |>
  arrange(change) |>
  get_summary_stats()

# decreasing species -- bayes
df |>
  filter(trend == 'decrease') |>
  dplyr::select(species, occupancy_2005, occupancy_2025, change, tau, sens) |>
  arrange(change) |>
  get_summary_stats()


df <- df %>%
  mutate(
    change = occupancy_2025 - occupancy_2005,
    trend = factor(trend, levels = c("decrease", "stable", "increase")),
    trend_order = as.integer(trend),
    species = reorder(species, trend_order * 1000 + change)
  )

trend_colours <- c(
  "decrease" = "#D55E00",
  "stable" = "#999999",
  "increase" = "#009E73"
)

# dumbbell (occupancy assoluta 2005 → 2025) ---
p_dumb <- df |>
  ggplot(aes(y = species)) +
  geom_segment(
    aes(x = occupancy_2005, xend = occupancy_2025, colour = trend),
    linewidth = 0.8,
    show.legend = FALSE
  ) +
  geom_point(
    aes(x = occupancy_2005, colour = trend),
    shape = 21,
    fill = "white",
    size = 3,
    stroke = 1.1
  ) +
  geom_point(
    aes(x = occupancy_2025, colour = trend),
    shape = 16,
    size = 3
  ) +
  scale_x_continuous(
    labels = scales::percent_format(accuracy = 1),
    limits = c(0, 1),
    expand = expansion(mult = c(0.02, 0.05))
  ) +
  scale_colour_manual(
    values = trend_colours,
    name = "Trend",
    labels = c("Decreasing", "Stable", "Increasing")
  ) +
  labs(x = "\nOccupancy (2005 - 2025)", y = NULL) +
  theme_bw(base_size = 13) +
  theme(
    legend.position = 'top',
    axis.text.y = element_text(face = "italic"),
    legend.title = element_text(size = 10),
    panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.4)
  )

# barre (delta occupancy)
p_bar <- df |>
  mutate(change = ifelse(trend == "stable", 0, change)) |>
  ggplot(aes(y = species, x = change, fill = trend)) +
  geom_col() +
  geom_vline(xintercept = 0, linewidth = 0.4, colour = "grey40") +
  scale_x_continuous(
    labels = scales::percent_format(accuracy = 1),
    expand = expansion(mult = c(0.05, 0.05))
  ) +
  scale_fill_manual(
    values = trend_colours,
    name = "Trend",
    labels = c("Decreasing", "Stable", "Increasing")
  ) +
  labs(x = "\n\u0394 Occupancy (2025 - 2005)", y = NULL) +
  theme_bw(base_size = 13) +
  theme(
    legend.position = 'none',
    legend.title = element_text(size = 10),
    panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.4),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

p_combined <- p_dumb + p_bar + plot_layout(widths = c(7, 3))

p_combined

ggsave(
  'output/figures/trendB_combined.pdf',
  plot = p_combined,
  width = 13,
  height = 11
)
