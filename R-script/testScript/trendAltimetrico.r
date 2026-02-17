# analisi trend altimetrico per le specie

glimpse(dfo)


dfo |>
  filter(code == 'Continental') |>
  group_by(cellcode, date_year) |>
  summarise(mean_tin = mean(tin_1, na.rm = T)) |>
  ggplot(aes(mean_tin, date_year, group = date_year)) +
  geom_violin(fill = 'green') +
  scale_x_sqrt()

dfo |>
  filter(code == 'Continental') |>
  group_by(species, date_year) |>
  summarise(mean_tin = mean(tin_1, na.rm = T)) |>
  ggplot(aes(mean_tin, species)) +
  geom_violin(fill = 'green') +
  scale_x_sqrt()


# Statistiche per anno
annual_summary <- dfo %>%
  filter(code == 'Continental') |>
  group_by(date_year) %>%
  summarise(
    mean_tin = mean(tin_1, na.rm = TRUE),
    median_tin = median(tin_1, na.rm = TRUE),
    sd_tin = sd(tin_1, na.rm = TRUE),
    n = n()
  )

ggplot(annual_summary, aes(x = date_year, y = mean_tin)) +
  geom_point() +
  geom_line() +
  geom_errorbar(
    aes(ymin = mean_tin - sd_tin, ymax = mean_tin + sd_tin),
    width = 0.2
  ) +
  geom_smooth(method = "lm", se = TRUE, color = "red") +
  labs(title = "Media annuale", x = "Anno", y = "tin_1") +
  theme_minimal()


# Test di Mann-Kendall per trend monotono
library(Kendall)

# Trend globale
mk_global <- MannKendall(annual_summary$mean_tin)
summary(mk_global)

# Per specie
mk_by_species <- dfo %>%
  filter(code == 'Continental') |>
  group_by(species) %>%
  summarise(
    tau = cor(date_year, tin_1, method = "kendall"),
    p_value = cor.test(date_year, tin_1, method = "kendall")$p.value,
    n = n()
  ) %>%
  arrange(desc(abs(tau)))

print(mk_by_species)

# Visualizza
ggplot(
  mk_by_species,
  aes(x = reorder(species, tau), y = tau, color = p_value < 0.05)
) +
  geom_point(size = 3) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() +
  scale_color_manual(
    values = c("grey", "red"),
    name = "p < 0.05"
  ) +
  labs(
    title = "",
    x = "Specie",
    y = "Tau di Kendall"
  ) +
  theme_minimal()


# glm
library(lme4)
library(lmerTest)

data_scaled <- dfo %>%
  filter(code == 'Continental') |>
  mutate(year_sc = as.numeric(scale(date_year)))

# Modello con intercetta random
m1 <- lmer(tin_1 ~ year_sc + (1 | species), data = data_scaled)

# Modello con slope random (permette trend diversi per specie)
m2 <- lmer(tin_1 ~ year_sc + (year_sc | species), data = data_scaled)

# Confronta modelli
anova(m1, m2)
# Se m2 è migliore (AIC più basso), le specie hanno trend diversi

# Risultati
summary(m2)

# Effetti fissi (trend medio)
fixef(m2)

# Effetti random (deviazione per specie)
ranef(m2)$species

# Coefficienti specifici per specie
coef(m2)$species

# Predizioni del modello
data_scaled$fitted <- predict(m2)

ggplot(data_scaled, aes(x = date_year, y = tin_1, color = species)) +
  geom_point(alpha = 0.3) +
  geom_line(aes(y = fitted), linewidth = 1) +
  facet_wrap(~species, scales = "free") +
  labs(title = "Trend per specie (modello mixed)", x = "Anno", y = "tin_1") +
  theme_minimal() +
  theme(legend.position = "none")

# Tabella coefficienti
tidy(m2, effects = "fixed", conf.int = TRUE)

# Trend per specie
species_trends <- coef(m2)$species %>%
  rownames_to_column("species") %>%
  arrange(desc(year_sc))

print(species_trends)

# Ordina per trend
species_coefs %>%
  rownames_to_column("species") %>%
  arrange(desc(year_sc))

# Identifica specie con shift forte (es. |year_sc| > 2)
species_coefs %>%
  rownames_to_column("species") %>%
  mutate(
    shift_direction = case_when(
      year_sc > 2 ~ "Shift verso l'alto (forte)",
      year_sc > 0.5 ~ "Shift verso l'alto (moderato)",
      year_sc < -2 ~ "Shift verso il basso (forte)",
      year_sc < -0.5 ~ "Shift verso il basso (moderato)",
      TRUE ~ "Stabile"
    )
  ) %>%
  arrange(desc(year_sc))
