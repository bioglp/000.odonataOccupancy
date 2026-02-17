library(dplyr)
library(tidyr)
library(unmarked)

set.seed(123)

# Simuliamo un dataset simile al tuo
df_occ <- expand.grid(
  cellcode = paste0("C", 1:100),     # 100 celle
  date_year = 2015:2020,             # 6 anni
  species = c("Rana_temporaria")     # una sola specie
) %>%
  mutate(
    n_occurrences = rpois(n(), lambda = runif(1, 0.2, 1.2))
  ) %>%
  mutate(presence = ifelse(n_occurrences > 0, 1, 0))

y <- df_occ %>%
  select(cellcode, date_year, presence) %>%
  pivot_wider(names_from = date_year, values_from = presence, values_fill = 0) %>%
  arrange(cellcode)

# Matrice di rilevamenti (celle × anni)
y_mat <- as.matrix(y[,-1])
rownames(y_mat) <- y$cellcode

head(y_mat)
umf <- unmarkedFrameOccu(y = y_mat)






