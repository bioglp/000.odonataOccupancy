set.seed(32)
library(dplyr)

set.seed(32)

df <- expand.grid(
  species = species_list,
  date_year = years,
  site_id = 1:max(effort_df$n_sites_censiti)
) %>%
  rowwise() %>%
  filter(site_id <= effort_df$n_sites_censiti[effort_df$date_year == date_year]) %>%
  ungroup() %>%
  mutate(
    prob = case_when(
      species == "SpecieA" ~ runif(n(), 0.3, 0.7),  # trend neutro
      species == "SpecieB" ~ runif(n(), 0.4, 0.9),  # trend positivo
      species == "SpecieC" ~ pmax(0.05, 0.8 - 0.07*(date_year - 2000))  # trend negativo
    ),
    observed = rbinom(n(), 1, prob)
  ) %>%
  filter(observed == 1) %>%
  select(-observed, -prob) %>%
  rename(cellcode = site_id)

# Risultato
head(df)

# Calcolo del numero di siti occupati dalle specie
species_sites <- df %>%
  group_by(species, date_year) %>%
  summarise(n_sites = n_distinct(cellcode), .groups = "drop")

species_sites %>% 
  ggplot(aes(x = date_year, y = n_sites)) +
  geom_smooth(method='lm',col='coral')+
  stat_regline_equation()+
  geom_line(size = 1) +
  geom_point() +
  facet_grid(~species)+
  labs(x = "", 
       y = "No. sites")+
  theme(plot.title = element_text(face = 'italic'))  

# Calcolo del numero di siti occupati dalle specie x effort
species_sites <- df %>%
  group_by(species, date_year) %>%
  summarise(n_sites = n_distinct(cellcode), .groups = "drop")

effort_df <- df %>% 
  group_by(date_year) %>%
  distinct(cellcode) %>% 
  count(date_year, name = 'n_sites_censiti')

species_sites <- species_sites %>%
  left_join(effort_df, by = "date_year")

results <- species_sites %>%
  group_by(species) %>%
  do({
    mod <- glm(n_sites ~ date_year, family = poisson,
               offset = log(n_sites_censiti), data = .)
    tidy(mod)  # restituisce coefficienti, SE, p-value
  }) %>%
  ungroup()

species_trends <- results %>%
  filter(term == "date_year") %>%
  mutate(trend = case_when(
    estimate > 0 ~ "increase",
    estimate < 0 ~ "decrease",
    TRUE ~ "stable"
  ))
species_trends

predicted_sites <- species_sites %>%
  group_by(species) %>%
  do({
    mod <- glm(n_sites ~ date_year, family = poisson,
               offset = log(n_sites_censiti), data = .)
    data.frame(date_year = .$date_year,
               pred = predict(mod, type = "response"))
  }) %>%
  ungroup()

# Grafico linee
predicted_sites %>% 
  ggplot(aes(x = date_year, y = pred)) +
  geom_area(fill = 'coral', size = 1) +   # linea dei predetti dal GLM
  geom_point(aes(y = n_sites), data = species_sites, size = 2) +  # punti osservati reali
  labs(x = "Anno", y = "Numero di siti (corretto per sforzo)") +
  facet_wrap(~species) +
  theme(plot.title = element_text(face = 'italic'))
