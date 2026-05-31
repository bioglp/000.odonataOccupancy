# UPDATE 2026-05-31

# libraries ----

library(pacman)
p_load(
  tidyverse,
  gt,
  skimr,
  ggpubr,
  rstatix,
  openxlsx,
  scales,
  patchwork,
  parallel
)
theme_set(theme_bw(base_size = 14))


# import data ----
year_start <- 2004
dfo <- read_csv('data/raw/datasetSegnalazioni_R20251106.csv')
dfo <- dfo |>
  rename_with(tolower) |>
  filter(date_year >= year_start) |>
  mutate(total_count = ifelse(is.na(total_count), 1, total_count))
dim(dfo)

dfo0 <- dfo |>
  filter(
    date_year < 2026,
    total_count == 0,
    code == 'Continental',
    !is.na(cellcode)
  )
table(dfo0$latin_taxon)
dim(dfo0)

dfo <- dfo |>
  filter(
    date_year < 2026,
    total_count != 0,
    code == 'Continental',
    !is.na(cellcode)
  )
dfo <- dfo |>
  rename(species = latin_taxon)
glimpse(dfo) #View(dfo)

dfo <- dfo |>
  select(date_year, species, cellcode, cellcode20, tin_1)

# iNaturalist
iNat <- read_csv('data/raw/iNaturalist_20241123.csv')
iNat <- iNat |>
  select(
    date_year = DATE_YEAR,
    species = scientific_name,
    cellcode,
    cellcode20,
    tin_1 = TIN_11
  )

dfo <- rbind(dfo, iNat)
dim(dfo)

# rinominare specie ----
# dfo |> distinct(species) |> arrange() |> pull()
dfo$species[dfo$species == 'Sympetrum fonscolombei'] <- 'Sympetrum fonscolombii'
dfo$species[
  dfo$species == 'Calopteryx splendens caprai'
] <- 'Calopteryx splendens'
dfo$species[
  dfo$species == 'Calopteryx splendens splendens'
] <- 'Calopteryx splendens'
dfo$species[dfo$species == 'Calopteryx virgo padana'] <- 'Calopteryx virgo'
dfo$species[
  dfo$species == 'Calopteryx virgo meridionalis'
] <- 'Calopteryx virgo'
dfo$species[dfo$species == 'Calopteryx virgo virgo'] <- 'Calopteryx virgo'
dfo$species[dfo$species == 'Gomphus flavipes'] <- 'Stylurus flavipes'
dfo$species[dfo$species == 'Lestes viridis'] <- 'Chalcolestes viridis'
dfo$species[dfo$species == 'Lestes virens vestalis'] <- 'Lestes virens'
dfo$species[
  dfo$species == 'Onychogomphus forcipatus unguiculatus'
] <- 'Onychogomphus forcipatus'
dfo$species[
  dfo$species == 'Onychogomphus forcipatus forcipatus'
] <- 'Onychogomphus forcipatus'
dfo$species[dfo$species == 'Oxygastra curtisi'] <- 'Oxygastra curtisii'
dfo$species[
  dfo$species == 'Cordulegaster bidentata bidentata'
] <- 'Cordulegaster bidentata'
dfo$species[
  dfo$species == 'Orthetrum coerulescens coerulescens'
] <- 'Orthetrum coerulescens'

dfo <- dfo |> filter(!str_detect(species, "sp\\.$|sp\\s|sp$"))
dfo <- dfo |> filter(species != 'No odonata')

dfo |>
  filter(species == 'Selysiothemis nigra') |>
  count(species) |>
  arrange(desc(n))

dim(dfo)
unique(dfo$species) # 69 specie totali


# controlli su celle ----

# totale celle
ncell <- dfo |>
  distinct(cellcode) |>
  count() |>
  pull()

# celle monitorate x anno
cells <- dfo |>
  group_by(date_year) |>
  distinct(cellcode)

# cells |>
#   count(date_year) |>
#   mutate(p = n / ncell) |>
#   ggplot(aes(date_year, p)) +
#   geom_area(fill = 'forestgreen', alpha = 0.6) +
#   labs(y = 'no. cells\n', x = '') +
#   scale_y_continuous(labels = percent_format())

# aggiunta dei codici celle ----
