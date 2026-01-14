dir.create('input', showWarnings = FALSE)
#dir.create('output', showWarnings = FALSE)
#dir.create('plots', showWarnings = FALSE)
dir.create('R-script', showWarnings = FALSE)
library(rtrim)
library(pacman)
p_load(tidyverse, gt, skimr, ggpubr, rstatix, openxlsx, scales, patchwork)
theme_set(theme_bw(base_size = 14))


year_start <- 2005

dfo <- read.xlsx('input/datasetSegnalazioni_R20251106.xlsx')
dfo <- dfo |> rename_with(tolower)
dfo <- dfo |>
  filter(
    date_year >= year_start,
    date_year < 2026,
    code == 'Continental',
    !is.na(cellcode20)
  )
dfo <- dfo |>
  rename(species = latin_taxon)
glimpse(dfo) #View(dfo)

# rinominare specie
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
