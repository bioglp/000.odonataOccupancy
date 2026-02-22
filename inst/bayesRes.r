library(tidyverse)
db <- readRDS('output_bayes/risultati_Trithemis annulata.rds')
db$summary$statistics

est <- db$summary$statistics[
  grep("mean_psi", rownames(db$summary$statistics)),
]

anno <- 2005:2025
est <- cbind(est, year = anno)

est |>
  ggplot(aes(year, Mean)) +
  geom_line() +
  ylim(0, 1)

# Estraiamo i campioni di tutte le catene per il trend
beta_samples <- as.matrix(db$samples)[, "beta_year_psi"]

# Quanti campioni sono sotto lo zero?
prob_declino <- mean(beta_samples < 0)

# Visualizziamo il risultato
cat("Probabilità di declino dell'areale:", round(prob_declino * 100, 2), "%\n")

db$summary
