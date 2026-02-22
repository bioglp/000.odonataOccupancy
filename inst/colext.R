library(unmarked)

data(crossbill)
colnames(crossbill)

DATE <- as.matrix(crossbill[,32:58])
y.cross <- as.matrix(crossbill[,5:31])
y.cross[is.na(DATE) != is.na(y.cross)] <- NA

sd.DATE <- sd(c(DATE), na.rm=TRUE)
mean.DATE <- mean(DATE, na.rm=TRUE)
DATE <- (DATE - mean.DATE) / sd.DATE

years <- as.character(1999:2007)
years <- matrix(years, nrow(crossbill), 9, byrow=TRUE)
umf <- unmarkedMultFrame(y=y.cross,
                         siteCovs=crossbill[,2:3], 
                         yearlySiteCovs=list(year=years),
                         obsCovs=list(date=DATE),
                         numPrimary=9)

# Models
# A model with constant parameters
fm0 <- colext(~1, ~1, ~1, ~1, umf)

# Like fm0, but with year-dependent detection
fm1 <- colext(~1, ~1, ~1, ~year, umf)

# Like fm0, but with year-dependent colonization and extinction
fm2 <- colext(~1, ~year-1, ~year-1, ~1, umf)

# A fully time-dependent model
fm3 <- colext(~1, ~year-1, ~year-1, ~year, umf)

# Like fm3 with forest-dependence of 1st-year occupancy
fm4 <- colext(~forest, ~year-1, ~year-1, ~year, umf)

# Like fm4 with date- and year-dependence of detection
fm5 <- colext(~forest, ~year-1, ~year-1, ~year + date + I(date^2),
              umf, starts=c(coef(fm4), 0, 0))

# Same as fm5, but with detection in addition depending on forest cover
fm6 <- colext(~forest, ~year-1, ~year-1, ~year + date + I(date^2) +
                forest, umf)

models <- fitList('psi(.)gam(.)eps(.)p(.)'    = fm0,
                  'psi(.)gam(.)eps(.)p(Y)'    = fm1,
                  'psi(.)gam(Y)eps(Y)p(.)'    = fm2,
                  'psi(.)gam(Y)eps(Y)p(Y)'    = fm3,
                  'psi(F)gam(Y)eps(Y)p(Y)'    = fm4,
                  'psi(F)gam(Y)eps(Y)p(YD2)'  = fm5,
                  'psi(F)gam(Y)eps(Y)p(YD2F)' = fm6)
ms <- modSel(models)
ms

op <- par(mfrow=c(1,2), mai=c(0.8,0.8,0.1,0.1))

nd <- data.frame(forest=seq(0, 100, length=50))
E.psi <- predict(fm6, type="psi", newdata=nd, appendData=TRUE)

with(E.psi, {
  plot(forest, Predicted, ylim=c(0,1), type="l",
       xlab="Percent cover of forest",
       ylab=expression(hat(psi)), cex.lab=0.8, cex.axis=0.8)
  lines(forest, Predicted+1.96*SE, col=gray(0.7))
  lines(forest, Predicted-1.96*SE, col=gray(0.7))
})

nd <- data.frame(date=seq(-2, 2, length=50),
                 year=factor("2005", levels=c(unique(years))),
                 forest=50)
E.p <- predict(fm6, type="det", newdata=nd, appendData=TRUE)
E.p$dateOrig <- E.p$date*sd.DATE + mean.DATE

with(E.p, {
  plot(dateOrig, Predicted, ylim=c(0,1), type="l",
       xlab="Julian date", ylab=expression( italic(p) ),
       cex.lab=0.8, cex.axis=0.8)
  lines(dateOrig, Predicted+1.96*SE, col=gray(0.7))
  lines(dateOrig, Predicted-1.96*SE, col=gray(0.7))
})

#-----------------

library(unmarked)

# Dati dalla tua tabella (escludendo la colonna 'siteX')
y_data <- matrix(c(
  1, 1, 0,
  0, 0, NA,
  NA, 1, 0,
  NA, 0, 1,
  1, 1, 1
), nrow = 5, byrow = TRUE, 
dimnames = list(paste0("site", 1:5), 
                paste0("y", 2023:2025)))

# Creazione dell'oggetto unmarkedFrameColExt
# Assumiamo che non ci siano covariate del sito (siteCovs) o delle stagioni (yearlySiteCovs)
# n-obs è il numero di rilevamenti per stagione, che in questo caso è 1 (una colonna per anno)
umf_colext <- unmarkedMultFrame(y = y_data, numPrimary = 3)

# Visualizza il sommario
summary(umf_colext)

# Modello:
# 1. Probabilità di Occupazione Iniziale (psi)
# 2. Probabilità di Colonizzazione (gamma)
# 3. Probabilità di Estinzione (epsilon)
# 4. Probabilità di Rilevamento (p)
# ~1~1~1~1 significa che tutti i parametri sono modellati come costanti (intercetta)
# Colext formula: ~cov_p ~cov_gamma ~cov_epsilon ~cov_psi
model_colext <- colext(psiformula = ~1, gammaformula = ~1,
                       epsilonformula = ~1, pformula = ~1,
                       data = umf_colext)

# Stima dei parametri
print(model_colext)

# ---------

# --- 1. Creazione dei DataFrame Iniziali ---

# Tabella Presenze (Species Detections = 1)
presenze <- tibble(
  site = c("site1", "site5", "site1", "site3", "site5", "site4", "site5"),
  year = c(2023, 2023, 2024, 2024, 2024, 2025, 2025),
  # Aggiungiamo la colonna "rilevamento" che sarà 1
  rilevamento = 6
)

# Tabella Siti Visitati (Survey Effort = 0 o 1)
visitati <- tibble(
  site = c("site1", "site2", "site5", "site1", "site2", "site3", "site4", "site5", "site1", "site3", "site4", "site5"),
  year = c(2023, 2023, 2023, 2024, 2024, 2024, 2024, 2024, 2025, 2025, 2025, 2025)
)

# --- 2. Preparazione della Matrice di Rilevamento (y) ---

dati_occupancy_long <- visitati %>%
  left_join(presenze, by = c("site", "year")) %>%
  
  # b) Sostituisci NA con 0 (i siti visitati in cui non c'è stata presenza)
  mutate(y = replace_na(rilevamento, 0)) %>%
  
  # Tieni solo le colonne necessarie
  select(site, year, y)


# b) Converti il formato da "long" a "wide" (formato Matrice y)
# Questo creerà automaticamente NA per le combinazioni site/year non visitate
matrice_y_df <- dati_occupancy_long %>%
  pivot_wider(
    names_from = year,
    values_from = y,
    # Riempie i valori mancanti (siti non visitati in un dato anno) con NA
    values_fill = NA
  ) %>%
  # Rimuovi la colonna site e converti in matrice R
  select(-site) 

# Converti il DataFrame in una matrice, essenziale per unmarked
matrice_y <- as.matrix(matrice_y_df)
# Assegna i nomi delle righe per chiarezza (opzionale)
rownames(matrice_y) <- matrice_y_df$site