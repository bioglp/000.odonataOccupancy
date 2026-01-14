library(dplyr)
library(unmarked)

df_sp <- dfp %>% filter(species == "Aeshna affinis")

df_sp <- df_sp %>%
  mutate(site = paste(cellcode20, date_year, sep = "_"))

y <- matrix(df_sp$total_count, ncol = 1)
rownames(y) <- df_sp$site

siteCovs <- data.frame(year = df_sp$date_year)

umf <- unmarkedFramePCount(y = y, siteCovs = siteCovs)

fm <- pcount(~ 1 ~ year, data = umf, K = K)
summary(fm)

newdat <- data.frame(year = sort(unique(df_sp$date_year)))

pred <- predict(fm, type = "state", newdata = newdat)
pred <- cbind(newdat,pred)

pred %>% 
  ggplot(aes(x = year, y = Predicted)) +
  geom_line(size = 1.2, color = "darkgreen") +
  # geom_ribbon(aes(ymin = lower, ymax = upper), 
  #             fill = "darkgreen", alpha = 0.25) +
  labs(
    x = "",
    y = "Abundance (λ)",
    title = paste("Trend for",name)
  )
