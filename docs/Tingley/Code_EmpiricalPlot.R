#######################################################################
## The following code was used to plot
## the empirical data presented in:
## Tingley et al. 2020. Multi-species occupancy models as robust 
## estimators of community richness. Methods in Ecology and Evolution.
########################################################################

########### DESCRIPTION
# Comparison of empirical vs estimated
###########

## Load Libraries
library(R2jags)
library(SpadeR)
library(viridis)

## Load data
folder<-getwd()
load("CommunityResults_2010.Rdata")
r2010 <- model
load("CommunityResults_2011.Rdata")
r2011 <- model
load("CommunityResults_2012.Rdata")
r2012 <- model
load("CommunityResults_2013.Rdata")
r2013 <- model
load("CommunityResults_2014.Rdata")
r2014 <- model
load("CommunityResults_2015.Rdata")
r2015 <- model
load("CommunityResults_2016.Rdata")
r2016 <- model
load("CommunityResults_2017.Rdata")
r2017 <- model
load("CommunityResults_2018.Rdata")
r2018 <- model
rm(model)
load("Data_EmpiricalRaw.Rdata")
attach(export.ms)

######
r2010 <- data.frame(rbind(r2010[[1]], r2010[[2]], r2010[[3]]))
r2011 <- data.frame(rbind(r2011[[1]], r2011[[2]], r2011[[3]]))
r2012 <- data.frame(rbind(r2012[[1]], r2012[[2]], r2012[[3]]))
r2013 <- data.frame(rbind(r2013[[1]], r2013[[2]], r2013[[3]]))
r2014 <- data.frame(rbind(r2014[[1]], r2014[[2]], r2014[[3]]))
r2015 <- data.frame(rbind(r2015[[1]], r2015[[2]], r2015[[3]]))
r2016 <- data.frame(rbind(r2016[[1]], r2016[[2]], r2016[[3]]))
r2017 <- data.frame(rbind(r2017[[1]], r2017[[2]], r2017[[3]]))
r2018 <- data.frame(rbind(r2018[[1]], r2018[[2]], r2018[[3]]))
r.all <- rbind(r2010, r2011, r2012, r2013, r2014, r2015, r2016, r2017, r2018)
r.all$Year <- rep(years, each = dim(r2018)[1])
n.zero <- 80
true.species <- sum(apply(apply(y, c(3,4), max, na.rm = T), 2, max))

##### Chao estimators
chao.est <- data.frame(chao = 1, s.e. = 1, l95 = 1, u95 = 1)
for(i in 1:n.years) {
	zY <- apply(y[,,i,], c(1,3), max, na.rm = T)
	zY <- zY[zY[,1] != "-Inf",] # this fixes the error warnings
	chao.est[i, ] <- matrix(ChaoSpecies(t(zY), datatype = "incidence_raw", k = 20, conf = 0.95)$Species_table[4, ], ncol = 4)
}
chao.est$Year <- years

###### Table 2: Bias, Accuracy, Precision
msom.est <- data.frame(mean = tapply(r.all$omega, INDEX = r.all$Year, FUN = mean)*(n.zero + n.species))
msom.est$sd <- tapply((r.all$omega*(n.zero + n.species)), INDEX = r.all$Year, FUN = sd)
table2 <- cbind(msom.est, chao.est[, 1:2])
table2$bias.msom <- table2$mean - true.species
table2$bias.chao <- table2$chao - true.species
table2$accur.msom <- (table2$bias.msom)^2
table2$accur.chao <- (table2$bias.chao)^2
table2$precis.msom <- (table2$mean - mean(table2$mean))^2
table2$precis.chao <- (table2$chao - mean(table2$chao))^2
write.csv(table2, file = "Table2.csv")

##### FIGURE 3 ######
accum <- t(apply(y, MARGIN = c(3,4), FUN = max, na.rm = T))

pdf(file = "Figure3.pdf", width = 8.5, height = 4)
par(mfrow = c(1, 2), mar=c(3,3,2.2,0.4), mgp = c(1.8,0.3,0), tcl = -0.2, cex.axis = 0.9, las = 1)

boxplot(colSums(accum), 
				chao.est$chao,
				c(mean(r2010$omega)*(n.zero + n.species),
								 mean(r2011$omega)*(n.zero + n.species),
								 mean(r2012$omega)*(n.zero + n.species),
								 mean(r2013$omega)*(n.zero + n.species),
								 mean(r2014$omega)*(n.zero + n.species),
								 mean(r2015$omega)*(n.zero + n.species),
								 mean(r2016$omega)*(n.zero + n.species),
								 mean(r2017$omega)*(n.zero + n.species),
				        mean(r2018$omega)*(n.zero + n.species)),
				names = c("Observed", "Chao", "MSOM"),
				boxwex = 0.5,
				staplewex = 0,
				lty = 1,
				lwd = 2,
				col = "#00000060",
				ylab = "Community Size Estimate (N)", ylim = c(70, 200))
abline(h = true.species, col = "red", lwd = 1, lty = 2)
legend("topright", bty = "n", col = c("red"), lwd = 1, lty = 2, cex = 1.3, legend = c("Truth"))
title(main = "(a)", adj = 0.02, line = 0.5, cex.main = 1.5)
##
plot(0,0, type = "n", ylim = c(70, 200), xlim = c(2009.5, 2018.5), xlab = "Year", ylab = "Community Size Estimate (N)")
abline(h = true.species, col = "red", lty = 2)
for(i in 1:n.years) {
  lines(c(years[i] - 0.15, years[i] - 0.15), y = quantile(r.all$omega[r.all$Year == years[i]]*(n.zero + n.species), probs = c(0.025, 0.975)))
  lines(c(years[i] - 0.15, years[i] - 0.15), y = c(mean(r.all$omega[r.all$Year == years[i]]*(n.zero + n.species)) + sd(r.all$omega[r.all$Year == years[i]]*(n.zero + n.species)), mean(r.all$omega[r.all$Year == years[i]]*(n.zero + n.species)) - sd(r.all$omega[r.all$Year == years[i]]*(n.zero + n.species))), lwd = 4)
  points(years[i] - 0.15, mean(r.all$omega[r.all$Year == years[i]]*(n.zero + n.species)), pch = 16, cex = 2)
  lines(c(years[i] + 0.15, years[i] + 0.15), c(chao.est$l95[i], chao.est$u95[i]))
  lines(c(years[i] + 0.15, years[i] + 0.15), c(chao.est$chao[i] + chao.est$s.e.[i], chao.est$chao[i] - chao.est$s.e.[i]), lwd = 4, col = "gray")
  points(years[i] + 0.15, chao.est$chao[i], pch = 16, col = "gray", cex = 2)
}
legend("topright", bty = "n", pch = 16, col = c("black", "gray"), lwd = 4, cex = 1.0, legend = c("MSOM", "Chao"))
title(main = "(b)", adj = 0.02, line = 0.5, cex.main = 1.5)
dev.off()

##### FIGURE 4 ######
pdf(file = "Figure4.pdf", width = 8.5, height = 4)
par(mfrow = c(1, 2), mar=c(3,3,2.2,0.4), mgp = c(1.8,0.3,0), tcl = -0.2, cex.axis = 0.9, las = 1)

plot(plogis(r.all$mu.b0), r.all$omega*(n.zero + n.species), pch = 16, col = viridis(n = 9, alpha = 0.2)[r.all$Year-2009],
     xlab = expression(paste("Mean probability of occupancy (",italic(mu[psi]), ")")), ylab = "Community Size Estimate (N)")
points(x = plogis(tapply(r.all$mu.b0, r.all$Year, mean)),
       y = tapply(r.all$omega, r.all$Year, mean)*(n.zero + n.species),
       pch = 9)
text(x = plogis(tapply(r.all$mu.b0, r.all$Year, mean))[-c(2,6)],
     y = as.vector(tapply(r.all$omega, r.all$Year, mean)*(n.zero + n.species))[-c(2,6)],
     labels = years[-c(2,6)],
     pos = 4, cex = 0.6)
text(x = plogis(tapply(r.all$mu.b0, r.all$Year, mean))[c(2,6)],
     y = as.vector(tapply(r.all$omega, r.all$Year, mean)*(n.zero + n.species))[c(2,6)],
     labels = years[c(2,6)],
     pos = 3, cex = 0.6)
title(main = "(a)", adj = 0.02, line = 0.5, cex.main = 1.5)
# Detect
plot(plogis(r.all$mu.a0), r.all$omega*(n.zero + n.species), pch = 16, col = viridis(n = 9, alpha = 0.2)[r.all$Year-2009],
     xlab = expression(paste("Mean probability of detection ",(mu[p]))), ylab = "Community Size Estimate (N)", xlim = c(0.3, 0.72))
points(x = plogis(tapply(r.all$mu.a0, r.all$Year, mean)),
       y = tapply(r.all$omega, r.all$Year, mean)*(n.zero + n.species),
       pch = 9)
text(x = plogis(tapply(r.all$mu.a0, r.all$Year, mean)),
     y = as.vector(tapply(r.all$omega, r.all$Year, mean)*(n.zero + n.species)),
     labels = years,
     pos = 3, cex = 0.6)
title(main = "(b)", adj = 0.02, line = 0.5, cex.main = 1.5)
dev.off()