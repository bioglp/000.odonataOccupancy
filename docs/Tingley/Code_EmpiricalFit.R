#######################################################################
## The following code was used to analyze
## the empirical data presented in:
## Tingley et al. 2020. Multi-species occupancy models as robust 
## estimators of community richness. Methods in Ecology and Evolution.
########################################################################

########### DESCRIPTION
# Run JAGS model and extract z-matrix
###########

## Load Libraries
library(R2jags)
library(rjags)
library(foreign)
library(dclone)     

## Load data
folder<-getwd()
load("Data_EmpiricalRaw.Rdata") 
attach(export.ms)
setwd(folder)

############ Pick year and limit data
year <- 2018
year.keep <- which(is.na(y[, 1, which(years == year), 1]) == F)
visits.keep <- which(colSums(is.na(y[year.keep, , which(years == year), 1])) < length(year.keep))
  
######## Augment data
y.year <- y[year.keep, visits.keep, which(years == year), ]
n.zero <- 80 # arbitrarily large number of all-zero species
h.zero <- is.na(y.year[,,1])*1
h.zero[, ] <- 0
# Augment the data
y.aug <- array(c(as.vector(y.year), rep(h.zero, times = n.zero)), 
               dim = dim(y.year) + c(0, 0, n.zero))
  
############# Prepare data for JAGS model

#### JAGS data
jags.data <- list(
  y = y.aug,
  survey = ef,
  date = day.of.year[year.keep, which(years == year)],
  time = time.of.day[year.keep, which(years == year)],
  elev = elev[year.keep],
  elev2 = elev2[year.keep],
  cc = cc[year.keep],
  precc = precc[year.keep],
  fa = as.vector(scale(fire.age[year.keep, which(years == year)])),
  n.sites = dim(y.aug)[1],
  n.visits = dim(y.aug)[2],
  n.species = dim(y.year)[3],
  n.zeroes = n.zero
  )

z.naive <- apply(jags.data$y, MARGIN = c(1,3), max, na.rm=T) 
w.naive <- apply(z.naive, 2, max)
omega.naive <- sum(w.naive) / (n.species + n.zero)

############### JAGS Model specification ############

community <- function() {         
  ######################
  ## Define hyper-parameters
  ######################  
  # Occupancy (mean and variance for each)
  mu.b0 ~ dnorm(0, 0.001)
  mu.b.elev ~ dnorm(0, 0.001)
  mu.b.elev2 ~ dnorm(0, 0.001)
  mu.b.burn ~ dnorm(0, 0.001)
  mu.b.canopy ~ dnorm(0, 0.001)
  mu.b.fa ~ dnorm(0, 0.001)
  tau.b0 ~ dgamma(0.1,0.1)
  tau.b.elev ~ dgamma(0.1,0.1)
  tau.b.elev2 ~ dgamma(0.1,0.1)
  tau.b.burn ~ dgamma(0.1,0.1)
  tau.b.canopy ~ dgamma(0.1,0.1)
  tau.b.fa ~ dgamma(0.1,0.1)
  # Detection (mean and variance for each)
  mu.a0 ~ dnorm(0, 0.001)
  mu.a.survey ~ dnorm(0, 0.001)
  mu.a.date ~ dnorm(0, 0.001)
  mu.a.time ~ dnorm(0, 0.001)
  tau.a0 ~ dgamma(0.1,0.1)
  tau.a.date ~ dgamma(0.1,0.1)
  tau.a.time ~ dgamma(0.1,0.1)
  tau.a.survey ~ dgamma(0.1,0.1)
  # Prior for the proportion of unobserved species that exist in community
  omega ~ dunif(0, 1)
  
  ######################
  ## Draw coefficients for each species from hyper-parameters
  ######################
  # Start loop over all species
  for(i in 1:(n.species + n.zeroes)) {
    # coefficients for each species
    b0[i] ~ dnorm(mu.b0, tau.b0)
    b.elev[i] ~ dnorm(mu.b.elev, tau.b.elev)
    b.elev2[i] ~ dnorm(mu.b.elev2, tau.b.elev2)
    b.burn[i] ~ dnorm(mu.b.burn, tau.b.burn)
    b.canopy[i] ~ dnorm(mu.b.canopy, tau.b.canopy)
    b.fa[i] ~ dnorm(mu.b.fa, tau.b.fa)
    a0[i] ~ dnorm(mu.a0, tau.a0)
    a.survey[i] ~ dnorm(mu.a.survey, tau.a.survey)
    a.date[i] ~ dnorm(mu.a.date, tau.a.date)
    a.time[i] ~ dnorm(mu.a.time, tau.a.time)
    # Create variable which is whether a species is a true member of the regional community
    w[i] ~ dbern(omega)
    
    #########################
    ## Define the likelihood
    #########################
    ##loop through all of the sites (N)
    for(j in 1:n.sites) {
      ## State process (i.e., the true occupancy of the site by species i)
      logit(psi[j, i]) <- b0[i] + b.elev[i]*elev[j] + b.elev2[i]*elev2[j] + b.burn[i]*cc[j] + b.canopy[i]*precc[j] + b.fa[i]*fa[j]
      mu.psi[j, i] <- psi[j, i] * w[i]
      z[j, i] ~ dbern(mu.psi[j, i])
      ## Observation process
      ##loop through the repeat surveys (T)
      for(k in 1:n.visits) {
        logit(p[j, k, i]) <- a0[i] + a.survey[i]*survey[k] + a.date[i]*date[j] + a.time[i]*time[j]
        mu[j, k, i] <- z[j, i] * p[j, k, i]
        y[j, k, i] ~ dbern(mu[j, k, i])  
      } ##close k loop
    } ##close the j loop
  } ##close the i loop
}


###################----------Running JAGS----------###################
inits <- function() {list(z = z.naive,
                          w = w.naive,
                          omega = omega.naive)}
parms <- c("mu.b0", 
           "mu.a0",
           "mu.a.survey",
           "mu.a.date",
           "mu.a.time",
           "omega")

setwd(folder)
nc <- 3
n.adapt <- 1000
n.burn <- 10000
n.iter <- 30000
thin <- 30

start.time<-Sys.time()
cl <- makePSOCKcluster(nc)   
tmp <- clusterEvalQ(cl, library(dclone))  
parLoadModule(cl, "glm")  
parListModules(cl)  
model <- jags.parfit(cl, jags.data, params = parms, community, inits=inits, n.chains=nc,
                              n.adapt=n.adapt, n.update = n.burn, thin = thin, n.iter = n.iter)
stopCluster(cl)   
end.time=Sys.time()
elapsed.time = difftime(end.time, start.time, units='mins')
elapsed.time #about 1000 minutes per model

## save
setwd(folder)
save(model, file=paste("CommunityResults_",year,".Rdata", sep = ""))


