#######################################################################
## The following code was used to simulate, analyze, and plot
## the simulation data presented in:
## Tingley et al. 2020. Multi-species occupancy models as robust 
## estimators of community richness. Methods in Ecology and Evolution.
########################################################################

library(R2jags)
library(SpadeR)
library(doParallel)
library(foreach)

######################
## Simulate the data
######################
  #################################################
  ## Define all the parameters for the simulations
  #################################################
    N = c(20, 125) ##number of species
    n.reps = c(2, 7) ##number of replicate surveys
    n.sites = c(60, 30) ##number of sampling locations
    mu.psi = c(qlogis(0.1), qlogis(0.7)) ##mean species occupancy probability (converted to logit scale)
    sd.psi = c(0.2, 1) ##SD of species occupancy probabilities (on logit scale)
    mu.p = c(qlogis(0.3), qlogis(0.7)) ##mean species detection probability (converted to logit scale)
    sd.p = c(0.2, 1) ##SD of species detection probabilities (on logit scale)

  #############################################################
  ## Create a table of parameters for each type of simulation
  #############################################################
    n.models = 7
    models = data.frame(mu.psi = mu.psi[c(2,2,1,2,2,1,1)],
                        sd.psi = sd.psi[c(1,2,1,1,1,2,1)],
                        mu.p = mu.p[c(2,2,2,2,1,1,1)],
                        sd.p = sd.p[c(1,1,1,2,1,2,1)],
                        N = c(rep(N[1], 2*n.models), rep(N[2], 2*n.models)),
                        n.reps = c(rep(n.reps[1], n.models),
                                   rep(n.reps[2], n.models),
                                   rep(n.reps[1], n.models),
                                   rep(n.reps[2], n.models)),
                        n.sites = c(rep(n.sites[1], n.models * 4),
                                    rep(n.sites[2], n.models * 4)))
    models$model = rep(1:7)

  ###################################
  ## A function to simulate the data
  ###################################
    sim.data = function(mu.psi, sd.psi, mu.p, sd.p, N, n.reps, n.sites){
      ##assign psi and p for each species
        psi = plogis(rnorm(N, mu.psi, sd.psi))
        p = plogis(rnorm(N, mu.p, sd.p))

      ##assign the number of truely occupied sites
        z = matrix(rbinom(n = n.sites * N,
                          size = 1,
                          prob = psi),
                   nrow = n.sites,
                   ncol = N,
                   byrow = T)

      ##assign the number of replicates surveys where sp. detected
        y = matrix(rbinom(n = n.sites * N,
                          size = n.reps,
                          prob = c(t(t(z)*p))),
                   nrow = n.sites,
                   ncol = N)

      ##get the latent number of species in the community
        sim.N = sum(apply(z, 2, sum)>0)
        sim.obs = sum(apply(y, 2, sum)>0)

      ##save all the data in a list
        d = list(z, y, psi, p, N, sim.N, sim.obs, n.sites, n.reps)
        names(d) = c("z", "y", "psi", "p", "N", "sim.N", "sim.obs", "n.sites", "n.reps")
        return(d)
    }

  ######################
  ## Simulate the data
  ######################
    sim.results = replicate(100,
                            mapply(sim.data,
                                   models$mu.psi,
                                   models$sd.psi,
                                   models$mu.p,
                                   models$sd.p,
                                   models$N,
                                   models$n.reps,
                                   models$n.sites))

    ##save the results
      save(sim.results, models,
           file = "sim_data.R")

##############################
## Analyze the simulated data
##############################
  ############################################################
  ## MSOM Function
  ## d is a simulated data object from the sim.data function
  ## m is the model number in the table models created above 
  ############################################################
    occ = function(d, m){
      ##augment the data
        y.aug = cbind(d$y, matrix(0, nrow = d$n.sites, ncol = d$sim.obs))

      ##bundle the data
        jags.data = list(N.aug = ncol(y.aug),
                         y = y.aug,
                         n.sites = d$n.sites,
                         n.reps = d$n.reps)

      ##create the initial values
        z.naive = (y.aug>0)*1
        w.naive = as.numeric(apply(z.naive, 2, sum)>0)
        omega.naive = mean(w.naive)

        inits = function() {list(z = z.naive,
                                 w = w.naive,
                                 omega = omega.naive)}

      ##specify parameters to monitor
        parms = c("mu.psi",
                  "sd.psi",
                  "mu.p",
                  "sd.p",
                  "N")

      ##run the model
        model = jags(data = jags.data,
                     model.file = "Simulation_JAGS_model.txt",
                     inits = inits,
                     parameters.to.save = parms,
                     n.burnin = 1000,
                     n.iter = 20000,
                     n.thin = 25,
                     n.chains = 3)

      ##save the results
        results = list(m = m,
                       model$BUGSoutput$sims.list$N,
                       model$BUGSoutput$sims.list$mu.p,
                       model$BUGSoutput$sims.list$mu.psi,
                       model$BUGSoutput$sims.list$sd.p,
                       model$BUGSoutput$sims.list$sd.psi,
                       model$BUGSoutput$summary)
        names(results) = c("model" ,"N", "mu.p", "mu.psi", "sd.p", "sd.psi", "summary")
        return(results)
    }

  ###################################
  ## Analyze all the data using MSOM
  ###################################
    ##load the simulated data
      load("sim_data.R")

    ##set up the parallel processing
      cl = makeCluster(6)
      registerDoParallel(cl)

    ##loop through the replicates
      occ.results = foreach(i = 1:56, .packages = c("R2jags")) %dopar% {
        apply(sim.results[, i, 1:dim(sim.results)[3]], 2, occ, m = i)
      }

    stopCluster(cl)

    save(occ.results,
         file = "OccupancyResults.RData")

  ###########################################################################################
  ## iChao2 Function
  ## d is a simulated data object from the sim.data function
  ## pooled.reps is used to either random select a single replicate survey (pooled.reps = F)
  ## or to combine the observations across the replicate surveys (pooled.reps = T)  
  ############################################################################################
    chao = function(d, pooled.reps = T){
      ##change the y matrix to presence/absence
      if(pooled.reps){
        y01 = d$y > 0
      }else{
        ##create a blank matrix of the right size
          y01 = matrix(0, nrow = nrow(d$y), ncol = ncol(d$y))

        ##loop through the y matrix
          for(i in 1:nrow(d$y)){
            for(j in 1:ncol(d$y)){
              ##create a vector with a simulated detection history for the number of reps
                v = c(rep(1, d$y[i,j]), 
                      rep(0, d$n.reps - d$y[i,j]))
  
              ##randomly select a sampling occassion from the detection history
              ##which determines whether species was detected during
              ##a single survey
                y01[i,j] = sample(v, 1)
  
            }##close species loop
          }##close sites loop
      }##close else

      ##run the ChaoSpecies function from the SpadeR package
      ##use tryCatch to catch errors caused when there are no
      ##singletons or doubletons (see Methods)
        out = tryCatch({
                  ChaoSpecies(t(y01), datatype = "incidence_raw", k = 10)
                },
                error = function(cond){
                  return(NA)
                })

      ##if ChaoSpecies did not throw an error, save the ChaoSpecies results
      ##otherwise save NA
        if(!is.na(out)){
          out$Species_table
        }else{
          NA
        }
        
    }##close function

  ######################
  ## Analyze using Chao
  ######################
    ##load simulation results
      load("sim_data.RData")

    ##analyze pooling all the reps
      chao.results.pooled.reps = lapply(1:56, function(i){
                                    t(apply(sim.results[, i, 1:dim(sim.results)[3]], 2, chao))
                                  })

    ##analyze selecting one rep
      chao.results.unpooled.reps = lapply(1:56, function(i){
                                    t(apply(sim.results[, i, 1:dim(sim.results)[3]], 2, chao, pooled.reps = F))
                                   })

    ##save the chao results
      diff = setdiff(ls(), c("chao.results.pooled.reps", "chao.results.unpooled.reps"))
      rm(list = diff)
      rm(diff)
      save.image("ChaoResults_SpadePackage.RData")

###################################
## Summarize and plot the results
###################################
  ##load the simulation data
    load("sim_data.RData")
      
  ##load the chao results
    load("ChaoResults_SpadePackage.RData")
  
  ##load all the  occupancy results
    load("OccupancyResults.RData")

  ############################################################
  ## Plot bias, precision, accuracy, and coverage probability
  ############################################################
    #################################
    ## Save the estimates in a table
    #################################
      ##occupancy
        occ.est = sapply(1:56, function(m){ 
          sapply(1:100, function(i){
            mean(occ.results[[m]][[i]]$N)
          })
        })
        
      ##calculate the chao estimates for the pooled results
        chao.pooled.est = sapply(1:56, function(m){
          sapply(1:100, function(i){
            if(is.data.frame(chao.results.pooled.reps[[m]][[i]])){
              chao.results.pooled.reps[[m]][[i]][4,"Estimate"]
            }else{
              sim.results[,m,i]$sim.N
            }
          })
        })
      
      ##unpooled chao  
        chao.unpooled.est = sapply(1:56, function(m){
          sapply(1:100, function(i){
            if(is.data.frame(chao.results.unpooled.reps[[m]][[i]])){
              chao.results.unpooled.reps[[m]][[i]][4,"Estimate"]     
            }else{
              sim.results[,m,i]$sim.N
            }
          })
        })
        
    #######################
    ## Calculate the error
    #######################
      occ.error = sapply(1:56, function(m){ 
                    sapply(1:100, function(i){
                      (mean(occ.results[[m]][[i]]$N) - models[m, "N"])
                    })
                  })
      
      pooled.error = sapply(1:56, function(m){
                      sapply(1:100, function(i){
                        if(is.data.frame(chao.results.pooled.reps[[m]][[i]])){
                          chao.results.pooled.reps[[m]][[i]][4,"Estimate"] -  models[m,"N"]     
                        }else{
                          sim.results[,m,i]$sim.N - models[m, "N"]
                        }
                      })
                    })
        
      unpooled.error = sapply(1:56, function(m){
                        sapply(1:100, function(i){
                          if(is.data.frame(chao.results.unpooled.reps[[m]][[i]])){
                            chao.results.unpooled.reps[[m]][[i]][4,"Estimate"] -  models[m,"N"]     
                          }else{
                            sim.results[,m,i]$sim.N - models[m, "N"]
                          }
                        })
                      })
        
    ###########################################################
    ## Save bias, precision, and accuracy for the 56 scenarios
    ###########################################################
      ##create a table to store the results
        model.stats = models
      
      ##save the mean estimate
        model.stats$occ.mean = apply(occ.est, 2, mean)
        model.stats$chao.p.mean = apply(chao.pooled.est, 2, mean)
        model.stats$chao.up.mean = apply(chao.unpooled.est, 2, mean)
      
      ##save the scaled bias
        model.stats$occ.bias = apply(occ.error, 2, mean)/model.stats$N
        model.stats$chao.p.bias = apply(pooled.error, 2, mean)/model.stats$N
        model.stats$chao.up.bias = apply(unpooled.error, 2, mean)/model.stats$N
      
      ##save the scaled precision
        model.stats$occ.pre = apply(occ.est, 2, FUN = function(x){sd(x)/mean(x)})
        model.stats$chao.p.pre = apply(chao.pooled.est, 2, FUN = function(x){sd(x)/mean(x)})
        model.stats$chao.up.pre = apply(chao.unpooled.est, 2, FUN = function(x){sd(x)/mean(x)})
      
      ##save the accuracy
        model.stats$occ.acc = apply(occ.error, 2, FUN = function(x){sqrt(mean(x^2))})/model.stats$N
        model.stats$chao.p.acc = apply(pooled.error, 2, FUN = function(x){sqrt(mean(x^2))})/model.stats$N
        model.stats$chao.up.acc = apply(unpooled.error, 2, FUN = function(x){sqrt(mean(x^2))})/model.stats$N
        
    ###################################################
    ## Save coverage probabililty for the 56 scenarios
    ###################################################
      ##occupancy
        model.stats$occ.coverage = sapply(1:56, function(m){
          mean(sapply(1:100, function(i){
            models[m, "N"] >= round(quantile(occ.results[[m]][[i]]$N, 0.025)) &
              models[m, "N"] <= round(quantile(occ.results[[m]][[i]]$N, 0.975))
          }))
        })
        
        
      ##pooled Chao
        model.stats$chao.p.coverage = sapply(1:56, function(m){
          mean(sapply(1:100, function(i){
            if(is.data.frame(chao.results.pooled.reps[[m]][[i]])){
              models[m, "N"] >= round(chao.results.pooled.reps[[m]][[i]][4,"95%Lower"]) &
                models[m, "N"] <= round(chao.results.pooled.reps[[m]][[i]][4,"95%Upper"])
            }else{
              NA
            }
          }), na.rm = T)
        })
        
      ##unpooled results
        model.stats$chao.up.coverage =sapply(1:56, function(m){
          mean(sapply(1:100, function(i){
            if(is.data.frame(chao.results.unpooled.reps[[m]][[i]])){
              models[m, "N"] >= round(chao.results.unpooled.reps[[m]][[i]][4,"95%Lower"]) &
                models[m, "N"] <= round(chao.results.unpooled.reps[[m]][[i]][4,"95%Upper"])
            }else{
              NA
            }
          }), na.rm = T)
        })
        
    #################################################################
    ## Remove scenarios where most of the Chao estimates failed
    ## because there were no singletons or doubletons (see Methods)
    #################################################################
      ##set the minimum number of successful estimates to keep the scenario
        min.to.keep = 10
          
      ##Pooled
        chao.pooled.NumEstimates = sapply(1:56, function(m){
          sum(sapply(1:100, function(i){
            is.data.frame(chao.results.pooled.reps[[m]][[i]])
          }))
        })
          
        chao.pooled.keep = as.numeric(chao.pooled.NumEstimates >= min.to.keep)
        chao.pooled.keep[chao.pooled.keep == 0] = NA
        
        model.stats$chao.p.coverage = model.stats$chao.p.coverage * chao.pooled.keep  
        
      ##Unpooled
        chao.unpooled.NumEstimates = sapply(1:56, function(m){
          sum(sapply(1:100, function(i){
            is.data.frame(chao.results.unpooled.reps[[m]][[i]])
          }))
        })
        
        chao.unpooled.keep = as.numeric(chao.unpooled.NumEstimates >= min.to.keep)
        chao.unpooled.keep[chao.unpooled.keep == 0] = NA
  
        model.stats$chao.up.coverage = model.stats$chao.up.coverage * chao.unpooled.keep

    ##########################################
    ## Plot each metric for the scenario with
    ## 60 sites and 125 species 
    ##########################################
      ##define the models to keep
        n.sites = 60
        N = 125
        var = "coverage" ##the variable to plot including coverage, bias, pre, acc
      
      ##get the model numbers associated with the parameters defined above
        m.2 = which(model.stats$N == N & model.stats$n.reps == 2 & model.stats$n.sites == n.sites)
        m.7 = which(model.stats$N == N & model.stats$n.reps == 7 & model.stats$n.sites == n.sites)
        
      ##combine all the results
        d.plot = cbind(model.stats[m.2, paste0("chao.up.", var)],
                       model.stats[m.2, paste0("chao.p.", var)],
                       model.stats[m.7, paste0("chao.p.", var)],
                       model.stats[m.2, paste0("occ.", var)],
                       model.stats[m.7, paste0("occ.", var)])
        
      ##plot
        par(mar = c(1.1, 2.1, 0.5, 0.1))
        
        ylim = c(min(d.plot) - abs(0.1*min(d.plot)),
                 max(d.plot) + 0.1*max(d.plot))
        
        bp = barplot(d.plot,
                     #ylim = ylim, ##all variable except coverage
                     ylim = c(0, 1), #coverage
                     beside = T,
                     col = c(rep("grey90", 21), rep("grey50", 14)),
                     axes = F)
        box()
        axis(side = 2)
        axis(side = 1, labels = F, at = bp)
        abline(h = 0.95, lty = "dashed") #coverage
        #abline(h = 0) #bias
        abline(v = 24.5)
        
    ######################################
    ## Plot bias, precision and accuracy
    ## for all scenarios together
    ######################################
      layout(matrix(1:4, nrow = 2))
      
      N = rep(c(20, 125), each = 2)
      n.sites = c(60, 30, 60, 30)
      var = "acc" #"bias" #"pre" "acc" "coverage"
      col = c("grey100", "grey80", "grey60", "grey30", "grey10")
      
      par(mar = c(1.1, 2.1, 0.5, 0.1))
      
      for(i in 1:4){
        d.filt = model.stats[model.stats$N == N[i] &
                               model.stats$n.sites == n.sites[i], ]
        d.stat = rbind(d.filt[d.filt$n.reps == 2, paste0("chao.up.", var)],
                       d.filt[d.filt$n.reps == 2, paste0("chao.p.", var)],
                       d.filt[d.filt$n.reps == 7, paste0("chao.p.", var)],
                       d.filt[d.filt$n.reps == 2, paste0("occ.", var)],
                       d.filt[d.filt$n.reps == 7, paste0("occ.", var)])
        barplot(d.stat,
                beside = T,
                col = col,
                #ylim = c(-0.15, 0.35)#bias
                #ylim = c(0, 0.52)#precision
                ylim = c(0, 0.52)
                
        )
        box()
      }
    
    ###############################################
    ## Plot coverage probability for all scenarios
    ###############################################
      ##make a plot function
        barplot.ci = function(N, n.sites, mar = c(2.1, 2.1, 2.1, 2.1)){
          ##get the model numbers
            m.2 = which(models$N == N & models$n.reps == 2 & models$n.sites == n.sites)
            m.7 = which(models$N == N & models$n.reps == 7 & models$n.sites == n.sites)
            
          ##combine all the results
            d = cbind(model.stats$chao.up.coverage[m.2],
                      model.stats$chao.p.coverage[m.2],
                      model.stats$chao.p.coverage[m.7],
                      model.stats$occ.coverage[m.2],
                      model.stats$occ.coverage[m.7])
            
          ##plot
            par(mar = mar)
            bp = barplot(d,
                         ylim = c(0,1.09),
                         cex.axis = 1.5,
                         beside = T,
                         col = c(rep("grey90", 21), rep("grey50", 14)),
                         axes = F)
            box()
            axis(side = 2, cex.axis = 1.25)
            axis(side = 1, labels = F, at = bp)
            abline(h = 0.95, lty = "dashed")
            abline(v = 24.5)
        }
      
      ##plot all the scenarios
        layout(matrix(1:4, nrow = 4, byrow = T))
        
        barplot.ci(N = 125, n.sites = 60, mar = c(1.1, 2.1, 2.5, 1.1))
        barplot.ci(N = 125, n.sites = 30, mar = c(2.5, 2.1, 1.1, 1.1))
        
        barplot.ci(N = 20, n.sites = 60, mar = c(1.1, 2.1, 2.5, 1.1))
        barplot.ci(N = 20, n.sites = 30, mar = c(2.5, 2.1, 1.1, 1.1))

  ###############################################
  ## Compare uncertatinty in the point estimates 
  ###############################################
    ##get the SE for every scenario
      occ.se = sapply(1:56, function(m){
        sapply(1:100, function(i){
          sd(occ.results[[m]][[i]]$N)
        })
      })
        
      chao.pooled.se = sapply(1:56, function(m){
        sapply(1:100, function(i){
          if(is.data.frame(chao.results.pooled.reps[[m]][[i]])){
            chao.results.pooled.reps[[m]][[i]][4, "   s.e."]     
          }else{
            NA
          }
        })
      })
        
      chao.unpooled.se = sapply(1:56, function(m){
        sapply(1:100, function(i){
          if(is.data.frame(chao.results.unpooled.reps[[m]][[i]])){
            chao.results.unpooled.reps[[m]][[i]][4, "   s.e."]     
          }else{
            NA
          }
        })
      })
        
    ##get the means
      mean(unlist(occ.se), na.rm = T)
      mean(unlist(chao.pooled.se), na.rm = T)
      mean(unlist(chao.unpooled.se), na.rm = T)
        
    ##t tests
      t.test(c(unlist(occ.se)),
             c(unlist(chao.pooled.se)),
             paired = T)
        
      t.test(c(unlist(occ.se)),
             c(unlist(chao.unpooled.se)),
             paired = T)
        
########################################
## Analyze addtional values of mu.psi
########################################
  ######################
  ## Simulate the data
  ######################
    ##parameter values
      N = 50
      n.reps = 3
      n.sites = 50
      sd.psi = 0.5
      mu.p = qlogis(0.5)
      sd.p = 0.5
      mu.psi = qlogis(seq(0.1, 0.9, by = 0.1))

    ##slimulate the data
      sim.results.2 = replicate(100,
                              sapply(mu.psi, sim.data,
                                     sd.psi = sd.psi,
                                     mu.p = mu.p,
                                     sd.p = sd.p,
                                     N = N,
                                     n.reps = n.reps,
                                     n.sites = n.sites))

    ##########################
    ## Analyze with Occupancy
    ##########################
      ##set up the parallel processing
        cl = makeCluster(6)
        registerDoParallel(cl)

      ##loop through the replicates
        occ.results.2 = foreach(i = 1:9, .packages = c("R2jags")) %dopar% {
          apply(sim.results.2[, i, 1:dim(sim.results.2)[3]], 2, occ, m = i)
        }

        stopCluster(cl)

        save(occ.results.2, sim.results.2,
             file = "OccupancyResults_muPsiTest.RData")

  #####################
  ## Analyze with Chao
  #####################
    ##analyze pooling all the reps
      chao.results.pooled.reps.2 = lapply(1:9, function(i){
        t(apply(sim.results.2[, i, 1:dim(sim.results.2)[3]], 2, chao))
      })
        
    ##analyze selecting one rep
      chao.results.unpooled.reps.2 = lapply(1:9, function(i){
        t(apply(sim.results.2[, i, 1:dim(sim.results.2)[3]], 2, chao, pooled.reps = F))
      })
    
  ##########################
  ## Summarize the results
  ##########################
    #################################
    ## Save the estimates in a table
    #################################
      ##occupancy
        occ.est.2 = sapply(1:9, function(m){ 
          sapply(1:100, function(i){
            mean(occ.results.2[[m]][[i]]$N)
          })
        })
      
      ##calculate the chao estimates for the pooled results
        chao.pooled.est.2 = sapply(1:9, function(m){
          sapply(1:100, function(i){
            if(is.data.frame(chao.results.pooled.reps.2[[m]][[i]])){
              chao.results.pooled.reps.2[[m]][[i]][4,"Estimate"]
            }else{
              sim.results.2[,m,i]$sim.N
            }
          })
        })
        
      ##unpooled chao  
        chao.unpooled.est.2 = sapply(1:9, function(m){
          sapply(1:100, function(i){
            if(is.data.frame(chao.results.unpooled.reps.2[[m]][[i]])){
              chao.results.unpooled.reps.2[[m]][[i]][4,"Estimate"]     
            }else{
              sim.results.2[,m,i]$sim.N
            }
          })
        })
      
    #######################
    ## Calculate the error
    #######################
      occ.error.2 = sapply(1:9, function(m){ 
        sapply(1:100, function(i){
          (mean(occ.results.2[[m]][[i]]$N) - N)
        })
      })
      
      pooled.error.2 = sapply(1:9, function(m){
        sapply(1:100, function(i){
          if(is.data.frame(chao.results.pooled.reps.2[[m]][[i]])){
            chao.results.pooled.reps.2[[m]][[i]][4,"Estimate"] -  N     
          }else{
            sim.results.2[,m,i]$sim.N - N
          }
        })
      })
      
      unpooled.error.2 = sapply(1:9, function(m){
        sapply(1:100, function(i){
          if(is.data.frame(chao.results.unpooled.reps.2[[m]][[i]])){
            chao.results.unpooled.reps.2[[m]][[i]][4,"Estimate"] -  N     
          }else{
            sim.results.2[,m,i]$sim.N - N
          }
        })
      })
      
    ###########################################################
    ## Save bias, precision, and accuracy for the 9 scenarios
    ###########################################################
      ##create a table to store the results
        model.stats.2 = data.frame(model = 1:9,
                                   mu.psi = plogis(mu.psi))
      
      ##save the mean estimate
        model.stats.2$occ.mean = apply(occ.est.2, 2, mean)
        model.stats.2$chao.p.mean = apply(chao.pooled.est.2, 2, mean)
        model.stats.2$chao.up.mean = apply(chao.unpooled.est.2, 2, mean)
      
      ##save the scaled bias
        model.stats.2$occ.bias = apply(occ.error.2, 2, mean)/N
        model.stats.2$chao.p.bias = apply(pooled.error.2, 2, mean)/N
        model.stats.2$chao.up.bias = apply(unpooled.error.2, 2, mean)/N
        
      ##save the scaled precision
        model.stats.2$occ.pre = apply(occ.est.2, 2, FUN = function(x){sd(x)/mean(x)})
        model.stats.2$chao.p.pre = apply(chao.pooled.est.2, 2, FUN = function(x){sd(x)/mean(x)})
        model.stats.2$chao.up.pre = apply(chao.unpooled.est.2, 2, FUN = function(x){sd(x)/mean(x)})
      
      ##save the accuracy
        model.stats.2$occ.acc = apply(occ.error.2, 2, FUN = function(x){sqrt(mean(x^2))})/N
        model.stats.2$chao.p.acc = apply(pooled.error.2, 2, FUN = function(x){sqrt(mean(x^2))})/N
        model.stats.2$chao.up.acc = apply(unpooled.error.2, 2, FUN = function(x){sqrt(mean(x^2))})/N
      
    ###################################################
    ## Save coverage probabililty for the 9 scenarios
    ###################################################
      ##occupancy
        model.stats.2$occ.coverage = sapply(1:9, function(m){
          mean(sapply(1:100, function(i){
            N >= round(quantile(occ.results.2[[m]][[i]]$N, 0.025)) &
            N <= round(quantile(occ.results.2[[m]][[i]]$N, 0.975))
          }))
        })
      
      
      ##pooled Chao
        model.stats.2$chao.p.coverage = sapply(1:9, function(m){
          mean(sapply(1:100, function(i){
            if(is.data.frame(chao.results.pooled.reps.2[[m]][[i]])){
              N >= round(chao.results.pooled.reps.2[[m]][[i]][4,"95%Lower"]) &
              N <= round(chao.results.pooled.reps.2[[m]][[i]][4,"95%Upper"])
            }else{
              NA
            }
          }), na.rm = T)
        })
      
      ##unpooled results
        model.stats.2$chao.up.coverage =sapply(1:9, function(m){
          mean(sapply(1:100, function(i){
            if(is.data.frame(chao.results.unpooled.reps.2[[m]][[i]])){
              N >= round(chao.results.unpooled.reps.2[[m]][[i]][4,"95%Lower"]) &
              N <= round(chao.results.unpooled.reps.2[[m]][[i]][4,"95%Upper"])
            }else{
              NA
            }
          }), na.rm = T)
        })
        
    #############################################################
    ## Remove scenarios where most of the Chao estimates failed
    ## because there were no singletons or doubletons
    #############################################################
      ##set the minimum number of successful estimates to keep the scenario
        min.to.keep = 10
      
      ##Pooled
        chao.pooled.NumEstimates.2 = sapply(1:9, function(m){
          sum(sapply(1:100, function(i){
            is.data.frame(chao.results.pooled.reps.2[[m]][[i]])
          }))
        })
      
        chao.pooled.keep.2 = as.numeric(chao.pooled.NumEstimates.2 >= min.to.keep)
        chao.pooled.keep.2[chao.pooled.keep.2 == 0] = NA
        
        model.stats.2$chao.p.coverage = model.stats.2$chao.p.coverage * chao.pooled.keep.2  
        
      ##Unpooled
        chao.unpooled.NumEstimates.2 = sapply(1:9, function(m){
          sum(sapply(1:100, function(i){
            is.data.frame(chao.results.unpooled.reps.2[[m]][[i]])
          }))
        })
        
        chao.unpooled.keep.2 = as.numeric(chao.unpooled.NumEstimates.2 >= min.to.keep)
        chao.unpooled.keep.2[chao.unpooled.keep.2 == 0] = NA
        
        model.stats.2$chao.up.coverage = model.stats.2$chao.up.coverage * chao.unpooled.keep.2
        
    ###############################################
    ## Plot bias, precision, accuracy, and coverage
    ###############################################
      var = "bias" ##change to plot different variable: bias, acc, pre, coverage
        
      ##combine all the results
        d.plot = cbind(model.stats.2[, paste0("chao.up.", var)],
                       model.stats.2[, paste0("chao.p.", var)],
                       model.stats.2[, paste0("occ.", var)])
        
      ##plot
        dev.off()
        par(mar = c(1.1, 2.1, 0.5, 0.1))
        
        ylim = c(min(d.plot) - abs(0.1*min(d.plot)),
                 max(d.plot) + 0.1*max(d.plot))
        
        bp = barplot(d.plot,
                     #ylim = ylim, ##all variables except coverage
                     #ylim = c(0,1), ##coverage
                     ylim = c(-0.001, 0.033),
                     beside = T,
                     col = c(rep("grey90", 18), rep("grey50", 9)),
                     axes = F)
        box()
        axis(side = 2)
        axis(side = 1, labels = F, at = bp)
        abline(h = 0.95, lty = "dashed") #coverage
        abline(h = 0) #bias
        abline(v = 20.5)