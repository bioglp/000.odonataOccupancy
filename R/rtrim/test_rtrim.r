library(rtrim)
data(skylark2)
summary(skylark2)

# idx <- which(names(skylark2)=="year")      # rename year->season
# names(skylark2)[idx] <- "season"
# count_summary(skylark2, year_col="season") # show that it works

z1 <- trim(count ~ site + season, 
  data=skylark2, 
  model=3, serialcor=TRUE, overdisp=TRUE)
summary(z1)

plot(overall(z1))

z2 <- trim(count ~ site + year + habitat, 
  data=skylark2, 
  model=3, serialcor=TRUE, overdisp=TRUE)
summary(z2)

z3 <- trim(count ~ site + year + habitat, 
          data=skylark2, model=2, changepoints="all",
           serialcor=TRUE, overdisp=TRUE)
summary(z3)
wald(z3)

z4 <- trim(count ~ site + year + habitat, 
    data=skylark2, model=2, changepoints="all",
           stepwise=TRUE, serialcor=TRUE, overdisp=TRUE)
summary(z4)

z6 <- trim(count ~ site + year + habitat, data=skylark2, 
          model=2, changepoints="auto",
          serialcor=TRUE, overdisp=TRUE, weights="weight")
idx = index(z6, "fitted", covars=TRUE)
plot(idx)
