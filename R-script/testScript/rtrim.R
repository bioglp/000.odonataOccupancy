library(rtrim)

data(skylark2) # use extended version of the Skylark dataset
head(skylark2)

summary(skylark2)
count_summary(skylark2, year_col="season") # show that it works

z1 <- trim(count ~ site + year, data=skylark2, 
           model=3, serialcor=TRUE, overdisp=TRUE)
summary(z1)

plot(overall(z1))

# covariate
z2 <- trim(count ~ site + year + habitat, 
           data=skylark2, model=3, 
           serialcor=TRUE, overdisp=TRUE)
summary(z2)

plot(index(z2,which="fitted",covars=TRUE))
