##################################################
### SET GLOBAL VARIABLES IN THIS FILE ############

#the lower limit of each age bracket of the age structured model
min_age<-c(0, seq(from=5, to=75, by=10))

#the years over which the 1950 DAW paper reports data
daw_years<-seq(from=1860, to=1940, by=10)

# beta curve parameters
# gamma: is the asymptote to the beta curve or the floor.
# alpha: is the year at which beta = beta_o, the y intercept.
# alpha=1860/1850 etc
gamma<-0.0001
alpha<- 1860