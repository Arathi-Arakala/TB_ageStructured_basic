####################################################################
## Create a basic 9 age class, 5 compartment TB model. 
## The model has S, L_a, L_b, I and R compartments.
## This model will be fitted to the Daw 1950 data 
## The data gives us the mortality due to TB in the 9 age classes
## Fitting will be done by allowing 3 parameters to vary: 
## (1) "mult" which defines a mixing matrix M (defines the rate of mixing between age classes) and
## (2) beta_o and tau, which define the contact rate beta modeled as a decaying exponential.
####################################################################


thiswd<-getwd()
source("mcmc_fitting_3params.R")
source("setGlobalVar.R") # this file sets the values for global variables


#retrieve the DAW table
daw<-getDAWtable()


# run the mcmc method to fit the mortality from the model 
# at each year n age class to the DAW table
# to fit to each year seprately, the second argument must be TRUE
# to fit to all the DAW data, the second argument must be FALSE
starttime<-timestamp()
M<-master_mcmc_runner_3params(5, FALSE) #second argument determines fitting to all DAW data points.
endtime<-timestamp()

# The best parameter is the one with the highest log likelihood
index<-which.max(M$log_lh)
M_best<-M[index,]
mult<-M_best$mult
beta_o<-M_best$beta_o
tau<-M_best$tau
#c(mult, beta_o, tau)

#plot the best beta curve and save the plot as a .pdf in "Plots" folder
plotBeta_decayExp(beta_o, tau,gamma,1)


#plot all the beta curves from accepted values
M_accepted<-M[which(M$accepted==1),]


# plot all the 90 points. 
# for windows OS, replace quartz() with windows()
for(t in 1:length(daw_years)){
  t_year<-daw_years[t]
  row_index<-which(daw_years==t_year)
  
  op<-ageStrTBModel_mortality_3params(mult, beta_o, tau, gamma, alpha, t_year)
  
  if(t==1){
    quartz()
    par(mfrow=c(1,1), oma=c(0,0,2,0))
    yrange<-range(c(daw$t1Total/2, op))
    plot(min_age, daw$t1Total[row_index,]/2, type='b', pch=16, lty=1,col=t, lwd=2, xlab="Age Groups", ylab="Annual death rates per million", main="Parameter fitting to ALL DAW data", sub="source:Daw 1950", ylim=yrange)
    lines(min_age, op, type='b', pch=16, lty=2, lwd=2, col=t)
    legend(x=40, y=5000, legend=c("DAW data", "Model output"), lty=c(1,2), lwd=2)
  }
  
  if(t>2){
    lines(min_age, daw$t1Total[row_index,]/2,type='b', pch=16, lty=1,col=t, lwd=2)
    lines(min_age, op, type='b', pch=16, lty=2, lwd=2, col=t)
  }
  
}

# plot the accepted beta_o and tau values.
quartz()
par(mfrow=c(2,1),oma=c(0,0,2,0))
plot(1:dim(M)[1], M$beta_o, type='p', pch=(M$accepted+16), col=3 )
plot(1:dim(M)[1], M$tau, type='p', pch=(M$accepted+16), col=4)
M$log_lh
M$accepted


##### Now that we got the best fit, let us use these parameters to see the disease dynamics
beta_o<-M_best$beta_o
tau<-M_best$tau
mult<-M_best$mult

# run to eqbm for the given beta_o n tau at the starting year 1860
beta1=getBeta_decayExp(beta_o, tau, gamma, alpha-1850)
#mult<-40
parameters<-getParameters(9, mult, beta1)
EqbmOutput<-IsEqbm(parameters)
eqbm<-EqbmOutput$eqbm

#use eqbm as starting conditions for ode run
nstart<-as.numeric(eqbm[-1])
# get the parameters for the ode model
parameters_timevarybeta<-getParameters_timevarybeta(9, mult, beta_o, tau, gamma, alpha)

#Run the TB model till 1940
# we use a model with a closed population to start with
# the mixing matrix is defined inside the model.
t_year<-1940
n_years<-t_year-1850 # start year
time<-seq(from=1, to=365*n_years, by=100 ) #timestep output monthly to speed up
output<-as.data.frame(ode(nstart,time,TBmodel_9ageclasses_timevarybeta_closedpopln,parameters_timevarybeta))

##### SOME PLOTS #######

# plot time series output by age compartment
# plotOutputByAge(output)
# 
# plot the age distribution for the different disease classes 
# plotDiseaseClassbyAge(output)
#
# This is the change in total population over time. 
# This should be flat for a closed population model
# plotTotalPopln(output)

#plot the distribution of disease classes in each age group as stacked bar plots
plotCompartmentalDisbnByAgeGroup(output, t_year)