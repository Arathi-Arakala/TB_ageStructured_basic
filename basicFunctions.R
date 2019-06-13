#### This file will contain a list of commonly used functions ####
##################################################################

source("modelDefine.R")

### 1. Function to display the table in DAW paper #####
## Excerpt from DAW1950
## "The data used are set out in Table I and consist of the 
## mortality rates in England and Wales for all forms of tuberculosis. 
## For the years 1851-1920, the rates were taken or derived from 
## Table 12 of the Registrar-General’sDecennial Supplement, 1921, Part III. 
## For 1921-40, the rates have been calculated from the numbers of deaths 
## and mean populations given in the annual volumes of the 
## Registrar-General's StatisticRaelview."


getDAWtable<-function(){
  table1Male<-data.frame(age0=c(6323,6018,5798,5004,4347,3129,1942,1067,612), age5=c(1166,967,828,727,615,552,545,321,184),
                         age15=c(3400, 3157,2493,1976,1641,1353, 1299,1083, 817), age25=c(4163, 4206, 3785, 3164, 2541,2158, 1840,1413, 1009),
                         age35=c(4119, 4244, 4198, 3685, 3251, 2622, 2204, 1645, 1134), age45=c(3957, 3969,3928,3611, 3296,2934,2335,1729,1390),
                         age55=c(3479, 3433, 3285, 3027, 2768, 2574, 2135, 1437, 1287), age65=c(2573, 2174, 2025, 1913, 1706, 1686, 1390, 972, 808),
                         age75=c(1061, 740, 650, 732,629, 668, 585, 398, 353))
  table1Female<-data.frame(age0=c(5232, 4917, 4663, 3987, 3516, 2636, 1619, 881, 522), age5=c(1388, 1109, 956, 949, 780, 704, 682, 408, 214), 
                           age15=c(4079, 3689, 2907, 2267, 1670, 1338, 1470, 1353, 1061), age25=c(4690, 4482, 3631, 2932, 2086, 1651, 1484, 1230, 929),
                           age35=c(4293, 3988, 3475, 2846, 2264, 1710, 1401, 936, 639), age45=c(3236, 2954, 2535, 2146, 1753, 1449, 1156, 729, 481), 
                           age55=c(2523, 2178, 1866, 1597, 1344, 1186, 943, 617, 423), age65=c(1783, 1354, 1193, 1058, 906, 894, 750, 511, 358),
                           age75=c(834, 528, 452, 452, 427, 494, 437, 326, 236))
  table1Total<-table1Male+table1Female
  
  list(t1Male=table1Male, t1Female=table1Female, t1Total=table1Total)
}
# table<-table1Male
# plotTablePerCY(table, "male")
# plotTableMeanBirthYear(table, "male")

############################################################################

getEqbmMortalitybyAgeClass<-function(output_pre, mu_I){
  eqbm<-output_pre[dim(output_pre)[1], ]
  
  #calculate total population in each age group
  N_tot<-rep(0, times=9)
  for(i in 2:10){
    N_tot[i-1]<-sum(eqbm[c(i,i+9, i+18, i+27, i+36)])
  }
  
  
  #columns 29:37 are the Infected classes of the 9 age classes
  modelMortality<-(eqbm[29:37]/N_tot)*mu_I*10^6
  modelMortality
}

ageStrTBModel_mortality<-function(mult, beta1){
  nstart=c(rep(3*10^6, times=9), rep(0, times=9), rep(0, times=9), c(0,1000, rep(0, times=7)), rep(0, times=9))
  # mult<-40
  # beta1=3*24/365
  parameters<-getParameters(9, mult, beta1)
  time<-seq(from=1, to=365*20, by=1 )
  
  output_pre<-as.data.frame(ode(nstart,time,TBmodel_9ageclasses,parameters))
  eqbm<-output_pre[dim(output_pre)[1], ]
  modelOP<-as.numeric(getEqbmMortalitybyAgeClass(output_pre, parameters$mu_I))
  modelOP*365
}
# mult=40
# beta1<-3*24/365
# ageStrTBModel_mortality(mult, beta1)
# 

# beta_o and tau are parameters of the time varying beta function, diff to tau_in and tau_out
#This function returns the age specific mortality for a given parameter set and goven year t_year
ageStrTBModel_mortality_3params<-function(mult, beta_o, tau, gamma, alpha, t_year){
 
  # run to eqbm for the given beta_o n tau
  beta1=getBeta_decayExp(beta_o, tau, gamma, t_year-1850)
  #mult<-40
  parameters<-getParameters(9, mult, beta1)
  EqbmOutput<-IsEqbm(parameters)
  eqbm<-EqbmOutput$eqbm
  
  #use eqbm as starting conditions for ode run
  nstart<-as.numeric(eqbm[-1])
  parameters_timevarybeta<-getParameters_timevarybeta(9, mult, beta_o, tau, gamma, alpha)
  
  
  n_years<-t_year-1850 # start year
  time<-seq(from=1, to=365*n_years, by=100 ) #timestep output monthly to speed up
  #run ode and get age specific mortality for the chosen t_year
  #output_pre<-as.data.frame(ode(nstart,time,TBmodel_9ageclasses_timevarybeta,parameters_timevarybeta))
  output_pre<-as.data.frame(ode(nstart,time,TBmodel_9ageclasses_timevarybeta_closedpopln,parameters_timevarybeta))
  
  modelOP<-as.numeric(getEqbmMortalitybyAgeClass(output_pre, parameters_timevarybeta$mu_I))
  modelOP*365
}

# beta_o and tau are parameters of the time varying beta function, diff to tau_in and tau_out
#This function returns the age specific mortality for a given parameter set 
#output is for all years.
#output should be a matrix of mortalities over Daw years versus age
ageStrTBModel_mortality_3params_allDawYears<-function(mult, beta_o, tau, gamma, alpha){

    # run to eqbm for the given beta_o n tau at the starting year 1860
  beta1=getBeta_decayExp(beta_o, tau, gamma, alpha-1850)
  #mult<-40
  parameters<-getParameters(9, mult, beta1)
  EqbmOutput<-IsEqbm(parameters)
  eqbm<-EqbmOutput$eqbm
  
  #use eqbm as starting conditions for ode run
  nstart<-as.numeric(eqbm[-1])
  parameters_timevarybeta<-getParameters_timevarybeta(9, mult, beta_o, tau, gamma, alpha)
  
  #rows are daw years, 1860 to 1940
  #columns are 9 age classes
  modelOP_all<-matrix(0, nrow=9, ncol=9)
  daw_years<-seq(from=1860, to=1940, by=10)
  
  for(t in 1:length(daw_years) ){
    t_year<-daw_years[t]
    n_years<-t_year-1850 # start year
    time<-seq(from=1, to=365*n_years, by=100 ) #timestep output monthly to speed up
    #run ode and get age specific mortality for the chosen t_year
    #output_pre<-as.data.frame(ode(nstart,time,TBmodel_9ageclasses_timevarybeta,parameters_timevarybeta))
    output_pre<-as.data.frame(ode(nstart,time,TBmodel_9ageclasses_timevarybeta_closedpopln,parameters_timevarybeta))
    
    modelOP<-as.numeric(getEqbmMortalitybyAgeClass(output_pre, parameters_timevarybeta$mu_I))
   
    modelOP_all[t,]<-modelOP*365

  }
  modelOP_all
 
}
#   
## This function returns the value of beta given the starting value and the decay
## It defines beta as a decaying exponential function.
## t will be a day from 1850 to 1940, day1 of 1850 is t=0, and increases daily from there.
## gamma is the value that beta will asymtotically converge to
getBeta_decayExp<-function(beta_o, tau,gamma, t){
#beta_o<-0.1001
#tau<-0.01
  t_values<-seq(from=1850, to=1940, by=10)
  t_values_timeshift_daily<-(t_values-1850)*365
  beta<-(beta_o * exp( - tau * ( t_values_timeshift_daily))) + gamma
  beta_func<-approxfun(t_values_timeshift_daily, beta, method="linear")
  return(beta_func(t))
}


#here alpha1 refers to the intensity of re-infection after first infection. Don't confuse with alpha here, which is asssociated with the time varying beta function
ageStrTBModel_mortality_4params<-function(mult, beta_o, tau, alpha1, gamma, alpha, eqbm, t_year){
  # nstart=c(rep(3*10^6, times=9), rep(0, times=9), rep(0, times=9), c(0,1000, rep(0, times=7)), rep(0, times=9))
  nstart<-as.numeric(eqbm[-1])
  # mult<-40
  # beta1=3*24/365
  # t_year can be any year, to compare with DAW table it could be 1860 (column 1) to 1940(column 9), in steps of 10.
  # beta1<-getBeta_decayExp(beta_o, tau, t_year)
  # parameters<-getParameters(9, mult, beta1)
  #parameters_timevarybeta<-getParameters_timevarybeta(9, mult, beta_o, tau, gamma, alpha)
  parameters_timevarybeta<-getParameters_timevarybeta_alpha(9, mult, beta_o, tau, alpha1, gamma, alpha) #alpha1 is the alpha we refer to
  
  n_years<-t_year-1860
  time<-seq(from=1, to=365*n_years, by=1 )
  
  output_pre<-as.data.frame(ode(nstart,time,TBmodel_9ageclasses_timevarybeta,parameters_timevarybeta))
  modelOP<-as.numeric(getEqbmMortalitybyAgeClass(output_pre, parameters_timevarybeta$mu_I))
  modelOP*365
}
#   
## This function returns the value of beta given the starting value and the decay
# ## It defines beta as a decaying exponential function.
# ## t will be a year from 1850 to 1940
# getBeta_decayExp<-function(beta_o, tau, t){
#   #beta_o<-0.1001
#   #tau<-0.01
#   t_values<-seq(from=1850, to=1940, by=10)
#   beta<-beta_o * exp( - tau * ( t_values-1860))
#   beta_func<-approxfun(t_values-1860, beta, method="linear")
#   return(beta_func(t))
#   
#   
#   quartz()
#   plot(t_values-1860, beta_func(t_values-1860), type='b', lwd=2, pch=16, lty=1, main="time varying beta")
#   points(t, beta_func(t), pch=18, col="red")
# }
# 
# This function checks if we have reached equilibrium. 
# t_period is the time period over which we check the population change, expressed in days.
# change_limit is the largest change in population over the time period that
# can be considered as eqbm. Expressed as a percentage.

checkEqbm<-function(output_pre, t_period, change_limit){
  isEqbm_flag<-FALSE
  t_max<-output_pre[dim(output_pre)[1],1]
  output_allCompartments<-matrix(0, nrow = dim(output_pre)[1], ncol = 5)
  CompIndex<- c(2,11,20,29,38)
  for(c in 1:length(CompIndex))
    output_allCompartments[,c]<-as.matrix(rowSums(output_pre[,c:(c+8)]) )
  
  total_pop<-rowSums(output_allCompartments)
  
  Infected<-output_allCompartments[,4]/total_pop
    
  percentage_change<-100*( abs(Infected[dim(output_pre)[1]] - Infected[ dim(output_pre)[1] - t_period])/Infected[ dim(output_pre)[1] - t_period])
  
  if(percentage_change<change_limit)
    isEqbm_flag <- TRUE
  
  isEqbm_flag  

}


# function that runs an ode model to equilibrium.
IsEqbm<-function(parameters){
  year_step<-5
  nstart=c(rep(3*10^6, times=9), rep(0, times=9), rep(0, times=9), c(0,1000, rep(0, times=7)), rep(0, times=9))
  time<-seq(from=1, to=365*year_step, by=1 )
  
  ### Run the model to equilibrium by running in steps of "year_step" years
  ### year_step is 5, running in steps of 5 years
  ### t_period to check for eqbm is 365 days
  ### change_limit as tol for change is 1%
  t_period<-365
  change_limit<-1
  eqbmFlag<-FALSE
  eqbm<-numeric()
  t_last<-numeric()
  output_pre_all<-numeric()
  repeat{
    output_pre<-as.data.frame(ode(nstart,time,TBmodel_9ageclasses,parameters))
    eqbmFlag<-checkEqbm(output_pre, t_period, change_limit)
    eqbm<-output_pre[dim(output_pre)[1],]
    t_last<-as.numeric(eqbm[1])
    output_pre_all<-rbind(output_pre_all, output_pre)
    tail(output_pre_all[,1])
    
    time<-seq(from=t_last+1, to=t_last+(365*year_step), by=1)
    nstart<-as.numeric(eqbm[-1])
    if(eqbmFlag==TRUE) 
      break
  }
  
  
 list(eqbm=eqbm, output_pre_all=output_pre_all)
}



