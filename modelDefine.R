#########################################################
# This file defines the ode model and the parameters
#########################################################
source("plotting.R")

getParameters_timevarybeta<-function(Nclasses, mult, beta_o, tau, gamma, alpha){
  # function that returns the model parameters
  # The model has 9 age classes - Daw 1950 paper
  # The model has demography parameters based on geographic location
  # It has epi parameters based on geography and age
  # LOCATION IS ENGLAND FROM 1850 TO 1940 - DAW 1950 PAPER
  
  #define the time varying beta
  # beta_o and tau are the y intercept and the slope of decay
  # gamma is the floor
  # alpha is the year at which beta = beta_o
  # alpha=1860/1850 etc
  # gamma = 0.01
  t_values<-seq(from=1850, to=2020, by=10)
  t_values_timeshift_daily<-(t_values-1850)*365
  beta<-(beta_o * exp( - tau * ( t_values_timeshift_daily))) + gamma
  beta_func<-approxfun(t_values_timeshift_daily, beta, method="linear")
  # quartz()
  # plot(t_values_timeshift_daily, beta_func(t_values_timeshift_daily), type='b', lwd=1, pch=16, ylim=c(0,beta_o+gamma+0.01) )
  
  
  
  # t_values<-seq(from=1850, to=2020, by=10)
  # t_values<-t_values-alpha
  # t_values<-t_values*365
  # beta<-( beta_o * exp( - tau * ( t_values))  ) + gamma
  # beta_func<-approxfun(t_values, beta, method="linear") # valid time value inputs: from -10*365 to 160*365
  # quartz()
  # plot(t_values, beta_func(t_values), type='b', lwd=1, pch=16, ylim=c(0,beta_o+gamma+0.01) )
  # abline(h=gamma, lty=2, lwd=1)
  parameters<-list(
    ### demography parameters ####
    #################################################
    pi=0.0193/365, #birth rate into England population in 1850 *<>
    #mu = rep((1/70)/365, Nclasses), # standard death rate, in England in 1850, if average life expectancy is 70 *<>
    mu = c( rep((1/70)/365, Nclasses-1),(10/70)/365 ),
    tau_in = c(0,1/(5*365),rep(1/(10*365), Nclasses-1) ), # ageing parameters
    tau_out = c(1/(5*365), rep(1/(10*365), Nclasses-1), 0),# ageing parameters
    
    ### disease parameters that vary with geography ###
    ####################################################
    mu_I = 0.5/(3*365), # daily death rate due to disease in India, WHO world TB report
    #mu_T = (0.15*2)/365, # death rate during treatment, half of mu_I
    #delta = 0.59*73/(0.41*365*70*3), # based on CDR 2015, glbal TB report for India, probability of treatment=0.59 , will be typically time dependent
    #delta =0,
    # delta=0*73/(0.41*365*70*3),# prob of treatment=0, no treatment case
    
    ### disease based parameters ####################
    #################################################
    #epsilon = rep(0.0011, times=3), #rate of progression to active" disease from early latency per day, 
    epsilon = c(6.6E-03, 2.7E-03, rep(2.7E-04, times=Nclasses-2)), #from Romain's paper "optimally......."
    #kappa = rep(0.01, times=3), # rate of progression from early to late latency
    kappa = c(1.2E-02 ,1.2E-02 , rep(5.4E-03, times=Nclasses-2) ), # rate of progression from early to late latency
    #nu = rep(5.5e-6, times=3) , #rate of progression to active disease from late latency per day.
    nu = c( 1.9E-11,6.4E-06, rep(3.3E-06, times=Nclasses-2) ), #rate of progression to active disease from late latency per day.
    #iota = c(1,1,1) , # same for all age classes
    iota = c(0, 0, rep(1, times=Nclasses-2)  ), # change in infectiousness compared to an adult for the 3 age classes
    
    gamma = 0.5/(3*365), # rate of spontaneous cure, Dowdy et al; 
    #o=0, # fraction of treated individuals contributing to the transmission rate
    phi= (0.74*2)/365, # rate of spontaneous recovery, per year. is 74% from golabl TB report for India
    #omega = (0.11*2)/365, # probability of defaulting treatment
    beta=10*24/365, #contact rate of infection
    beta_func=beta_func, #the beta_func will create the correct beta_mx inside the ode model
    mult=mult, # the multiplier for infection within younger ageclass.
    chi = 0.21, # fractional reduction in the force of infection corresponding to return from late to early latency
    alpha = 0.21, # fractional reduction in force of infection due to recovery from first infection
    rho = 0.3 # proportion of infected people who are infectious. avg over years, see Notifications table , WHO TB report.
  )
  
}


# this is the TB model with births at the rate pi
TBmodel_9ageclasses_timevarybeta<-function(t,n,parameters){
  
  with(as.list(parameters),{
    #read in parameters
    pi<-parameters$pi
    mu<-parameters$mu
    tau_in<-parameters$tau_in
    tau_out<-parameters$tau_out
    mu_I<-parameters$mu_I
    epsilon<-parameters$epsilon
    kappa<-parameters$kappa
    nu<-parameters$nu
    gamma<-parameters$gamma
    phi<-parameters$phi
    beta<-parameters$beta
    beta_func<-parameters$beta_func
    chi<-parameters$chi
    alpha<-parameters$alpha
    iota<-parameters$iota
    rho<-parameters$rho
    
    
    # read in intial values of state variables
    S_u<-n[1:9] #number of people unvaccinated
    L_a<-n[10:18] #number in early latency
    L_b<-n[19:27] #number in late latency
    I<-n[28:36] # number infected
    S_r<-n[37:45] # number of people recovered after treatment
    
    N<-rep(0, times=9) #total population in each age class
    N_tot<-sum(S_u)+sum(L_a)+sum(L_b)+sum(I)+sum(S_r)
    dS_udt<-rep(0, times=9)
    dL_adt<-rep(0, times=9)
    dL_bdt<-rep(0, times=9)
    dIdt<-rep(0, times=9)
    dS_rdt<-rep(0, times=9)
    #calculate the new lamda for age structured model
    #lamda=rho*(beta* (sum(I)+(o*sum(T_r)))/ N_tot ) # force of infection, homogenous mixing
    
    # get the correct time specific beta and set up beta_mx
    #in beta_func time is in days, in the ode time step is in days, (AA: testing-divide by 365 will get it into the right unit?)
    beta1<-beta_func(t) #/365
    beta_mx=matrix(c(mult*beta1, mult*beta1, rep(beta1, times=7), mult*beta1, mult*beta1, rep(beta1, times=7), beta1, beta1, rep(mult*beta1, times=7), beta1, beta1, rep(mult*beta1, times=7), beta1, beta1, rep(mult*beta1, times=7), beta1, beta1, rep(mult*beta1, times=7), beta1, beta1, rep(mult*beta1, times=7), beta1, beta1, rep(mult*beta1, times=7), beta1, beta1, rep(mult*beta1, times=7)), nrow = 9, ncol = 9)
    
    
    lamda<-0
    for(a in 1:9){
      N[a]=S_u[a]+L_a[a]+L_b[a]+I[a]+S_r[a] #total population in age class a
      incr=iota[a]*rho*beta*(I[a])/N_tot
      lamda=lamda+incr
    }
    
    # calculate lamda as a vector of length 9, one for each class
    lamda_vec<-rep(0, times=9)
    for(l in 1:9){
      lamda_vec[l]<-0
      for(a in 1:9){
        #N[a]=S_v[a]+S_u[a]+L_a[a]+L_b[a]+I[a]+T_r[a]+S_r[a] #total population in age class a
        incr=iota[a]*rho*beta_mx[l,a]*(I[a])/N_tot
        lamda_vec[l]=lamda_vec[l]+incr
      }
    }
    
    
   
    # using lamda_vec
    for(a in 1:9){
      if(a==1){
        dS_udt[a]=(pi*N_tot) - (mu[a]*S_u[a]) - lamda_vec[a]*S_u[a] - ((tau_out[a])*S_u[a])
        dL_adt[a]=lamda_vec[a]*S_u[a] + (alpha*lamda_vec[a]*S_r[a]) + (chi*lamda_vec[a]*L_b[a]) - (mu[a]+epsilon[a]+kappa[a])*L_a[a] - ((tau_out[a])*L_a[a])
        dL_bdt[a]= kappa[a]*L_a[a] + gamma*I[a] -(mu[a]+nu[a]+(chi*lamda_vec[a]) )*L_b[a] - ((tau_out[a])*L_b[a])
        dIdt[a]= nu[a]*L_b[a] + epsilon[a]*L_a[a] - (mu_I+gamma+phi) * I[a] - ((tau_out[a])*I[a])
        dS_rdt[a]= phi*I[a]- ( (alpha*lamda_vec[a]) + mu[a]) * S_r[a] - ((tau_out[a])*S_r[a])
      }
      if(a%in% 2:8){
        dS_udt[a]= - (mu[a]*S_u[a]) - lamda_vec[a]*S_u[a] +(tau_in[a])*S_u[a-1] - (tau_out[a])*S_u[a]
        dL_adt[a]=lamda_vec[a]*S_u[a] + (alpha*lamda_vec[a]*S_r[a]) + (chi*lamda_vec[a]*L_b[a]) - (mu[a]+epsilon[a]+kappa[a])*L_a[a] + (tau_in[a])*L_a[a-1] - (tau_out[a])*L_a[a]
        dL_bdt[a]= kappa[a]*L_a[a] + gamma*I[a] -(mu[a]+nu[a]+(chi*lamda_vec[a]) )*L_b[a] + (tau_in[a])*L_b[a-1] - (tau_out[a])*L_b[a]
        dIdt[a]= nu[a]*L_b[a] + epsilon[a]*L_a[a] - (mu_I+gamma+phi) * I[a] + (tau_in[a])*I[a-1] - (tau_out[a])*I[a]
        dS_rdt[a]= phi*I[a] - ( (alpha*lamda_vec[a]) + mu[a]) * S_r[a] + (tau_in[a])*S_r[a-1] - (tau_out[a])*S_r[a]
      }
      
      if(a==9){
        dS_udt[a]= - (mu[a]*S_u[a]) - lamda_vec[a]*S_u[a] + (tau_in[a])*S_u[a-1]
        dL_adt[a]= lamda_vec[a]*S_u[a] + (alpha*lamda_vec[a]*S_r[a]) + (chi*lamda_vec[a]*L_b[a]) - (mu[a]+epsilon[a]+kappa[a])*L_a[a] + (tau_in[a])*L_a[a-1]
        dL_bdt[a]= kappa[a]*L_a[a] + gamma*I[a] -(mu[a]+nu[a]+(chi*lamda_vec[a]) )*L_b[a] + (tau_in[a])*L_b[a-1]
        dIdt[a]= nu[a]*L_b[a] + epsilon[a]*L_a[a] - (mu_I+gamma+phi) * I[a] + (tau_in[a])*I[a-1]
        dS_rdt[a]= phi*I[a] - ( (alpha*lamda_vec[a]) + mu[a]) * S_r[a] + (tau_in[a])*S_r[a-1]
        
      }
      
    }
    
    return(list(c(dS_udt, dL_adt,dL_bdt, dIdt,dS_rdt )))
  })
}

# TB model with closed popln, make all the deaths go in as new births.
TBmodel_9ageclasses_timevarybeta_closedpopln<-function(t,n,parameters){
  
  with(as.list(parameters),{
    #read in parameters
    pi<-parameters$pi
    mu<-parameters$mu
    tau_in<-parameters$tau_in
    tau_out<-parameters$tau_out
    mu_I<-parameters$mu_I
    epsilon<-parameters$epsilon
    kappa<-parameters$kappa
    nu<-parameters$nu
    gamma<-parameters$gamma
    phi<-parameters$phi
    beta<-parameters$beta
    beta_func<-parameters$beta_func
    chi<-parameters$chi
    alpha<-parameters$alpha
    iota<-parameters$iota
    rho<-parameters$rho
    mult<-parameters$mult
    
    
    # read in intial values of state variables
    S_u<-n[1:9] #number of people unvaccinated
    L_a<-n[10:18] #number in early latency
    L_b<-n[19:27] #number in late latency
    I<-n[28:36] # number infected
    S_r<-n[37:45] # number of people recovered after treatment
    
    N<-rep(0, times=9) #total population in each age class
    N_tot<-sum(S_u)+sum(L_a)+sum(L_b)+sum(I)+sum(S_r)
    dS_udt<-rep(0, times=9)
    dL_adt<-rep(0, times=9)
    dL_bdt<-rep(0, times=9)
    dIdt<-rep(0, times=9)
    dS_rdt<-rep(0, times=9)
    #calculate the new lamda for age structured model
    #lamda=rho*(beta* (sum(I)+(o*sum(T_r)))/ N_tot ) # force of infection, homogenous mixing
    
    # get the correct time specific beta and set up beta_mx
    #in beta_func time is in days, in the ode time step is in days, (AA: testing-divide by 365 will get it into the right unit?)
    beta1<-beta_func(t) #/365
    #beta_mx=matrix(c(mult*beta1, mult*beta1, rep(beta1, times=7), mult*beta1, mult*beta1, rep(beta1, times=7), beta1, beta1, rep(mult*beta1, times=7), beta1, beta1, rep(mult*beta1, times=7), beta1, beta1, rep(mult*beta1, times=7), beta1, beta1, rep(mult*beta1, times=7), beta1, beta1, rep(mult*beta1, times=7), beta1, beta1, rep(mult*beta1, times=7), beta1, beta1, rep(mult*beta1, times=7)), nrow = 9, ncol = 9, byrow = TRUE)
    # beta_mx=matrix(c(mult*beta1, mult*beta1, rep(beta1, times=7), 
    #                  mult*beta1, mult*beta1, rep(beta1, times=7), 
    #                  beta1, beta1, rep(mult*beta1, times=3),rep(beta1, times=4), 
    #                  beta1, beta1, rep(mult*beta1, times=3),rep(beta1, times=4), 
    #                  beta1, beta1, rep(mult*beta1, times=3),rep(beta1, times=4), 
    #                  beta1, beta1, rep(mult*beta1, times=3),rep(beta1, times=4), 
    #                  rep(beta1, times=9), 
    #                  rep(beta1, times=9), 
    #                  rep(beta1, times=9)), 
    #                nrow = 9, ncol = 9, byrow = TRUE)
    
    #staggered drop for older age group interactions
    beta_mx=matrix(c(mult*beta1, mult*beta1, rep(beta1, times=7), 
                     mult*beta1, mult*beta1, rep(beta1, times=7), 
                     beta1, beta1, rep(mult*beta1, times=4),rep(beta1, times=3), 
                     beta1, beta1, rep(mult*beta1, times=4),rep(beta1, times=3), 
                     beta1, beta1, rep(mult*beta1, times=4),rep(beta1, times=3), 
                     beta1, beta1, rep(mult*beta1, times=4),rep(beta1, times=3), 
                     rep(beta1, times=2), rep( (mult/2)*beta1, times=7),
                     rep(beta1, times=2), rep((mult/4)*beta1, times=7),
                     rep(beta1, times=2), rep((mult/8)*beta1, times=7) ),
                   nrow = 9, ncol = 9, byrow = TRUE)
    
    
    lamda<-0
    for(a in 1:9){
      N[a]=S_u[a]+L_a[a]+L_b[a]+I[a]+S_r[a] #total population in age class a
      incr=iota[a]*rho*beta*(I[a])/N_tot
      lamda=lamda+incr
    }
    
    # calculate lamda as a vector of length 9, one for each class
    lamda_vec<-rep(0, times=9)
    for(l in 1:9){
      lamda_vec[l]<-0
      for(a in 1:9){
        #N[a]=S_v[a]+S_u[a]+L_a[a]+L_b[a]+I[a]+T_r[a]+S_r[a] #total population in age class a
        incr=iota[a]*rho*beta_mx[l,a]*(I[a])/N_tot
        lamda_vec[l]=lamda_vec[l]+incr
      }
    }
    
    #total death using N, when age soecific death rates occur
    tot_death<-0
    for(a in 1:9){
      tot_death<-tot_death+(mu[a]*N[a])
    }
    tot_death<-tot_death+(sum(I)*mu_I)
    
    # using lamda_vec
    for(a in 1:9){
      if(a==1){
        #AA: the deaths are all coming in as births. first term is possible as all death rates are equal.
        ##AA: (mu[a]*N_tot) + (sum(I)*mu_I) -this is tot_death when all death rates are same
        dS_udt[a]=( tot_death ) - (mu[a]*S_u[a]) - lamda_vec[a]*S_u[a] - ((tau_out[a])*S_u[a])
        dL_adt[a]=lamda_vec[a]*S_u[a] + (alpha*lamda_vec[a]*S_r[a]) + (chi*lamda_vec[a]*L_b[a]) - (mu[a]+epsilon[a]+kappa[a])*L_a[a] - ((tau_out[a])*L_a[a])
        dL_bdt[a]= kappa[a]*L_a[a] + gamma*I[a] -(mu[a]+nu[a]+(chi*lamda_vec[a]) )*L_b[a] - ((tau_out[a])*L_b[a])
        dIdt[a]= nu[a]*L_b[a] + epsilon[a]*L_a[a] - (mu_I+gamma+phi+mu[a]) * I[a] - ((tau_out[a])*I[a])
        dS_rdt[a]= phi*I[a]- ( (alpha*lamda_vec[a]) + mu[a]) * S_r[a] - ((tau_out[a])*S_r[a])
      }
      if(a%in% 2:8){
        dS_udt[a]= - (mu[a]*S_u[a]) - lamda_vec[a]*S_u[a] +(tau_in[a])*S_u[a-1] - (tau_out[a])*S_u[a]
        dL_adt[a]=lamda_vec[a]*S_u[a] + (alpha*lamda_vec[a]*S_r[a]) + (chi*lamda_vec[a]*L_b[a]) - (mu[a]+epsilon[a]+kappa[a])*L_a[a] + (tau_in[a])*L_a[a-1] - (tau_out[a])*L_a[a]
        dL_bdt[a]= kappa[a]*L_a[a] + gamma*I[a] -(mu[a]+nu[a]+(chi*lamda_vec[a]) )*L_b[a] + (tau_in[a])*L_b[a-1] - (tau_out[a])*L_b[a]
        dIdt[a]= nu[a]*L_b[a] + epsilon[a]*L_a[a] - (mu_I+gamma+phi+mu[a]) * I[a] + (tau_in[a])*I[a-1] - (tau_out[a])*I[a]
        dS_rdt[a]= phi*I[a] - ( (alpha*lamda_vec[a]) + mu[a]) * S_r[a] + (tau_in[a])*S_r[a-1] - (tau_out[a])*S_r[a]
      }
      
      if(a==9){
        dS_udt[a]= - (mu[a]*S_u[a]) - lamda_vec[a]*S_u[a] + (tau_in[a])*S_u[a-1]
        dL_adt[a]= lamda_vec[a]*S_u[a] + (alpha*lamda_vec[a]*S_r[a]) + (chi*lamda_vec[a]*L_b[a]) - (mu[a]+epsilon[a]+kappa[a])*L_a[a] + (tau_in[a])*L_a[a-1]
        dL_bdt[a]= kappa[a]*L_a[a] + gamma*I[a] -(mu[a]+nu[a]+(chi*lamda_vec[a]) )*L_b[a] + (tau_in[a])*L_b[a-1]
        dIdt[a]= nu[a]*L_b[a] + epsilon[a]*L_a[a] - (mu_I+gamma+phi+mu[a]) * I[a] + (tau_in[a])*I[a-1]
        dS_rdt[a]= phi*I[a] - ( (alpha*lamda_vec[a]) + mu[a]) * S_r[a] + (tau_in[a])*S_r[a-1]
        
      }
      
    }
    
    return(list(c(dS_udt, dL_adt,dL_bdt, dIdt,dS_rdt )))
  })
}

