require(deSolve)
require(graphics)
require(ggplot2)
require(xlsx)

source("basicFunctions.R")

daw<-getDAWtable()
my_data<-daw$t1Total[1,]/2

metropolis_hastings <- function(log_lh_func, log_priors_func, proposal_func, init_params, n_accepted, max_iterations=1e4){
  # log_lh_func: a function returning the log-likelihood value. This is log( P(data | theta) )
  # priors_func: a function returning the joint prior log-likelihood of the parameters. This is  log( P(theta) )
  # proposal_func (or jumping function): a funciton returning a new parameter set, starting from another parameter set. P(theta' | theta)
  # init_params: a list containing the parameters. e.g list(beta=0.5, gamma=2.6)
  # n_accepted: number of accepted runs
  # max_iterations: maximum number of iterations allowed
  
  # preliminary check: max_iterations has to be >= n_accepted
  stopifnot(max_iterations>=n_accepted)
  
  # prepare storage for the results. columns will be named "log_lh", "accepted", "param1", "param2", ...
  results = data.frame(log_lh=double(), accepted=integer())
  for (param_name in names(init_params)){
    results[[param_name]] = double()
  }
  
  # initialise counters
  count_iterations = 0
  count_accepted = 0
  
  # calculate the likelihoods of the initial parameter set
  current_params = init_params  # current_params is the last accepted set of paramaters 
  current_log_lh = log_lh_func(init_params)
  current_log_prior = log_priors_func(init_params)
  
  while (count_accepted < n_accepted){
    
    if (count_iterations>=max_iterations){
      print("The maximum number of iterations has been reached. The simulation has been aborted.")  
      break
    }
    
    # Generate a new candidate parameter set
    proposed_params = proposal_func(current_params)
    
    # Evaluate the likelihood of the new parameter set
    proposed_log_lh = log_lh_func(proposed_params)
    proposed_log_prior = log_priors_func(proposed_params)
    
    # Acceptance or rejection?
    accepted = 0
    log_proba_of_acceptance = proposed_log_prior + proposed_log_lh - (current_log_prior + current_log_lh)  # we could have stored (current_log_prior + current_log_lh) in a variable
    proba_of_acceptance = exp(log_proba_of_acceptance)  # transform to actual proba
    if (proba_of_acceptance>=1){  # the proposed parameter set is "better" than the current one
      accepted = 1
    }else{
      accepted = rbinom(n=1, size=1, prob=proba_of_acceptance)
    }
    
    # storage
    new_row = list(log_lh=proposed_log_lh + proposed_log_prior, accepted=accepted)
    for (param_name in names(proposed_params)){
      new_row[[param_name]] = proposed_params[[param_name]]
    }
    results = rbind(results, new_row)
    
    # If the run is accepted, we update the relevant variables
    if (accepted==1){ 
      current_params = proposed_params
      current_log_prior = proposed_log_prior
      current_log_lh = proposed_log_lh
      count_accepted = count_accepted + 1
    }
    
    count_iterations = count_iterations + 1
    print(c(count_iterations, count_accepted))
  }
  
  # post-simulation calculation of acceptance ratio. 
  acceptance_ratio = count_accepted/count_iterations
  print(paste('Acceptance ratio:', acceptance_ratio, sep=' '))
  
  return(results)
  
}

my_log_lh_func <- function(params){
  # params is a list of parameters and associated values.
  # return the likelihood value associated with the parameter set.
  
  # run the TB model with the input parameters, output will be mortality over 9 age groups
  tb_results = ageStrTBModel_mortality(params$mult, params$beta1)
  
  # Our likelihood is obtained by multiplying Gaussian elements centered on the model estimate for each datapoint.
  # We assume that the standard deviation is 0.05. This means that 95% of the Gaussian dentisty sits within an 
  # interval of width 0.1 (2*sd).
  sd = 500
  overall_log_lh = 0
  for (i in 1:length(my_data)){  # for each value of the daw data
    model_output=tb_results[i]
    single_log_lh = dnorm(x=as.numeric(my_data[i]) , mean=model_output ,sd = sd, log = TRUE)
    overall_log_lh = overall_log_lh + single_log_lh
  }
  
  return(overall_log_lh)
}
#my_log_lh_func(params)

my_log_priors_func <- function(params){
  # params is a list of parameters and associated values
  # this funciton returns the joint prior distribution (actually log version)
  joint_log_prior = 0
  for (param_name in names(params)){
    if (param_name == "mult"){  # flat prior
      y = dunif(x=params[[param_name]], min = 35, max = 45,log = TRUE)
    }else if(param_name == "beta1"){  # beta distribution on interval [0,1]
      y = dunif(x=params[[param_name]], min = 24/365, max = 5*24/365, log=TRUE)
    }
    joint_log_prior = joint_log_prior + y
  }
  return(joint_log_prior)
}
#my_log_priors_func(params)

my_proposal_func <- function(params){
  # params is a list of parameters and associated values
  # this function returns another list of parameters
  # We use normal distributions to generate new parameter values
  
  # standard deviations
  sd=list(mult=2.5, beta1=0.02)
  
  new_params = list()
  for (param_name in names(params)){
    new_value = -1
    while (new_value < 0){  # we want the parameter value to be positive
      new_value = params[[param_name]] + rnorm(n=1, mean=0, sd=sd[[param_name]])
    }
    new_params[[param_name]] = new_value
  }
  
  return(new_params)
}

master_mcmc_runner <- function(n_accepted=5){
  init_params = list(mult = 40, beta1 = 3*24/365)
  M = metropolis_hastings(log_lh_func = my_log_lh_func, log_priors_func = my_log_priors_func, proposal_func = my_proposal_func, init_params = init_params, n_accepted = n_accepted)
  return(M)
}


