#3.3 A) SOLVE THE SRLM BASED ON THE INITIAL CONDITIONS BETWEEN NOW AND TAU
#run the SRLM for the specified time tau, using initial conditions
#new.ICs and relevant parameters
#patches for which the "config=0" are left out
#making sure "off" patches don't colonize or get colonized in the rhs of the ODE system

setNrun_SRLMODE<-function(parameter.values, f.names, ICs, t, tau, sim.data, configs.data){
  
  n.patches<-parameter.values[1]
  config<-parameter.values[(4+n.patches*4):(3+n.patches*5)]
  
  #3.3 A1) CREATE A NAMED LIST OF PARAMETERS FOR THE SRLM ODEs
  #creating a list of the parameters we calculated to be passed to our ODE function
  ############################################
  parameters<-name_SRLMODEparams(parameter.values, f.names)
  ###########################################
  
  #3.3 A2) Setting our initial values for the SRLM ODE
  ###########################################
  IC.names<-rep(paste0("p",1:n.patches))
  IC.values<-ICs
  p<-setNames(IC.values,IC.names)
  ###########################################
  
  #3.3 A3) Make sure tau is not too small or this will create an error in the ODE solver
  ###########################################
  if(tau<1){ #if tau is less than one
    step.length<-tau/2 #let the step length for the SRLM be half of tau, to ensure the DE's are 
    #solved properly on this interval
  }else{step.length<-1} #otherwise just use a step length of 1
  ###########################################
  
  #3.3 A4) Solve the system of ODE's describing the SRLM for a time interval between now and tau time
  # from now
  ###########################################
  output<-ode(SRLM.ODE, parameters,times=seq(0,tau,step.length),y=ICs)
  ###########################################
  
  #3.3 A5) append the metapopulation dynamics over the time spent in the current habitat configuration
  # to the simulated metapopulation data and what configuration the habitat was in during this time to 
  # the habitat configuration data
  
  # record how much time elapsed and the habitat configuration it was spent in 
  ###########################################
  output<-data.frame(output)
  output$time<-output$time+t
  sim.data<-rbind(sim.data, output)
  configs.data<-rbind(configs.data, c(t, config))
  ###########################################

  return(list(sim.data, configs.data))
}
  

  