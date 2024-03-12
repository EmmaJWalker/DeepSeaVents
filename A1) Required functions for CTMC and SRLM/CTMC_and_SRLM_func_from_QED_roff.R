CTMC.SRLM<-function(total.t, landscape, e.rate, c.rate, alpha, gamma, epsilon, 
                    self.rec, r.on, r.off, initial.p,
                    initial.config){
  
  #1) SETTING UP ALL OUR PARAMETERS FOR THE SRLM's ODEs
  #########################################################################################
  n.patches<-length(landscape$patch.ID)
  x.coord<-landscape$x.coord
  y.coord<-landscape$y.coord
  areas<-landscape$areas
  #calculate the extinction rate for each patch
  extinction.rates<-e.rate/areas 
  #calculate delta
  delta<-(e.rate+r.off)/c.rate
  # construct the distance matrix
  dist.mat<-matrix(rep(0,n.patches*n.patches),n.patches,n.patches)
  for (i in 1:n.patches){
    for (j in 1:n.patches){
      dist.mat[i,j]<-sqrt((x.coord[i]-x.coord[j])^2+(y.coord[i]-y.coord[j])^2)
    }
  }
  #CHECK: dist.mat
  #Getting the f(dij) values to plug into the rhs
  f.vals<-disp.kernel(dist.mat,alpha,gamma,epsilon,self.rec,n.patches)
  #CHECK: f.vals
  f.vals<-as.vector(f.vals)
  f.names<-matrix(rep(NA,n.patches*n.patches),n.patches,n.patches)
  for (i in 1:n.patches){
    for (j in 1:n.patches){
      f.names[i,j]<-paste0("f.vals",i,j)
    }
  }
  f.names<-as.vector(f.names)
  #########################################################################################
  
  #2) SET THE INITIAL CONDITIONS
  #########################################################################################
  config<-initial.config
  new.ICs<-initial.p
  #########################################################################################
  
  #3) BEGIN THE SIMULATION
  #########################################################################################
  t<-0
  tau<-0
  tau.times<-0 #no initial time until transition
  ###############################################
  ##IF WANT A TRANSITION LIMIT USE THE FOLLOWING:
  ##trans.limit<-total.t
  #while ((sum(config)>0) & (length(tau.times<=trans.limit))) {
  ###############################################
  ###############################################
  #IF WANT A TIME LIMIT USE THE FOLLOWING:
  time.limit<-total.t #you can remove the time.limit if you want
  while ((sum(config)>0) & (t<=time.limit)){#(t+tau<=time.limit)) { #while at least one patch remains on 
    #and we are within the time limit
  ###############################################
    
    #3.1 SET THE RATES OF TRANSITION FOR EACH PATCH
    ###############################################
    rate.vec<-rep(NA, n.patches) #construct a vector of transition rates 
    #depending on initial configuration
    rate.vec[config==0]<-r.on #if the patch is off set it to turn on at rate r.on
    rate.vec[config==1]<-r.off #if the patch is on set it to turn off at rate r.off
    ###############################################
    
    #3.2 CALCULATE THE TIME TO THE NEXT TRANSITION
    #select the time until transition (exponential, with mean 1/rate.sum)
    #####################################################################
    rate.sum<-sum(rate.vec)
    tau<-rexp(1, rate.sum)
    #tau<-1000 #for test purposes ONLY!!!
    #####################################################################
    
    #3.3 RECORD THE TIME OF EACH TRANSITION
    #######################################
    if (sum(tau.times) == 0){
      tau.times<-tau} else{
        tau.times<-c(tau.times, tau)}
    #######################################
    
    #3.3 A) SOLVE THE SRLM BASED ON THE INITIAL CONDITIONS BETWEEN NOW AND TAU
    #run the SRLM for the specified time tau, using initial conditions
    #new.ICs and relevant parameters
    #patches for which the "config=0" are left out
    #making sure "off" patches don't colonize or get colonized in the rhs of the ODE system
    #########################################################################################
    
    #3.3 A1) CREATE A NAMED LIST OF PARAMETERS FOR THE SRLM ODEs
    #creating a list of the parameters we calculated to be passed to our ODE function
    ############################################
    parameter.names<-c("n.patches", "c.rate", "self.rec",
                       rep(paste0("extinction.rates",1:n.patches)),
                       rep(paste0("x.coord",1:n.patches)),
                       rep(paste0("y.coord",1:n.patches)),
                       rep(paste0("areas",1:n.patches)),
                       rep(paste0("config",1:n.patches))
                       ,f.names)
    parameter.values<-c(n.patches,c.rate,self.rec,extinction.rates,x.coord,y.coord,areas,config
                        ,f.vals)
    parameters<-list(c(setNames(parameter.values, parameter.names)))
    ###########################################
    
    #3.3 A2) Setting our initial values for the SRLM ODE
    ###########################################
    IC.names<-rep(paste0("p",1:n.patches))
    IC.values<-new.ICs
    p<-setNames(IC.values,IC.names)
    ###########################################
    
    #3.3 A3) set time interval and step size to solve over**
    ##########################################FASTER
    #if(tau<1){ #if tau is less than one
    #  step.length<-tau/2 #let the step length for the SRLM be half of tau, 
    #}else{step.length<-1} #otherwise just use a step length of 1
    endt<-tau 
    if (t+tau>time.limit){ #if the time till the next transition will exceed our simulation time limit
      endt<-time.limit-t #solve the model only up until that time limit
    }
    #if(endt<(1*(10^(-10)))){
    #  end.t<-(1*(10^(-10)))#cap the smallest tau can be to prevent this breaking the ODE solver
    #}
    if(endt<1){ #if tau is less than one
      step.length<-endt/2 #let the step length for the SRLM be half of tau, 
    }else{step.length<-1} #otherwise just use a step length of 1
    
    ##########################################FASTER
    #3.3 A4) Solve the SRLM ODEs over the time interval between transitions
    #######################################################################################################
    output<-ode(y=new.ICs, parms=parameters, times=seq(0, endt, step.length), func=SRLM.ODE,
                rootfun=at.equilibrium)
    #######################################################################################################
    
    
    #3.3 A5) Check if equilibrium was reached before the end of the time interval between 
    #transitions was reached and append the equilibrium occupancy as the occupancy at the end 
    #of that time interval
    #########################################
    #tail(output, n = 2) #look at the end of the simulation
    #calculate the final change in patch occupancy
    x1<-output[length(output[,1]),-1]-output[length(output[,1])-1,-1] 
    if(any(is.na(x1))){ #if an NA value occured in solving the ODE
      t<-time.limit #end the simulation immediately!
    } else if (Norm(x1)<=(1e-4)){ #if the last change in occupancy was less than 1e-4
      #the SRLM reached equilibrium before endt, ending the solver early to save on computation
      #so we shall append the output with that same equilibrium occupancy at endt
      output<-rbind(output, c(endt, output[length(output[,1]),-1]))
    }
    ######################################
    
    #3.3 A6) Shift the time by the already accumulated time (t) and append to past values
    #############################################################################
    output<-data.frame(output)
    output$time<-output$time+t
    if (t==0){sim.data<-output
    configs.data<-c(t, config)}
    if (t!=0){sim.data<-rbind(sim.data, output)
    configs.data<-rbind(configs.data, c(t, config))}
    ############################################################################
    
    #3.4 UPDATE THE TIME ELAPSED IN SIMULATION
    ##########################################
    t<-t+tau
    ##########################################
    
    #3.5 MAKE THE TRANSITION
    ######################################################################################################
    cum.rate<-cumsum(rate.vec/rate.sum) #calculate the cumulative rate of transition
    #define intervals between 0 and 1 by the proportion of time each patch should transition
    upper<-cum.rate
    lower<-c(0,cum.rate[-n.patches])
    unif.val<-runif(1, min=0, max=1) #pick a random number between 0 and 1
    if(config[unif.val<=upper & unif.val>lower]==1){ #if the patch lying on the selected interval 
      #within which the random value fell, is on
      config[unif.val<=upper & unif.val>lower]<-0 #turn it off
    } else {config[unif.val<=upper & unif.val>lower]<-1} #otherwise, turn it on
    #print(config)
    ######################################################################################################
    
    #3.5 A) set the new initial conditions
    ######################################################################################################
    new.ICs<-as.numeric(as.character(output[length(output[,1]),-1])) #to be the end-values from the SRLM
    new.ICs[config==0]<-0 #but zero if a patch is now off
    ######################################################################################################
    
    #for testing only
    print(t)
    print(config)
    
  }
  end.t<-t #time to absorption or simulation end
  #sim.data<-rbind(rep(end.t, ncol(sim.data))) #to recording end.t to know if hit timelimit or absorption in the sim.data
  return(list(sim.data=sim.data, configs.data=configs.data))
}
