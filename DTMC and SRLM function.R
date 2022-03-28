DTMC.SRLM<-function(total.t, landscape, e.rate, c.rate, alpha, gamma, epsilon, 
                    self.rec, r.on, r.off, P.mat){
  
  n<-n.patches
  
  #1) SETTING UP ALL OUR PARAMETERS FOR THE SRLM's ODEs
  n.patches<-length(landscape$patch.ID)
  x.coord<-landscape$x.coord
  y.coord<-landscape$y.coord
  areas<-landscape$areas
  #calculate the extinction rate for each patch
  extinction.rates<-e.rate/areas 
  #calculate delta
  delta<-e.rate/c.rate
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
  
  #2) SET THE INITIAL CONDITIONS
  #initial habitat configuration of patches on/off
  config<-rbinom(n.patches, 1, 0.5)
  #initial occupancy levels
  new.ICs<-0.5*config
  #configs.data<-config
  
  ########################################################################
  ########################################################################
  #STUFF SO WE DON"T HAVE TO FLIP PATCHES ON AND OFF ALL THE TIME AND CAN 
  #JUST PICK CONFIGURATIONS IN OUR SIMULATIONS
  #########################################################################
  #1. START IN THE INITIAL CONFIGURATION
  #since we'll just keep track of the state associated with each configuration
  #rather than updating the configuration over and over again in the simulation
  #let's convert the configuration in to the row/column index that corresponds to that state
  #Because of how P.mat is ordered, with the first half of the matrix rows 
  #corresponding to the first patch in the config being off, the first half 
  #of those corresponding to the 2nd patch being off as well, and so on...
  #the following algorithm should provide us with which row corresponds to 
  #our current configuration
  for (i in 1:n){
    if (i==1) {
      if (config[n+1-i]==0){x<-2^0} else {x<-2^1}
    } else {
      if (i>1) {
        if (config[n+1-i]==0){x<-x+0} else {x<-x+2^(i-1)}
      }
    }
  }
  #1. create a dataframe with all possible configurations
  #so we can quickly convert state x to a config throughout the simulation
  #and save on storing configs for every timestep
  #ordered according to the sequence they occur as rows in our transition matrix
  #NOTE: updated to store these as sparse matrices to save memory :)
  all.configs <- Matrix(0, nrow=2^n, ncol=n, sparse=TRUE) #sparse matrix version
  #all.configs<-matrix(rep(NA, (2^n)*n), nrow=2^n, ncol=n) #no sparse matrix version
  for (i in 1:n){
    #if (i==1){BB<-c(0,1) #no sparse vector version
    if (i==1){BB<-as(c(0,1), "sparseVector") #sparse vector version
    }else{BB<-c(BB[1:(2^(i-2))],BB,BB[((2^(i-2))+1):(2^(i-1))])}
    all.configs[,(n-i+1)]<-t(rep(BB,2^(n-i)))}
  #######################################################################
  ########################################################################
  
  #3) BEGIN THE SIMULATION
  t<-0
  time.limit<-total.t #you can remove the time.limit if you want
  while (x>1 & t<time.limit) { #while at least one patch remains on
    
    ################################################################
    #Find the time to the next transition the good old fashioned way:
    #i.e. simulate the DTMC until a transition occurs
    tau<-0
    x.new<-sample((1:2^n), size=1, prob=P.mat[x,])
    while (x==x.new){ #while we are stuck in one state
      tau<-tau+1 #and keeping track of how much time is spent in the same state
      x.new<-sample((1:2^n), size=1, prob=P.mat[x,]) #simulate the transition of the DTMC
    } #when we finally make a transition to a new state, BEFORE WE MAKE THAT TRANSITION,
    #let's simulate the SRLM over the time we spent in that state...
    
    #set the initial config for the SRLM to the previous state
    config<-all.configs[x,]
    
    #3.3 A) SOLVE THE SRLM BASED ON THE INITIAL CONDITIONS BETWEEN NOW AND TAU
    #run the SRLM for the specified time tau, using initial conditions
    #new.ICs and relevant parameters
    #patches for which the "config=0" are left out
    #making sure "off" patches don't colonize or get colonized in the rhs of the ODE system
    
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
    
    #3.3 A3) Make sure tau is not <1 or this will create an error in the ODE solver
    if (tau<=1){ #if tau is less than one it will be rounded to 0 which will produce an error in ode
      tau<-1 #so if tau is less than 1, let it be rounded up to 1 
    }
    
    #3.3 A4) SOLVE AND FORMAT OUTPUT INTO A DATA TABLE
    output<-ode(SRLM.ODE, parameters,times=seq(0,round(tau, 0),1),y=new.ICs)
    #shift the time by the already accumulated time (t) and append to past values
    output<-data.frame(output)
    output$time<-output$time+t
    if (t==0){sim.data<-output
    configs.data<-c(t, config)}
    if (t!=0){sim.data<-rbind(sim.data, output)
    configs.data<-rbind(configs.data, c(t, config))}
    
    #NOW we can switch to the next configuration
    x<-x.new
    
    #3.2 UPDATE THE TIME ELAPSED IN SIMULATION
    t<-t+tau
    
  } #end of while loop
  
  end.t<-t #time to absorption or simulation end
  return(list(sim.data=sim.data, configs.data=configs.data))
}