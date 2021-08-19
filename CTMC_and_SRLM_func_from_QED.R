CTMC.SRLM<-function(total.t, landscape, e.rate, c.rate, alpha, gamma, epsilon, 
                    self.rec, r.on, r.off){
  
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
  #randomly set initial habitat configuration of patches on/off
  #config<-rbinom(n.patches, 1, 0.5)
  #start from QED
  config<-
  #initial occupancy levels (all on patches fully occupied)
  new.ICs<-1*config
  #configs.data<-config
  
  #3) BEGIN THE SIMULATION
  t<-0
  tau.times<-0 #no initial time until transition
  time.limit<-total.t #you can remove the time.limit if you want
  while (sum(config)>0 & t<time.limit) { #while at least one patch remains on
    
    #3.1 SET THE RATES OF TRANSITION FOR EACH PATCH
    rate.vec<-rep(NA, n.patches) #construct a vector of transition rates 
    #depending on initial configuration
    rate.vec[config==0]<-r.on #if the patch is off set it to turn on at rate r.on
    rate.vec[config==1]<-r.off #if the patch is on set it to turn off at rate r.off
    
    #3.2 CALCULATE THE TIME TO THE NEXT TRANSITION
    #select the time until transition (exponential, with mean 1/rate.sum)
    rate.sum<-sum(rate.vec)
    tau<-rexp(1, 1/rate.sum)
    #tau<-1000 #for test purposes ONLY!!!
    
    #3.3 RECORD THE TIME OF EACH TRANSITION
    if (sum(tau.times) == 0){
      tau.times<-tau
    } else{
      tau.times<-c(tau.times, tau)
    }
    
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
    #if (tau<=1){ #if tau is less than one it will be rounded to 0 which will produce an error in ode
    #  tau<-1 #so if tau is less than 1, let it be rounded up to 1 
    #}
    
    #3.3 A4) SOLVE AND FORMAT OUTPUT INTO A DATA TABLE
    #must check order of magnitude of tau and set step length accordingly
    ###########################################TOO SLOW
    #if(tau<=10){
    #  mag<-ceiling(log10(tau))
    #} else { mag<-floor(log10(tau))}
    #step.length<-1*(10^(mag-1)) #set step length to be an order of magnitude less than tau
    ##########################################TOO SLOW
    ##########################################FASTER
    if(tau<1){ #if tau is less than one
      step.length<-tau/2 #let the step length for the SRLM be half of tau, to ensure the DE's are 
      #solved properly on this interval
    }else{step.length<-1} #otherwise just use a step length of 1
    ##########################################FASTER
    output<-ode(SRLM.ODE, parameters,times=seq(0,tau,step.length),y=new.ICs)
    
    #shift the time by the already accumulated time (t) and append to past values
    output<-data.frame(output)
    output$time<-output$time+t
    if (t==0){sim.data<-output
    configs.data<-c(t, config)}
    if (t!=0){sim.data<-rbind(sim.data, output)
    configs.data<-rbind(configs.data, c(t, config))}
    #pick what configuration to switch to
    
    #3.4 UPDATE THE TIME ELAPSED IN SIMULATION
    t<-t+tau
    
    #3.5 MAKE THE TRANSITION
    cum.rate<-cumsum(rate.vec/rate.sum) #calculate the cumulative rate of transition
    #define intervals between 0 and 1 by the proportion of time each patch should transition
    upper<-cum.rate
    lower<-c(0,cum.rate[-n.patches])
    unif.val<-runif(1, min=0, max=1) #pick a random number between 0 and 1
    if(config[unif.val<=upper & unif.val>lower]==1){ #if the patch lying on the selected interval 
      #within which the random value fell, is on
      config[unif.val<=upper & unif.val>lower]<-0 #turn it off
    } else {config[unif.val<=upper & unif.val>lower]<-1} #otherwise, turn it on
    print(config)
    
    #3.5 A) set the new initial conditions
    new.ICs<-as.numeric(as.character(output[length(output[,1]),-1])) #to be the end-values from the SRLM
    new.ICs[config==0]<-0 #but zero if a patch is now off
    
    print(t)
    print(config)
  }
  end.t<-t #time to absorption or simulation end
  return(list(sim.data=sim.data, configs.data=configs.data))
}

#Check:
landscape<-create.landscape(n.patches=10, landscape.type="linear", landscape.limit=100, 
                            patch.distribution="uniform", areas.distribution="uniform", areas.limit=1, 
                            clustering.iters=0)
print(ggplot(landscape, aes(x.coord,y.coord)) + theme_classic() + geom_point(aes(size = areas)))
n.patches<-length(landscape$patch.ID) #for future use
#
test<-CTMC.SRLM(total.t<-2000, landscape=landscape, e.rate<-0.1, c.rate<-0.2, alpha<-10, gamma<-0, 
                epsilon<-n.patches, self.rec<-1, r.on<-1/100, r.off<-1/100000)
test #~1 min for 15 timesteps -> 2.2 hrs :/
