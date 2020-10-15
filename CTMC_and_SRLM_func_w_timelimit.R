#This function simulates the combined CTMC+SRLM for the dynamics habitat metapopulation
#The CTMC portion is equivalent to the "CTMC_simulation" script, and the "SRLM_ODE" function
#is called on to run the ODEs in between the habitat configuration changes
CTMC.SRLM<-function(total.t, landscape, e.rate, c.rate, alpha, gamma, epsilon, 
                    self.rec, r1, r2){
  #setting up all our parameters
  ##########################################
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
  #################################################
  
  #################################################
  #initial habitat configuration of patches on/off
  config<-rbinom(n.patches, 1, 0.5)
  #initial occupancy levels
  new.ICs<-0.5*config
  #################################################
  
  t<-0
  time.limit<-total.t #you can remove the time.limit if you want
  while (sum(config)>0 & t<time.limit) { #while at least one patch remains on
    rate.vec<-rep(NA, n.patches) #construct a vector of transition rates 
    #depending on initial configuration
    for (i in 1:n.patches){ #for each patch
      if(config[i]==0){ #if the patch is off
        rate.vec[i]<-r1} # set it to turn from off to on at rate r1
      else{ #otherwise
        rate.vec[i]<-r2}} #set it to turn from on to off at rate r2
    
    #select the time until transition (exponential, with mean 1/run.sum)
    #here equivalently specified by an exponential distribution with rate=run.sum
    run.sum<-sum(rate.vec)
    tau<-rexp(1, rate=run.sum)
    print(tau)
    #tau<-1000 #for test purposes ONLY!!!
    #run the SRLM for the specified time tau, using initial conditions
    #new.ICs and relevant parameters
    #patches for which the "config=0" are left out
    #making sure "off" patches don't colonize or get colonized in the rhs of the ODE system
    if (t+tau<time.limit){ #if the time it will take to run to tau is within the time limit run 
      #the SRLM up to tau 
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
      
      #setting our initial values for the SRLM ODE
      ###########################################
      IC.names<-rep(paste0("p",1:n.patches))
      IC.values<-new.ICs
      p<-setNames(IC.values,IC.names)
      ###########################################
      if (tau<=1){ #if tau is less than one it will be rounded to 0 which will produce an error in ode
        tau<-1 #if tau is less than 1, let it be rounded up to 1 
      }
      output<-ode(SRLM.ODE, parameters,times=seq(0,round(tau, 0),1),y=new.ICs)
      #shift the time by the already accumulated time (t) and append to past values
      output<-data.frame(output)
      output$time<-output$time+t
      if (t==0){sim.data<-output}
      if (t!=0){sim.data<-rbind(sim.data, output)}
      #pick what configuration to switch to
      #use a uniform random draw to choose the next state based on rates of transition
      unif.val<-runif(1, min=0, max=1) #pick a random number between 0 and 1
      cum.rate<-cumsum(rate.vec/run.sum) #calculate the cumulative rate of transition
      #next state is determined by comparing the unif.value to the cumulative rate
      for (i in 1:n.patches){ #for each patch
        if (unif.val<=cum.rate[i]) { #if the unif.value <= the cumulative rate
          if (config[i]==0){ #if the patch was off
            config[i]<-1} #turn it on
          else { #otherwise
            config[i]<-0} #turn it off
        }
      }
      #set the new initial conditions
      new.ICs<-as.numeric(as.character(output[length(output$time),-1])) #to be the end-values from the SRLM
      new.ICs[config==0]<-0 #but zero if a patch is now off
      t<-t+tau #update the timestep to be where the simulation ended
      print(t)
      print(config)
    } else { #otherwise just do the same but up to the time limit
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
      
      #setting our initial values for the SRLM ODE
      ###########################################
      IC.names<-rep(paste0("p",1:n.patches))
      IC.values<-new.ICs
      p<-setNames(IC.values,IC.names)
      ###########################################
      output<-ode(SRLM.ODE, parameters,times=seq(0,time.limit,1),y=new.ICs)
      #shift the time by the already accumulated time (t) and append to past values
      output<-data.frame(output)
      output$time<-output$time+t
      if (t==0){sim.data<-output}
      if (t!=0){sim.data<-rbind(sim.data, output)}
      #pick what configuration to switch to
      #use a uniform random draw to choose the next state based on rates of transition
      unif.val<-runif(1, min=0, max=1) #pick a random number between 0 and 1
      cum.rate<-cumsum(rate.vec/run.sum) #calculate the cumulative rate of transition
      #next state is determined by comparing the unif.value to the cumulative rate
      for (i in 1:n.patches){ #for each patch
        if (unif.val<=cum.rate[i]) { #if the unif.value <= the cumulative rate
          if (config[i]==0){ #if the patch was off
            config[i]<-1} #turn it on
          else { #otherwise
            config[i]<-0} #turn it off
        }
      }
      #set the new initial conditions
      new.ICs<-as.numeric(as.character(output[length(output$time),-1])) #to be the end-values from the SRLM
      new.ICs[config==0]<-0 #but zero if a patch is now off
      t<-time.limit #update the timestep to be where the simulation ended
      print(t)
      print(config)
    }
    
  }
  end.t<-t #time to absorption or simulation end
  return(sim.data)
}
