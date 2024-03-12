CTMC.sim<-function(timesteps, n.patches, r.on, r.off, initial.config) {
  
  #1. SET UP THE PATCHES IN THE INITIAL CONFIGURATION
  config<-initial.config #setting up patches in the initial configuration
  tau.times<-0 #no initial time until transition
  
  #2. END IF WE BEGAN IN THE ABSORBING STATE
  if (sum(config)==0){ #if we are already in the absorbing state
    end.t<-0 #time to absorption is 0
    
    #3. OTHERWISE, PROVIDED WE'RE STILL NOT IN THE ABSORBING STATE...  
  } else { #otherwise
    for (t in 0:timesteps) { #comment out this line and end bracket if don't want to set a time limit
      #***actually can't do that in R
      while (sum(config)>0) { #while at least one patch remains on
        
        #3.1 SET THE RATES OF TRANSITION FOR EACH PATCH
        rate.vec<-rep(NA, n.patches) #construct a vector of transition rates 
        #depending on initial configuration
        rate.vec[config==0]<-r.on #if the patch is off set it to turn on at rate r.on
        rate.vec[config==1]<-r.off #if the patch is on set it to turn off at rate r.off
        
        #3.2 CALCULATE THE TIME TO THE NEXT TRANSITION
        #select the time until transition (exponential, with mean 1/rate.sum)
        rate.sum<-sum(rate.vec)
        tau<-rexp(1, rate.sum) #tau = time until next transition
        
        #3.3 RECORD THE TIME OF EACH TRANSITION
        #############################
        if (sum(tau.times) == 0){
          tau.times<-tau
        } else{
          tau.times<-c(tau.times, tau)
        }
        ##############################
        
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
        
        end.t<-t #time to absorption or simulation end
      } #end of while loop
    } #end of for loop
  } #end of else bracket
  return(list(end.t=end.t, tau.times=tau.times))
}

