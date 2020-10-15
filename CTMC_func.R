#This script simulates the CTMC portion of the model, 
#either for a certain amount of time (t.duration),
#or until it reaches the absorbing state.
#Each configuration is represented by a binary vector, 
#with the first element (starting from the left) representing whether site 1 is on/off
#the seccond whether site 2 is on/off, ...etc.

#NOTE: this function is actually very fast provided the expected time to absorption is low!

##################################################################################################
#SIMULATION using Gillespie algorithm to simulate Markov chain
CTMC.sim<-function(timesteps, n.patches, r.on, r.off, initial.config) {
  config<-initial.config
  tau.times<-0 #no initial time until transition
  
  if (sum(config)==0){ #if we are already in the absorbing state
    end.t<-0 #time to absorption is 0
  } else { #otherwise
    for (t in 0:timesteps) { #comment out this line and end bracket if don't want to set a time limit
      #***actually can't do that in R
      while (sum(config)>0) { #while at least one patch remains on
        rate.vec<-rep(NA, n.patches) #construct a vector of transition rates 
        #depending on initial configuration
        for (i in 1:n.patches){ #for each patch
          if(config[i]==0){ #if the patch is off
            rate.vec[i]<-r.on} # set it to turn from off to on at rate r.on
          else{ #otherwise
            rate.vec[i]<-r.off}} #set it to turn from on to off at rate r.off
        
        #select the time until transition (exponential, with mean 1/rate.sum)
        rate.sum<-sum(rate.vec)
        tau<-rexp(1, 1/rate.sum) #tau = time until next transition #********************double check if this should be run.sum or 1/run.sum
        
        #############################
        if (sum(tau.times) == 0){
          tau.times<-tau
        } else{
          tau.times<-c(tau.times, tau)
        }
        ##############################
        t<-t+tau
        
        unif.val<-runif(1, min=0, max=1) #pick a random number between 0 and 1
            if ((unif.val>=r.off/(rate.sum))){ #if the random number is greater than or equal to 
              #the porportion of times a patch should turn off
              if (any(config==0)){ #if any patches are off
                off.patches<-which(config==0) #find which patches are off
                x<-sample(off.patches, 1)#choose a random patch that's off 
                config[x]<-1# and turn it on 
              }
              }else{#otherwise,
                if (any(config==1)){ #if any or all patches are on
                  on.patches<-which(config==1) #find which patches are on
                  x<-sample(on.patches, 1) #choose a random patch
                  config[x]<-0 #and turn it off
                }
            }
        print(config)
        #cum.rate<-cumsum(rate.vec/rate.sum) #calculate the cumulative rate of transition
    
        #next state is determined by comparing the unif.value to the cumulative rate:
        #for (i in 1:n.patches){ #for each patch
        #  if (unif.val<=cum.rate[i]) { #if the unif.value <= the cumulative rate
        #   if (config[i]==0){ #if the patch was off
        #      config[i]<-1} #turn it on
        #    else { #otherwise
        #      config[i]<-0} #turn it off
        #  }}
        end.t<-t #time to absorption or simulation end
      } #end of while loop
    } #end of for loop
  } #end of else bracket
  return(list(end.t, tau.times))
}
#TEST

#INITIALIZING PARAMETERS
#timesteps<-1000000 #total number of timesteps to run the model (optional, 
#if unspecified will run until absorption)
#n.patches<-2 #number of patches in the landscape
#r.on<-100 #rate at which patches turn from off to on
#r.off<-1 #rate at which patches turn from on to off
#initial.config<-c(1,1)#rbinom(n.patches, 1, 0.5) #initial habitat configuration of patches on and off
#drawn from the binomial distribution with 0.5 probability of success getting 1 (or 0)

#out<-CTMC.sim(1000, 2, r.on=100, r.off=1, initial.config = c(1,1))
#out
#warnings are fine because they are indeed telling me that the part of the function they correspond 
#to are doing what I want
#warnings()
#out[[1]]
#out[[2]]

