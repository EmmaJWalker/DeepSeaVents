##########################################################################################################
# DESCRIPTION:
##########################################################################################################

# This function makes a transition from one habitat configuration to another based on the probability of 
# transitioning to that configuration from the current configuration given by the generator matrix

# INPUTS:
# config: a binary vector indicating which patches are present or absent in a landscape
# r.on: the rate of patch gain
# r.off: the rate of patch loss

# OUTPUTS:
# config: a binary vector indicating which patches are present or absent in the new habitat configuration 
# in a landscape

# REQUIRED PACKAGES:
# none

# REQUIRED FUNCTIONS:
# none
##########################################################################################################
transition<-function(config, r.on, r.off){
  n.patches<-length(config)
  rate.vec<-rep(NA, n.patches) #construct a vector to hold transition rates
  rate.vec[config==0]<-r.on #if the patch is off set it to turn on at rate r.on
  rate.vec[config==1]<-r.off #if the patch is on set it to turn off at rate r.off
  rate.sum<-sum(rate.vec)
  
  cum.rate<-cumsum(rate.vec/rate.sum) #calculate the cumulative rate of transition
  #define intervals between 0 and 1 by the proportion of time each patch should transition
  upper<-cum.rate
  lower<-c(0,cum.rate[-n.patches])
  unif.val<-runif(1, min=0, max=1) #pick a random number between 0 and 1
  if(config[unif.val<=upper & unif.val>lower]==1){ #if the patch lying on the selected interval 
    #within which the random value fell, is on
    config[unif.val<=upper & unif.val>lower]<-0 #turn it off
  } else {config[unif.val<=upper & unif.val>lower]<-1} #otherwise, turn it on
  
  return(config)
}
##########################################################################################################
