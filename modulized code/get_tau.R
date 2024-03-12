##########################################################################################################
# DESCRIPTION:
##########################################################################################################

# This function draws a value for when the system will leave a given habitat configuration (i.e. how long
# until the system will on average jump to a new configuration from the current configuration)

# INPUTS:
# config: a binary vector indicating which patches are present or absent in a landscape
# r.on: the rate of patch gain
# r.off: the rate of patch loss

# OUTPUTS:
# tau: a numeric value of how much time will pass on average before the habitat configuration will change

# REQUIRED PACKAGES:
# none

# REQUIRED FUNCTIONS:
# none
##########################################################################################################
get_tau<-function(config, r.on, r.off){
  n.patches<-length(config)
  rate.vec<-rep(NA, n.patches) #construct a vector to hold transition rates
  rate.vec[config==0]<-r.on #if the patch is off set it to turn on at rate r.on
  rate.vec[config==1]<-r.off #if the patch is on set it to turn off at rate r.off
  
  #select the time until transition (exponential, with mean 1/rate.sum)
  rate.sum<-sum(rate.vec)
  tau<-rexp(1, rate.sum)
  return(tau)
}







