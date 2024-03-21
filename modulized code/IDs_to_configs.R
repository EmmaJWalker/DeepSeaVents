##########################################################################################################
# DESCRIPTION:
##########################################################################################################

# This function convert identifying numbers for configurations into the configurations (binary vectors 
# indicating which patches are on and off) they represent. It uses an algorithm detailed in the 
# Configuration IDs document.

# INPUTS:
# x: an positive integer identifying a unique configuration of patches on and off for a given total 
# number of patches the landscape could have present (on)
# n.patches: how many patches total are in the system (i.e. the total number of patches the landscape 
# could have present)

# OUTPUTS:
# config: a binary vector indicating which patches are on or off as ordered by there patchIDs

# REQUIRED PACKAGES:
# none

# REQUIRED FUNCTIONS:
# none

####################################################################################
# FUNCTION CODE:
####################################################################################
IDsToConfigs<-function(config_ID, n.patches){
  x<-config_ID
  n.on.vec<-1:n.patches
  z<-factorial(n.patches)/(factorial(n.on.vec)*factorial(n.patches-n.on.vec)) 
  #number of different configurations possible with any given number of patches on
  z<-c(1, z) #add in the all off config
  #splitting these possibilities in half to utilize the symmetry of the system to rule out 
  #half the possibilities right from the get go
  if(even(n.patches)){#need to know if even or odd
    half.z<-z[(n.patches/2)+1]/2
    half.z<-c(z[1:(n.patches/2)], half.z)
  }else{half.z<-z[1:ceiling(n.patches/2)]}
  half.possible<-(2^n.patches)/2 #index at which the configs become an inverted mirror of 
  #the other configs
  if(x>half.possible){
    x<-(2^n.patches)-x
    invert<-T 
  }else{invert<-F}  
  if(x==0){#if x is zero it's the all on state
    config<-rep(0,n.patches)
  }else{
    if(x==1){ #if x is 1, the first patch is off
      config<-rep(0,n.patches)
      config[1]<-1
    }else{
      if(invert==T){x<-x+1} ##############**Because of the assymetry of excluding the all off state**
      config<-rep(0, n.patches)#vector storing the config (aka which patches are off)
      #let them all be on to start and we'll find which are off
      b<-which(x-cumsum(half.z)<=0)
      n.off<-b[1]-1 #tells us how many patches are off in the configuration 
      #now want to know, given how many patches are off,(therefore
      prev.x<-cumsum(half.z[1:n.off]) #gives the number of eliminated configs with fewer patches)
      prev.x<-prev.x[n.off]
      if(n.off==1){ #if only one patch is off
        config[x-1]<-1 #which patch is off is simply given by x-1
      }else{
        #now...given whichever may be the first off patch,
        #how many ways can we arrange the other n.off-1 patches in the remaining 
        #n-1, n-2, ... n-(n-n.off-1) locations?
        #n choose k, n locations, k are off
        #and likewise for the 2nd, 3rd etc... off patches until we have placed them all?
        locs<-n.patches-1:(n.patches-(n.off-1))
        prev.x2<-0
        r<-n.patches
        for (i in 1:(n.off-1)){#for each off patch
          z2<-factorial(locs)/(factorial(n.off-i)*factorial(locs-(n.off-i))) 
          #calculate how many ways can we arrange the other n.off-1 patches in the remaining 
          #n-1, n-2, ... n-(n-n.off-1) locations? (n choose k, n locations, k are off)
          b<-which(x-prev.x-prev.x2-cumsum(z2)<=0)
          off<-n.patches-locs[1]+b[1]-1 #tells us which is off patch i
          config[off]<-1 #set the patch identified as off to be off
          r<-n.patches-off #update the number of remaining number patches that could be on or off
          locs<-r-1:(r-(n.off-i-1)) #set locs to be the remaining locations off patches could be in
          z2<-c(0,z2)
          prev.x2<-prev.x2+sum(z2[1:b[1]])
        }
        #for the remaining last patch
        locs<-1:r
        b<-which(x-prev.x-prev.x2-locs<=0)
        off<-n.patches-r+b[1]#tells us which is off
        config[off]<-1 #set the patch identified as off to be off
      }
    }
  }
  if(invert==T){
    c<-as.logical(config)
    config[c]<-0
    config[!c]<-1
  }
  return(config)
}

####################################################################################
# QUICK CHECK: (uncomment and run to check)
####################################################################################
#rm(list=ls()) #to ensure a clean environment
## now run the FUNCTION CODE!!!
#IDsToConfigs(config_ID=31, n.patches=5)
#IDsToConfigs(config_ID=32512, n.patches=15)
#IDsToConfigs(config_ID=32767, n.patches=15)
####################################################################################
