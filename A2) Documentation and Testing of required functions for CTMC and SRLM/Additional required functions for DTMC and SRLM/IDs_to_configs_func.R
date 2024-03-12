#CONVERTING THESE CHOSEN CONFIGURATION IDs INTO THEIR CORRESPONDING CONFIGURATIONS:
IDs.to.configs<-function(x, n.patches){
  x<-25
  n.patches<-6
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
    invert<-T #since I want to use an sparse identity matrix to keep track of which patches are off
    #in the following algorithm to convert config indices to their configurations based on how patches are off
    #the configuration needs to be inverted to give it with 0's for off patches and 1's for on patches as
    #has been convention in the rest of my code
  }else{invert<-F} #however, if the config index is greater than half possible configs
  #there is no need to invert because these configurations are an inverted mirror of the first half anyways 
  config<-rep(0, n.patches)#vector storing the config (aka which patches are off)
  b<-which(x-cumsum(half.z)<=0)
  n.off<-b[1]-1 #tells us how many patches are off in the configuration 
  r.locs<-n.patches #how many different locations any remaining off patches could be in
  prev.x<-1
  for(i in 1:n.off){#now need to find the location of each off patch
    i<-1
    I<-Diagonal(r.locs) #identity matrix indicating where a single 
    #patch could be turned off in the remaining possible locations
    if(i==n.off){#if there's only 1 more off patch to place,
      #then it must just be in one of the remaining locations so
      r<-which(x-prev.x-1:r.locs==0)
      config[n.patches-r.locs+r]<-1
    }else{
      if(i==1){prev.x<-prev.x+n.patches}
      ##r.locs<-r.locs-1 #since an off patch must be in one of those locations, 
      #the remaining locations any of the other patches could be in is decreased by one 
      ##n<-r.locs:1 #locations the next off patch could be in
      ##**************************
      n.locs<-(n.patches-1):1 #how many different locations the 1st off patch could be in
      n<-rep(NA, length(n.locs))
      for (j in 1:length(n.locs)){# and within each of those locations:
        n.locs<-c(n.locs, n.locs[-1])
        n[j]<-sum(n.locs) #the number of locations the next off patch could be in
      }
      ##i.locs<-(n.patches-1):n.off #how many different locations the 1st off patch could be in 
      #given how many off patches there are
      ##n.locs<-(n.patches-n.off):1 #how many different locations the last patch could be in given 
      #how many off patches their are
      ##n<-rep(NA, length(n.locs))
      ##  for (j in 1:length(n.locs)){# and within each of those locations:
      ##    n.locs<-c(n.locs, n.locs[-1])
      ##    n[j]<-sum(n.locs) #the number of configurations
      ##  }
      ###*************************
      r<-which(x-prev.x-cumsum(n)<=0) #*gives location of the purple off patch, want pink off patch*
      #which is the first location and off patch must be in for this number of patches to be off?
      c<-r[1] #c is where a patch must be off
      config[(r.locs+1):1]<-rev(I[,c]) #gives the config updated with that patch off
      r.locs<-n.patches-c #update the total number of locations the remaining off patches could be in
      #given what we now know
      prev.x<-prev.x+sum(n[-r]) #update the index configs we know x is not
    }
  }
  if(invert==F){
    c<-as.logical(config)
    config[c]<-0
    config[!c]<-1
  }
  return(config)
}

#CHECK:
IDs.to.configs(x=32512, n.patches=15)
