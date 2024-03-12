#Finds the equilibrium metapopulation size for each configuration
###############################################################################
QED.metrics<-function(s.size, landscape, e.rate, c.rate, gamma, epsilon, self.rec, alpha, r.on, r.off){
  
  n.patches<-length(landscape$patch.ID)
  
  #get the QED of the reduced CTMC
  QED.output<-quasi.eq(n.patches, r.on, r.off)
  QED.dist<-QED.output[[1]] #QED distribution
  # use this to get the probabilities of just the subsampled configurations
  #converting the QED of the reduced CTMC to the QED of the full CTMC:
  n.on.vec<-1:n.patches
  z<-factorial(n.patches)/(factorial(n.on.vec)*factorial(n.patches-n.on.vec)) #calculating N_choose_k
  #aka number of different configurations possible with that number of patches on
  QED.p<-QED.dist/z #converting the probabilities n patches are on at QED to 
  #probabilities of being in each possible configuration with n patches on
  
  #calculate the extinction rate for each patch
  #extinction.rates<-e.rate/areas #if want this to be patch dependent
  #calculate delta
  delta<-e.rate/c.rate
  
  #USE ONLY A SUBSAMPLE OF THE POSSIBLE CONFIGURATIONS:
  z.tracker<-z
  configs.sample<-rep(NA, s.size)
  c.tracker<-data.frame(matrix(rep(NA, s.size*2), ncol=2, nrow=s.size))
  for (i in 1:s.size){
    n.on.sample <- sample(n.on.vec[z.tracker>0], 1, prob=QED.dist[z.tracker>0]/sum(QED.dist[z.tracker>0])) #sub sample how many patches could 
    #possibly be on with the probability that number of patches in on at QED
    prior.configs<-cumsum(z[1:(n.on.sample)])
    p<-(1:z[n.on.sample]) #index of possible configurations with that number of patches on
      q<-filter(c.tracker, X1==n.on.sample)
      if(any(q$X2)){
        p<-p[!p %in% q$X2]
        #p<-p[-which(p==q$X2)] #taking out any combinations of n and k we have already used
        #from the possibilities to sample
      }
      if(length(p)==1){#if only one possible config with that many patches on remains,
        k<-p #then choose that config
      }else{ #otherwise,
        k<-sample(p, 1)  # within the
        #number of configurations with that number of patches on,
        #randomly choose a configuration with that number of patches on
      }
      configs.sample[i] <- prior.configs[n.on.sample-1] + k + 1 #get the index of that unique config out 
      #of all the possible configs
    z.tracker[n.on.sample]<-z.tracker[n.on.sample]-1 #now one less configuration with that many patches on is possible
    c.tracker$X1[i]<-n.on.sample #keep track of combinations of n's and
    c.tracker$X2[i]<-k #k's chosen in a data frame
  }
  
  #CONVERTING THESE CHOSEN CONFIGURATION IDs INTO THEIR CORRESPONDING CONFIGURATIONS:
  IDs.to.configs<-function(x, n.patches){
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
        r.locs<-r.locs-1 #since an off patch must be in one of those locations, 
        #the remaining locations any of the other patches could be in is decreased by one 
        n<-r.locs:1 #locations the other off patch could be in
        r<-which(x-prev.x-cumsum(n)<=0)
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
  }
  
    
  #all.configs<-matrix(rep(NA, (n.patches*s.size)), nrow=s.size, ncol=n.patches)
  #for (i in 1:length(configs.sample)){
  #  all.configs[i,]<-BigintToBin(configs.sample[i], digits=(n.patches))
  #}
  #reordering all.configs by the number of patches on
  #all.configs<-data.frame(all.configs)
  #all.configs$n.on<-rowSums(all.configs)
  #all.configs<-all.configs[order(all.configs$n.on)]
  #all.configs<-all.configs[,-(n.patches+1)]
  
  #########################################
  #IN PARALLEL:
  library("doParallel") 
  cl <- makeCluster(detectCores(), type='PSOCK')
  registerDoParallel(cl)
  
  each.config.calcs<-function(i, n.patches, s.size, landscape, e.rate, c.rate, gamma, epsilon, self.rec, alpha){
    
    i=i
    n.patches=n.patches
    s.size=s.size
    landscape=landscape
    e.rate=e.rate
    c.rate=c.rate
    gamma=gamma
    epsilon=epsilon
    self.rec=self.rec
    alpha=alpha
    
    
    #LOADING CODE SPECIFIC PACKAGES:
    library("R.utils") #for int to binary
    library("pracma") #for matrix math
    library("R.utils") #for integer to binary conversion
    library("dplyr") #for various stuff
    library("gdata") #to get upper triangle of a matrix
    library("lava") #to get the reverse diagonal of a matrix
    library("Matrix") #for sparse matrices
    library("reticulate") #for something...
    library("gtools") #for na.replace
    #LOADING REQUIRED FUNCTIONS
    setwd("C:/Users/Administrator/Desktop/DEC 2022")
    source("dispersal_kernel_func.r") #disp.kernel()
    source("lambda_M_func.r") #get.lambda.M()
    source("Pstar_func.r") #pstar.function()
    
    # foreach (1:45){
    v<-floor((s.size)/48) #define the set of configs to be run on core i
    #if (i==46){}else{ #core number 46 doesn't work it seems
    if (i==48 & ((s.size)/48)>floor((s.size)/48)){ #if it's core number 48
      v<-v+(s.size)-(floor((s.size)/48)*48) #give the remaining configs to it 
      sub.i.configs<-c(((s.size)-v+1):(s.size))
    }else{
      sub.i.configs<-c((((i-1)*v)+1):(i*v))
    }
    i.metrics<-data.frame(matrix(data=rep(NA, v*(5+n.patches)), nrow=v, ncol=(5+n.patches)))
    for (k in sub.i.configs) { #for each config core i...
      config<-all.configs[k,] #selects a config k...
        #calculates the persistence capacity in that configuration
        landscape.config<-landscape[config!=0,]
        lambda.M<-get.lambda.M(landscape.config, alpha, gamma, epsilon, self.rec, e.rate, c.rate)
        lambda.M.delta<-get.lambda.M(landscape.config, alpha, gamma, epsilon, self.rec, e.rate, c.rate)/delta
        if (lambda.M.delta < delta){
          pstar.configs<-rep(0, n.patches)
        } #if below the persistence threshold, the metapop is expected extinct at equilibrium 
        else { #otherwise
          #calculates the equilibrium patch occupancy of each patch in config k
          pstar.configs<-rep(NA, n.patches)
          pstar.configs[config!=0]<-pstar.function(landscape.config, alpha, e.rate, c.rate, self.rec, gamma, epsilon, iterations=1000)
          pstar.configs<-na.replace(pstar.configs, 0)
          
        }
        metapop.size<-sum(pstar.configs) #calculate total metapop size expected at eq
        n.on<-sum(config)
        i.metrics[k,]<-data.frame(lambda.M, lambda.M.delta, metapop.size, n.on, t(pstar.configs))
    }
    i.metrics<-na.omit(i.metrics)
    return(i.metrics)
  } 
  
  config.calcs<-foreach(i=1:48, .combine=rbind) %dopar% each.config.calcs(i=(i), n.patches=n.patches, s.size=s.size,
                                                                     landscape=landscape, e.rate=e.rate, c.rate=c.rate, 
                                                                     gamma=gamma, epsilon=epsilon, self.rec=self.rec, 
                                                                     alpha=alpha)
  # turn parallel processing off and run sequentially again:
  registerDoSEQ()
  #stopCluster()
  #stopImplicitCluster()
  
  #END PARALLEL:
  ######################################
  lambda.M<-config.calcs[,1]
  lambda.M.delta<-config.calcs[,2]
  metapop.size<-config.calcs[,3]
  n.on<-config.calcs[,4]
  pstar.configs<-config.calcs[,(5:(4+n.patches))]
  
  
  n.on<-rowSums(all.configs) #get how many patches on in each of our configs
  QED<-QED.p[n.on] #extracting the probability of each of our configurations given by the number of patches on in it
  
  
  #get the time-averaged metapopulation sizes @QED
  ar.pstar.QED<-dot(QED,metapop.size) #this is the expected metapop size in each configuration x the porportion of time spent in it at QED and summed (therefore the arithmentic mean)
  #should actually be the geometric mean
  geom.pstar.QED<-exp(dot(QED,log(metapop.size)))
  lm.QED.delta<-exp(dot(QED,log(lambda.M.delta)))
  lm.QED<-exp(dot(QED,log(lambda.M)))
  
  #get expected number of patches in habitat @QED
  n.on.in.each.config<-rowSums(all.configs)
  #geometric mean habitat expected across configurations and time spent in them:
  exp.n.QED<-exp(dot(QED,log(n.on.in.each.config)))
  
  print(QED)
  print(n.on.in.each.config)
  print(lambda.M)
  print(lambda.M.delta)
  print(lm.QED)
  print(lm.QED.delta)
  
  return(list(lm.QED.delta, lm.QED, geom.pstar.QED, ar.pstar.QED, exp.n.QED, QED, all.configs, pstar.configs))
}