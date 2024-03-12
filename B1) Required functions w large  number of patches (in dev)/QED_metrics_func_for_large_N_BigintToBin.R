#Finds the equilibrium metapopulation size for each configuration
###############################################################################
QED.metrics<-function(s.size, landscape, e.rate, c.rate, gamma, epsilon, self.rec, alpha, r.on, r.off){
  
  n.patches<-length(landscape$patch.ID)
  
  #calculate the extinction rate for each patch
  #extinction.rates<-e.rate/areas #if want this to be patch dependent
  #calculate delta
  delta<-e.rate/c.rate
  
  #USE ONLY A SUBSAMPLE OF THE POSSIBLE CONFIGURATIONS:
  configs.sample <- floor(runif(s.size, min=2, max=(2^n.patches)))
  all.configs<-matrix(rep(NA, (n.patches*s.size)), nrow=s.size, ncol=n.patches)
  for (i in 1:length(configs.sample)){
    all.configs[i,]<-BigintToBin(configs.sample[i], digits=(n.patches))
  }
  
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