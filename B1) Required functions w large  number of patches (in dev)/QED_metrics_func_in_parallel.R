#Finds the equilibrium metapopulation size for each configuration
###############################################################################
QED.metrics<-function(landscape, e.rate, c.rate, gamma, epsilon, self.rec, alpha){
  
  n.patches<-length(landscape$patch.ID)
  
  #calculate the extinction rate for each patch
  #extinction.rates<-e.rate/areas #if want this to be patch dependent
  #calculate delta
  delta<-e.rate/c.rate
  
  #1. create a dataframe with all possible configurations
  #ordered according to the sequence they occur as rows in our transition matrix
  #NOTE: updated to store these as sparse matrices to save memory :)
  n<-n.patches
  all.configs <- Matrix(0, nrow=2^n, ncol=n, sparse=TRUE) #sparse matrix version
  #all.configs<-matrix(rep(NA, (2^n)*n), nrow=2^n, ncol=n) #no sparse matrix version
  for (i in 1:n){
    #if (i==1){BB<-c(0,1) #no sparse vector version
    if (i==1){BB<-as(c(0,1), "sparseVector") #sparse vector version
    }else{BB<-c(BB[1:(2^(i-2))],BB,BB[((2^(i-2))+1):(2^(i-1))])}
    all.configs[,(n-i+1)]<-t(rep(BB,2^(n-i)))}
  
  #########################################
  #IN PARALLEL:
  library("doParallel") 
  cl <- makeCluster(detectCores(), type='PSOCK')
  registerDoParallel(cl)
  
  each.config.calcs<-function(i, n.patches, landscape, e.rate, c.rate, gamma, epsilon, self.rec, alpha){
    
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
    v<-floor((2^n.patches)/48) #define the set of configs to be run on core i
    #if (i==46){}else{ #core number 46 doesn't work it seems
    if (i==48 & ((2^n.patches)/48)>floor((2^n.patches)/48)){ #if it's core number 48
      v<-v+(2^n.patches)-(floor((2^n.patches)/48)*48) #give the remaining configs to it 
      sub.i.configs<-c(((2^n.patches)-v+1):(2^n.patches))
    }else{
      sub.i.configs<-c((((i-1)*v)+1):(i*v))
    }
    i.metrics<-data.frame(matrix(data=rep(NA, v*(5+n.patches)), nrow=v, ncol=(5+n.patches)))
    for (k in sub.i.configs) { #for each config core i...
    #  if (k>(2^n.patches)){ #checks if k is greater than the number of configs it isn't a config so do nothing
    #  }else{ #otherwise,
        config<-all.configs[k,] #selects a config k...
        if (k==1){ #checks if it's the absorbing state...where
          #the first row would correspond with the absorbing state 
          #therefore, the persistence capacity and metapop size for this config should be 0
          lambda.M<-0
          lambda.M.delta<-0
          pstar.configs<-rep(0, n.patches)
          metapop.size<-0
          n.on<-0
          i.metrics[k,]<-data.frame(k, lambda.M, lambda.M.delta, t(pstar.configs), metapop.size, n.on)
          #i.metrics<-data.frame(k, lambda.M, lambda.M.delta, t(pstar.configs), metapop.size, n.on)
        } else { #otherwise
          #calculates the persistence capacity in that configuration
          landscape.config<-landscape[config!=0,]
          lambda.M<-get.lambda.M(landscape.config, alpha, gamma, epsilon, self.rec, e.rate, c.rate)
          lambda.M.delta<-get.lambda.M(landscape.config, alpha, gamma, epsilon, self.rec, e.rate, c.rate)/delta
          if (lambda.M.delta < delta){
            #metapop.size[i]<-0
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
          print(k)
          i.metrics[k,]<-data.frame(k, lambda.M, lambda.M.delta, metapop.size, n.on, t(pstar.configs))
          #combines these all in a data frame of all the calculations done on core i (lazy and inefficient way...)
          #if (k>45){
          #  i.metrics<-data.frame(k, lambda.M, lambda.M.delta, t(pstar.configs), metapop.size, n.on) 
          #}else{i.metrics<-rbind(i.metrics, data.frame(k, lambda.M, lambda.M.delta, t(pstar.configs), metapop.size, n.on))}
        }
      }
   # }
   # }
    i.metrics<-na.omit(i.metrics) #to remove extra empty rows of data since 2^n.patches may not be a multiple of 48
    return(i.metrics)
  } 
  
  results<-foreach(i=1:48, .combine=rbind) %dopar% each.config.calcs(i=(i), n.patches=n.patches,
                                                                     landscape=landscape, e.rate=e.rate, c.rate=c.rate, 
                                                                     gamma=gamma, epsilon=epsilon, self.rec=self.rec, 
                                                                     alpha=alpha)
  # turn parallel processing off and run sequentially again:
  registerDoSEQ()
  #stopCluster()
  #stopImplicitCluster()
  
  #END PARALLEL:
  
  ######################################
  return(results)
}