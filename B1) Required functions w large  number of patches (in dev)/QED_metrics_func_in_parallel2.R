#Finds the equilibrium metapopulation size for each configuration
###############################################################################
QED.metrics<-function(c.sample.size, landscape, e.rate, c.rate, gamma, epsilon, self.rec, alpha){
  
  n.patches<-length(landscape$patch.ID)
  
  #calculate the extinction rate for each patch
  #extinction.rates<-e.rate/areas #if want this to be patch dependent
  #calculate delta
  delta<-e.rate/c.rate
  
  ##1. create a dataframe with all possible configurations
  ##ordered according to the sequence they occur as rows in our transition matrix
  ##NOTE: updated to store these as sparse matrices to save memory :)
  #n<-n.patches
  #all.configs <- Matrix(0, nrow=2^n, ncol=n, sparse=TRUE) #sparse matrix version
  ##all.configs<-matrix(rep(NA, (2^n)*n), nrow=2^n, ncol=n) #no sparse matrix version
  #for (i in 1:n){
  #  #if (i==1){BB<-c(0,1) #no sparse vector version
  #  if (i==1){BB<-as(c(0,1), "sparseVector") #sparse vector version
  #  }else{BB<-c(BB[1:(2^(i-2))],BB,BB[((2^(i-2))+1):(2^(i-1))])}
  #  all.configs[,(n-i+1)]<-t(rep(BB,2^(n-i)))}
  
  
  #USE ONLY A SUBSAMPLE OF THE POSSIBLE CONFIGURATIONS:
  library(R.utils)
  configs.sample<-sample((2:(2^n.patches)), size=c.sample.size, replace=F)
  configs.sample<-intToBin(configs.sample)
  all.configs<-matrix(rep(NA, (n.patches*c.sample.size)), nrow=c.sample.size, ncol=n.patches)
  for (i in 1:length(configs.sample)){
    x<-strsplit(configs.sample[i], "")
    all.configs[i,]<-ifelse(x[[1]]=="1",1,0)
  }
  
  #########################################
  #IN PARALLEL:
  library("doParallel") 
  cl <- makeCluster(detectCores(), type='PSOCK')
  registerDoParallel(cl)
  
  each.config.calcs<-function(i, n.patches){
    
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
    getwd()
    setwd("I:/one rep 30 patches")  
    
    # foreach (1:45){
    sub.i.configs<-c((((i-1)*48)+1):(i*48)) #define the set of configs to be run on core i
    i.metrics<-data.frame(matrix(data=rep(NA, 48*(5+n.patches)), nrow=48, ncol=(5+n.patches)))
    for (k in sub.i.configs) { #for each config core i...
      if (k>(2^n.patches)){ #checks if k is greater than the number of configs it isn't a config so do nothing
      }else{ #otherwise,
        config<-all.configs[k,] #selects a config k...
        if (k==1){ #checks if it's the absorbing state...where
          #the first row would correspond with the absorbing state 
          #therefore, the persistence capacity and metapop size for this config should be 0
          lambda.M<-0
          lambda.M.delta<-0
          pstar.configs<-rep(0, n.patches)
          metapop.size<-0
          n.on<-0
          i.metrics[(k-(48*(i-1))),]<-data.frame(k, lambda.M, lambda.M.delta, t(pstar.configs), metapop.size, n.on)
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
          i.metrics[(k-(48*(i-1))),]<-data.frame(k, lambda.M, lambda.M.delta, metapop.size, n.on, t(pstar.configs))
          #combines these all in a data frame of all the calculations done on core i (lazy and inefficient way...)
          #if (k>45){
          #  i.metrics<-data.frame(k, lambda.M, lambda.M.delta, t(pstar.configs), metapop.size, n.on) 
          #}else{i.metrics<-rbind(i.metrics, data.frame(k, lambda.M, lambda.M.delta, t(pstar.configs), metapop.size, n.on))}
        }
      }
    }
    i.metrics<-na.omit(i.metrics) #to remove extra empty rows of data since 2^n.patches may not be a multiple of 48
    return(i.metrics)
  } 
  
  config.calcs<-foreach(i=1:48, .combine=rbind) %dopar% each.config.calcs(i=(i), n.patches=n.patches)
  # turn parallel processing off and run sequentially again:
  registerDoSEQ()
  #stopCluster()
  #stopImplicitCluster()
  
  #END PARALLEL:
  ######################################
  
  metapop.size<-rowSums(config.calcs$pstar.configs) #calculate total metapop size expected at eq in each config
  
  #get the QED of the reduced CTMC
  QED.output<-quasi.eq(n.patches, r.on, r.off)
  QED<-QED.output[[1]] #QED distribution
  # use this to get the probabilities of just the subsampled configurations
  #converting the QED of the reduced CTMC to the QED of the full CTMC:
  n.on.vec<-n.on.vec[-1] #can't have none on at QED by definition
  n.on.vec<-n.on.vec[-(n.patches+1)] #can't have more than n.patches on either (this was just an extra calc for the sake of making the matrices)
  x<-factorial(n.patches)/(factorial(n.on.vec)*factorial(n.patches-n.on.vec)) #calculating N_choose_k
  #aka number of different configurations possible with that number of patches on
  QED.dist/x #converting the QED to the QED for each configuration rather than for the number of patches on
  
  n.on<-all.configs
  
  #get the time-averaged metapopulation sizes @QED
  ar.pstar.QED<-dot(QED,metapop.size[2:2^n.patches]) #this is the expected metapop size in a given configuration x the porportion of time spent in it at QED and summed (therefore the arithmentic mean)
  #should actually be the geometric mean
  geom.pstar.QED<-exp(dot(QED,log(metapop.size[2:2^n.patches])))
  lm.QED.delta<-exp(dot(QED,log(lambda.M.delta[2:2^n.patches])))
  lm.QED<-exp(dot(QED,log(lambda.M[2:2^n.patches])))
  
  #get expected number of patches in habitat @QED
  n.on.in.each.config<-rowSums(all.configs)
  #geometric mean habitat expected across configurations and time spent in them:
  exp.n.QED<-exp(dot(QED,log(n.on.in.each.config[2:2^n.patches])))
  
  all.configs<-as.matrix(all.configs)
  pstar.configs<-as.matrix(pstar.configs)
  
  return(list(lm.QED.delta, lm.QED, geom.pstar.QED, ar.pstar.QED, exp.n.QED, QED, all.configs, pstar.configs))
}