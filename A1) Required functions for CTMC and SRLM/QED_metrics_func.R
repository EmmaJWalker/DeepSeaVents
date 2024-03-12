#Finds the equilibrium metapopulation size for each configuration
###############################################################################
QED.metrics<-function(landscape, e.rate, c.rate, gamma, epsilon, self.rec, alpha, r.on, r.off){
  
  n.patches<-length(landscape$patch.ID)
  
  #calculate the extinction rate for each patch
  #extinction.rates<-e.rate/areas #if want this to be patch dependent
  #calculate delta
  delta<-e.rate/c.rate
  
  lambda.M<-rep(0,2^n.patches)
  lambda.M.delta<-rep(0,2^n.patches)
  pstar.configs<-matrix(rep(0,2^(n.patches+1)), ncol=(n.patches), nrow=(2^n.patches))
  metapop.size<-rep(0,2^n.patches)
  
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
  
  for(i in 1:nrow(all.configs)){
    if (i==1){
      #the first row would correspond with the absorbing state 
      #therefore, the persistence capacity and metapop size for this config should be 0
      lambda.M[i]<-0
      #metapop.size[i]<-0
      pstar.configs[i,]<-rep(0, n.patches)
    } else { #otherwise
      #calculate the persistence capacity in that configuration
      landscape.config<-landscape[all.configs[i,]!=0,]
      lambda.M[i]<-get.lambda.M(landscape.config, alpha, gamma, epsilon, self.rec, e.rate, c.rate)
      lambda.M.delta[i]<-get.lambda.M(landscape.config, alpha, gamma, epsilon, self.rec, e.rate, c.rate)/delta
      if (lambda.M.delta[i] < delta){
        #metapop.size[i]<-0
        pstar.configs[i,]<-rep(0, n.patches)
        } #if below the persistence threshold, the metapop is expected extinct at equilibrium 
      else { #otherwise
        #calculate the equilibrium patch occupancy of each patch in each config
        pstar.configs[i,all.configs[i,]!=0]<-pstar.function(landscape.config, alpha, e.rate, c.rate, self.rec, gamma, epsilon, iterations=1000)
        }
    }
    }
  
  metapop.size<-rowSums(pstar.configs) #calculate total metapop size expected at eq in each config
  
  #get the time-averaged metapopulation sizes @QED
  QED.output<-quasi.eq(n.patches, r.on, r.off)
  QED<-QED.output[[1]] #QED distribution
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
