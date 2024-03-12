# This function computes expectations of metapopulation size and persistence, and number of habitat patches at QED by
# taking the landscape capacity (lambda.M), persistence capcity (lambda.M.delta), equilibrium occupancy (metapop.size), 
# and number of habitat patches present within habitat configurations and computing weighted averages based on the
# probabilities with which these habitat configurations occur at QED.

# INPUTS:
# config.metrics:a dataframe as output by the where each row contains:
#                 - lambda.M: the landscape capacity 
#                 - lambda.M.delta: the persistence capacity 
#                 - metapop.size: total equilibrium occupancy 
#                 - n.on: total.number of patches on in the habitat configuration
#                 - p.star: equilibrium occupancy of each patch in each subsequent column
#                 within a given configuration

# OUTPUTS:
# a list of objects containing:
#                 - lm.QED.delta: the geometric mean persistence capacity at QED
#                 - lm.QED: the geometric mean landscape capacity at QED
#                 - geom.pstar.QED: geometric mean equilibrium occupancy of the landscape by the metapopulation at QED
#                 - ar.pstar.QED: arithmetic mean equilibrium occupancy of the landscape by the metapopulation at QED
#                 - exp.n.QED: expected mean number of patches present at QED
#                 - pstar.configs: a dataframe of the equilibrium occupancies of patches within each configuration 
#                   supplied to this function 

# REQUIRED PACKAGES:
# none 

# REQUIRED FUNCTIONS:
# none
##########################################################################################################
get_QEDMetapopMetrics<-function(config.metrics, QED){
  
  lambda.M<-config.metrics[,1]
  lambda.M.delta<-config.metrics[,2]
  metapop.size<-config.metrics[,3]
  n.on<-config.metrics[,4]
  pstar.configs<-config.metrics[,-(1:4)]
  
  #get the time-averaged metapopulation sizes @QED
  ar.pstar.QED<-dot(QED,metapop.size)
  #this is the expected metapop size in each configuration x the porportion of time spent in it at QED and summed (therefore the arithmentic mean)
  #should actually be the geometric mean
  geom.pstar.QED<-exp(dot(QED,log(metapop.size)))
  lm.QED.delta<-exp(dot(QED,log(lambda.M.delta)))
  lm.QED<-exp(dot(QED,log(lambda.M)))
  #geometric mean habitat expected across configurations and time spent in them:
  exp.n.QED<-exp(dot(QED,log(n.on)))
  
  return(list(lm.QED.delta, lm.QED, geom.pstar.QED, ar.pstar.QED, exp.n.QED, pstar.configs))
}
##########################################################################################################
# QUICK CHECK: (uncomment and run to check)
##########################################################################################################
#rm(list=ls()) #to ensure a clean environment
## now run the FUNCTION CODE!!!
#library("Matrix") #to get sparse matrices
#library("gdata") #to get upper and lower triangles of matrices
#library("lava") #to get anti-diagonal of a matrix using revdiag
#library("gtools") #to get even() function
#library("tidyverse") #to get filter function
#library("doParallel") #to run code in parallel
#setwd("C:/Users/Administrator/Desktop/JAN 2024/modulized code")
#source("make_Gmat.r")
#source("Gmat_to_Csubmat.r")
#source("get_QED.r")
#source("make_landscape.r")
#source("get_distmat.r")
#source("get_dispkernel.r")
#source("get_lambdaM.r")
#source("get_pstar.r")
#source("get_MetapopMetrics.r")
#source("make_allConfigs.r")
#source("get_config_metrics.r")
#n.patches<-10
#landscape<-make_landscape(n.patches=n.patches, landscape.type="line", landscape.size=100, 
#                          p.dist="uniform", p.size.dist="uniform", max.p.size=1, 
#                          clust.its=1000)
#G.mat<-make_Gmat(n.patches=n.patches,r.on=0.1,r.off=0.1)
#C.submat<-Gmat_to_Csubmat(G.mat)
#QED.output<-get_QED(C.submat)
#QED<-QED.output[[1]]
#QED<-QED/sum(QED) #standardizing to account for rounding errors
#configs<-make_allConfigs(n.patches=n.patches)
#configs<-configs[-1,] #remove the all off state
#config.metrics<-matrix(rep(NA, (n.patches+4)*nrow(configs)), nrow=nrow(configs), ncol=(n.patches+4))
#for (i in 1:nrow(configs)){
#  config.metrics[i,]<-vec(get_config_metrics(config=configs[i,], landscape=landscape, e.rate=0.1, c.rate=0.2, gamma=0, 
#                     self.rec=1, alpha=0.1, iterations=1000))
#}
#QEDmetapop.metrics<-get_QEDMetapopMetrics(config.metrics, QED)
#QEDmetapop.metrics
##########################################################################################################











