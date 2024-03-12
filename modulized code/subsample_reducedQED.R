# subsample configs w probability QED (calculated from reduced CTMC)

##########################################################################################################
# DESCRIPTION:
##########################################################################################################

# This function constructs the configurations (binary vectors indicating which patches are on and off) for 
# just a subset of all the configurations possible for a given number of patches chosen according to the
# probability with which these configurations occur at Q.E.D. It uses an algorithm detailed in the 
# Configuration IDs document.

# INPUTS:
# QED: the quasi-equilibrium distribution of states with a given number of habitat patches on
# n.patches: how many patches total are in the system
# s.size: how many configurations to choose (max all possible configurations 2^n.patches-1 because the all 
#         off state is excluded)

# OUTPUTS:
# configs.sample: a dataframe of binary vectors indicating which patches are on and off in which each
#                 row constitutes a configuration of patches on (1) or off (0) and each column corresponds
#                 to a given patch

# REQUIRED PACKAGES:
# tidyverse: to get filter function
# oddeven: #to get even() function

# REQUIRED FUNCTIONS:
# IDsToConfigs.r: 
# reduced_to_full_QED: 

####################################################################################
# FUNCTION CODE:
####################################################################################
subsample_reducedQED<-function(QED, n.patches, s.size){
  
  #STEP 1: convert the QED of the reduced CTMC to the QED of the full CTMC:
  QED.dist<-reduced_to_full_QED(QED=QED, n.patches=n.patches)
  
  #STEP 2: use the following algorithm to choose configurations as ordered by their configuration IDs for the full CTMC
  # based on the number of patches on in them according to their probability at QED (see Configuration ID document for details)
  n.on.vec<-1:n.patches
  z<-factorial(n.patches)/(factorial(n.on.vec)*factorial(n.patches-n.on.vec)) #calculating N_choose_k
  z.tracker<-z
  #aka number of different configurations possible with that number of patches on
  config.IDsample<-rep(NA, s.size)
  c.tracker<-data.frame(matrix(rep(NA, s.size*2), ncol=2, nrow=s.size))
  for (i in 1:s.size){
    if(sum(z.tracker>0)==1){ #if only one possibility for the number of on patches remains
      n.on.sample<-n.on.vec[z.tracker>0] #set that to be the number of patches on
    }else{ #otherwise,
      n.on.sample <- sample(n.on.vec[z.tracker>0], 1, prob=QED.dist[z.tracker>0]/sum(QED.dist[z.tracker>0])) #sub sample how many patches could 
      #possibly be on with the probability that number of patches in on at QED
    }
    prior.configs<-cumsum(z[1:(n.on.sample)])
    ########################################### REQUIRES TOO MUCH MEMORY (requires vector size (2^n.patches)/n.patches):
    #### NOTE: this method is better for moderate numbers of patches (e.g. up to 30-40) when s.size is large relative to 
    #### the number of configs possible since the probability of choosing a configuration with the same number of patches 
    #### on will often be high, and this method doesn't require repeated sampling until a unique index of a configuration 
    #### with that number of patches on is chosen #####NAH! THE SECCOND METHOD IS FASTER EVEN THEN!!!
    ####p<-(1:z[n.on.sample]) #index of possible configurations with that number of patches on
    ####q<-filter(c.tracker, X1==n.on.sample)
    ####if(any(q$X2)){
    ####  p<-p[!p %in% q$X2]
    ####  #p<-p[-which(p==q$X2)] #taking out any combinations of n and k we have already used
    ####  #from the possibilities to sample
    ####}
    ####if(length(p)==1){#if only one possible config with that many patches on remains,
    ####  k<-p #then choose that config
    ####}else{ #otherwise,
    ####  k<-sample(p, 1)  # within the
    ####  #number of configurations with that number of patches on,
    ####  #randomly choose a configuration with that number of patches on
    #####}
    ############################################
    ############################################ MUCH MORE MEMORY EFFICIENT!:
    ##### NOTE: this method is better for very large numbers of patches or when s.size is small 
    ##### relative to the number of possible configs because the probability of choosing the same
    ##### configuration with the same number of patches is on is likely almost always very low
    
    k<-round(runif(n=1, min=1, max=z[n.on.sample])) #randomly choose the index of a config with a given number of patches on
    #by randomly choosing a number between 1 and the number of configs possible with that number of patches on
    if (all(is.na(c.tracker$X2[c.tracker$X1==n.on.sample]))){ # if no config with the given number of patches on has been 
      #chosen before, there is no need to pick a new k, do nothing
    } else { #otherwise,
      # check if the number of unique configs chosen including the one chosen with a given number of patches on 
      # is less than the number of configs chosen including k, if not then k must have been chosen before and 
      # keep picking a new k until that is not the case
      c.used<-c.tracker[!is.na(c.tracker$X1),]
      while (length(unique(c(c.used$X2[c.used$X1==n.on.sample],k)))<length(c(c.used$X2[c.used$X1==n.on.sample],k))){ 
        k<-round(runif(1, 1, z[n.on.sample])) #keep choosing new k's till k is unique
      }
    }
    ###########################################  
    if(n.on.sample==1){
      config.IDsample[i]<-k+1
    }else{
      config.IDsample[i] <- prior.configs[n.on.sample-1] + k + 1 #get the index of that unique config out 
      #of all the possible configs
    }
    z.tracker[n.on.sample]<-z.tracker[n.on.sample]-1 #now one less configuration with that many patches on is possible
    c.tracker$X1[i]<-n.on.sample #keep track of combinations of n's and
    c.tracker$X2[i]<-k #k's chosen in a data frame
  }
  
  #STEP 3: Convert the chosen configuration IDs to the configurations of on and off patches they represent
  configs.sample<-matrix(rep(NA, (n.patches*s.size)), nrow=s.size, ncol=n.patches)
  for (i in 1:length(config.IDsample)){
    configs.sample[i,]<-IDsToConfigs(x=config.IDsample[i], n.patches=n.patches)
  }
  
  return(configs.sample)
  ###return(config.IDsample)
}

####################################################################################
# QUICK CHECK: (uncomment and run to check)
####################################################################################
#rm(list=ls()) #to ensure a clean environment
## now run the FUNCTION CODE!!!
#library("Matrix") #to get sparse matrices
#library("gdata") #to get upper and lower triangles of matrices
#library("lava") #to get anti-diagonal of a matrix using revdiag
#library("gtools") #to get even() function
#library("tidyverse") #to get filter function
#setwd("C:/Users/Administrator/Desktop/JAN 2024/modulized code")
#source("make_reduced_Gmat.r")
#source("Gmat_to_Csubmat.r")
#source("get_QED.r")
#source("reduced_to_full_QED.r")
#source("IDs_to_configs.r")
#G.mat<-make_reduced_Gmat(n.patches=50,r.on=0.1,r.off=0.1)
#C.submat<-Gmat_to_Csubmat(G.mat)
#output<-get_QED(C.submat)
#QED<-output[[1]]
#Sys.time()
#configs.sample<-subsample_reducedQED(QED=QED, n.patches=50, s.size=1000)
#Sys.time()
#configs.sample
##########################################################################################################






