##########################################################################################################
# DESCRIPTION:
##########################################################################################################

# Because patches randomly and independently are gained and lost, the probabilities of habitat configurations
# with the same number of habitat patches present are equal at QED. Therefore, the QED of the reduced CTMC
# describing only the number of habitat patches on converts to the QED of the full CTMC describing every 
# possible habitat configuration by dividing those probabilities by the number of those habitat configurations 
# with the same number of habitat patches on.
# 
# This function uses this property to covert the QED probabilities of states with a given number of patches 
# on to the QED probabilities of each habitat configuration with that number of habitat patches on. 

# INPUTS:
# QED: a vector is length n.patches where each element gives the probability of being in a state with a given 
# number of habitat patches on
# n.patches: how many patches total may be in the landscape

# OUTPUTS:
# QED: a vector is length n.patches where each element gives the probability of being in a
# specific habitat configuration at QED given the number of habitat patches present in that configuration

# REQUIRED PACKAGES:
# none

# REQUIRED FUNCTIONS:
# none

####################################################################################
# FUNCTION CODE:
####################################################################################
reduced_to_full_QED<-function(QED, n.patches){
  n.on.vec<-1:n.patches
  z<-factorial(n.patches)/(factorial(n.on.vec)*factorial(n.patches-n.on.vec)) #calculating N_choose_k
  #aka number of different configurations possible with that number of patches on
  QED<-QED/z #converting the probabilities n patches are on at QED to 
  #probabilities of being in each possible configuration with n patches on
  
  return(QED)
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
#setwd("C:/Users/Administrator/Desktop/JAN 2024/modulized code")
#source("make_reduced_Gmat.r")
#source("Gmat_to_Csubmat.r")
#source("get_QED.r")
#source("IDs_to_configs.r")
#G.mat<-make_reduced_Gmat(n.patches=30,r.on=0.1,r.off=0.1)
#C.submat<-Gmat_to_Csubmat(G.mat)
#output<-get_QED(C.submat)
#QED<-output[[1]]
#QED<-reduced_to_full_QED(QED=QED, n.patches=30)
#QED
####################################################################################