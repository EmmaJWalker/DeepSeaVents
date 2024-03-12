##########################################################################################################
# DESCRIPTION:
##########################################################################################################

# This function calculates the average time to absorption (all-off state) from each state (organized by 
# row) using C.submat.

# Mathematically, this is given by summing across rows of the fundamental matrix and multiplying by the 
# trace of (I-C.submat) (i.e. the rate of leaving states, where I is an identity matrix)

# Note: This requires calculating the fundamental matrix which for large numbers of patches is a 
# computationally costly endeavor (should you wish to calculate this and the ratio of means distribution
# many times over it would be better to first calculate F.mat outside this function and then use it as an
# argument to run both the Fmat_to_RMD function and this function (adding it as an argument and 
# commenting out step 1). I have not bothered to do this here as I only intend to calculate the ratio
# of means distribution for a few test cases just for checking coding and convergence to Q.E.D.

# INPUTS:
# C.submat: the submatrix of transition rates among the transient (non-absorbing states) of the habitat CTMC
# optionally add F.mat (the fundamental matrix) to cut out step 1.

# OUTPUTS:
# T.absorp: the average time to absorption (all-off state) from each state

# REQUIRED PACKAGES:
# Matrix: for working with sparse matrices

####################################################################################
# FUNCTION CODE:
####################################################################################
get_Tabsorp<-function(C.submat){
  
  #1. calculate the fundamental matrix
  F.mat<--Matrix::solve(C.submat) #solve gives the inverse of a matrix
  
  #2. calculate the time to absorption
  c<-sum(-Matrix::diag(C.submat)) #the trace is just the sum of the diagonal elements
  T.absorp<-rowSums(F.mat)*c
  
  return(T.absorp)
}

####################################################################################
# QUICK CHECK: (uncomment and run to check)
####################################################################################
#rm(list=ls()) #to ensure a clean environment
## now run the FUNCTION CODE!!!
#library("Matrix") #to get sparse matrices
#library("gdata") #to get upper and lower triangles of matrices
#library("lava") #to get anti-diagonal of a matrix using revdiag
#setwd("C:/Users/Administrator/Desktop/JAN 2024/modulized code")
#source("make_Gmat.r")
#source("Gmat_to_Csubmat.r")
#G.mat<-make_Gmat(n.patches=2,r.on=0.1,r.off=0.001)
#C.submat<-Gmat_to_Csubmat(G.mat)
#T.absorp<-get_Tabsorp(C.submat)
#T.absorp
##########################################################################################################

