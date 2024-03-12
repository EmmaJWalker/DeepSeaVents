##########################################################################################################
# DESCRIPTION:
##########################################################################################################

# This function uses the fundamental matrix (F.mat) for the CTMC to calculate the ratio of means 
# distribution (RMD), aka the frequency of each state (given by each column) provided the system began 
# in a given state (given by each row) before absorption (ending in  the all-off state).

# The ratio of means distribution converges to the Q.E.D. provided time to absorption is long enough
# and thus, can be used as a way to check the CTMC converges to the Q.E.D. sufficiently quickly relative 
# to the absorbing state, such that the Q.E.D. describes the frequency of each state regardless of
# what state the CTMC may have began in.

# INPUTS:
# F.mat: the submatrix of transition rates among the transient (non-absorbing states) of the habitat CTMC

# OUTPUTS:
# RMD: the ratio of means distribution

# REQUIRED PACKAGES:
# none

####################################################################################
# FUNCTION CODE:
####################################################################################
Fmat_to_RMD<-function(F.mat){
RMD<-F.mat/rowSums(F.mat)
  return(RMD)
}
####################################################################################
# QUICK CHECK: (uncomment and run to check)
####################################################################################
#rm(list=ls()) #to ensure a clean environment
## now run the FUNCTION CODE!!!
#library("Matrix") #to get sparse matrices
#library("gdata") #to get upper and lower triangles of matrices
#library("lava") #to get antidiagonal of a matrix using revdiag
#setwd("C:/Users/Administrator/Desktop/JAN 2024/modulized code")
#source("make_Gmat.r")
#source("Gmat_to_Csubmat.r")
#source("Csubmat_to_Fmat.r")
#G.mat<-make_Gmat(n.patches=2,r.on=0.1,r.off=0.001)
#C.submat<-Gmat_to_Csubmat(G.mat)
#F.mat<-Csubmat_to_Fmat(C.submat)
#RMD<-Fmat_to_RMD(F.mat)
#RMD
##########################################################################################################
