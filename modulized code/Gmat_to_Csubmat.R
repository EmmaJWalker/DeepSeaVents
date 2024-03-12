##########################################################################################################
# DESCRIPTION:
##########################################################################################################

# This function simply removes transitions into and out of the absorbing (all-off) state the CTMC's 
# generator matrix (G.mat) to provide the submatrix (C.mat), describing only the transition rates 
# among the transient (non-absorbing states). Works the same regardless of whether G.mat is describing
# transitions among states corresponding to the number of patches on or among specific configurations
# of patches on and off since the first row and column always corresponds to the absorbing (all-off) 
# state.

# INPUTS:
# G.mat: the generator matrix for the habitat CTMC

# OUTPUTS:
# C.submat: the submatrix of transition rates among the transient (nonabsorbing states) of the habitat CTMC

# REQUIRED PACKAGES:
# none

####################################################################################
# FUNCTION CODE:
####################################################################################

Gmat_to_Csubmat<-function(G.mat){
  C.submat<-G.mat[-1,-1]
  return(C.submat)
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
#G.mat<-make_Gmat(n.patches=2,r.on=0.1,r.off=0.001)
#G.mat
#C.submat<-Gmat_to_Csubmat(G.mat)
#C.submat
##########################################################################################################
