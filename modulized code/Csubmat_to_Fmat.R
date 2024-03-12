##########################################################################################################
# DESCRIPTION:
##########################################################################################################

# This function uses the generator matrix (G.mat) for the CTMC to calculate the Fundamental Matrix (F.mat),
# mathematically defined as the negative inverse of C.mat

# INPUTS:
# C.submat: the submatrix of transition rates among the transient (nonabsorbing states) of the habitat CTMC

# OUTPUTS:
# F.mat: the Fundamental Matrix

# REQUIRED PACKAGES:
# Matrix: for working with sparse matrices

####################################################################################
# FUNCTION CODE:
####################################################################################
Csubmat_to_Fmat<-function(C.submat){
F.mat<--Matrix::solve(C.submat) #solve gives the inverse of a matrix 
return(F.mat)
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
#F.mat<-Csubmat_to_Fmat(C.submat)
#F.mat
##########################################################################################################
