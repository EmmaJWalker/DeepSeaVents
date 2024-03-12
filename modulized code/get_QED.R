##########################################################################################################
# DESCRIPTION:
##########################################################################################################

# This function calculates the quasi-equilibrium distribution (Q.E.D. or QED) which provides the frequency
# of each of the transient states while not in the absorbing (all-off) state, regardless of the initial
# state provided convergence to the Q.E.D. is fast relative to the absorbing (all-off state). It also 
# calculates a rule of thumb to check if this is the case.

# INPUTS:
# C.submat: the submatrix of transition rates among the transient (non-absorbing states) of the habitat CTMC

# OUTPUTS:
# QED: the quasi-equilibrium distribution (Q.E.D.)
# rho1vrho2x2: rule of thumb threshold for when the QED well describes the frequency of states regardless
#              of the initial state (when <1)

# REQUIRED PACKAGES:
# Matrix: for working with sparse matrices

####################################################################################
# FUNCTION CODE:
####################################################################################
get_QED<-function(C.submat){

#1. calculate the eigenvalues and eigenvectors (right and left) and organize these
E<-eigen(C.submat) #calculate the right eigenvectors and eigenvalues 
EL <- eigen(Matrix::t(C.submat)) #calculate the left eigenvectors   
C1<-E$vectors #F1 is the set of right eigenectors
C2<-diag(E$values) #F2 is diagonal matrix of the eigenvalues
C3 <- EL$vectors #F3 is the set of left eigenvectors

#2. Calculate the QED given by the left eigenvector corresponding to the eigenvalue 
#,with maximal real part, of subgenerator matrix C
rho.1<-max(Re(E$values)) #rho.1 is the eigenvalue with max real part
rho.2<-max(Re(E$values[(E$values!=rho.1)])) #rho.2 is the next largest eigenvalue
max.index<-which.max(diag(C2)) #finds the location of rho.1
QED<-abs(C3[,max.index]) #take the absolute values of the left eigenvector corresponding to rho.1
#QED<-QED/sum(QED) #gives relative frequencies *******************************moved this to outside the function
#*****************************************************************************so that this step can be performed after
#*****************************************************************************converting from reduced to full QED

#3. Calculate a rule of thumb to check the rate of convergence to QED vs. absorption
rho1vrho2x2<-2*rho.1/rho.2 #compares the rate of convergence to QED versus to absorption
#the QED provides a good estimate of dynamics when this is <1, poor when >1

return(list(QED, rho1vrho2x2))
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
#QED<-get_QED(C.submat)
#QED
#output<-get_QED(C.submat)
#QED<-output[[1]]
#QED<-QED/sum(QED) #standardize to account for any rounding errors (so total probability=1)
#rho1vrho2x2<-output[[2]]
##########################################################################################################









