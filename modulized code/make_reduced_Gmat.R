##########################################################################################################
# DESCRIPTION:
##########################################################################################################

# This function constructs the generator matrix (G.mat) for the reduced CTMC, describing the transition rates 
# from habitat configurations with some number of patches on (corresponding to the row in G.mat)
# to some those with some other number of patches on (corresponding to the column in G.mat - see the Configuration IDs 
# document for details). 

# It does so by taking advantage of the toeplitz structure of the generator matrix to save on
# computation and memory - see the Configuration IDs document for details

# NOTE: the following warning messages are suppressed in the function code because they 
# inconsequential and to be ignored: 
# Warning messages:
# 1: In x[upper.tri(x, diag = diag)] <- value :
#  number of items to replace is not a multiple of replacement length
# 2: In ret[rev(lower.tri(x, diag = diag))] <- value :
#  number of items to replace is not a multiple of replacement length

# INPUTS:
# n.patches: the number of habitat patches that could be on or off in the system
# r.on: the rate at which patches turn on
# r.off: the rate at which patches turn off

# OUTPUTS:
# G.mat: the generator matrix for the habitat CTMC

# REQUIRED PACKAGES:
# Matrix: to get sparse matrix functions
# gdata: to get upper and lower triangles of matrices
# lava: to get the anti-diagonal of a matrix using revdiag

##########################################################################################################
# FUNCTION CODE:
##########################################################################################################

make_reduced_Gmat<-function(n.patches, r.on, r.off){
  
  #1. Build the toeplitz matrix
  n.on.vec<-seq(0,(n.patches+1),1) #provides vector of states indicating how many patches are on
  up.diag<-(n.patches-n.on.vec)*r.on #calculates the values for the upper diagonal of the matrix
  up.diag<-up.diag[-(n.patches+1)] #removing the last since that calculation goes beyond the range necessary
  lo.diag<-n.on.vec*r.off #calculates the values for the lower diagonal of the matrix
  lo.diag<-lo.diag[-1] #removing the first since that calculation goes beyond the range necessary
  A<-matrix(0,(n.patches+1),(n.patches+1)) #make an empty matrix of correct dimensions
  M<-A #make a copy to insert the upper diagonal we want 
  N<-A #make a copy to insert the lower diagonal we want
  #put these in as the diagonals to these respective matrices
  diag(M)<-up.diag 
  diag(N)<-lo.diag
  #and use them to populate our matrix appropiately
  suppressWarnings(upperTriangle(A, diag=F, byrow=F)<-upperTriangle(M, diag=T, byrow=F))
  suppressWarnings(lowerTriangle(A, diag=F, byrow=T)<-lowerTriangle(N, diag=T, byrow=T))
  #ignore the warning messages: it by default trims the matrices as I want them to be trimmed such that the 
  #appropriate values populate the upper and lower triangles
  
  #2. Since a CTMC, ensure that the rows sum to 0 with the negative sum on the diagonal
  diag(A)<--rowSums(A) 
  
  #3. Set the all-off state to be absorbing
  A[1,]<-0 #since it's an absorbing matrix the first row corresponding to all patches off has no rates of transition out of this state
  
  #4. Make it a sparse matrix (for consistency)
  A <- Matrix(A, sparse=TRUE) 
  
  #A now provides us with the infinitessimal generator matrix!
  G.mat<-(A)
  
  return(G.mat)
}

##########################################################################################################
# QUICK CHECK: (uncomment and run to check)
##########################################################################################################
#rm(list=ls()) #to ensure a clean environment
## now run the FUNCTION CODE!!!
#library("Matrix") #to get sparse matrices
#library("gdata") #to get upper and lower triangles of matrices
#library("lava") #to get antidiagonal of a matrix using revdiag
#G.mat<-make_reduced_Gmat(n.patches=5,r.on=0.1,r.off=0.001)
#G.mat

##########################################################################################################

