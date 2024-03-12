##########################################################################################################
# DESCRIPTION:
##########################################################################################################

# This function constructs the generator matrix (G.mat) for the CTMC, describing the transition rates 
# from each habitat configuration (with a configuration ID corresponding to its row in G.mat)
# to another (with a configuration ID corresponding to its column in G.mat - see the Configuration IDs 
# document for details). 

# It does so by taking advantage of the block toeplitz structure of the generator matrix to save on
# computation and memory - see the Configuration IDs document for details

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

####################################################################################
# FUNCTION CODE:
####################################################################################

make_Gmat<-function(n.patches, r.on, r.off){
  
  #1. Build the base block toeplitz matrix
  #corresponding with the 4 base transitions that can occur
  #1) a patch that was off can turn on
  #2) a patch that was on can turn off
  #3) a patch that was on can stay on
  #4) a patch that was off can stay off
  A<-matrix(0,4,4) #make a 4 by 4 matrix
  #where the upper triangle of transitions is to at least one patch being on
  upperTriangle(A, diag=FALSE, byrow=FALSE) <- r.on 
  #and the lower is to at least patches being off
  lowerTriangle(A, diag=FALSE, byrow=FALSE) <- r.off 
  revdiag(A)<-0   #since both patches can't go off or on at the same time
  A <- Matrix(A, sparse=TRUE)
  
  #2. Build the generator matrix from the base blocks
  if (n.patches==2){} else {
    #if there are only 2 patches the generator matrix is just composed of this base 
    #toeplitz matrix A, but for each additional patch, 4 possible transitions x each 
    #of the previously possible transitions excluding that single patch are added):
      #starting with that base toeplitz matrix for two patches
      #the offset diagonal toeplitz matrices iteratively scaling as patches are added 
      #iteratively as follows:
    for (n in 2:(n.patches-1)){ #for each patch excluding the added patch
      #a sparse diagonal matrix of transitions to all other states from that patch 
      #being on is added
      r.on.diag<-bandSparse(2^(n),2^(n),0,list(rep(r.on, 2^n+1))) 
      #a sparse diagonal matrix of transitions to all other states from that patch 
      #being off is added
      r.off.diag<-bandSparse(2^(n),2^(n),0,list(rep(r.off, 2^n+1))) 
      #then make a toeplitz matrix of those 3 toeplitz matrices by binding them 
      #appropriately together
      A1<-cbind(A, r.on.diag) 
      A2<-cbind(r.off.diag, A)
      A<-rbind(A1, A2)
    }
  }
    
  #3. Since a CTMC, ensure that the rows sum to 0 with the negative sum on the diagonal
  diag(A)<--Matrix::rowSums(A)
    
  #4. Set the all-off state to be absorbing
  A[1,]<-0 
    
  #A now provides us with the infinitessimal generator matrix!
  G.mat<-(A)
  return(G.mat)
}

####################################################################################
# QUICK CHECK: (uncomment and run to check)
####################################################################################
#rm(list=ls()) #to ensure a clean environment
## now run the FUNCTION CODE!!!
#library("Matrix") #to get sparse matrices
#library("gdata") #to get upper and lower triangles of matrices
#library("lava") #to get antidiagonal of a matrix using revdiag
#G.mat<-make_Gmat(n.patches=5,r.on=0.1,r.off=0.001)
#G.mat
##########################################################################################################

