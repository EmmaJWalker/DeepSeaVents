#rm(list=ls()) #clear the workspac

quasi.eq<-function(n.patches, r.on, r.off){
  #library("R.utils") #needed for creating a binary number
  
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
  upperTriangle(A, diag=F, byrow=F)<-upperTriangle(M, diag=T, byrow=F) 
  lowerTriangle(A, diag=F, byrow=T)<-lowerTriangle(N, diag=T, byrow=T)
  #ignore the warning messages: it by default trims the matrices as I want them to be trimmed such that the 
  #appropiate values populate the upper and lower triangles
  diag(A)<--rowSums(A) #since this is a CTMC (therefore transitions are instantaneous) rows must sum to 0 and therefore transitions out of static states happen at a rate equal to - the sum across rows
  A[1,]<-0 #since it's an absorbing matrix the first row corresponding to all patches off has no rates of transition out of this state
  A <- Matrix(A, sparse=TRUE) #just making it a sparse matrix for consistency with the code that uses the full CTMC
  
  #this provides us with the infinitessimal generator matrix A 
  G.mat<-(A)
  #G.mat
  
  #C.submat is the submatrix corresponding to the transient states (non-absorbing states)
  C.submat<-G.mat[2:(n.patches+1),2:(n.patches+1)]
  #sparse_C <- as(C.submat, "sparseMatrix")
  #print(C.submat)
  
  #the negetive inverse of C gives us the fundamental matrix N=(I-Q)^-1, 
  #which gives the ratio of means distribution
  fundamental.matrix<--Matrix::solve(C.submat) #solve gives the inverse of a matrix *not inv*
  
  #from C we can easily calculate the quasi equilibrium distribution 
  #and time to absorption using the eigenvectors and eigenvalues of the matrix
  
  #first calculate the eigenvalues and eigenvectors (right and left) and organize these in
  E<-eigen(C.submat) #calculate the right eigenvectors and eigenvalues 
  EL <- eigen(Matrix::t(C.submat)) #calculate the left eigenvectors   
  C1<-E$vectors #F1 is the set of right eigenectors
  C2<-diag(E$values) #F2 is diagonal matrix of the eigenvalues
  C3 <- EL$vectors #F3 is the set of left eigenvectors
  
  #the Quasi-stationary distribution is the left eigenvector corresponding to the eigenvalue 
  #,with maximal real part, of subgenerator matrix C
  rho.1<-max(Re(E$values)) #rho.1 is the eigenvalue with max real part
  rho.2<-max(Re(E$values[(E$values!=rho.1)])) #rho.2 is the next largest eigenvalue
  max.index<-which.max(diag(C2)) #finds the location of rho.1
  QED.dist<-abs(C3[,max.index]) #take the absolute values of the left eigenvector corresponding to rho.1
  QED.dist<-QED.dist/sum(QED.dist) #gives relative frequencies
  
  ##CONVERT EACH STATE "X" TO THE NUMBER OF PATCHES OFF IN THE CONFIGURATION IT REPRESENTS:
  ##1. build a vector with the number of patches off in all possible configurations
  ##ordered according to the sequence they occur as rows in our full CTMC transition matrix
  #for (i in 1:n.patches){
  #  if (i==1){X<-c(1,0)
  #  }else{X<-c(X+1,X)}}
  #X<-X[-1] #removing the number off in the all off state
  ##then using this to obtain the full CTMC's QED for each configuration
  #  QED.dist<-QED.dist[n.patches-X]
  

  #the time to absorption given we begin in a given state is given by summing across rows of
  #the fundamental matrix (a short cut to do this is to multiply the matrix by 
  #a column vector of 1's) and multiplying by the trace of (I-Q)
  c<-sum(-Matrix::diag(C.submat)) #the trace is just the sum of the diagonal elements
  T.absorp<-fundamental.matrix%*%rep(1,n.patches)*c
  
  rate.comparison<-2*rho.1/rho.2 #compares the rate of convergence to QED versus to absorption
  #the QED provides a good estimate of dynamics when this is <1, poor when >1
  
  ratio.of.means.dist<-matrix(rep(NA,n.patches),n.patches,n.patches)
  for (i in 1:n.patches){
    ratio.of.means.dist[i,]<-fundamental.matrix[i,]/sum(fundamental.matrix[i,])}
  #note: the ratio of means distribution is for the reduced CTMC not the full CTMC here
  
  #print(G.mat)
  return(list(QED.dist=QED.dist, T.absorp=T.absorp, rate.comparison=rate.comparison, ratio.of.means.dist=ratio.of.means.dist))
}