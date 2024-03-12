#rm(list=ls()) #clear the workspace

#The first step is to create binary vectors i and j for each configuration
#done by expressing each value from 0 to 2^(N-1) in binary, 
#and a zero in the n'th decimal place of the value [i.e., starting from the right]
#indicates the nth patch is off, while a one indicates the nth patch is on
quasi.eq<-function(n.patches, r.on, r.off){
  #library("R.utils") #needed for creating a binary number
  
  #REPLACING Austin's Algorithm for creating the generator matrix to optimize speed and memory usage and allow us to calculate the QED for larger numbers of patches. Basically, The generator matrix has a Block Toeplitz Structure and therefore can be built iteratively from a series of characteristic Toeplitz matrices. This saves time calculating each possible state transition and the rate at which each transition occurs and instead contructs the generator matrix according to the pattern by which it grows for each patch added to the system.
  
  #there are 4 base transitions that can occur
  #1. a patch that was off can turn on
  #2. a patch that was on can turn off
  #3. a patch that was on can stay on
  #4. a patch that was off can stay off
  #this creates the first base toeplitz matrix 
  if (n.patches==2){
    A<-matrix(0,4,4) # a 4 by 4
    upperTriangle(A, diag=FALSE, byrow=FALSE) <- r.on #where the upper triangle of transitions is to at least one patch being on
    lowerTriangle(A, diag=FALSE, byrow=FALSE) <- r.off #and the lower is to at least patches being off
    revdiag(A)<-0 #both patches can't go off or on at the same time
    diag(A)<--rowSums(A) #since this is a CTMC (therefore transitions are instantaneous) rows must sum to 0 and therefore transitions out of static states happen at a rate equal to - the sum across rows
    A[1,]<-0 #since it's an absorbing matrix the first row corresponding to all patches off has no rates of transition out of this state
    A <- Matrix(A, sparse=TRUE) #setting our matrix to be sparse
  } else { #if there are only 2 patches the generator matrix is just composed of this base toeplitz matrix A but for each additional patch 4 possible transitions x each of the previously possible transitions excluding that single patch are added)
    A<-matrix(0,4,4) #thus, we start with that base toeplitz matrix for two patches
    upperTriangle(A, diag=FALSE, byrow=FALSE) <- r.on
    lowerTriangle(A, diag=FALSE, byrow=FALSE) <- r.off
    revdiag(A)<-0
    A <- Matrix(A, sparse=TRUE)
    #the offset diagonal toeplitz matrices iteratively scaling as patches are added iteratively as follows
    for (n in 2:(n.patches-1)){ #for each patch excluding the added patch
      r.on.diag<-bandSparse(2^(n),2^(n),0,list(rep(r.on, 2^n+1))) #a sparse diagonal matrix of transitions to all other states from that patch being on is added
      r.off.diag<-bandSparse(2^(n),2^(n),0,list(rep(r.off, 2^n+1))) #a sparse diagonal matrix of transitions to all other states from that patch being off is added
      #Make a toeplitz matrix of those 3 toeplitz matrices by binding them appropiately together
      A1<-cbind(A, r.on.diag) 
      A2<-cbind(r.off.diag, A)
      A<-rbind(A1, A2)
    }
    diag(A)<--Matrix::rowSums(A) #ensure again that the rows sum to 0 with the negative sum on the diagonal
    A[1,]<-0 #again, since it's an absorbing matrix the first row corresponding to all patches off has no rates of transition out of this state
  }
  #this provides us with the infinitessimal generator matrix A 
  G.mat<-(A)
  #G.mat
  
  #C.submat is the submatrix corresponding to the transient states (i.e. on to off or off to on 
  #(not staying on or turning off))
  C.submat<-G.mat[2:2^(n.patches),2:2^(n.patches)]
  #sparse_C <- as(C.submat, "sparseMatrix")
  #print(C.submat)
  
  #the negetive inverse of C gives us the fundamental matrix N=(I-Q)^-1, 
  #which gives the ratio of means distribution
  fundamental.matrix<--Matrix::solve(C.submat) #solve gives the inverse of a matrix 
  
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
  
  #the time to absorption given we begin in a given state is given by summing across rows of
  #the fundamental matrix (a short cut to do this is to multiply the matrix by 
  #a column vector of 1's) and multiplying by the trace of (I-Q)
  c<-sum(-Matrix::diag(C.submat)) #the trace is just the sum of the diagonal elements
  T.absorp<-fundamental.matrix%*%rep(1,2^n.patches-1)*c
  
  rate.comparison<-2*rho.1/rho.2 #compares the rate of convergence to QED versus to absorption
  #the QED provides a good estimate of dynamics when this is <1, poor when >1
  
  ratio.of.means.dist<-matrix(rep(NA,(2^n.patches-1)^2),2^n.patches-1,2^n.patches-1)
  for (i in 1:2^n.patches-1){
    ratio.of.means.dist[i,]<-fundamental.matrix[i,]/sum(fundamental.matrix[i,])}
  
  #print(G.mat)
  return(list(QED.dist=QED.dist, T.absorp=T.absorp, rate.comparison=rate.comparison, ratio.of.means.dist=ratio.of.means.dist))
}