#NOTE: same as QED calculator
#This function takes as its arguments the number of patches (N), 
#and the rates of patches turning from off to on (r.on) and on to off (r.off).
#It constructs the transition matrix (A) that describes the transition rates from each 
#configuration of on/off to each other.
#It then takes the submatrix of transient states (A_submat) and calculates the quasi-equilibrium
#distribution (dist) and expected time to absorption into the all the all-zero state (t_absorp)
#using the dominant eigenvector/eigenvalue pair from A_submat

#The first step is to create binary vectors i and j for each configuration
#done by expressing each value from 0 to 2^(N-1) in binary, 
#and a zero in the n'th decimal place of the value [i.e., starting from the right]
#indicates the nth patch is off, while a one indicates the nth patch is on
quasi.eq<-function(n.patches, r.on, r.off){
  library("R.utils") #needed for creating a binary number
  A<-matrix(rep(0,2^n.patches*2^n.patches),2^n.patches,2^n.patches)
  
  for (i in 0:(2^n.patches-1)){ # for each posible state transition from off to on
    for (j in 0:(2^n.patches-1)){ #and for each possible state transition from on to off
      bin.val.i<-intToBin(i) #obtain the binary representation of i 
      bin.val.i<-strsplit(bin.val.i,"")
      bin.val.j<-intToBin(j) #obtain the binary representation of j
      bin.val.j<-strsplit(bin.val.j,"")
      states.i<-rep(0,n.patches) #create a vector to hold all possible future states of each patch
      states.j<-rep(0,n.patches) #create a vector to hold all possible previous states of each patch
      for (k in 1:length(bin.val.i[[1]])){ #for each binary digit in i
        states.i[length(bin.val.i[[1]])-k+1]<-as.numeric(as.character(bin.val.i[[1]][k]))}
      #set the possible future states for each patch to be every possible binary state (on or off) 
      for (k in 1:length(bin.val.j[[1]])){ #for each binary digit in j
        states.j[length(bin.val.j[[1]])-k+1]<-as.numeric(as.character(bin.val.j[[1]][k]))}
      #set the possible previous states for each patch to be every possible binary state (on or off)
      
      rate.vec.i<-rep(0,n.patches) #create a vector to hold the trans rates of patches turning on/off
      if (sum(states.i!=0)){ #provided we're not in the absorbing state when all patches are off
        #(in which case we never move out of it at any rate so rate.vec remains 0)
        for (q in 1:n.patches){ #for each patch
          if (states.i[q]==0){rate.vec.i[q]<-r.on} #if it is off, set it to turn on at rate r.on
          else{rate.vec.i[q]<-r.off} #otherwise set it to turn off at rate r.off
        }
      }
      # Hamming distance between the two vectors representing the on/off
      #determines if patches are set to turn on or off, or both (/neither) within a timestep
      dist.H<-sum(abs(states.i-states.j)) #using the difference (distance) between the two
      #if the state vectors are the same (i.e. a patch that's off will remain off and visa versa), 
      if (dist.H == 0) {A[i+1,j+1]<--sum(rate.vec.i)} #the transition rate is the negetive sum of 
      #all the other transition rates
      #where as if the state vectors differ by one spot, (i.e. a patch that's on will turn off and visa versa)
      if (dist.H == 1){A[i+1,j+1]<-rate.vec.i[states.i!=states.j]}#then the transition rate is given 
      #by whatever the differing spot is (i.e., from on to off or vice versa)
      #where as if the vectors are different in two or more spots, (i.e. a patch that's on will go both on and off)
      if (dist.H > 1){A[i+1,j+1]<-0} #then the infinitesimal transition rate is 0 
    }
  }
  #this provides us with the infinitessimal generator matrix A 
  G.mat<-(A)
  
  #C.submat is the submatrix corresponding to the transient states (i.e. on to off or off to on 
  #(not staying on or turning off))
  C.submat<-G.mat[2:2^(n.patches),2:2^(n.patches)]
  #the negetive inverse of C gives us the fundamental matrix N=(I-Q)^-1, 
  #which gives the ratio of means distribution
  fundamental.matrix<--inv(C.submat)
  
  #from C we can easily calculate the quasi equilibrium distribution 
  #and time to absorption using the eigenvectors and eigenvalues of the matrix
  
  #first calculate the eigenvalues and eigenvectors (right and left) and organize these in
  E<-eigen(C.submat) #calculate the right eigenvectors and eigenvalues 
  EL <- eigen(t(C.submat)) #calculate the left eigenvectors   
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
  c<-sum(-diag(C.submat)) #the trace is just the sum of the diagonal elements
  T.absorp<-fundamental.matrix%*%rep(1,2^n.patches-1)*c
  
  rate.comparison<-2*rho.1/rho.2 #compares the rate of convergence to QED versus to absorption
  #the QED provides a good estimate of dynamics when this is <1, poor when >1
  
  ratio.of.means.dist<-matrix(rep(NA,(2^n.patches-1)^2),2^n.patches-1,2^n.patches-1)
  for (i in 1:2^n.patches-1){
    ratio.of.means.dist[i,]<-fundamental.matrix[i,]/sum(fundamental.matrix[i,])}
  
  print(G.mat)
  return(list(QED.dist, T.absorp, rate.comparison, ratio.of.means.dist))
}
#TEST 
#quasi.eq(n.patches=2, r.on=100, r.off=1)

