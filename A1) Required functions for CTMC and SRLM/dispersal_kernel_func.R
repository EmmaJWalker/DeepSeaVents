disp.kernel<-function(mat, alpha, gamma, epsilon, self.rec, n.patches){
  k<-matrix(rep(0,n.patches*n.patches),n.patches,n.patches)
  
  for (i in 1:n.patches){ #for each patch
    for(j in 1:n.patches){ #to each other patch
      if((i<j) && (abs(i-j)<=epsilon)){ #if patch j is "upstream" of patch i and reachable
        #by dispersers from patch i
        k[i,j]<-(1-gamma)*exp(-alpha*mat[i,j])} #dispersal decreases according to the strength
      #to which upstream dispersal is limited, and further decreases exponentially with 
      #1/(avg.dispersal distance) (alpha) and the distance between the patches
      else if ((i>j) && (abs(i-j)<=epsilon)){ #if patch i is "downstream" of patch i and reachable
        #by dispersers from patch j
        k[i,j]<-exp(-alpha*mat[i,j])} #dispersal decreases exponentially with 
      #1/(avg. dispersal distance) (alpha) and the distance between the patches
      else if ((i==j) && (self.rec==1)){ #otherwise if the other patch is itself and self recruitment
        #is allowed
        k[i,j]<-1} #dispersal to itself occurs
      else{ #otherwise,
        k[i,j]<-0 #no dispersal between the two patches or itself occurs
      }
    }
  }
  return(k) #provide the appropriate dispersal kernel
}