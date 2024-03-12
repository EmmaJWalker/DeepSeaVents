#source("dispersal_kernel_func.r")

#CALCULATES lambda.M
#################################################################################################################
#calculates lambda.M for a given landscape (landscape), species specific average dispersal distance (a), 
#and species specific ratio of extinction to colonization rate (delta)

#landscape = a dataframe created by the create landscape function which specifies patch locations, areas, 
#and and an interpatch distance matrix
#a = 1/(average dispersal distance) for a species
#delta = the ratio of the extinction to colonization rate of a species

#NOTE: delta has been pulled out of this calculation
##################################################################################################################
#F2: DETERMINE lambda.M FUNCTION ***With delta (otherwise exactly the same)
##################################################################################################################
get.lambda.M<-function(r.off, landscape, alpha, gamma, epsilon, self.rec, e.rate, c.rate){
  n.patches<-length(landscape$patch.ID)
  x.coord<-landscape$x.coord
  y.coord<-landscape$y.coord
  AREA<-landscape$areas
  delta<-(e.rate+r.off)/c.rate
  
  # construct the distance matrix
  dist.mat<-matrix(rep(0,n.patches*n.patches),n.patches,n.patches)
  for (i in 1:n.patches){
    for (j in 1:n.patches){
      dist.mat[i,j]<-sqrt((x.coord[i]-x.coord[j])^2+(y.coord[i]-y.coord[j])^2)
    }
  }
  #Getting the f(dij) values to plug into the rhs
  DISPKERNEL<-disp.kernel(dist.mat,alpha,gamma,epsilon,self.rec,n.patches)
  
  #STEP 1: SET UP OF MATRIX M
  M<-matrix(rep(NA, n.patches*n.patches), n.patches,n.patches)
  for (i in 1:n.patches){ #for each patch
    for (j in 1:n.patches){ #to each other patch (including itself)
          M[i,j]<-AREA[j]*DISPKERNEL[i,j]
    }
  }

  #FIND lambda M (The Metapop persistence capacity)
  eigen.M<-eigen(M) #determine eigenvalues and eigenvectors for matrix M
  lambda.M<-eigen.M$values[1] 
  #lambda.M is the leading eigenvalue of matrix M which is given as the first output value
  lambda.M<-lambda.M#/delta #<-we pull delta out to be able to scale Lm by it to compare simulations across landscapes
  return(lambda.M)
}
#END lambda.M FUNCTION
#################################################################################################################################
#CHECK
#lambda.M<-get.lambda.M(landscape=landscape, alpha=0.1, e.rate=1, c.rate=1, self.rec=1, gamma=0, epsilon=n.patches)
#lambda.M
#20/lambda.M
#delta<-20/(lambda.M=get.lambda.M(landscape=landscape, alpha=0.1, e.rate=0.1, c.rate=4, self.rec=1, gamma=0, epsilon=n.patches))
#delta
#e.rate<-1/delta
#c.rate<-1
#lambda.M<-get.lambda.M(landscape=landscape, alpha=0.1, e.rate=1/delta, c.rate=1, self.rec=1, gamma=0, epsilon=n.patches)
#lambda.M*delta
