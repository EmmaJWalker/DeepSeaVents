##########################################################################################################
# DESCRIPTION:
##########################################################################################################

# This function calculates the landscape or metapopulation capacity (lambda.M, whether e.rate/c.rate=1 or 
# not respectively) within a given configuration of habitat

# INPUTS:
# landscape: a dataframe containing patch.ID, p.sizes, coordinates, distance.matrix
#            -patch.ID is a unique integer ID number for each patch
#            -p.sizes provides the size of each patch
#            -y.coord provides vertical (e.g. latitudinal) patch locations on a "map"
#            -x.coord provides horizontal (e.g. longitudinal) patch locations on a "map" (all 0 if 1D landscape)
#            -subsequent columns contain interpatch distances forming a distance (or adjacency) matrix
# e.rate: positive numeric value indicating the within patch population extinction rate
# c.rate: positive numeric value indicating colonization rate when dispersers arrive at a patch
# disp.kernel: a dataframe of values weighted relative to a species ability to move from one patch 
#              location to another (i.e. effective adjacency of patches)

# OUTPUTS:
# lambda.M: the landscape or metapopulation capacity (lambda.M, whether e.rate/c.rate=1 or 
# not respectively) in the input landscape

# REQUIRED PACKAGES:
# none

#NOTE: delta can be factored out of this calculation by letting e.rate/c.rate = 1
####################################################################################
# FUNCTION CODE:
####################################################################################
get_lambdaM<-function(landscape, e.rate, c.rate, disp.kernel){
  n.patches<-length(landscape$patch.ID)
  delta<-e.rate/c.rate
  SIZE<-landscape$p.sizes
  DISPKERNEL<-disp.kernel
  
  #STEP 1: SET UP OF MATRIX M
  M<-matrix(rep(NA, n.patches*n.patches), n.patches,n.patches)
  for (i in 1:n.patches){ #for each patch
    for (j in 1:n.patches){ #to each other patch (including itself)
      M[i,j]<-SIZE[j]*DISPKERNEL[i,j]
    }
  }
  
  #STEP 2: Calculate lambda M 
  eigen.M<-eigen(M) #determine eigenvalues and eigenvectors for matrix M
  lambda.M<-eigen.M$values[1] 
  #lambda.M is the leading eigenvalue of matrix M which is given as the first output value
  lambda.M<-lambda.M#/delta #<-we pull delta out to be able to scale Lm by it to compare simulations across landscapes
  return(lambda.M)
}

####################################################################################
# QUICK CHECK: (uncomment and run to check)
####################################################################################
#rm(list=ls()) #to ensure a clean environment
## now run the FUNCTION CODE!!!
#library("ggplot2")
#landscape<-make_landscape(n.patches=4, landscape.type="plane", landscape.size=20, 
#                 p.dist="clustered", p.size.dist="random", max.p.size=2, 
#                 clust.its=10)
#head(landscape)
#print(ggplot(landscape, aes(x.coord,y.coord)) + theme_classic() + geom_point(aes(size = p.sizes)))
#dist.mat<-get_distmat(landscape)
#disp.kernel<-get_dispkernel(dist.mat, alpha=0.1, gamma=0.5, self.rec=1)
#lambda.M<-get_lambdaM(landscape, e.rate=0.1, c.rate=0.2, disp.kernel)
#lambda.M
####################################################################################





