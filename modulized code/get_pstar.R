##########################################################################################################
# DESCRIPTION:
##########################################################################################################

# This function calculates the equilibrium occupancy (probabilities with which patches occupied at 
# equilibrium) within a given configuration of habitat

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
# iterations: number of iterations for which the iterative function f for finding pstar is iterated 
#             -greater iterations = greater accuracy, less iterations = lower accuracy

# OUTPUTS:
# p.star: the expected equilibrium occupancy of each patch of a metapopulation in the input landscape

# REQUIRED PACKAGES:
# none

#NOTE: delta can be factored out of this calculation by letting e.rate/c.rate = 1
####################################################################################
# FUNCTION CODE:
####################################################################################
get_pstar<-function(landscape, e.rate, c.rate, disp.kernel, iterations){
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
  
  #STEP 2: ITERATING FUNCTION TO FIND P*
  p<-rep(NA, n.patches*iterations);dim(p)<-c(iterations,n.patches)
  p[1,1:n.patches]<-rep(0.1,n.patches)
  for (t in 1:(iterations-1)){ #iterate the following for a given number of iterations
    P<-p[t,] #use p of previous iteration
    g<-(M%*%P) #SRLM g function from Ovaskainen and Hanski, 2001
    f<-(g/(g+delta)) #SRLM f function from Ovaskainen and Hanski, 2001, 
    #which when iterated many times will give p*
    p[(t+1),]<-t(f)} 
  #set p for the next iteration = to the output value of p for this iteration
  p.star<-p[iterations,] #p* = the final iterations value of p after many iterations
  
  return(p.star)
}

####################################################################################
# QUICK CHECK: (uncomment and run to check)
####################################################################################
#rm(list=ls()) #to ensure a clean environment
## now run the FUNCTION CODE!!!
#library("ggplot2")
#setwd("C:/Users/Administrator/Desktop/JAN 2024/modulized code")
#source("make_landscape.r")
#source("get_distmat.r")
#source("get_dispkernel.r")
#landscape<-make_landscape(n.patches=4, landscape.type="plane", landscape.size=20, 
#                 p.dist="clustered", p.size.dist="random", max.p.size=2, 
#                 clust.its=10)
#head(landscape)
#print(ggplot(landscape, aes(x.coord,y.coord)) + theme_classic() + geom_point(aes(size = p.sizes)))
#dist.mat<-get_distmat(landscape)
#disp.kernel<-get_dispkernel(dist.mat, alpha=0.1, gamma=0.5, self.rec=1)
#p.star<-get_pstar(landscape, e.rate=0.1, c.rate=0.2, disp.kernel, iterations=10)
#p.star
####################################################################################