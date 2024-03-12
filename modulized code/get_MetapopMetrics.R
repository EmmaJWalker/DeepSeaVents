# INPUTS:
# landscape: a dataframe containing patch.ID, p.sizes, coordinates, distance.matrix
#            -patch.ID is a unique integer ID number for each patch
#            -p.sizes provides the size of each patch
#            -y.coord provides vertical (e.g. latitudinal) patch locations on a "map"
#            -x.coord provides horizontal (e.g. longitudinal) patch locations on a "map" (all 0 if 1D landscape)
#            -subsequent columns contain interpatch distances forming a distance (or adjacency) matrix
# e.rate: positive numeric value indicating the within patch population extinction rate
# c.rate: positive numeric value indicating colonization rate when dispersers arrive at a patch
# alpha: a positive numeric value indicating 1/(avg. dispersal distance of a species)
# gamma: a positive numeric value <1 indicating the strength to which upstream (where upstream are those with j>i
#        as indexed by patchID when "landscape" is ordered by patchID) dispersal may be limited (0 no limitation)
# self.rec: a positive numeric value indicating the extent to which disperses from a given patch may resettle on 
#        the same patch (0 makes resettlement impossible)
# iterations: number of iterations for which the iterative function f for finding pstar is iterated 
#             -greater iterations = greater accuracy, less iterations = lower accuracy


# OUTPUTS:
# a list of objects containing:
#                 - lambda.M.delta: the persistence capacity of a metapopulation in a landscape
#                 - lambda.M: the landscape capacity 
#                 - pstar: a vector containing the equilibrium occupancy of each patch within a landscape
#                 - metapop.size: the total equilibrium occupancy of a landscape (aka how many patches occupied in it)

# REQUIRED PACKAGES:
# none 

# REQUIRED FUNCTIONS:
# get_lambdaM
# get_pstar
# get_distmat
# get_dispkernel
##########################################################################################################
get_MetapopMetrics<-function(landscape, e.rate, c.rate, alpha, gamma, self.rec, iterations){
  delta<-e.rate/c.rate
  dist.mat<-get_distmat(landscape=landscape)
  disp.kernel<-get_dispkernel(dist.mat=dist.mat, alpha=alpha, gamma=gamma, self.rec=self.rec)
  lambda.M<-get_lambdaM(landscape=landscape, e.rate=e.rate, c.rate=c.rate, disp.kernel=disp.kernel)
  lambda.M.delta<-get_lambdaM(landscape=landscape, e.rate=e.rate, c.rate=c.rate, disp.kernel=disp.kernel)/delta
  if (lambda.M.delta < delta){
    pstar<-0
  } else {
    p.star<-get_pstar(landscape=landscape, e.rate=e.rate, c.rate=c.rate, disp.kernel=disp.kernel, iterations=iterations)
  }
  metapop.size<-sum(p.star)
  
  return(list(lambda.M.delta, lambda.M, metapop.size,  p.star))
}
####################################################################################
# QUICK CHECK: (uncomment and run to check)
####################################################################################
#rm(list=ls()) #to ensure a clean environment
## now run the FUNCTION CODE!!!
#setwd("C:/Users/Administrator/Desktop/JAN 2024/modulized code")
#source("make_landscape.r")
#source("get_distmat.r")
#source("get_dispkernel.r")
#source("get_lambdaM.r")
#source("get_pstar.r")
#landscape<-make_landscape(n.patches=4, landscape.type="plane", landscape.size=20, 
#                 p.dist="clustered", p.size.dist="random", max.p.size=2, 
#                 clust.its=10)
#metapop.metrics<-get_MetapopMetrics(landscape=landscape, e.rate=0.1, c.rate=0.2,
#                                    alpha=0.1, gamma=0, self.rec=1, iterations=1000)
#metapop.metrics
####################################################################################
