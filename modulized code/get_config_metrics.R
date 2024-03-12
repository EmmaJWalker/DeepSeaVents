##########################################################################################################
# DESCRIPTION:
##########################################################################################################

# This function calculates the landscape capacity (lambda.M), persistence capacity (lambda.M.delta),
# equilibrium occupancy of each patch (p.star), total equilibrium occupancy across patches in the 
# landscape (metapop.size), and the number of patches on (n.on) within a single habitat configuration 
# in a dynamic landscape

# INPUTS:
# config: a binary vector indicating with patches are present (1) or absent (0) within the landscape 
# landscape: a dataframe containing patch.ID, p.sizes, coordinates, distance.matrix
#            -patch.ID is a unique integer ID number for each patch
#            -p.sizes provides the size of each patch
#            -y.coord provides vertical (e.g. latitudinal) patch locations on a "map"
#            -x.coord provides horizontal (e.g. longitudinal) patch locations on a "map" (all 0 if 1D landscape)
#            -subsequent columns contain interpatch distances forming a distance (or adjacency) matrix
# e.rate: positive numeric value indicating the within patch population extinction rate
# c.rate: positive numeric value indicating colonization rate when dispersers arrive at a patch
# gamma: a positive numeric value <1 indicating the strength to which upstream (where upstream are those with j>i
#        as indexed by patchID when "landscape" is ordered by patchID) dispersal may be limited (0 no limitation)
# self.rec: a positive numeric value indicating the extent to which disperses from a given patch may resettle on 
#        the same patch (0 makes resettlement impossible)
# iterations: number of iterations for which the iterative function f for finding pstar is iterated 
#             -greater iterations = greater accuracy, less iterations = lower accuracy

# OUTPUTS:
# config.metrics: a single row of a dataframe containing:
#                 - lambda.M: the landscape capacity 
#                 - lambda.M.delta: the persistence capacity 
#                 - metapop.size: total equilibrium occupancy 
#                 - n.on: total.number of patches on in the habitat configuration
#                 - p.star: equilibrium occupancy of each patch in each subsequent column

# REQUIRED PACKAGES:
# none

# REQUIRED FUNCTIONS:
# get_distmat
# get_dispkernel
# get_lambdaM
# get_pstar

#NOTE: delta can be factored out of this calculation by letting e.rate/c.rate = 1
####################################################################################
# FUNCTION CODE:
####################################################################################
get_config_metrics<-function(config, landscape, e.rate, c.rate, gamma, 
                                        self.rec, alpha, iterations){
  n.patches<-length(config)
  delta<-e.rate/c.rate
  landscape.config<-landscape[config!=0,] #for a landscape comprised of only the on patches
  dist.mat<-get_distmat(landscape)
  disp.kernel<-get_dispkernel(dist.mat, alpha=alpha, gamma=gamma, self.rec=self.rec)
  lambda.M<-get_lambdaM(landscape=landscape.config, e.rate=e.rate, c.rate=c.rate, disp.kernel=disp.kernel)
  lambda.M.delta<-get_lambdaM(landscape=landscape.config, e.rate=e.rate, c.rate=c.rate, disp.kernel=disp.kernel)/delta
  # now can save on calculation since
  # if below the persistence threshold, the metapop is expected extinct at equilibrium
  if (lambda.M.delta < delta){
    pstar.configs<-rep(0, n.patches) # so p.star is just 0 in all patches
  } else { #otherwise
    #calculates the equilibrium patch occupancy of each on patch in config k
    pstar.configs<-rep(NA, n.patches)
    pstar.configs[config!=0]<-get_pstar(landscape=landscape.config, e.rate=e.rate, c.rate=c.rate, disp.kernel=disp.kernel, iterations=iterations)
    pstar.configs<-na.replace(pstar.configs, 0)    
    }
  metapop.size<-sum(pstar.configs) #calculate total metapop size expected at eq
  n.on<-sum(config)
  config.metrics<-data.frame(lambda.M.delta, lambda.M, metapop.size, n.on, t(pstar.configs))
  config.metrics<-na.omit(config.metrics) 
  #rename(config.metrics, lambda.M=X1, lambda.M.delta=X2, metapop.size=X3, n.on=X4)
  return(config.metrics)
}
####################################################################################
# QUICK CHECK: (uncomment and run to check)
####################################################################################
#rm(list=ls()) #to ensure a clean environment
## now run the FUNCTION CODE!!!
#source("make_landscape.r")
#source("get_distmat.r")
#source("get_dispkernel.r")
#source("get_lambdaM.r")
#source("get_pstar.r")
#config<-c(1,1,0,1,0)
#landscape<-make_landscape(n.patches=5, landscape.type="plane", landscape.size=20, 
#                 p.dist="clustered", p.size.dist="random", max.p.size=1, 
#                 clust.its=10)
#config.metrics<-get_config_metrics(config, landscape, c.rate=0.2, e.rate=0.1, gamma=0, 
#                   self.rec=1, alpha=0.1, iterations=1000)
#config.metrics
####################################################################################
