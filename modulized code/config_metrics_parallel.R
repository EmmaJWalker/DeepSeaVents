##########################################################################################################
# DESCRIPTION:
##########################################################################################################

# This function calculates the landscape capacity (lambda.M), persistence capacity (lambda.M.delta),
# equilibrium occupancy of each patch (p.star), total equilibrium occupancy across patches in the 
# landscape (metapop.size), and the number of patches on (n.on) and places these in a dataframe
# for some subset of all possible habitat configurations of a dynamic landscape assigned to be 
# performed by a single core, when calculating these metrics over some number of cores (N) in parallel

# INPUTS:
# i: an index indicating which core on which to perform the calculations
# N: the number of cores over which the calculations are being performed in parallel over
# configs: a dataframe where each row is a binary vector indicating witch patches are present (1) 
#          or absent (0) within the landscape, where each column is indexed by each patch.ID
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
# configs.metrics:a dataframe where each row contains:
#                 - lambda.M: the landscape capacity 
#                 - lambda.M.delta: the persistence capacity 
#                 - metapop.size: total equilibrium occupancy 
#                 - n.on: total.number of patches on in the habitat configuration
#                 - p.star: equilibrium occupancy of each patch in each subsequent column
#                 within a given configuration

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
config_metrics_parallel<-function(i, N, configs, landscape, e.rate, c.rate, gamma, 
                                self.rec, alpha, iterations){
  
  #i=i
  #landscape=landscape
  #e.rate=e.rate
  #c.rate=c.rate
  #gamma=gamma
  #self.rec=self.rec
  #alpha=alpha
  #configs=configs
  #N=N
  
  #LOADING CODE SPECIFIC PACKAGES ONTO CORE i:
  library("R.utils") #for int to binary
  library("pracma") #for matrix math
  library("R.utils") #for integer to binary conversion
  library("dplyr") #for various stuff
  library("gdata") #to get upper triangle of a matrix
  library("lava") #to get the reverse diagonal of a matrix
  library("Matrix") #for sparse matrices
  library("reticulate") #for something...
  library("gtools") #for na.replace
  #LOADING REQUIRED FUNCTIONS
  setwd("C:/Users/Administrator/Desktop/JAN 2024/modulized code")
  source("get_distmat.r")
  source("get_dispkernel.r")
  source("get_lambdaM.r")
  source("get_pstar.r")
  source("split_jobs.r")
  source("IDs_to_configs.r")
  source("get_config_metrics.r")
  
  # 1 config == 1 job
  configs.for.i<-configs[split_jobs(J=nrow(configs), N=N, i=i),]
  tot.configs.i<-nrow(configs.for.i)
  n.patches<-ncol(configs.for.i)
  configs.metrics<-data.frame(matrix(data=rep(NA, tot.configs.i*(5+n.patches)), nrow=tot.configs.i, ncol=(5+n.patches)))
  for (k in 1:tot.configs.i) { #for each config assigned to core i...
    # calculate the metapopulation metrics within that config
    config<-configs.for.i[k,]
    config.metrics<-get_config_metrics(config=config, landscape=landscape, e.rate=e.rate, c.rate=c.rate, gamma=gamma, 
                                       self.rec=self.rec, alpha=alpha, iterations=iterations)
    # and add to a dataframe of those metrics 
    configs.metrics[k,]<-config.metrics
  }
  configs.metrics<-na.omit(configs.metrics)
  return(configs.metrics)
}
####################################################################################
# QUICK CHECK: (uncomment and run to check)
####################################################################################
# check by running the _________ script.