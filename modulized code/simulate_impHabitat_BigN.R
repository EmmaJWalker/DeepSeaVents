##########################################################################################################
# DESCRIPTION:
##########################################################################################################

# This function simulates the dynamics of a metapopulation in a landscape of impermanent habitat patches
# gained and lost according to a CTMC describing how the landscape's configurations of available habitat 
# patches change over time.

# INPUTS:
# limit: the total number of transitions between habitat configuration overwhich the dynamics of the
# limit.type:
# system will be simulated
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
# r.on: rate of patch gain
# r.off: rate of patch loss
# QED: the quasi-equilibrium distribution of states where a given number of habitat configurations are present in the
# landscape

# OUTPUTS:
# sim.data:a dataframe providing the occupancy of each patch in the landscape over the increments of time on which the
#          SRLM's system of ODE's is solved
# configs.data: a dataframe providing the configuration of the landscape (with which patches are on or off indicated by
#               1's and 0's respectively) at each time point in which that configuration occurred

# REQUIRED PACKAGES:
# deSolve: for ode solver
# gtools: for even function

# REQUIRED FUNCTIONS:
#transition
#setNrun_SRLMODE
#transition
#SRLMODE
#at_equilibrium_rootfunc
#get_distmat
#get_dispkernel
#get_pstar
#name_SRLMODEparams
#subsample_reducedQED
#reduced_to_full_QED
#IDs_to_configs: IDsToConfigs
#get_tau
#get_Mij

################################################################################
# FUNCTION CODE:
################################################################################
simulate_impHabitat_BigN<-function(limit, limit.type="time", n.patches, r.on, r.off, QED){
  #*****************************************************************************
  # SET THE INITIAL CONDITIONS FOR THE SYSTEM:
  #*****************************************************************************
  # starting config chosen w pr(config @ QED)
  config<-vec(subsample_reducedQED(QED=QED, n.patches=n.patches, s.size=1))
  t<-0 #start at time 0
  tau.times<-0 #no initial time until transition
  trans<-0 #no transitions have occured initially
  configs.data<-rep(NA, n.patches+1)
  #*****************************************************************************
  # PROVIDE THE STOPPING CONDITION TO KEEP SIMULATIONS FROM TAKING TOO LONG
  #*****************************************************************************
  if (limit.type=="trans"){
    dont.stop<-(sum(config) > 0 & trans < limit)
  } else if (limit.type=="time"){
    dont.stop<-(sum(config) > 0 & t < limit)
  }
  #*****************************************************************************
  # BEGIN THE SIMULATION:
  #*****************************************************************************
  while (dont.stop==T) {
    #while at least one patch remains on and the stopping condition is not met:
    # 1) RECORD THE CONFIGURATION THE HABITAT IS IN AT THIS TIME
    configs.data<-rbind(configs.data, c(t, config))
    # 2) DETERIME WILL BE SPENT IN IT:
    tau<-get_tau(config=config, r.on=r.on, r.off=r.off)
    # 3) UPDATE THE TIME ELAPSED IN SIMULATION
    t<-t+tau
    # 4) TRANSITION INTO A NEW CONFIGURATION:
    config<-transition(config=config, r.on=r.on, r.off=r.off)
    # 5) UPDATE A RECORD OF HOW MANY TRANSITIONS HAVE OCCURED
    trans<-trans + 1 #since we simulate on intervals of transition this 
    #increments by 1 each time we run through the algorithm
  }
  #*****************************************************************************
  return(configs.data=configs.data)
}
################################################################################
# QUICK CHECK: (uncomment and run to check)
################################################################################
#rm(list=ls()) #to ensure a clean environment
#library("pryr") #to get memory used function
#mem_used()
## now run the FUNCTION CODE!!!
#library("Matrix") #to get sparse matrices
#library("gdata") #to get upper and lower triangles of matrices
#library("lava") #to get anti-diagonal of a matrix using revdiag
#library("deSolve")
#setwd("/Users/abuga/Dropbox/Mac/Desktop/git/DeepSeaVents/modulized code")
#source("make_landscape.r")
#source("get_distmat.r")
#source("get_dispkernel.r")
#source("get_lambdaM.r")
#source("get_pstar.r")
#source("transition.r")
#source("name_SRLMODEparams.r")
#source("SRLMODE.r")
#source("at_equilibrium_rootfunc.r")
#source("get_distmat.r")
#source("make_reduced_Gmat.r")
#source("Gmat_to_Csubmat.r")
#source("get_QED.r")
#source("subsample_reducedQED.r")
#source("reduced_to_full_QED.r")
#source("IDs_to_configs.r")
#source("get_tau.r")
#source("name_Mij.r")
#n.patches<-50
#r.on<-0.1
#r.off<-0.1
#G.mat<-make_reduced_Gmat(n.patches=n.patches,r.on=r.on,r.off=r.off)
#C.submat<-Gmat_to_Csubmat(G.mat)
#output<-get_QED(C.submat)
#QED<-output[[1]]
#landscape<-make_landscape(n.patches=n.patches, landscape.type="line", 
#                          landscape.size=100, p.dist="uniform", 
#                          p.size.dist="uniform", max.p.size=1, clust.its=0)
#Sys.time()
#output<-simulate_Metapop_impHabitat_BigN(limit=1000, limit.type="time", landscape=landscape, 
#                                         e.rate=0.1, c.rate=0.2, alpha=0.1, 
#                                         gamma=0, self.rec=1, r.on=r.on, 
#                                         r.off=r.off, QED=QED)
#Sys.time()
#tail(output$sim.data)
#tail(output$configs.data)
## takes 7 min 50 patches, r.on=0.1, r.off=0.01
#nrow(output$configs)
















