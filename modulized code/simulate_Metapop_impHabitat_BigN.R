##########################################################################################################
# DESCRIPTION:
##########################################################################################################

# This function simulates the dynamics of a metapopulation in a landscape of impermanent habitat patches
# gained and lost according to a CTMC describing how the landscape's configurations of available habitat 
# patches change over time.

# INPUTS:
# total.trans: the total number of transitions between habitat configuration overwhich the dynamics of the
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
# deSolve

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

####################################################################################
# FUNCTION CODE:
####################################################################################
simulate_Metapop_impHabitat_BigN<-function(limit, landscape, e.rate, c.rate, alpha, gamma, 
                                  self.rec, r.on, r.off, QED){
  
  #1) SETTING UP ALL OUR PARAMETER VALUES FOR THE SRLM's ODEs
  #################################################################################################
  n.patches<-length(landscape$patch.ID)
  p.sizes<-landscape$p.sizes
  x.coord<-landscape$x.coord
  y.coord<-landscape$y.coord
  ##calculate the extinction rate for each patch (if areas are all 1 the extinction rate is patch independent)
  extinction.rates<-e.rate/p.sizes
  delta<-e.rate/c.rate
  dist.mat<-get_distmat(landscape=landscape)
  disp.kernel<-get_dispkernel(dist.mat=dist.mat, alpha=alpha, gamma=gamma, self.rec=self.rec)  
  #################################################################################################
  
  #2) SET THE INITIAL CONDITIONS FOR THE SYSTEM:
  #################################################################################################
      # starting config chosen w pr(config @ QED)
      config<-vec(subsample_reducedQED(QED=QED, n.patches=n.patches, s.size=1))
      # starting at the equilibrium occupancy of that habitat configuration ########CHECK!!!!!!!!!!
      ICs<-config
      ICs[config!=0]<-get_pstar(landscape=landscape[config!=0,], e.rate=e.rate, c.rate=c.rate, 
                                disp.kernel=disp.kernel, iterations=1000)
      
      t<-0
      tau.times<-0 #no initial time until transition
      trans<-0 #no transitions have occured initially
      sim.data<-rep(NA, n.patches+1) #########################CHECK!!!!!
      configs.data<-rep(NA, n.patches+1) #######################CHECK!!!!!
      
      #3.3 A1) CREATE A NAMED LIST OF PARAMETERS FOR THE SRLM ODEs
      #note: while obviously the values of some of these parameters will change with configuration
      #changes, the names of these parameters don't so we just set these up with the initial values 
      #just to start *********************************
      ############################################
      parameter.values<-c(n.patches,c.rate,self.rec,extinction.rates,x.coord,y.coord,p.sizes,config
                          ,as.vector(disp.kernel))
      parameters<-name_SRLMODEparams(parameter.values)
      ###########################################
      
     
  # 3) PROVIDE A STOPPING CONDITION TO KEEP SIMULATIONS FROM TAKING TO LONG
  #################################################################################################
      #transition.limit<-limit #setting a limit to the number of transitions we sampl
      time.limit<-limit
  #################################################################################################
  
  #4) BEGIN THE SIMULATION
  #################################################################################################
  #while (sum(config) > 0 & trans < transition.limit) { 
  while (sum(config) > 0 & t < time.limit) {
    #while at least one patch remains on and we are within the transition limit
    
    # 4.1) GET HOW LONG WILL BE SPENT IN THE CURRENT HABITAT CONFIGURATION
    tau<-get_tau(config=config, r.on=r.on, r.off=r.off)
    
    # 4.2) SET UP AND RUN THE SRLM AND RECORD THE HABITAT AND METAPOPULATION DYNAMICS BETWEEN t AND tau
    parameter.values<-c(n.patches,c.rate,self.rec,extinction.rates,x.coord,y.coord,p.sizes,config
                        ,as.vector(disp.kernel))
    output<-setNrun_SRLMODE(parameter.values=parameter.values, f.names=f.names, 
                            ICs=ICs, t=t, tau=tau, 
                    sim.data=sim.data, configs.data=configs.data)
    sim.data<-output[[1]]
    configs.data<-output[[2]]
    # NOTE: probably what I am going to need to do to prevent the machine running out of memory is to 
    # save the metapopulation dynamics over each interval tau as separate files
    
    # 4.3) STOP THE SIMULATION IF THE METAPOPULATION WENT EXTINCT
    if (sum(sim.data[nrow(sim.data),-1]) <= 0){ break }
    
    # 4.4) UPDATE THE TIME ELAPSED IN SIMULATION
    t<-t+tau
    
    # 4.5) MAKE THE TRANSITION & AND UPDATE INITIAL CONDITIONS WITH THE CURRENT OCCUPANCIES OF 
    # PATCHES THAT REMAINED ON:
    config<-transition(config=config, r.on=r.on, r.off=r.off)
    ICs<-as.numeric(as.character(sim.data[length(sim.data[,1]),-1])) #to be the end-values from the SRLM
    ICs[config==0]<-0 #but zero if a patch is now off
    
    # 4.6) RECORD THE TIME OF EACH TRANSITION & HOW MANY TRANSITIONS HAVE OCCURED
    # THIS WILL STOP THE SIMULATION IF THE TRANSITION LIMIT HAS BEEN EXCEEDED
    if (sum(tau.times) == 0){
      tau.times<-tau
    } else {
      tau.times<-c(tau.times, tau)
    }
    trans<-trans + 1 #since we simulate on intervals of transition this 
    #increments by 1 each time we run through the algorithm

  }
  end.t<-t #time to absorption or simulation end
  return(list(sim.data=sim.data, configs.data=configs.data))
}
####################################################################################
# QUICK CHECK: (uncomment and run to check)
####################################################################################
rm(list=ls()) #to ensure a clean environment
#library("pryr") #to get memory used function
#mem_used()
## now run the FUNCTION CODE!!!
library("Matrix") #to get sparse matrices
library("gdata") #to get upper and lower triangles of matrices
library("lava") #to get anti-diagonal of a matrix using revdiag
library("deSolve")
setwd("C:/Users/Administrator/Desktop/JAN 2024/modulized code")
source("make_landscape.r")
source("get_distmat.r")
source("get_dispkernel.r")
source("get_lambdaM.r")
source("get_pstar.r")
source("transition.r")
source("setNrun_SRLMODE.r")
source("name_SRLMODEparams.r")
source("SRLMODE.r")
source("at_equilibrium_rootfunc.r")
source("get_distmat.r")
source("make_reduced_Gmat.r")
source("Gmat_to_Csubmat.r")
source("get_QED.r")
source("subsample_reducedQED.r")
source("reduced_to_full_QED.r")
source("IDs_to_configs.r")
source("get_tau.r")
n.patches<-50
r.on<-0.1
r.off<-0.1
G.mat<-make_reduced_Gmat(n.patches=n.patches,r.on=r.on,r.off=r.off)
C.submat<-Gmat_to_Csubmat(G.mat)
output<-get_QED(C.submat)
QED<-output[[1]]
landscape<-make_landscape(n.patches=n.patches, landscape.type="line", landscape.size=100, 
                 p.dist="uniform", p.size.dist="uniform", max.p.size=1, 
                 clust.its=0)
Sys.time()
output<-simulate_Metapop_impHabitat_BigN(limit=1000, landscape=landscape, e.rate=0.1, c.rate=0.2, alpha=0.1, gamma=0, 
                                           self.rec=1, r.on=r.on, r.off=r.off, QED=QED)
Sys.time()
tail(output$sim.data)
tail(output$configs.data)
# takes 7 min 50 patches, r.on=0.1, r.off=0.01
nrow(output$configs)



  
  
  
  
  
  
  
  
  
  
  
  
  
