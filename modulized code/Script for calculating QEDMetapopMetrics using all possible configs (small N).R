##########################################################################################################
# DESCRIPTION:
##########################################################################################################

# Source this script to calculate the persistence capacity (lambda.M.delta), landscape capacity (lambda.M),
# geometric mean equilibrium occupancy of the landscape by the metapopulation at QED (geom.pstar.QED), 
# arithmetic mean equilibrium occupancy of the landscape by the metapopulation at QED (ar.pstar.QED),
# expected mean number of patches present at QED (exp.n.QED), and a dataframe of the equilibrium 
# occupancies of patches within every habitat configuration possible.

# Because the number of habitat configurations possible becomes huge (size 2^n.patches) for large numbers of 
# patches (n.patches), this script can only be run for landscapes with a small number of patches.
##########################################################################################################
#rm(list=ls()) #to ensure a clean environment

# LOADING REQUIRED PACKAGES
library("Matrix") #to get sparse matrices
library("gdata") #to get upper and lower triangles of matrices
library("lava") #to get anti-diagonal of a matrix using revdiag
library("gtools") #to get even() function
library("tidyverse") #to get filter function
library("doParallel") #to run code in parallel
# SET THE WORKING DIRECTORY TO WHERE THE REQUIRED FUNCTIONS ARE STORED
#setwd("C:/Users/Administrator/Desktop/JAN 2024/modulized code")
# SOURCING THE REQUIRED FUNCTIONS
source("make_Gmat.r")
source("Gmat_to_Csubmat.r")
source("get_QED.r")
source("make_landscape.r")
source("get_distmat.r")
source("get_dispkernel.r")
source("get_lambdaM.r")
source("get_pstar.r")
source("get_MetapopMetrics.r")
source("make_allConfigs.r")
source("get_config_metrics.r")
source("get_QEDMetapopMetrics.r")

#SET UP PARAMETERS FOR DESIRED METAPOPULATION AND DYNAMIC LANDSCAPE
####################################################################
#n.patches<-10 # integer equal to the total number of patches in the system
#landscape.type<-"line" # set to "line" or "plane" indicating whether to have a 1D or 2D patch arrangement
#landscape.size<-100 # numerical value indicating max y (and x coordinates if 2D)
#p.dist<-"uniform" # set to "clustered", "random", "uniform" for how patches should be distributed in the landscape
#p.size.dist<-"uniform" # set to "lognormal", "random" or "uniform" indicating how patch sizes are to be distributed
#max.p.size<-1 # numerical value indicating largest patch size
#clust.its<-1000 # integer equal to the number of times to iterate the algorithm increasing clustering or uniformity for
##                 the distribution of patches in the landscape
#r.on<-0.1 # the rate at which patch gain
#r.off<-0.01 # the rate at which patch loss
#N<-48 # integer number of computational cores available
#e.rate<-0.1 # within patch population extinction rate
#c.rate<-0.2 # colonization rate of patches when a disperse arrives at a patch
#alpha<-0.1 # 1/(the avg. dispersal distance of a species)
#iterations<-1000 # integer value indicating the number of times to perform the interative function for finding pstar
##                  (larger => greater accuracy)
#############################################################################
#MAKING A LANDSCAPE (SET LANDSCAPE PARAMETERS AS DESIRED)
landscape<-make_landscape(n.patches=n.patches, landscape.type=landscape.type, landscape.size=landscape.size, 
                          p.dist=p.dist, p.size.dist=p.size.dist, max.p.size=max.p.size, 
                          clust.its=clust.its)
#MAKING THE GENERATOR MATRIX OF A CTMC WHERE PATCHES ARE LOST AT RATE r.on AND GAINED AT RATE r.off
#AND CALCULATING THE QED OF HABITAT CONFIGURATIONS
G.mat<-make_Gmat(n.patches=n.patches,r.on=r.on,r.off=r.off)
C.submat<-Gmat_to_Csubmat(G.mat)
QED.output<-get_QED(C.submat)
QED<-QED.output[[1]]
QED<-QED/sum(QED) #standardizing to account for rounding errors
#OBTAINING ALL POSSIBLE HABITAT CONFIGURATIONS 
configs<-make_allConfigs(n.patches=n.patches)
configs<-configs[-1,] #remove the all off state
#RUNNING ALL THE WITHIN CONFIG METAPOPULATION METRIC CALCULATIONS IN SERIAL 
#(ALTERNATIVELY COULD BE DONE IN PARALLEL FOR SPEED)
#############################################################################
config.metrics<-matrix(rep(NA, (n.patches+4)*nrow(configs)), nrow=nrow(configs), ncol=(n.patches+4))
for (i in 1:nrow(configs)){
  config.metrics[i,]<-vec(get_config_metrics(config=configs[i,], landscape=landscape, e.rate=0.1, c.rate=0.2, gamma=0, 
                     self.rec=1, alpha=0.1, iterations=1000))
}
#############################################################################
#CALCULATING ESTIMATES OF THE QED METAPOPULATION METRICS USING THE SUBSET OF ALL POSSIBLE HABITAT CONFIGURATIONS
QEDmetapop.metrics<-get_QEDMetapopMetrics(config.metrics, QED)
QEDmetapop.metrics
