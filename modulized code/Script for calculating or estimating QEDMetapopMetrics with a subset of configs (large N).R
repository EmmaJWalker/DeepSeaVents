##########################################################################################################
# DESCRIPTION:
##########################################################################################################

# Source this script to calculate the persistence capacity (lambda.M.delta), landscape capacity (lambda.M),
# geometric mean equilibrium occupancy of the landscape by the metapopulation at QED (geom.pstar.QED), 
# arithmetic mean equilibrium occupancy of the landscape by the metapopulation at QED (ar.pstar.QED),
# expected mean number of patches present at QED (exp.n.QED), and a dataframe of the equilibrium 
# occupancies of patches within each configuration used in the calculation of these metrics. It does so for
# a metapopulation (as specified by the metapop parameters *see below) in a dynamic landscape of impermanent 
# patches with a given QED of habitat configurations but only for a given set of those habitat configurations
# (with how many to use specified by the parameter s.size). This set may include up to all possible 
# 2^n.patches configurations (when computationally feasible) or a subset of these configurations chosen by 
# their probability of occurrence at QED to provide estimates of these metrics in landscapes containing larger 
# numbers of patches (when using all 2^n.patches configurations is computationally infeasible). 

# To increase computational efficiency within habitat configuration calculations are performed in parallel 
# over N cores (where N is the number of cores available). However, computation of these metrics may still 
# be limited by the amount of memory available even using a relatively small number of habitat configurations
# possible in landscapes containing larger numbers of patches. This is because the algorithm used to convert
# ID #'s associated with each configuration possible in a landscape still requires storage of vectors at least 
# as long as the number of configurations possible with the same number of patches. This is still much less 
# than 2^n.patches but these still become big quickly as the number of patches in landscapes increases.
# e.g. Using a server with ~ 50 cores and lots of memory it is very feasible to perform these computations with 
# landscapes of ~30 patches, but 50 patches becomes computationally impossible. Using a larger percentage of all
# the configurations possible in a landscape increases the accuracy of these metrics but increases the time 
# required to compute these metrics based on the number of cores available to you. e.g. If you have ~50 cores
# available for these calculations then calculations for ~50x more configurations can be used and computed in the 
# same time it would take to compute these metrics using the same number of configurations on 1 core.

##########################################################################################################
#rm(list=ls()) #to ensure a clean environment

#LOADING REQUIRED PACKAGES
library("Matrix") #to get sparse matrices
library("gdata") #to get upper and lower triangles of matrices
library("lava") #to get anti-diagonal of a matrix using revdiag
library("gtools") #to get even() function
library("tidyverse") #to get filter function
library("doParallel") #to run code in parallel
#SET WORKING DIRECTORY TO WHERE THE REQUIRED FUNCTIONS ARE
#setwd("C:/Users/Administrator/Desktop/JAN 2024/modulized code")
#LOADING REQUIRED FUNCTIONS
source("make_reduced_Gmat.r")
source("Gmat_to_Csubmat.r")
source("get_QED.r")
source("IDs_to_configs.r")
source("reduced_to_full_QED.r")
source("subsample_reducedQED.r")
source("make_landscape.r")
source("config_metrics_parallel.r")
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
#s.size<-100 # integer number of possible configurations to use (max 2^n.patches)
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
G.mat<-make_reduced_Gmat(n.patches=nrow(landscape),r.on=r.on,r.off=r.off)
C.submat<-Gmat_to_Csubmat(G.mat)
QED.output<-get_QED(C.submat)
QED<-QED.output[[1]]
QED<-QED/sum(QED) #standardizing to account for rounding errors
#OBTAINING A SUBSET OF ALL POSSIBLE HABITAT CONFIGURATIONS CHOSEN WITH THEIR PROBABILITY OF OCCURENCE AT QED
configs.subset<-subsample_reducedQED(QED=QED, n.patches=n.patches, s.size=s.size)
QED<-reduced_to_full_QED(QED=QED, n.patches=n.patches)
#RUNNING ALL THE WITHIN CONFIG METAPOPULATION METRIC CALCULATIONS IN PARALLEL:
##############################################################################
cl <- makeCluster(detectCores(), type='PSOCK')
registerDoParallel(cl)
config.metrics<-foreach(i=1:N, 
                      .combine=rbind) %dopar% 
  config_metrics_parallel(i=(i), N=N, configs=configs.subset,
                        landscape=landscape, e.rate=e.rate, c.rate=c.rate, 
                        gamma=gamma, self.rec=self.rec, 
                        alpha=alpha, iterations=iterations)
registerDoSEQ()
#END PARALLEL:
#############################################################################
n.on<-config.metrics[,4]
QED<-QED[n.on]
#CALCULATING ESTIMATES OF THE QED METAPOPULATION METRICS USING THE SUBSET OF ALL POSSIBLE HABITAT CONFIGURATIONS
QEDmetapop.metrics<-get_QEDMetapopMetrics(config.metrics=config.metrics, QED=QED)
QEDmetapop.metrics



