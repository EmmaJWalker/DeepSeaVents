# Set parameters script

rm(list=ls()) #to ensure a clean environment

#SET WORKING DIRECTORY TO WHERE THE REQUIRED FUNCTIONS ARE
setwd("C:/Users/Administrator/Desktop/JAN 2024/modulized code")

#SET UP PARAMETERS FOR DESIRED METAPOPULATION AND DYNAMIC LANDSCAPE
####################################################################
n.patches<-10 # integer equal to the total number of patches in the system
landscape.type<-"line" # set to "line" or "plane" indicating whether to have a 1D or 2D patch arrangement
landscape.size<-100 # numerical value indicating max y (and x coordinates if 2D)
p.dist<-"uniform" # set to "clustered", "random", "uniform" for how patches should be distributed in the landscape
p.size.dist<-"uniform" # set to "lognormal", "random" or "uniform" indicating how patch sizes are to be distributed
max.p.size<-1 # numerical value indicating largest patch size
clust.its<-1000 # integer equal to the number of times to iterate the algorithm increasing clustering or uniformity for
#                 the distribution of patches in the landscape
r.on<-0.1 # the rate at which patch gain
r.off<-0.01 # the rate at which patch loss
s.size<-100 # integer number of possible configurations to use (max 2^n.patches)
N<-48 # integer number of computational cores available
e.rate<-0.1 # within patch population extinction rate
c.rate<-0.2 # colonization rate of patches when a disperse arrives at a patch
alpha<-0.1 # 1/(the avg. dispersal distance of a species)
iterations<-1000 # integer value indicating the number of times to perform the interative function for finding pstar
#                  (larger => greater accuracy)
#############################################################################

#source("Script for calculating or estimating QEDMetapopMetrics with a subset of configs (large N).r")
#annoyingly I can't seem to source and run scripts running things in parallel...

source("Script for calculating QEDMetapopMetrics using all possible configs (small N).r")
QEDmetapop.metrics
