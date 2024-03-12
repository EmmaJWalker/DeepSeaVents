rm(list=ls()) #clear the workspace

#LOADING CODE SPECIFIC PACKAGES:
library("R.utils") #for int to binary
library("ggplot2") #for plotting
library("pracma") #for matrix math
library("deSolve") #for ode solving
library("RColorBrewer") #for plotting colours
library("colorspace") #for plotting colours
library("ggthemes") #for plotting themes
library("R.utils") #for integer to binary conversion
library("dplyr") #for various stuff
library("gdata") #to get upper triangle of a matrix
library("lava") #to get the reverse diagonal of a matrix
library("Matrix") #for sparse matrices
library("reticulate") #for something...
library("gtools") #for na.replace
#LOADING REQUIRED FUNCTIONS
setwd("C:/Users/Administrator/Desktop/DEC 2022")
source("BigintToBin_func.r") #create.landscape()
source("create_landscape_func.r") #create.landscape()
source("CTMC_and_SRLM_func_from_QED.r") #CTMC.SRLM()
source("dispersal_kernel_func.r") #disp.kernel()
source("lambda_M_func.r") #get.lambda.M()
source("QED_metrics_func_for_large_N_BigintToBin.r") #QED.metrics()
source("dewoody_func.r") #dewoody()
source("Pstar_func.r") #pstar.function()
source("reduced_quasi_eq_sparse_func2.r") #quasi.eq()
source("SRLMODE_func.r") #SRLM.ODE()
source("at_equilibrium_rootfunc.r") #at.equilibrium()
#source("single_parallel_batch_func.r") #single.batch()
getwd()
setwd("I:/one rep 30 patches")  


#SETTING EXPERIMENT PARAMETERS:
r.offs<-c(1/100000, 1/10000,1/1000,1/100,1/10)
r.ons<-c(1/100,1/10,1)
alphas<-c(1/100,1/10,10) #avg. disp = 1/10 interpatch distance, the interpatch distance and global disp
n.patches<-10
landscape.type<-"linear"
landscape.limit<-100
patch.distribution<-"uniform"
areas.distribution<-"uniform"
areas.limit<-1
clustering.iters<-0

#CREATING EXPERIMENTAL LANDSCAPES
landscape<-create.landscape(n.patches=n.patches, landscape.type=landscape.type, landscape.limit=landscape.limit, 
                            patch.distribution=patch.distribution, areas.distribution=areas.distribution, areas.limit=areas.limit, 
                            clustering.iters=clustering.iters)
#print(ggplot(landscape, aes(x.coord,y.coord)) + theme_classic() + geom_point(aes(size = areas))) #plotting
n.patches<-length(landscape$patch.ID) #for future use

Sys.time()
results<-QED.metrics(sample=1000, landscape=landscape, e.rate=0.1, c.rate=0.2, gamma=0, epsilon=n.patches, self.rec=1, alpha=1/100, r.on=1/100, r.off=1/10000)
Sys.time()
results

#WORKS FAST EVEN FOR 100 PATCHES NOW!!! WOOOOOOOOOOOOO!!!!!!