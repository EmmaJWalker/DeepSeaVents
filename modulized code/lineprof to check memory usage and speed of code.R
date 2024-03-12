#devtools::install_github("hadley/lineprof")
library("lineprof")

source("simulate_Metapop_impHabitat_BigN.r")

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
n.patches<-20
r.on<-0.1
r.off<-0.01
G.mat<-make_reduced_Gmat(n.patches=n.patches,r.on=r.on,r.off=r.off)
C.submat<-Gmat_to_Csubmat(G.mat)
output<-get_QED(C.submat)
QED<-output[[1]]
landscape<-make_landscape(n.patches=n.patches, landscape.type="line", landscape.size=100, 
                          p.dist="uniform", p.size.dist="uniform", max.p.size=1, 
                          clust.its=0)

prof <- lineprof(simulate_Metapop_impHabitat_BigN(total.trans=10, landscape=landscape, e.rate=0.1, c.rate=0.2, alpha=0.1, gamma=0, 
                                                  self.rec=1, r.on=r.on, r.off=r.off, QED=QED))
shine(prof)