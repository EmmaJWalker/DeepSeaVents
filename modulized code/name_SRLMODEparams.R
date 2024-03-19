##########################################################################################################
# DESCRIPTION:
##########################################################################################################

# This function creates a named list of the parameters to be supplied to 
# the SRLMODE function because R’s ode solvers require all parameters to be 
# supplied in a single named list of with unique names for every parameter use 
# in an equation (aka if a parameter is a vector or a matrix every single 
# element of that vector or matrix must have a unique parameter name! how 
# annoying!).

# INPUTS:
# parameter.values: a concatenated list of containing all the objects defining 
#                   all the parameter values to be used to construct the system 
#                   of differential equations for the SRLM: n.patches, c.rate, 
#                   self.rec, extinction.rates, x.coord, y.coord, areas, config, 
#                   f.vals

# OUTPUTS: a list containing:
# parameters: a named list where each entry contains the unique name of a 
#             parameter and its value, for every single value needed to 
#             construct the system of differential equations for the SRLM

# REQUIRED PACKAGES:
# none

# REQUIRED FUNCTIONS: 
# name_Mij 

####################################################################################
# FUNCTION CODE:
####################################################################################
name_SRLMODEparams<-function(parameter.values){
  n.patches<-parameter.values[1]

  # make up names to each element of the dispersal kernel's matrix
    f.names<-as.vector(name_Mij(mat=matrix(parameter.values[-(1:(3+5*n.patches))], 
                                       nrow=n.patches, ncol=n.patches),
                                mat.name="f"))
  
  parameter.names<-c("n.patches", "c.rate", "self.rec",
                     rep(paste0("extinction.rates",1:n.patches)),
                     rep(paste0("x.coord",1:n.patches)),
                     rep(paste0("y.coord",1:n.patches)),
                     rep(paste0("areas",1:n.patches)),
                     rep(paste0("config",1:n.patches))
                     ,f.names)
  parameters<-list(c(setNames(parameter.values, parameter.names)))

  return(parameters)
}
####################################################################################
# QUICK CHECK: (uncomment and run to check)
####################################################################################
#rm(list=ls()) #to ensure a clean environment
## now run the FUNCTION CODE!!!
#setwd("/Users/abuga/Dropbox/Mac/Desktop/git/DeepSeaVents/modulized code")
#source("make_landscape.r")
#source("get_distmat.r")
#source("get_dispkernel.r")
#source("name_Mij.r")
#e.rate<-0.1
#c.rate<-0.2
#alpha<-0.01
#gamma<-0
#self.rec<-1
#landscape<-make_landscape(n.patches=4, landscape.type="line", landscape.size=20, 
#                 p.dist="uniform", p.size.dist="uniform", max.p.size=1, 
#                 clust.its=0)
#n.patches<-length(landscape$patch.ID)
#x.coord<-landscape$x.coord
#y.coord<-landscape$y.coord
#areas<-landscape$p.sizes
##calculate the extinction rate for each patch if area dependent
#extinction.rates<-e.rate/areas 
## construct the distance matrix
#dist.mat<-get_distmat(landscape = landscape)
##Getting the dispersal kernel f(dij)
#f.vals<-get_dispkernel(dist.mat=dist.mat,alpha=alpha,gamma=gamma,self.rec=self.rec)
#f.vals<-as.vector(f.vals)
#config<-c(1,0,1,1) #making up a configuration
#parameter.values<-c(n.patches,c.rate,self.rec,extinction.rates,x.coord,y.coord,areas,config
#                    ,f.vals)
#name_SRLMODEparams(parameter.values = parameter.values)
##########################################################################################################

