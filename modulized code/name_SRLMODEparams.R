##########################################################################################################
# DESCRIPTION:
##########################################################################################################

# This function both creates 1) a named list of the parameters to be supplied to 
# the SRLMODE function because R’s ode solvers require all parameters to be 
# supplied in a single named list of with unique names for every parameter use 
# in an equation (aka if a parameter is a vector or a matrix every single 
# element of that vector or matrix must have a unique parameter name! how 
# annoying!) and 2) the names it should to use for the interpatch distances if 
# these do not already exist and haven't been supplied to this function.

# INPUTS:
# parameter.values: a concatenated list of containing all the objects defining 
#                   all the parameter values to be used to construct the system 
#                   of differential equations for the SRLM: n.patches, c.rate, 
#                   self.rec, extinction.rates, x.coord, y.coord, areas, config, 
#                   f.vals
# f.names: a vector of character values providing a "name" for each interpatch 
#          distance, if not supplied to this function it will generate these 
#          itself

# OUTPUTS: a list containing:
# parameters: a named list where each entry contains the unique name of a 
#             parameter and its value, for every single value needed to 
#             construct the system of differential equations for the SRLM
# f.names: a vector of character values providing a "name" for each interpatch 
#          distance

# REQUIRED PACKAGES:
# Matrix: for working with sparse matrices

# REQUIRED FUNCTIONS: 
# none


name_SRLMODEparams<-function(parameter.values, f.names){
  n.patches<-parameter.values[1]

  # to avoid generating the parameter names for all the interpatch distances 
  # more than once, we only do this if they haven't already been set
  if (is.na(f.names)==T){ #
    f.names<-matrix(rep(NA,n.patches*n.patches),n.patches,n.patches)
    for (i in 1:n.patches){
      for (j in 1:n.patches){
        f.names[i,j]<-paste0("f.vals",i,j)
      }
    }
    f.names<-as.vector(f.names)
  }
  
  parameter.names<-c("n.patches", "c.rate", "self.rec",
                     rep(paste0("extinction.rates",1:n.patches)),
                     rep(paste0("x.coord",1:n.patches)),
                     rep(paste0("y.coord",1:n.patches)),
                     rep(paste0("areas",1:n.patches)),
                     rep(paste0("config",1:n.patches))
                     ,f.names)
  parameters<-list(c(setNames(parameter.values, parameter.names)))

  return(list(parameters, f.names))
}

####################################################################################
# QUICK CHECK: (uncomment and run to check)
####################################################################################
#rm(list=ls()) #to ensure a clean environment
## now run the FUNCTION CODE!!!
source("make_landscape.r")
source("get_distmat.r")
source("get_dispkernel.r")
landscape<-make_landscape(n.patches=4, landscape.type="line", landscape.size=20, 
                 p.dist="uniform", p.size.dist="uniform", max.p.size=1, 
                 clust.its=0)
n.patches<-length(landscape$patch.ID)
x.coord<-landscape$x.coord
y.coord<-landscape$y.coord
areas<-landscape$areas
#calculate the extinction rate for each patch if area dependent
extinction.rates<-e.rate/areas 
# construct the distance matrix
dist.mat<-matrix(rep(0,n.patches*n.patches),n.patches,n.patches)
for (i in 1:n.patches){
  for (j in 1:n.patches){
    dist.mat[i,j]<-sqrt((x.coord[i]-x.coord[j])^2+(y.coord[i]-y.coord[j])^2)
  }
}
#Getting the dispersal kernel f(dij)
f.vals<-disp.kernel(dist.mat,alpha,gamma,epsilon,self.rec,n.patches)
f.vals<-as.vector(f.vals)
parameter.values<-c(n.patches,c.rate,self.rec,extinction.rates,x.coord,y.coord,areas,config
                    ,f.vals)