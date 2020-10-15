rm(list=ls()) #clear the workspace
library("pracma")
library("deSolve")
setwd("C:/Users/abuga/OneDrive/Desktop/HVMcode")
source("create_landscape_func.r")
source("SRLM_ODE_func.r")
source("at_equilibrium_rootfunc.r")
source("dispersal_kernel_func.r")
source("lambda_M_func.r")

#Initializing Parameters:
###################################################################################################
total.t<-5000 #total time to run the model (optional; could run to absorption)
landscape<-create.landscape(n.patches=4, landscape.type="2D", landscape.limit=20, 
                patch.distribution="clustered", areas.distribution="random", areas.limit=2, 
                clustering.iters=10)
n.patches<-length(landscape$patch.ID)
alpha<-0.3 #dispersal occurs across the whole landscape
e.rate<-0.01
c.rate<-4
self.rec<-1
gamma<-0
epsilon<-n.patches
#e.rate<-0.3
#c.rate<-0.5
delta<-e.rate/c.rate
#alpha<-0.1 #inverse of the mean dispersal distances (negetive exponential)
#network structure parameters:
#gamma<-0
#epsilon<-1
#self.rec<-1
#rates of vents turning off to on and on to off respectively
r1<-0.02
r2<-0.01
config<-rbinom(n.patches, 1, 0.5)
####################################################################################################

#calculations and set up to obtain all our parameters
####################################################################################################
#plot the habitat network and calculate necessary preliminary matrices/values
x.coord<-landscape$x.coord
y.coord<-landscape$y.coord
areas<-landscape$areas
#calculate E for each patch
extinction.rates<-e.rate/areas
# construct the distance matrix
dist.mat<-matrix(rep(0,n.patches*n.patches),n.patches,n.patches)
for (i in 1:n.patches){
  for (j in 1:n.patches){
    dist.mat[i,j]<-sqrt((x.coord[i]-x.coord[j])^2+(y.coord[i]-y.coord[j])^2)
  }
}
#Getting the f(dij) values to plug into the rhs
f.vals<-disp.kernel(dist.mat,alpha,gamma,epsilon,self.rec,n.patches)
f.vals<-as.vector(f.vals)
f.names<-matrix(rep(NA,n.patches*n.patches),n.patches,n.patches)
for (i in 1:n.patches){
  for (j in 1:n.patches){
    f.names[i,j]<-paste0("f.vals",i,j)
  }
}
f.names<-as.vector(f.names)
###################################################################################################

#creating a list of the parameters we calculated to be passed to our ODE function
###################################################################################################
names<-c("n.patches", "c.rate", "self.rec",
         rep(paste0("extinction.rates",1:n.patches)),
         rep(paste0("x.coord",1:n.patches)),
         rep(paste0("y.coord",1:n.patches)),
         rep(paste0("areas",1:n.patches)),
         rep(paste0("config",1:n.patches))
         ,f.names
         )
values<-c(n.patches,c.rate,self.rec,extinction.rates,x.coord,y.coord,areas,config
          ,f.vals
          )
parameters<-list(c(setNames(values, names)))
###################################################################################################

#setting our initial values
###################################################################################################
new.ICs<-0.5*config
IC.names<-rep(paste0("p",1:n.patches))
IC.values<-new.ICs
p<-setNames(IC.values,IC.names)
tau<-4000 #just for checking
###################################################################################################

#first let's scale our e.rate and c.rate to ensure a high enough persistence capacity the
#metapopulation has the opportunity to experience some growth
lambda.M<-get.lambda.M(landscape, alpha, gamma, epsilon, self.rec, e.rate, c.rate)
lambda.M
e.rate/c.rate
#CHECK against matlab
#just running the ODE
output<-ode(SRLM.ODE, parameters,times=seq(0,round(tau,0),1),y=p)
head(output)
#running the ODE with stopping function when the change in patch occupancy is less than 1e-4
#and therefore essentially at equilibrium
out <- ode(SRLM.ODE, y=p, parms=parameters, times =seq(0,5000,1), rootfun = at.equilibrium)
head(out)
#and to check the out
#tail(out, n = 2) #look at the end of the simulation
#x<-out[42,-1]-out[41,-1] #calculate the final change in patch occupancy
#Norm(x) #check that the norm is below 1e-4
plot(output[,1], rowSums(output[,-1]))
