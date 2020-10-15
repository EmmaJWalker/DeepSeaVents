rm(list=ls()) #clear the workspace
library("pracma")
library("deSolve")
setwd("C:/Users/abuga/OneDrive/Desktop/HVMcode")
source("create_landscape_func.r")
source("SRLM_ODE_func.r")
source("at_equilibrium_rootfunc.r")
source("dispersal_kernel_func.r")
source("CTMC_and_SRLM_func.r")
source("lambda_M_func.r")

#Initializing Parameters:
###################################################################################################
total.t<-200#10 #total time to run the model (optional; could run to absorption)
landscape<-create.landscape(n.patches=4, landscape.type="linear", landscape.limit=100, 
                            patch.distribution="uniform", areas.distribution="uniform", areas.limit=1, 
                            clustering.iters=0)
n.patches<-length(landscape$patch.ID)
e.rate<-0.1#0.3
c.rate<-4#0.5
alpha<-0.1 #inverse of the mean dispersal distances (negetive exponential)
#network structure parameters:
gamma<-0
epsilon<-n.patches #1
self.rec<-1
#rates of vents turning off to on and on to off respectively
r1<-0.02
r2<-0.01
####################################################################################################
#Initial.Lm<-3000
#first let's scale our e.rate and c.rate to ensure a high enough persistence capacity the
#metapopulation has the opportunity to experience some growth
#delta<-Initial.Lm/(get.lambda.M(landscape=landscape, alpha=alpha, e.rate=e.rate, c.rate=c.rate, self.rec=self.rec, gamma=gamma, epsilon=epsilon))
#delta
#e.rate<-1
#c.rate<-1/delta
#lambda.M<-get.lambda.M(landscape=landscape, alpha=alpha, e.rate=e.rate, c.rate=c.rate, self.rec=self.rec, gamma=gamma, epsilon=epsilon)*delta
#lambda.M


sim.data<-CTMC.SRLM(total.t, landscape, e.rate, c.rate, alpha, gamma, epsilon, self.rec, r1, r2)
  sim.data<-sim.data[[1]]
#weight the occupancy values by area and plot metapop size over time
area.weights<-landscape$areas/sum(landscape$areas)
#metapop.size<-sim.data[,-1]*area.weights *what Austin has, which gives the size of each patch 
#but we don't wanna plot this for all of them, I'm just going to plot the sum
metapop.size<-rowSums(sim.data[,-1]) #*area.weights)
time<-sim.data$time
#figure
plot(time, metapop.size)