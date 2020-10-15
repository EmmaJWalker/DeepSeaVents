rm(list=ls()) #clear the workspace
library("ggplot2")
library("pracma")
library("deSolve")
library("RColorBrewer")
library("colorspace")
library("wesanderson") #SUPER NICE COLOUR SCHEMES!
library("ggthemes")
setwd("C:/Users/abuga/Desktop/HVMcode")
source("create_landscape_func.r")
source("quasi_eq_func.r")
source("CTMC_func.r")
source("SRLM_ODE_func.r")
source("at_equilibrium_rootfunc.r")
source("dispersal_kernel_func.r")
source("CTMC_and_SRLM_func.r")
source("lambda_M_func.r")
source("metapop_size_at_QED.r")

n.patches<-2
landscape.type<-"linear"
landscape.limit<-10
patch.distribution<-"uniform"
areas.distribution<-"uniform"
areas.limit<-1
clustering.iters<-0

e.rate<-0.1
c.rate<-4
alpha<-0.3 
self.rec<-1
gamma<-0
epsilon<-n.patches

r1=0.001
r2=0.001

#Let's start with creating a uniform landscape of 4 patches arranged linearly:
###################################################################################################
landscape<-create.landscape(n.patches=n.patches, landscape.type=landscape.type, landscape.limit=landscape.limit, 
                            patch.distribution=patch.distribution, areas.distribution=areas.distribution,
                            areas.limit=areas.limit, clustering.iters=clustering.iters)
print(ggplot(landscape, aes(x.coord,y.coord)) + theme_classic() + geom_point(aes(size = areas)))
###################################################################################################

#Analytically calculating time to absorption and the Q.E.D.
###################################################################################################
QED<-quasi.eq(n.patches=n.patches, r.on=r1, r.off=r2)
dist<-QED[[1]] #this is the porportion of time spent between transitions at Q.E.D (or the proability 
# of any given configuration at Q.E.D.)
#eg. for a 2 patch system this would be the probability one patch, or the other, or both is on
dist
####################################################################################################









#simulating patch occupancy
###################################################################################################
sim.data<-CTMC.SRLM(total.t=3000, landscape=landscape, e.rate<-e.rate, c.rate<-c.rate, alpha=alpha,
                    gamma=gamma, epsilon=epsilon, self.rec=self.rec, r1=r1, r2=2)
config.data<-sim.data[[2]]
sim.data<-sim.data[[1]]
#weight the occupancy values by area and plot metapop size over time
area.weights<-landscape$areas/sum(landscape$areas)
#metapop.size<-sim.data[,-1]*area.weights *what Austin has, which gives the size of each patch 
#but we don't wanna plot this for all of them, I'm just going to plot the sum
metapop.size<-rowSums(sim.data[,-1]) #*area.weights) <-if we wanna add the area weights
metapop.size
time<-sim.data$time
#figure
plot(time, metapop.size)
#print(ggplot() + theme_classic()
#      + geom_line(aes(x = time, y = metapop.size),size=1, color="dodgerblue3") 
#      + labs(x = "time", y = "Metapopulation Size", 
#             title = "Combined CTMC and SRLM Simulation")
#      + theme(text = element_text(size=15)))


#now we can calculate
#average simulated metaop size over time calculated using trapezoidal method (As used by Austin)
#################################################################################################
avg.size=trapz(time,metapop.size)/time[length(time)]
avg.size
#################################################################################################

#versus the 
#weighted metapopulation size predicted at Q.E.D.
#################################################################################################
pred.avg.size<-metapop.size.at.QED(landscape=landscape, e.rate=e.rate, c.rate=c.rate, gamma=gamma, 
                                   epsilon=epsilon, self.rec=self.rec, alpha=alpha, r1=r1, r2=r2)
#################################################################################################
pred.avg.size

(3/2)-((3*e.rate)+(e.rate*(exp(-alpha*5)))/(2*c.rate*(1+exp(-alpha*5))))
(dist[3]*(1-e.rate/(c.rate*(1+exp(-alpha*5)))))+(dist[1]*2*(1-(e.rate/c.rate)))








config<-c(1,0)

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

#just running the ODE
output<-ode(SRLM.ODE, parameters,times=seq(0,round(tau,0),1),y=p)
head(output)
#running the ODE with stopping function when the change in patch occupancy is less than 1e-4
#and therefore essentially at equilibrium
out <- ode(SRLM.ODE, y=p, parms=parameters, times =seq(0,5000,1), rootfun = at.equilibrium)
head(out)

dist
