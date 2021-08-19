rm(list=ls()) #clear the workspace
library("ggplot2")
library("pracma")
library("deSolve")
library("RColorBrewer")
library("colorspace")
library("wesanderson") #SUPER NICE COLOUR SCHEMES!
library("ggthemes")
setwd("C:/Users/abuga/OneDrive/Desktop/DeepSeaVents-master")
source("quasi_eq_func_w_scipy.r")
source("CTMC_func.r")
source("SRLM_ODE_func2.r")
#source("SRLM_ODE_func.r")
#source("SRLM_ODE_func_nointeractions.r")
source("at_equilibrium_rootfunc.r")
source("dispersal_kernel_func.r")
source("CTMC_and_SRLM_func_w_transition_limit.r")
source("lambda_M_func.r")
source("lm and pstar at QED function.r")
source("create_landscape_func.r")
source("Pstar Function.r")
getwd()


#Let's start with creating a uniform landscape of 4 patches arranged linearly:
###################################################################################################
landscape<-create.landscape(n.patches=10, landscape.type="linear", landscape.limit=100, 
                 patch.distribution="uniform", areas.distribution="uniform", areas.limit=1, 
                 clustering.iters=0)
print(ggplot(landscape, aes(x.coord,y.coord)) + theme_classic() + geom_point(aes(size = areas)))
n.patches<-length(landscape$patch.ID) #for future use
###################################################################################################

#Analytically calculating time to absorption and the Q.E.D.
###################################################################################################
r.on=0.1
r.off=0.00001
QED<-quasi.eq(n.patches=n.patches, r.on=r.on, r.off=r.off)
QED
dist<-QED[[1]] #this is the porportion of time spent between transitions at Q.E.D (or the proability 
# of any given transition at Q.E.D.)
#eg. for a 2 patch system this would be the probability one patch, or the other, or both 
#switches either from on to off or visa versa
dist
T.absorp<-QED[[2]] #this is the expected time to absorption
T.absorp<-weighted.mean(T.absorp, dist) #expected time to absorption weighted by the QED
r.con<-QED[[3]]
r.con #QED supposedly only appropriate when r.con < 1
###################################################################################################

#Let's check for what ratios of r.on to r.off r.con < 1, to know the parameter space we can expect 
#the QED to be appropiate in
###################################################################################################
#note: takes < 10 min to run these calculations (if we wanted to not do these calculations 
#more than once in the whole project we could just save their output in files and read them in)
r.offs<-c(0.000001,0.00001,0.0001,0.001,0.01,0.1,1,10,100,1000,10000,100000)
for (i in 1:length(r.offs)){
  QED<-quasi.eq(n.patches=n.patches, r.on=1, r.off=r.offs[i])
  r.con<-QED[[3]]
  if (r.con < 1 ){
    print(r.offs[i])
  }
}
#ok so breaks at r.off 100x greater than r.on but let's just check out what that one looks like
#to see things just beyon the limits
QED<-quasi.eq(n.patches=n.patches, r.on=1, r.off=100)
QED[[4]] #so we can't even get a ratio of means distribution in this case (which would tell us 
#the probability we go to each state give we start in a particular state instead of the QED) so 
#I don't know how much we will care about this metric and might just want to take it out of our 
#calculations to save on computing power. but we are able to calculate everything else so could 
#otherwise explore this limit if desired

#but for now let's keep our analyses within the acceptable range for the QED
r.offs<-c(0.000001,0.00001,0.0001,0.001,0.01,0.1,1,10)
###################################################################################################

#Now we can compare this to the simulated time to absorption, given an initial configuration and a 
#ratio of r.on vs. r.off
###################################################################################################
QED<-quasi.eq(n.patches=n.patches, r.on=1, r.off=0.000001)
QED
dist<-QED[[1]] #this is the porportion of time spent between transitions at Q.E.D (or the proability 
# of any given transition at Q.E.D.)
#eg. for a 2 patch system this would be the probability one patch, or the other, or both 
#switches either from on to off or visa versa
dist
T.absorp<-QED[[2]] #this is the expected time to absorption
T.absorp<-weighted.mean(T.absorp, dist) #expected time to absorption weighted by the QED

#let's start with a random initial configuration
initial.config<-rbinom(n.patches, 1, 0.5)
T.sim<-CTMC.sim(timesteps=1000, n.patches=n.patches, r.on=1, r.off=0.000001, initial.config=initial.config)
T.sim$end.t
T.absorp
#NOTE: Okay, this can actually take a fairly long time to run a simulation even to just 1000 timesteps
#so we will probably want to do this in parallel when this is the case. For now since we're just testing 
#things let's pick a ratio of r.on to r.off for which T.absorp is relatively low just to speed things up
r.offs<-c(0.000001,0.00001,0.0001,0.001,0.01,0.1,1,10,100,1000,10000,100000)
for (i in 1:length(r.offs)){
  QED<-quasi.eq(n.patches=n.patches, r.on=1, r.off=r.offs[i])
  T.absorp<-QED[[2]]
  print(mean(T.absorp))
  if (mean(T.absorp) < 1000 ){
    print(r.offs[i])
  }
} #okay so actually for all of those we can calculate we are way over 1000 t so we can't really assess
#accuracy of T.absorp in estimating the simulated mean time to absorption without expending a huge 
#amount of computing power doing so. We don't actually need to verify this trusting that mathematicians
#that came up with this estimate are correct, it's just nice to be able to check our simulations and
#calculations match up with what they theoretically should be. Since T.absorp drops dramatically with
#fewer patches we can verify this in those cases and then reasonably trust things hold for more patches.

#we don't expect this to be exactly the same as analytically derived but to get an idea how 
#reasonable our estimate is we need to do multiple simulations
###################################################################################################

#so let's do a number of simulations
#to look at the distribution of simulated time to absorption values versus the expected analytically 
#derived time to absorption to see how they compare -we're not going to do this for 10 patches because 
#it will take too long
####################################################################################################
#runs<-2000
#times.to.absorption<-rep(NA, runs)
#mean.transition.times<-rep(NA, runs)
#for (i in 1:runs){
#  initial.config<-rbinom(n.patches, 1, 0.5)
#  print(initial.config)
#  output<-CTMC.sim(timesteps=1000, n.patches=n.patches, r.on=0.1, r.off=1, 
#                   initial.config=initial.config)
#  times.to.absorption[i]<-output[[1]]
#  mean.transition.times[i]<-output[[2]]
#}
#NOTE: this function is actually very fast provided the expected time to absorption is low!
#times.to.absorption
#hist(times.to.absorption, breaks = 100, main="Predicted vs Simulated Time to Absorption", xlab="Time to Absorption")
#mean<-mean(times.to.absorption[!is.na(times.to.absorption)])
#abline(v=T.absorp, col="red")
#abline(v=mean, lty=2)
#setwd("C:/Users/abuga/OneDrive/Desktop/HVMcode/preliminary plots")
#dev.copy(png,'Predicted vs Simulated Time to Absorption sample.png')
#dev.off()
#colours() #look up R's colours if you want to get fancy
####################################################################################################

#Now we can explore how well we analytically predict T.absorption as we vary the parameters
#begining with r1 vs r2 -let's not because this will take too long with 10 patches
####################################################################################################
#runs<-100 #higher this number drops the speed but not too badly
#r.ons<-c(0.01, 1, 10)
#r.offs<-c(0.01, 1, 10)
#patch.numbers<-c(2,4)#seq(2,6, 2) #NOTE: fast and under a min, until you add varying the patch numbers too
#get's slower the more patches you add for the low r1 and r2
#you can always comment out doing this if you just want to look at varying r1 and r2
#can always run in parallel on the server and then should take about a min
#time.limit<-100000000#0000 #higher this number doesn't actually affect the speed at all really
#for(n in 1:length(patch.numbers)){
#  n.patches<-patch.numbers[n]
  
#  for (j in 1:length(r.ons)){
#    for (k in 1:length(r.offs)){
#      r.on<-r.ons[j]
#      r.off<-r.offs[k]
#      QED<-quasi.eq(n.patches=n.patches, r.on=r.on, r.off=r.off)
#      dist<-QED[[1]]
#      T.absorp<-QED[[2]] #this is the expected time to absorption
#      T.absorp<-weighted.mean(T.absorp, dist)
#      print(T.absorp) #so we can check whether this is likely to take a really long time to run or not
      #and can abort if we don't want to wait so long
      
#      if (T.absorp<time.limit){ #if T.absorp falls within the timelimit we are willing to simulate for
        #simulate
#        print(T.absorp)
        
#        times.to.absorption<-rep(NA, runs)
#        mean.transition.times<-rep(NA, runs)
#        for (i in 1:runs){
#          initial.config<-rbinom(n.patches, 1, 0.5)
#          output<-CTMC.sim(timesteps=time.limit, n.patches=n.patches, r.on=r.on, r.off=r.off, 
#                           initial.config=initial.config)
#          times.to.absorption[i]<-output[[1]]
#          mean.transition.times[i]<-output[[2]]
#        }
#        hist(times.to.absorption, breaks = 100, main=paste0("n=",n.patches," r.on=",r.on," r.off=",r.off),
#             xlab="Time to Absorption")
#        mean.Tabsorp<-mean(times.to.absorption)
#        abline(v=T.absorp, col="red")
#        abline(v=mean, lty=2)
#        setwd("C:/Users/abuga/OneDrive/Desktop/HVMcode/preliminary plots/Accuracy histograms")
#        dev.copy(png, paste0("N=",n.patches, "r.on=",r.on, "r.off=",r.off, "Tabsorp comparison 4 test.png"))
#        dev.off()
        
#        mean.trans<-mean(mean.transition.times)
#      } else {
#        mean.Tabsorp<-NA
#      }
      
#      out<-data.frame(n.patches,r.on,r.off,T.absorp, mean.Tabsorp, mean.trans)
#      if (j==1 & k==1 & n==1){
#        T.absorp.data<-out
#      } else{
#        T.absorp.data<-rbind(T.absorp.data, out)
#      }
      #NOTE: this function is actually very fast!
#    }
#  }
#}
#setwd("C:/Users/abuga/OneDrive/Desktop/HVMcode/data files")
#write.csv("Tabsorp predicted vs simulated comparison data 4 test.csv")
#T.absorp.data<-read.csv("Tabsorp predicted vs simulated comparison data 4 test.csv")
#pal <- choose_palette()
#print(ggplot(data=T.absorp.data)
#      + theme_classic()
#      + geom_line(aes(x = log(r2/r1), y = T.absorp-mean.Tabsorp, group=n.patches, colour=n.patches),
#                  size=1)) 
#      #+ labs(x = "ln(r2/r1)", y = "Predicted Avg. T.absorption - Simulated Avg. T.absorption", 
#      #      title = "Comparison of Predicted versus Simulated Time to Absorption")
#      #+ scale_color_continuous("Number of Patches", high=pal(8), low=pal(1), guide="legend", 
#      #                         breaks=c(2,4,6), labels=c(2,4,6))
#      #+ theme(text = element_text(size=15))
#      #+ theme(legend.text=element_text(size=12)))
# *** It's not that it's actually bad with more patches, it's that we don't run till absorption in the simulation***
# so if we wanted to look at accuracy in this parameter space we would need to run the CTMC longer
#although, now I've subset to only look at the data where 
#print(ggplot(data=T.absorp.data[!is.na(T.absorp.data$mean.Tabsorp) & !T.absorp.data$T.absorp==-Inf,])
#      + theme_classic()
#      + geom_point(aes(x = log(r2/r1), y = mean.trans-(mean.trans/T.absorp), group=n.patches, colour=n.patches),
#                  size=2) 
#      + labs(x = "ln(r2/r1)", y = "Mean Transitional Time - (Mean Transitional Time/Time to Absorption)", 
#             title = "Accuracy")
#      + scale_color_continuous("Number of Patches", high=pal(8), low=pal(1), guide="legend", 
#                               breaks=c(2,4,6), labels=c(2,4,6))
#      + theme(text = element_text(size=15))
#      + theme(legend.text=element_text(size=12)))
#mean transition time-(mean transition time / time to absorption) vs r2/r1
#print(ggplot(data=T.absorp.data[!is.na(T.absorp.data$mean.Tabsorp) & !T.absorp.data$T.absorp==-Inf,])
#      + theme_classic()
#      + geom_point(aes(x = (r2/r1), y = mean.trans-(mean.trans/T.absorp), group=n.patches, colour=n.patches),
#                   size=2) 
#      + labs(x = "r2/r1", y = "Mean Transitional Time - (Mean Transitional Time/Time to Absorption)", 
#             title = "Accuracy")
#      + scale_color_continuous("Number of Patches", high=pal(8), low=pal(1), guide="legend", 
#                               breaks=c(2,4,6), labels=c(2,4,6))
#      + theme(text = element_text(size=15))
#      + theme(legend.text=element_text(size=12)))
#mean transition time-(mean transition time / time to absorption) vs. r2
#print(ggplot(data=T.absorp.data[!is.na(T.absorp.data$mean.Tabsorp) & !T.absorp.data$T.absorp==-Inf,])
#      + theme_classic()
#      + geom_point(aes(x = (r2), y = mean.trans-(mean.trans/T.absorp), group=n.patches, colour=n.patches),
#                   size=2) 
#      + labs(x = "r2", y = "Mean Transitional Time - (Mean Transitional Time/Time to Absorption)")
#      + scale_color_continuous("Number of Patches", high=pal(8), low=pal(1), guide="legend", 
#                               breaks=c(2,4,6), labels=c(2,4,6))
#      + theme(text = element_text(size=15))
#      + theme(legend.text=element_text(size=12)))
#mean transition time-(mean transition time / time to absorption) vs. r1
#print(ggplot(data=T.absorp.data[!is.na(T.absorp.data$mean.Tabsorp) & !T.absorp.data$T.absorp==-Inf,])
#      + theme_classic()
#      + geom_point(aes(x = (r1), y = mean.trans-(mean.trans/T.absorp), group=n.patches, colour=n.patches),
#                   size=2) 
#      + labs(x = "r1", y = "Mean Transitional Time - (Mean Transitional Time/Time to Absorption)")
#      + scale_color_continuous("Number of Patches", high=pal(8), low=pal(1), guide="legend", 
#                               breaks=c(2,4,6), labels=c(2,4,6))
#      + theme(text = element_text(size=15))
#      + theme(legend.text=element_text(size=12)))
####################################################################################################


#SKIPPING AHEAD TO HERE for 10 patches!
#now let's look at the simulation of the combined model
# let's just use r1=0.2 and r2=0.1 (r.off 0.5x r.on) to test this since that runs pretty quick
###################################################################################################
#let's create a landscape for this
landscape<-create.landscape(n.patches=10, landscape.type="linear", landscape.limit=100, 
                            patch.distribution="uniform", areas.distribution="uniform", areas.limit=1, 
                            clustering.iters=0)
print(ggplot(landscape, aes(x.coord,y.coord)) + theme_classic() + geom_point(aes(size = areas)))
n.patches<-length(landscape$patch.ID) #for future use

r.on<-1/100
r.off<-1/100000

e.rate<-1#0.1 #can pick just any number to start since we will just rescale this
c.rate<-1#4  #can pick just any number to start since we will just rescale this
alpha<-1/100 
self.rec<-1
gamma<-0
epsilon<-n.patches


#THIS IS A SLOW FUNCTION EVEN WHEN IT HAS THE QED TO START :/
#first let's scale our e.rate and c.rate to ensure a high enough persistence capacity the
#metapopulation has the opportunity to experience some growth
lambda.M<-get.lambda.M(landscape=landscape[1,], alpha=alpha, e.rate=e.rate, c.rate=c.rate, self.rec=self.rec, gamma=gamma, epsilon=epsilon)
landscape[1,] #landscape w only 1 patch
lambda.M
Initial.Lm<-100
#first let's scale our e.rate and c.rate to ensure a high enough persistence capacity the
#metapopulation has the opportunity to experience some growth
#delta<-Initial.Lm/(get.lambda.M(landscape=landscape, alpha=alpha, e.rate=e.rate, c.rate=c.rate, self.rec=self.rec, gamma=gamma, epsilon=epsilon))
#delta
#e.rate<-1/delta
#c.rate<-1
#lambda.M<-get.lambda.M(landscape=landscape, alpha=alpha, e.rate=e.rate, c.rate=c.rate, self.rec=self.rec, gamma=gamma, epsilon=epsilon)*delta
#lambda.M

#scaling lmQED rather than the max lm with all patches on:
lm.n.pstar<-lm.n.pstar.QED(landscape=landscape, e.rate=e.rate, c.rate=c.rate, gamma<-0, epsilon=n.patches, self.rec<-1, alpha=alpha, r1<-r.on, r2<-r.off)
exp.QED.size<-lm.n.pstar[[1]]
delta<-Initial.Lm/lm.n.pstar[[2]]
delta
e.rate<-1/delta
c.rate<-1
lm.n.pstar<-lm.n.pstar.QED(landscape=landscape, e.rate=e.rate, c.rate=c.rate, gamma<-0, epsilon=n.patches, self.rec<-1, alpha=alpha, r1<-r.on, r2<-r.off)
exp.QED.size<-lm.n.pstar[[1]]
exp.QED.size
exp.lm<-lm.n.pstar[[2]]*delta
exp.lm

#expected metapop size at QED
#lm.n.pstar<-lm.n.pstar.QED(landscape=landscape, e.rate=e.rate, c.rate=c.rate, gamma<-0, epsilon=n.patches, self.rec<-1, alpha=alpha, r1<-r.on, r2<-r.off)
#exp.QED.size<-lm.n.pstar[[1]]
#exp.QED.size
#exp.lm<-lm.n.pstar[[2]]*delta
#exp.lm

output<-CTMC.SRLM(total.trans=1000, landscape=landscape, e.rate<-e.rate, c.rate<-c.rate, alpha<-alpha,
                  gamma<-0, epsilon=n.patches, self.rec<-1, r1<-r.on, r2<-r.off) 
simu.data<-output[[1]]
configs.data<-output[[2]]
write.csv(simu.data, paste0("Simudata Combined CTMC and SRLM simulation sample lmQED100 d", delta, 
                            " ron", r.on, " roff", r.off, " a", alpha,  ".csv"))
write.csv(configs.data, paste0("configdata Combined CTMC and SRLM simulation sample lmQED100 d", delta, 
                               " ron", r.on, " roff", r.off, " a", alpha,  ".csv"))

#weight the occupancy values by area and plot metapop size over time
area.weights<-landscape$areas/sum(landscape$areas)
#metapop.size<-sim.data[,-1]*area.weights *what Austin has, which gives the size of each patch 
#but we don't wanna plot this for all of them, I'm just going to plot the sum
metapop.size<-rowSums(simu.data[,2:(n.patches+1)]) #*area.weights) <-if we wanna add the area weights
time<-simu.data$time
n.on<-rowSums(configs.data[,2:(n.patches+1)])
taus<-configs.data[,1]


#average simulated metaop size over time calculated using trapezoidal method (As used by Austin)
avg.size=trapz(time,metapop.size)/time[length(time)]
avg.size

arithmetic.mean.size<-lm.n.pstar[[3]]

#figure
#plot(time, metapop.size)
metapop.data<-data.frame(time, metapop.size)
habitat.data<-data.frame(taus,n.on)
print(ggplot() + theme_classic()
      + geom_rect(aes(xmin = taus, xmax = lead(taus), 
                      ymin = 0, ymax = n.on), 
                  fill = "dodgerblue")
      + geom_line(data = metapop.data, aes(x = time, y = metapop.size),size=1, color="yellowgreen") 
      + geom_step(data = habitat.data, aes(x = taus, y = n.on), color="dodgerblue", alpha=0.5)
      #+ geom_area(aes(x=taus, y = n.on), color="dodgerblue3", alpha=0.5)
      #+ geom_col(data = habitat.data, aes(x = taus, y = n.on), fill="dodgerblue3", alpha=0.5)
      + labs(x = "time", y = "Metapopulation Size", 
             title = "Combined CTMC and SRLM Simulation")
      + theme(text = element_text(size=15))
      + scale_y_continuous(sec.axis = sec_axis(~.*1, name = "Number of Habitable Patches"))
      + geom_hline(yintercept=exp.QED.size, linetype="dashed", color = "red")
      + geom_hline(yintercept=avg.size, linetype="dashed", color = "yellow")
      #+ geom_hline(yintercept=arithmetic.mean.size, linetype="dashed", color = "blue")
      )
#setwd("C:/Users/abuga/Documents/HVM rscripts/Emma HVM cleaned/Integrated code for clustered landscapes/preliminary plots")
dev.copy(png,'Combined CTMC and SRLM simulation sample lmQED100 d20.png')
dev.off()

output[[1]]
output[[2]]

habitat.data
metapop.data

plot(taus, n.on)
print(ggplot() + theme_classic()
      + geom_line(data= habitat.data, aes(x = taus, y = n.on)))


#attempting to plot with configurations
plot<-ggplot() + theme_classic() + 
  geom_line(aes(x = time, y = metapop.size),size=1, color="dodgerblue3") + 
  geom_col(aes(x = taus, y = n.on), fill="dodgerblue3", alpha=0.5) + 
  labs(x = "time", y = "Metapopulation Size", 
       title = "Combined CTMC and SRLM Simulation") + theme(text = element_text(size=15)) 
ypos<-rep(NA, length(taus))
size<-rep(NA, length(taus))
on<-rep(NA, length(taus))
for (i in 1:n.patches){
  ypos<-rep(landscape$y.coord[i]/10*2, length(taus))
  size<-rep(landscape$areas[i], length(taus))
  on<-as.factor(configs.data[,i+1])
  plot<-plot + geom_point(aes_string(x=taus, y=ypos[i], fill=factor(on[i])), size=size[i]*5)  
}
print(plot + scale_y_continuous(sec.axis = sec_axis(~.*1, name = "Number of Habitable Patches")) +
        scale_fill_manual(values = c("dodgerblue3", "black")))
#I just can't get the circles colours to alternate according to whether they are on or not

#################################################################################################

#now we can calculate
#average simulated metaop size over time calculated using trapezoidal method (As used by Austin)
#################################################################################################
avg.size=trapz(time,metapop.size)/time[length(time)]
avg.size
#################################################################################################

#versus the 
#weighted metapopulation size predicted at Q.E.D.
#################################################################################################
pred.avg.size<-metapop.size.at.QED(landscape=landscape, e.rate=e.rate, c.rate=c.rate, gamma<-0, 
                            epsilon=n.patches, self.rec<-1, alpha<-0.3, r1<-0.2, r2<-0.1)
pred.avg.size
#################################################################################################
#pred.avg.size is slightly higher 
# will need to perform multiple simulations to know how they compare

#now can perform multiple simulations to compare the average simulated metapop size versus the 
#weighted predicted metapop size at Q.E.D.
#################################################################################################
r1s<-c(0.01, 0.25, 0.5, 0.75, 1) #0.2
r2s<-c(0.01, 0.25, 0.5, 0.75, 1) #0.2
runs<-10 #100
time.limit<-1000
patch.numbers<-seq(2,12, 2) #OOOOOOOOH this is varying patch number....
#n.patches<-10
for(n in 1:length(patch.numbers)){
  n.patches<-patch.numbers[n]
  
  for (j in 1:length(r1s)){
    for (k in 1:length(r2s)){
      r1<-r1s[j]
      r2<-r2s[k]
      
      times.to.absorption<-rep(NA, runs)
      mean.transition.times<-rep(NA, runs)
      avg.sim.size<-rep(NA, runs)
      for (i in 1:runs){
        #get the times to absorption and mean transition times
        ##########################
        initial.config<-rbinom(n.patches, 1, 0.5)
        output<-CTMC.sim(timesteps=time.limit, n.patches=n.patches, r1=r1, r2=r2, 
                         initial.config=initial.config)
        times.to.absorption[i]<-output[[1]]
        mean.transition.times[i]<-output[[2]]
        ###########################
        #create a landscape
        #######################################
        landscape<-create.landscape(n.patches=n.patches, landscape.type="linear", landscape.limit=10, 
                                    patch.distribution="uniform", areas.distribution="uniform", areas.limit=1, 
                                    clustering.iters=0)
        n.patches<-length(landscape$patch.ID)
        #######################################
        
        #run our simulation of the CTMC and SRLM
        ########################################
        sim.data<-CTMC.SRLM(total.t=time.limit, landscape=landscape, e.rate=e.rate, c.rate=c.rate, alpha<-0.3,
                            gamma<-0, epsilon=n.patches, self.rec<-1, r1=r1, r2=r2)
        #weight the occupancy values by area and plot metapop size over time
        area.weights<-landscape$areas/sum(landscape$areas)
        metapop.size<-sim.data[,-1]*area.weights #*what Austin has, which gives the size of each patch 
        #but we don't wanna plot this for all of them, I'm just going to plot the sum
        metapop.size<-rowSums(sim.data[,-1]) #*area.weights) <-if we wanna add the area weights
        time<-sim.data$time
        T.eq<-time[length(time)]
        ########################################
        
        #now we can calculate
        #average simulated metaop size over time calculated using trapezoidal method (As used by Austin)
        ########################################
        avg.sim.size[i]<-trapz(time,metapop.size)/time[length(time)]
        ########################################
      }
      
      #versus the 
      #weighted metapopulation size predicted at Q.E.D.
      ########################################
      pred.avg.size<-metapop.size.at.QED(landscape=landscape, e.rate=e.rate, c.rate=c.rate, gamma<-0, 
                                         epsilon=n.patches, self.rec<-1, alpha<-0.3, r1=r1, r2=r2)
      ########################################
      QED<-quasi.eq(n.patches=n.patches, r1=r1, r2=r2)
      T.absorp<-QED[[2]] #this is the expected time to absorption
      T.QED<-QED[[3]]
      
      hist(avg.sim.size, breaks=20, main="Predicted vs Simulated Avg. Metapop Size", xlab = "Metapop Size")
      mean.sim.size<-mean(avg.sim.size)
      abline(v=pred.avg.size, col="red")
      abline(v=mean.sim.size, lty=2)
      #setwd("C:/Users/abuga/Documents/HVM rscripts/Emma HVM cleaned/Integrated code for clustered landscapes/preliminary plots/Accuracy histograms")
      dev.copy(png, paste0("N=",n.patches, "r1=",r1, "r2=",r2, "Metapop Size comparison.png"))
      dev.off()
      
      mean.transition.time<-mean(mean.transition.times)
      mean.sim.Tabsorp<-mean(times.to.absorption)
      
      out<-data.frame(n.patches,r1,r2,pred.avg.size, mean.sim.size, mean.transition.time, T.absorp, T.eq, T.QED)
      if (j==1 & k==1 & n==1){
        size.data<-out
      } else{
        size.data<-rbind(size.data, out)
      }
    }
  } #takes about 10 min for 4 patches to run!
}
#setwd("C:/Users/abuga/Documents/HVM rscripts/Emma HVM cleaned/Integrated code for clustered landscapes/data files")
write.csv(size.data, file="Metapop Size predicted vs simulated comparison data varying N.csv")

#Comparison of Predicted versus Simulated Metapop Size
#(transition time) / (time to absorption)
####################################################################################################
print(ggplot(data=size.data)
      + theme_classic()
      + geom_point(aes(x = mean.transition.time/T.absorp, y = (pred.avg.size-mean.sim.size)),
                  size=2) 
      + labs(x = "(transition time) / (time to absorption)", y = "(Predicted Avg. metapop Size at Q.E.D.) - (Simulated Avg. Metapop Size)", 
             title = "Comparison of Predicted versus Simulated Metapop Size")
      + theme(text = element_text(size=15)))
#setwd("C:/Users/abuga/Documents/HVM rscripts/Emma HVM cleaned/Integrated code for clustered landscapes/preliminary plots")
dev.copy(png,'Comparison of Predicted versus Simulated Metapop Size sample.png')
dev.off()
####################################################################################################

#Comparison of Predicted versus Simulated Metapop Size relative to total area
#(transition time) / (time to absorption) 
####################################################################################################
total.area<-2*n.patches
print(ggplot(data=size.data)
      + theme_classic()
      + geom_point(aes(x = mean.transition.time/T.absorp, y = (pred.avg.size-mean.sim.size)/total.area),
                   size=2) 
      + labs(x = "(transition time) / (time to absorption)", y = "((Predicted Avg. metapop Size at Q.E.D.) - (Simulated Avg. Metapop Size)) / (Total Area)", 
             title = "Comparison of Predicted versus Simulated Metapop Size")
      + theme(text = element_text(size=15)))
#setwd("C:/Users/abuga/Documents/HVM rscripts/Emma HVM cleaned/Integrated code for clustered landscapes/preliminary plots")
dev.copy(png,'Comparison of Predicted versus Simulated Metapop Size Relative to Total Area sample.png')
dev.off()
########################################################################################################################


#Comparison of Predicted versus Simulated Metapop Size relative to total area
#time till change in configuration (CTMC)
####################################################################################################
total.area<-2*n.patches
print(ggplot(data=size.data)
      + theme_classic()
      + geom_point(aes(x = log(mean.transition.time), y = (pred.avg.size-mean.sim.size)/total.area),
                   size=2) 
      + labs(x = "ln(time till change in configuration (CTMC))", y = "((Predicted Avg. metapop Size at Q.E.D.) - (Simulated Avg. Metapop Size)) / (Total Area)", 
             title = "Comparison of Predicted versus Simulated Metapop Size")
      + theme(text = element_text(size=15)))
#setwd("C:/Users/abuga/Documents/HVM rscripts/Emma HVM cleaned/Integrated code for clustered landscapes/preliminary plots")
dev.copy(png,'Predicted versus Simulated Metapop Size vs ln Tswitch sample.png')
dev.off()
#####################################################################################################

#Comparison of Predicted versus Simulated Metapop Size relative to total area
#time to absorption (CTMC)
####################################################################################################
total.area<-2*n.patches
print(ggplot(data=size.data)
      + theme_classic()
      + geom_point(aes(x = log(T.absorp), y = (pred.avg.size-mean.sim.size)/total.area),
                   size=2) 
      + labs(x = "ln(time to absorption (CTMC))", y = "((Predicted Avg. metapop Size at Q.E.D.) - (Simulated Avg. Metapop Size)) / (Total Area)", 
             title = "Comparison of Predicted versus Simulated Metapop Size")
      + theme(text = element_text(size=15)))
#setwd("C:/Users/abuga/Documents/HVM rscripts/Emma HVM cleaned/Integrated code for clustered landscapes/preliminary plots")
dev.copy(png,'Predicted versus Simulated Metapop Size vs ln Tabsorption sample.png')
dev.off()
####################################################################################################

#Comparison of Predicted versus Simulated Metapop Size relative to total area
#time to Eq (SRLM)
####################################################################################################
total.area<-2*n.patches
print(ggplot(data=size.data)
      + theme_classic()
      + geom_point(aes(x = log(T.eq), y = (pred.avg.size-mean.sim.size)/total.area),
                   size=2) 
      + labs(x = "ln(time to Eq (SRLM))", y = "((Predicted Avg. metapop Size at Q.E.D.) - (Simulated Avg. Metapop Size)) / (Total Area)", 
             title = "Comparison of Predicted versus Simulated Metapop Size")
      + theme(text = element_text(size=15)))
#setwd("C:/Users/abuga/Documents/HVM rscripts/Emma HVM cleaned/Integrated code for clustered landscapes/preliminary plots")
dev.copy(png,'Predicted versus Simulated Metapop Size vs ln Teq sample.png')
dev.off()
####################################################################################################

#Comparison of Predicted versus Simulated Metapop Size relative to total area
#time to Q.E.D. (CTMC)
####################################################################################################
total.area<-2*n.patches
print(ggplot(data=size.data)
      + theme_classic()
      + geom_point(aes(x = log(T.QED), y = (pred.avg.size-mean.sim.size)/total.area),
                   size=2) 
      + labs(x = "ln(time to Q.E.D. (CTMC))", y = "((Predicted Avg. metapop Size at Q.E.D.) - (Simulated Avg. Metapop Size)) / (Total Area)", 
             title = "Comparison of Predicted versus Simulated Metapop Size")
      + theme(text = element_text(size=15)))
#setwd("C:/Users/abuga/Documents/HVM rscripts/Emma HVM cleaned/Integrated code for clustered landscapes/preliminary plots")
dev.copy(png,'Predicted versus Simulated Metapop Size vs ln Tqed sample.png')
dev.off()
####################################################################################################
