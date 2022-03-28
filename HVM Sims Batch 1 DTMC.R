rm(list=ls()) #clear the workspace
library("ggplot2")
library("pracma")
library("deSolve")
library("RColorBrewer")
library("colorspace")
library("wesanderson") #SUPER NICE COLOUR SCHEMES!
library("ggthemes")
library("R.utils")
source("quasi_eq_func_w_scipy_and_DTMC.r")
source("DTMC_func.r")
source("SRLM_ODE_func2.r")
source("at_equilibrium_rootfunc.r")
source("dispersal_kernel_func.r")
source("DTMC_and_SRLM_func_w_transition_limit.r")
source("lambda_M_func.r")
source("lm and pstar at QED function.r")
source("create_landscape_func.r")
source("Pstar Function.r")
getwd()


r.offs<-c(1/100000,1/10000,1/1000,1/100,1/10)
r.ons<-c(1/100,1/10,1)
alphas<-c(10,1/10,1/100) #avg. disp = 1/10 interpatch distance, the interpatch distance and global disp

landscape<-create.landscape(n.patches=10, landscape.type="linear", landscape.limit=100, 
                            patch.distribution="uniform", areas.distribution="uniform", areas.limit=1, 
                            clustering.iters=0)
print(ggplot(landscape, aes(x.coord,y.coord)) + theme_classic() + geom_point(aes(size = areas)))
n.patches<-length(landscape$patch.ID) #for future use

for (i in 1:length(r.offs)){
  for (j in 1:length(r.ons)){
    for (k in 1:length(alphas)){
      
      r.off<-r.offs[i]
      r.on<-r.ons[j]
      
      e.rate<-1#0.1 #can pick just any number to start since we will just rescale this
      c.rate<-1#4  #can pick just any number to start since we will just rescale this
      alpha<-alphas[k] #global dispersal 1/landscape.limit
      self.rec<-1
      gamma<-0
      epsilon<-n.patches
      
      #first let's scale our e.rate and c.rate to ensure a high enough persistence capacity the
      #metapopulation has the opportunity to experience some growth
      Initial.Lm<-20
      
      
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
      
      stuff<-quasi.eq(n.patches, r.on, r.off)
      P.mat<-stuff$P.mat
      output<-DTMC.SRLM(total.trans=200, landscape=landscape, e.rate=e.rate, c.rate=c.rate, alpha=alpha,
                        gamma=0, epsilon=n.patches, self.rec=1, r.on=r.on, r.off=r.off, P.mat=P.mat) 
      simu.data<-output[[1]]
      configs.data<-output[[2]]
      write.csv(simu.data, paste0("Simudata Combined DTMC and SRLM simulation sample lmQED100 d", delta, 
                                  " ron", r.on, " roff", r.off, " a", alpha,  ".csv"))
      write.csv(configs.data, paste0("configdata Combined DTMC and SRLM simulation sample lmQED100 d", delta, 
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
      
      #figure
      #plot(time, metapop.size)
      print(ggplot() + theme_classic()
            + geom_line(aes(x = time, y = metapop.size),size=1, color="dodgerblue3") 
            + geom_col(aes(x = taus, y = n.on), fill="dodgerblue3", alpha=0.5)
            + labs(x = "time", y = "Metapopulation Size", 
                   title = "Combined DTMC and SRLM Simulation")
            + theme(text = element_text(size=15))
            + scale_y_continuous(sec.axis = sec_axis(~.*1, name = "Number of Habitable Patches"))
            + geom_hline(yintercept=exp.QED.size, linetype="dashed", color = "red")
            + geom_hline(yintercept=avg.size, linetype="dashed", color = "black"))
      #setwd("C:/Users/abuga/Documents/HVM rscripts/Emma HVM cleaned/Integrated code for clustered landscapes/preliminary plots")
      dev.copy(png, paste0("Combined DTMC and SRLM simulation sample lmQED100 d", delta, 
                           " ron", r.on, " roff", r.off, " a", alpha,  ".png"))
      dev.off()
    }
  }
}