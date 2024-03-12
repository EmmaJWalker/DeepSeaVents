rm(list=ls()) #clear the workspace

#Summary Data & Plots:
setwd("I:/one rep 30 patches")
library(tidyverse)
library(data.table) #gives the efficient rbindlist function
library(pracma) #gives trapz function
#library(ggplot2)
library("cowplot") #to get ggsave2 function
library(wesanderson)
library(ghibli)
wes_palettes$Darjeeling1
ghibli_palettes

#load in all the summary data files and put them into one data.frame
#file.names <- list.files(pattern="20 patch summary data lmQED20")
#sim.summaries <- lapply(file.names, read_csv)
#summary.data<-rbindlist(sim.summaries)
#write.csv(summary.data, "summary_data_1000reps.csv")
summary.data<-read_csv("summary_data_1000reps.csv")

head(summary.data)
plot.data<-summary.data

# get equilibrium occupancy expected using metric provided by dewoody et al. 2005
#######################################################################################################
setwd("C:/Users/Administrator/Desktop/DEC 2022")
source("dewoody_pstar_func.r")
source("create_landscape_func.r")
source("dispersal_kernel_func.r")
setwd("I:/one rep 30 patches")
n.patches<-20
landscape.type<-"linear"
landscape.limit<-100
patch.distribution<-"uniform"
areas.distribution<-"uniform"
areas.limit<-1
clustering.iters<-0
landscape<-create.landscape(n.patches=n.patches, landscape.type=landscape.type, 
                            landscape.limit=landscape.limit, 
                            patch.distribution=patch.distribution, 
                            areas.distribution=areas.distribution, 
                            areas.limit=areas.limit, 
                            clustering.iters=clustering.iters)
self.rec<-1
gamma<-0
epsilon<-n.patches
iterations<-1000
plot.data$dewoody.p.star<-rep(NA, nrow(plot.data))
for (i in 1:nrow(plot.data)){
  plot.data$dewoody.p.star[i]<-sum(dewoody.pstar(landscape=landscape, alpha=plot.data$alpha[i], e.rate=plot.data$delta[i], c.rate=1, 
                self.rec=self.rec, gamma=gamma, epsilon=epsilon, iterations=iterations, 
                r.on=plot.data$r.on[i], r.off=plot.data$r.off[i]))
}
#######################################################################################################

# data manipulation
########################################################################################################
plot.data$r.off<-as.factor(log10(plot.data$r.off))
plot.data$r.on<-as.factor(log10(plot.data$r.on))
plot.data$alpha<-factor(plot.data$alpha, 
                        levels=c("0.01", "0.1", "10"),
                        labels=c("Global Dispersal", "Stepping Stone", "1/10th Stepping Stone"))
plot.data$percent.occupied<-plot.data$avg.size/10*100
########################################################################################################

# grouped boxplot avg. metapop size
########################################################################################################
ggplot(plot.data) + 
  geom_boxplot(aes(x=r.on, y=dewoody.p.star/10*100, fill=r.off), col="blue", alpha=0.75) + 
  geom_boxplot(aes(x=r.on, y=exp.QED.size/10*100, fill=r.off), col="red", alpha=0.75) +
  geom_boxplot(aes(x=r.on, y=percent.occupied, fill=r.off), col="black", alpha=0.75) + 
  theme_classic() + scale_fill_manual(values=ghibli_palettes$MononokeMedium[-1]) +
  labs(title = "Expected versus Average Simulated Metapopulation Sizes Across Disturbance Regimes",
       subtitle = "(measured in % habitat occupied over 10,000 simulated timesteps)",
       y = "Average % Habitat Occupied", x = "Log Rate of Recovery") + 
  guides(fill=guide_legend("Log Rate of Disturbance")) +
  facet_grid(~ alpha)
ggsave2("20 patch Expected versus Average Simulated Metapopulation Sizes Across Disturbance Regimes.jpg", 
        height=4.5, width=10, units="in", dpi=800)
########################################################################################################

# grouped boxplot avg. number of habitat patches
########################################################################################################
ggplot(plot.data) + 
  #geom_boxplot(aes(x=r.on, y=exp.n.on, fill=r.off), col="red", alpha=0.75) +
  geom_boxplot(aes(x=r.on, y=avg.n.on, fill=r.off), col="black", alpha=0.9) + 
  theme_classic() + scale_fill_manual(values=ghibli_palettes$MarnieMedium2[-2]) +
  labs(title = "Expected versus Average Amount of Habitat Across Disturbance Regimes",
       subtitle = "(measured in number of habitatable patches over 10,000 simulated timesteps)",
       y = "Average Number of Habitatable Patches", x = "Log Rate of Recovery") + 
  guides(fill=guide_legend("Log Rate of Disturbance"))
########################################################################################################
ggsave2("20 patch Expected versus Average Amount of Habitat Across Disturbance Regimes.jpg", 
        height=4.5, width=10, units="in", dpi=800)

#RE-DEFINING ALL THE EXPERIMENT PARAMETERS THAT WERE RUN: 
###############################################################################
time.limit<-10000
r.offs<-c(1/100000, 1/10000,1/1000,1/100,1/10)
r.ons<-c(1/100,1/10,1)
alphas<-c(10,1/10,1/100) #avg. disp = 1/10 interpatch distance, the interpatch 
#distance and global disp
n.patches<-10
landscape.type<-"linear"
landscape.limit<-100
patch.distribution<-"uniform"
areas.distribution<-"uniform"
areas.limit<-1
clustering.iters<-0
landscape<-create.landscape(n.patches=n.patches, landscape.type=landscape.type, 
                            landscape.limit=landscape.limit, 
                            patch.distribution=patch.distribution, 
                            areas.distribution=areas.distribution, 
                            areas.limit=areas.limit, 
                            clustering.iters=clustering.iters)
n<-n.patches
###############################################################################


# preliminary plots to look at transient dynamics
########################################################################################################
setwd("C:/Users/Administrator/Desktop/DEC 2022")
source("PStar_func.r")
setwd("H:/1000 reps")
head(summary.data)

#pick a random iteration number (replicate)
#rx<-sample(summary.data$it.no, 1)
#rx
rx<-1 #cluster that generated the data
pcombos<-c(1:8,10:13,15:28,30:45) #parameter combos for 9,14,29 had no change in habitat over 10000 years!
rx<-1 #cluster that generated the data
pcombos<-c(9,14,29)
for (w in pcombos){
    w=1
    print(w)
    rx.summary.data<-summary.data%>%
      filter(it.no==rx)
    r.on<-rx.summary.data$r.on[w]
    r.off<-rx.summary.data$r.off[w]
    alpha<-rx.summary.data$alpha[w]
    delta<-rx.summary.data$delta[w]
    Initial.Lm<-rx.summary.data$Initial.Lm[w]
    exp.QED.size<-rx.summary.data$exp.QED.size[w]
    
    configs.data<-read.csv(paste0("configdata Combined CTMC and SRLM simulation sample lmQED",Initial.Lm,
                                  " ron",r.on," roff",r.off," a",alpha,"d",delta,"_",rx,".csv"))
    simu.data<-read.csv(paste0("Simudata Combined CTMC and SRLM simulation sample lmQED",Initial.Lm,
                               " ron",r.on," roff",r.off," a",alpha,"d",delta,"_",rx,".csv"))
    
    #if (length(t(as.matrix(configs.data==(2*n.patches+2))))){ #to be able to plot the cases where no habitat changes occured in 10000 years
    #  configs.data<-t(as.matrix(configs.data))
    #  configs.data<-configs.data[-1,]
    #  configs.data<-data.frame(rbind(configs.data, configs.data))
    #  configs.data[2,1]<-10000
    #  just.configs.data<-configs.data[,-1]
    #  simu.data<-simu.data[,-1]
    #}else{
      just.configs.data<-configs.data[,-(1:2)]
      configs.data<-configs.data[,-1]
      simu.data<-simu.data[,-1]
    #}
    #configs.data<-read.csv("configdata Combined CTMC and SRLM simulation sample lmQED20 ron1 roff0.001 a0.01d0.370311750341306_665.csv")
    #simu.data<-read.csv("Simudata Combined CTMC and SRLM simulation sample lmQED20 ron1 roff0.001 a0.01d0.370311750341306_665.csv")
    
    #Simulation Figure:
    ########################################################################################################
    ##weight the occupancy values by area and plot metapop size over time
    #area.weights<-landscape$areas/sum(landscape$areas)
    #metapop.size<-sim.data[,-1]*area.weights *what Austin has, which gives the size of each patch 
    #but we don't wanna plot this for all of them, I'm just going to plot the sum
    metapop.size<-rowSums(simu.data[,2:(n.patches+1)]) #*area.weights) <-if we wanna add the area weights
    time<-simu.data$time
    n.on<-rowSums(configs.data[,2:(n.patches+1)])
    taus<-configs.data[,1]
    if (taus[length(taus)]<time[length(time)]){
      #  #so the bar of the last habitat configuration actually shows up if it ended from all habitat turning off
      #n.on<-c(n.on, n.on[length(n.on)]) 
      taus<-c(taus, time[length(time)])
    }
    
    #average simulated metaop size over time calculated using trapezoidal method (As used by Austin)
    avg.size=trapz(time,metapop.size)/time[length(time)]
    avg.size
    
    #figure
    plot(time, metapop.size)
    metapop.data<-data.frame(time, metapop.size)
    tausmin<-taus[-length(taus)]
    tausmax<-lead(taus)
    tausmax<-tausmax[-length(tausmax)]
    habitat.data<-data.frame(tausmin,tausmax,n.on)
    print(ggplot() + theme_classic()
          + geom_rect(aes(xmin = tausmin, xmax = tausmax, 
                          ymin = 0, ymax = n.on), 
                      fill = ghibli_palettes$MarnieMedium2[7], alpha=0.5)
          + geom_step(data = metapop.data, aes(x = time, y = metapop.size),size=1, color=ghibli_palettes$MononokeMedium[3], alpha=0.7) 
          + labs(x = "time", y = "Metapopulation Size", 
                 title = "Combined CTMC and SRLM Simulation")
          + theme(text = element_text(size=15))
          + scale_y_continuous(sec.axis = sec_axis(~.*1, name = "Number of Habitable Patches"))
          + geom_hline(yintercept=exp.QED.size, linetype="solid", size=1, color = "red")
          + geom_hline(yintercept=avg.size, linetype="solid", size=1, color = ghibli_palettes$MononokeMedium[1])
          + xlim(0,taus[length(taus)])
          )
    
    ggsave2(paste0("Combined CTMC and SRLM simulation sample lmQED20", 
                   " ron", r.on, " roff", r.off, " a", alpha,  " d", delta, ".jpeg"), 
            height=5, width=10, units="in", dpi=800)
    
    #Transients Figure:
    ########################################################################################################
    
    e.rate<-delta
    c.rate<-1
    config.pstar<-rep(NA, nrow(configs.data))
    for (i in 1:nrow(configs.data)){
      landscape.config<-landscape[just.configs.data[i,]!=0,]
      config.pstar[i]<-sum(pstar.function(landscape.config, alpha, e.rate, c.rate, self.rec, gamma, epsilon, iterations=1000))
    }
    configs.data$config.pstar<-config.pstar
    y<-c()
    for (i in 1:nrow(configs.data)){
      start.t<-configs.data[i,1]
      end.t<-configs.data[i+1,1]
      if(i==1){
        xend.t<-10001 #arbitrary number that start.t can't initially take on
      }
      xend.t==start.t
      if(xend.t==start.t){
        start.t<-start.t+0.000001 #this prevents double counting when xend.t=start.t
      }
      x<-simu.data[start.t<=simu.data$time & simu.data$time<end.t,]
      xend.t<-x$time[length(x$time)]
      x$p.star<-rep(config.pstar[i], nrow(x))
      z<-length(x$p.star)
      y<-c(y, x$p.star)
    }
    #length(y)
    #nrow(simu.data)
    simu.data$pstar<-y
    simu.data$p<-rowSums(simu.data[,2:(n.patches+1)]) #*area.weights) <-if we wanna add area weights
    simu.data$transient<-simu.data$p/simu.data$pstar
    head(simu.data)
    
    print(ggplot(simu.data) + 
      #geom_point(aes(x=time, y=transient), pch=16, size=2, colour=ghibli_palettes$MononokeMedium[1], alpha=0.5) + 
      geom_step(aes(x=time, y=transient), size=1, colour=ghibli_palettes$MononokeMedium[1], alpha=0.5) +
      geom_hline(aes(yintercept=1), size=1, col="red") +
      ylim(0,max(simu.data$transient)+0.02) +
      theme_classic() +
      labs(y = "P/P* Within Each Habitat Configuration", x = "Time"))
    
    
    ggsave2(paste0("Transients sample lmQED20", 
                   " ron", r.on, " roff", r.off, " a", alpha,  " d", delta, ".jpeg"), 
            height=5, width=10, units="in", dpi=800)
    
    ########################################################################################################
    
    
    #library(TSA)
    #periodogram(simu.data$transient)
    
  
}
configs.data
simu.data
