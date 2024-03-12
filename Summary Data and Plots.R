rm(list=ls()) #clear the workspace

#Summary Data & Plots:
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
library("gtools") #to get even() function for IDsToConfigs
library("ivmte") #to get the magnitude() of a number
library("Bolstad2")
setwd("I:/600 reps 2023") #current figure sims norm(x1)<1e-4
#setwd("H:/2023_07_03_reps")  #new sims norm(x1)<1e-10
library(tidyverse)
library(data.table) #gives the efficient rbindlist function
library(pracma) #gives trapz function
#library(ggplot2)
library("cowplot") #to get ggsave2 function
library(wesanderson)
library(ghibli)
wes_palettes$Darjeeling2
ghibli_palettes

#######################################################################################################
#load in all the summary data files and put them into one data.frame
#takes about 30 min - 1 hr for ~10000 files
#file.names <- list.files(pattern="summary data lmQED20")
#sim.summaries <- lapply(file.names, read_csv)
#summary.data<-rbindlist(sim.summaries)
#write.csv(summary.data, "summary_data_600reps.csv")
#######################################################################################################

#####################################################################################################
#setwd("H:/2023_07_03_reps")  #new sims norm(x1)<1e-10
summary.data<-read_csv("summary_data_600reps.csv")
#q<-group_by(summary.data, r.on, r.off, alpha) %>%
#  tally(is.na(avg.size))
#q #NA's indicating a problem with the ODE solver only happen for r.on=0.01, r.off=0.1, and alpha=10
##Exclude this parameter combination from my results
#summary.data<-na.omit(summary.data)
####################################################################################################

#get number of replicates to check:
q<-group_by(summary.data, r.on, r.off, alpha)
tally(q)
#want only the first X many replicates for plots (since unequal numbers of replicates since I stopped my replicate sims early)
X<-200
r.offs<-c(1/100000,1/10000,1/1000,1/100,1/10)
r.ons<-c(1/100,1/10,1)
alphas<-c(1/100, 1/10,10)
for (j in 1:length(r.offs)){
  for (k in 1:length(r.ons)){
    for (l in 1:length(alphas)){
      sub.data<-filter(summary.data, r.on==r.ons[k] & r.off==r.offs[j] & alpha==alphas[l])
      sub.data<-sub.data[1:X,]
      if(j==1 & k==1 & l==1){summary.data2<-sub.data
      } else {summary.data2<-rbind(summary.data2, sub.data)}
    }
  }
}
head(summary.data2)
q<-group_by(summary.data2, r.on, r.off, alpha)
tally(q)
plot.data<-summary.data2

# get equilibrium occupancy expected using metric provided by dewoody et al. 2005
#######################################################################################################
setwd("C:/Users/Administrator/Desktop/DEC 2022")
source("dewoody_pstar_func.r")
source("create_landscape_func.r")
source("dispersal_kernel_func.r")
#source("QED_metrics_func.r")
source("QED_metrics_func_w_median_n.r")
source("dispersal_kernel_func.r") #disp.kernel()
source("lambda_M_func.r") #get.lambda.M()
source("Pstar_func.r") #pstar.function()
source("quasi_eq_sparse_func.r") #quasi.eq()
setwd("I:/600 reps 2023") #current figure sims norm(x1)<1e-4
#setwd("H:/2023_07_03_reps")  #new sims norm(x1)<1e-10
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
#to get expected number of patches on at QED in sims and add that to our data: 
plot.data$exp.n.on<-rep(NA, nrow(plot.data))
plot.data$T.absorp<-rep(NA, nrow(plot.data))
plot.data$rate.comparison<-rep(NA, nrow(plot.data))
r.offs<-c(1/100000,1/10000,1/1000,1/100,1/10)
r.ons<-c(1/100,1/10,1)
alphas<-c(0.01,0.1,10)
for (j in 1:length(r.offs)){
  for (k in 1:length(r.ons)){
    for (l in 1:length(alphas)){
      delta<-plot.data$delta[(plot.data$r.off==r.offs[j]&plot.data$r.on==r.ons[k]&plot.data$alpha==alphas[l])]
      delta<-delta[1]
      x<-QED.metrics(landscape=landscape, e.rate=delta, c.rate=1, gamma=gamma, epsilon=epsilon, self.rec=self.rec, alpha=alphas[l],r.on=r.ons[k], r.off=r.offs[j])
      s<-nrow(plot.data[(plot.data$r.off==r.offs[j]&plot.data$r.on==r.ons[k]&plot.data$alpha==alphas[l]),])
      plot.data$exp.QED.size[(plot.data$r.off==r.offs[j]&plot.data$r.on==r.ons[k]&plot.data$alpha==alphas[l])]<-rep(x[[3]],s)
      }
      y<-quasi.eq(n.patches=10, r.on=r.ons[k], r.off=r.offs[j])
      r<-nrow(plot.data[(plot.data$r.off==r.offs[j]&plot.data$r.on==r.ons[k]),])
      plot.data$exp.n.on[(plot.data$r.off==r.offs[j]&plot.data$r.on==r.ons[k])]<-rep(x[[5]],r)
      plot.data$T.absorp[(plot.data$r.off==r.offs[j]&plot.data$r.on==r.ons[k])]<-rep(y[[2]],r)
      plot.data$rate.comparison[(plot.data$r.off==r.offs[j]&plot.data$r.on==r.ons[k])]<-rep(y[[3]],r)
      plot.data$med.exp.n.on[(plot.data$r.off==r.offs[j]&plot.data$r.on==r.ons[k])]<-rep(x[[9]],r)
      plot.data$QED.T.absorp[(plot.data$r.off==r.offs[j]&plot.data$r.on==r.ons[k])]<-rep(x[[10]],r)
  }}
#takes awhile :/
#write.csv(plot.data, "plot_data_600reps2.csv")
plot.data<-read_csv("plot_data_600reps2.csv")
#######################################################################################################
#Looking at times to absorption
p<-plot.data%>%
  select(r.on, r.off, QED.T.absorp)%>%
  group_by(r.on, r.off)%>%
  reframe(unique(QED.T.absorp))
#write.csv(p, "Time to absorption table.csv")

# data manipulation
########################################################################################################
plot.data$r.off<-as.factor(log10(plot.data$r.off))
plot.data$r.on<-as.factor(log10(plot.data$r.on))
plot.data$alpha<-factor(plot.data$alpha, 
                        levels=c("0.01", "0.1", "10"),
                        labels=c("Global Dispersal", "Stepping Stone", "1/10th Stepping Stone"))
plot.data$percent.occupied<-plot.data$avg.size/10*100
########################################################################################################

#NEED TO RECALCULATE AVG.N.ON AT QED AND PLOT HABITAT FIGURE, also calculating variance in p versus p* in sims and p*
########################################################################################################
r.offs<-c(1/100000, 1/10000,1/1000,1/100,1/10)
r.ons<-c(1/100,1/10,1)
alphas<-c(10,1/10,1/100)
Initial.Lm<-20
#calculating the geometric avg. number of habitat patches on to compare against our expected geometric mean ammount of habitat
##need to recalculate avg.n.on in sims so doesn't include time spent in no-habitat state:
##reloading individual sims and recalculating the avg avg.n.on across sims
z<-1
new.avgs<-c()
geomavg.n.on<-c()
geom.avg.size<-c()
avg.size<-c()
trapz.avg.size<-c()
simp.avg.size<-c()
new.rons<-c()
new.roffs<-c()
new.alphas<-c()
p.star.vars<-c()
p.vars<-c()
for (j in 1:length(r.offs)){
  for (k in 1:length(r.ons)){
    for(l in 1:length(alphas)){
      #j<-1
      #k<-1
      #l<-1
      r.on<-r.ons[k] 
      r.off<-r.offs[j]
      alpha<-alphas[l]
      print(c(r.on, r.off, alpha))
      config.files <- list.files(pattern=paste0("configdata Combined CTMC and SRLM simulation sample lmQED", Initial.Lm, " ron", 
                                                r.on, " roff", r.off, " a", alpha,"*"))
      simu.files<-list.files(pattern=paste0("Simudata Combined CTMC and SRLM simulation sample lmQED", Initial.Lm,
                                            " ron",r.on," roff",r.off," a",alpha,"*"))
      configs.n.pstars.files<-list.files(pattern=paste0("configs.n.pstars",Initial.Lm,
                                                        " ron",r.on," roff",r.off," a",alpha,"*"))
      #n<-length(config.files)
      n<-200 #to use just the first 200 files
      for (i in 1:n) {
        #i<-1
        print(i)
        configs.data<-read.csv(config.files[i])
        #dealing with cases where no change in habitat occured so that it's in a dataframe with
        #the same format as in simulations where changes occured
        #########################################################################################
        if (ncol(configs.data)==2){
          configs.data<-t(configs.data)
          configs.data<-configs.data[-1,]
          for (N in 1:(n.patches+1)){
            if(N==1){configs.data2<-as.data.frame(configs.data[N], colnames=paste0("V", N))}else{
              configs.data2<-cbind(configs.data2,as.data.frame(configs.data[N]))
            }
            names(configs.data2)[N]<-paste0("V",N)
          }
          configs.data<-configs.data2
        } else{
          configs.data<-configs.data[,-1]
        }
        #########################################################################################
        simu.data<-read.csv(simu.files[i])
        configs.n.pstars<-read.csv(configs.n.pstars.files[i])
        #give the configs data and pstar data common names to merge by
        for (N in 1:n.patches){
          names(configs.n.pstars)[names(configs.n.pstars) == paste0("X",N)] <- paste0("n",N)
          names(configs.data)[names(configs.data) == paste0("V",N+1)] <- paste0("n",N)
          names(configs.n.pstars)[names(configs.n.pstars) == paste0("X",N,".1")] <- paste0("pstar",N)
        }
        pstars.data<-merge(configs.data[-1], configs.n.pstars, by=paste0("n",1:n.patches))
        pstar.sizes<-rowSums(pstars.data[,(n.patches+2):(n.patches*2+1)])
        if(length(pstar.sizes)==1){
          p.star.var<-0
        }else{
          p.star.var<-var(pstar.sizes)
        }
        metapop.sizes<-rowSums(simu.data[,3:n.patches+2])
        p.var<-var(metapop.sizes)
        #getting the avg. number of habitat patches on over time
        n.on<-rowSums(configs.data[,2:(n.patches+1)])
        metapop.size<-rowSums(simu.data[,3:(n.patches+2)])
        lead.time<-simu.data$time[-1]
        lag.time<-simu.data$time[-length(simu.data$time)]
        #arithmetic avg. size using the riemann sum formula (rectangular area)
        avg.size[z]<-sum(metapop.size[-length(metapop.size)]*(lead.time-lag.time))/lead.time[length(lead.time)]
        trapz.avg.size[z]<-trapz(simu.data$time, metapop.size)/lead.time[length(lead.time)]
        tiny.vec<-seq(0:0.0000000000000000000000000000000000000000000000001,length.out=length(metapop.size))
        unique.time<-simu.data$time+tiny.vec
        simp.avg.size[z]<-sintegral(unique.time, metapop.size, n.pts=length(metapop.size)*10)$int/(simu.data$time[length(simu.data$time)])
        if(any(metapop.size==0)){
          geom.avg.size[z]<-0
        }else{
          geom.avg.size[z]<-exp(sum(log(metapop.size[-length(metapop.size)])*(lead.time-lag.time))/lead.time[length(lead.time)])
        }
        ts<-configs.data[,1]
        if (length(ts)==1){ #if no habitat change occured
          #average amount of habitat = the amount of habitat
          avg.n.on<-n.on[1]
          geomavg.n.on[z]<-n.on[1]
        }else{
          lag.ts<-ts
          lead.ts<-c(ts[-1],simu.data$time[length(simu.data$time)])
          #reimann sum formula (rectangular area)
          avg.n.on<-sum(n.on*(lead.ts-lag.ts))/lead.ts[length(lead.ts)] #time averaged mean number of patches on in a simulation
          #using the reimann sum formula (rectangular areas under a function)
          #time averaged geometric mean:
          geomavg.n.on[z]<-exp(sum(log(n.on)*(lead.ts-lag.ts))/lead.ts[length(lead.ts)])
        }
        new.avgs[z]<-avg.n.on
        new.rons[z]<-r.on
        new.roffs[z]<-r.off
        new.alphas[z]<-alpha
        p.star.vars[z]<-p.star.var
        p.vars[z]<-p.var
        z<-z+1
      }
    }
  }
}
new.avgs<-data.frame(new.rons, new.roffs, geomavg.n.on, new.alphas, p.star.vars, p.vars, geom.avg.size, avg.size, trapz.avg.size, simp.avg.size)
#write.csv(new.avgs, "newavgs.csv")
new.avgs<-read.csv("newavgs.csv")
new.avgs$r.off<-as.factor(log10(new.avgs$new.roffs))
new.avgs$r.on<-as.factor(log10(new.avgs$new.rons))
new.avgs$alpha<-factor(new.avgs$new.alphas, 
                       levels=c("0.01", "0.1", "10"),
                       labels=c("Global Dispersal", "Stepping Stone", "1/10th Stepping Stone"))
table1p2<-new.avgs%>%
  group_by(r.off, r.on, alpha)%>%
  summarise(mean.pstar.vars=mean(p.star.vars), mean.p.vars=mean(p.vars), mean.size=mean(geom.avg.size))


pal<-c("#EC4E20", "#655A7C","#A4392A", "#AC8B41", "#E58507") #custom pallete "#462521","#D74E09", "#353531"
#"#016FB9"
################################################################################################ w/out p* dewoody
#p.data<-filter(plot.data, alpha=="Global Dispersal")
#ggplot() + 
#  geom_boxplot(data=p2.data, aes(x=r.on, y=geom.avg.size, fill=r.off, col=r.off), alpha=0.25) + 
#  #dummy line just to get the symbol for the legend:
#  geom_boxplot(data=p.data, aes(x=r.on, y=exp.QED.size, fill=r.off, col="p*_QED"), alpha=1) +
#  geom_line(data=p.data, aes(x=r.on, y=exp.QED.size, linetype="p*_QED"), alpha=0)+
#  scale_color_manual(labels=rev(c("","","", "","", "p*_QED")), breaks=rev(c("-5","-4","-3", "-2","-1", "p*_QED")), values = rev(c(pal, "black"))) +
#  theme_classic() + scale_fill_manual(values=pal) +
#  labs(y = "Avg. Occupancy", x = "Log Rate of Recovery") + 
#  guides(fill=guide_legend("Log Rate of Disturbance", override.aes = list(fill=pal, alpha=0.25, colour=pal)), 
#         linetype=guide_legend("Estimate", override.aes = list(linetype=1, alpha=1, colour=c("black"))),
#         colour="none")
#ggsave2("Expected versus GeomAverage Simulated Metapopulation Sizes Across Disturbance Regimes Global Dispersal relaxing rootfunc.jpg", 
#        #height=4.5, width=10, units="in", dpi=800)
#        height=4.5, width=10, units="in", dpi=800)
######################################################################################################## w p* dewoody
p.data<-filter(plot.data, alpha=="Global Dispersal")
p2.data<-new.avgs[new.avgs$new.alphas==1/100,]
ggplot(p.data) + 
  geom_boxplot(data=p2.data, aes(x=r.on, y=trapz.avg.size, fill=r.off, col=r.off), alpha=0.25) + 
  geom_boxplot(data=p.data, aes(x=r.on, y=dewoody.p.star, fill=r.off, col="p*_DeWoody"), alpha=1) + 
  #dummy line just to get the symbol for the legend:
  geom_line(data=p.data, aes(x=r.on,y=dewoody.p.star, linetype="p*_DeWoody"), alpha=0)+
  geom_boxplot(data=p.data, aes(x=r.on, y=exp.QED.size, fill=r.off, col="p*_QED"), alpha=1) +
  geom_line(data=p.data, aes(x=r.on, y=exp.QED.size, linetype="p*_QED"), alpha=0)+
  scale_color_manual(labels=rev(c("","","", "","","p*_DeWoody", "p*_QED")), breaks=rev(c("-5","-4","-3", "-2","-1","p*_DeWoody", "p*_QED")), values = rev(c(pal, "#016FB9", "black"))) +
  theme_classic() + scale_fill_manual(values=pal) +
  labs(y = "Avg. Occupancy", x = "Log Rate of Recovery") + 
guides(fill=guide_legend("Log Rate of Disturbance", override.aes = list(fill=pal, alpha=0.25, colour=pal)), 
       linetype=guide_legend("Estimate", override.aes = list(linetype=1, alpha=1, colour=c("#016FB9", "black"))),
       colour="none")
  #facet_grid(~ alpha)
ggsave2("Expected versus Average Simulated Metapopulation Sizes Across Disturbance Regimes Global Dispersal only.jpg", 
        #height=4.5, width=10, units="in", dpi=800)
        height=4.5, width=10, units="in", dpi=800)
########################################################################################################
p.data<-filter(plot.data, alpha=="Stepping Stone")
p2.data<-new.avgs[new.avgs$new.alphas==1/10,]
ggplot(p.data) + 
  geom_boxplot(data=p2.data, aes(x=r.on, y=trapz.avg.size, fill=r.off, col=r.off), alpha=0.25) + 
  geom_boxplot(data=p.data, aes(x=r.on, y=dewoody.p.star, fill=r.off, col="p*_DeWoody"), alpha=1) + 
  #dummy line just to get the symbol for the legend:
  geom_line(data=p.data, aes(x=r.on,y=dewoody.p.star, linetype="p*_DeWoody"), alpha=0)+
  geom_boxplot(data=p.data, aes(x=r.on, y=exp.QED.size, fill=r.off, col="p*_QED"), alpha=1) +
  geom_line(data=p.data, aes(x=r.on, y=exp.QED.size, linetype="p*_QED"), alpha=0)+
  scale_color_manual(labels=rev(c("","","", "","","p*_DeWoody", "p*_QED")), breaks=rev(c("-5","-4","-3", "-2","-1","p*_DeWoody", "p*_QED")), values = rev(c(pal, "#016FB9", "black"))) +
  theme_classic() + scale_fill_manual(values=pal) +
  labs(y = "Avg. Occupancy", x = "Log Rate of Recovery") + 
  guides(fill=guide_legend("Log Rate of Disturbance", override.aes = list(fill=pal, alpha=0.25, colour=pal)), 
         linetype=guide_legend("Estimate", override.aes = list(linetype=1, alpha=1, colour=c("#016FB9", "black"))),
         colour="none")
#facet_grid(~ alpha)
ggsave2("Expected versus Average Simulated Metapopulation Sizes Across Disturbance Regimes Stepping Stone only.jpg", 
        #height=4.5, width=10, units="in", dpi=800)
        height=4.5, width=10, units="in", dpi=800)
########################################################################################################
p.data<-filter(plot.data, alpha=="1/10th Stepping Stone")
p2.data<-new.avgs[new.avgs$new.alphas==10,]
ggplot(p.data) + 
  geom_boxplot(data=p2.data, aes(x=r.on, y=geom.avg.size, fill=r.off, col=r.off), alpha=0.25) + 
  geom_boxplot(data=p.data, aes(x=r.on, y=dewoody.p.star, fill=r.off, col="p*_DeWoody"), alpha=1) + 
  #dummy line just to get the symbol for the legend:
  geom_line(data=p.data, aes(x=r.on,y=dewoody.p.star, linetype="p*_DeWoody"), alpha=0)+
  geom_boxplot(data=p.data, aes(x=r.on, y=exp.QED.size, fill=r.off, col="p*_QED"), alpha=1) +
  geom_line(data=p.data, aes(x=r.on, y=exp.QED.size, linetype="p*_QED"), alpha=0)+
  scale_color_manual(labels=rev(c("","","", "","","p*_DeWoody", "p*_QED")), breaks=rev(c("-5","-4","-3", "-2","-1","p*_DeWoody", "p*_QED")), values = rev(c(pal, "#016FB9", "black"))) +
  theme_classic() + scale_fill_manual(values=pal) +
  labs(y = "Avg. Occupancy", x = "Log Rate of Recovery") + 
  guides(fill=guide_legend("Log Rate of Disturbance", override.aes = list(fill=pal, alpha=0.25, colour=pal)), 
         linetype=guide_legend("Estimate", override.aes = list(linetype=1, alpha=1, colour=c("#016FB9", "black"))),
         colour="none")
#facet_grid(~ alpha)
ggsave2("Expected versus Average Simulated Metapopulation Sizes Across Disturbance Regimes 1 10th Stepping Stone only.jpg", 
        #height=4.5, width=10, units="in", dpi=800)
        height=4.5, width=10, units="in", dpi=800)
########################################################################################################


# grouped boxplot avg. number of habitat patches
#pal<-c("#436969", "#444F30", "#858C6A", "#4F772D", "#90A955")#
#pal<-c("#4B371C", "#586523", "#A59E2E", "#648E07", "#7C985C")
pal<-c("#436969", "#4B371C", "#444F30", "#A59E2E", "#648E07")
#pal<-c("#70622D", "#2E4B6A","#585123", "#588157", "#A59E2E") #custom pallete "#042A2B","#3B252C", "#49496A", "#576CA8"
########################################################################################################
p2.data<-new.avgs[new.avgs$new.alphas==1/100,]
ggplot() + 
  geom_boxplot(data=p2.data, aes(x=r.on, y=geomavg.n.on, fill=r.off, col=r.off), alpha=0.25) +  
  geom_boxplot(data=p.data, aes(x=r.on, y=exp.n.on, fill=r.off, col="X"), alpha=0.25) +
  geom_line(data=p.data, aes(x=r.on,y=exp.n.on, linetype="Expected Avg. Number of Patches at QED"), alpha=0)+
  theme_classic() + scale_fill_manual(values=pal) + scale_color_manual(breaks=c("-5","-4","-3", "-2","-1","X"), values=c(pal,"black")) +
  labs(y = "Avg. Number of Habitatable Patches", x = "Log Rate of Recovery") + 
  guides(fill=guide_legend("Log Rate of Disturbance", override.aes = list(fill=pal, alpha=0.25, colour=pal)), 
         linetype=guide_legend("Estimate", override.aes = list(linetype=1, alpha=1, colour=c("black"))),
         colour="none")
########################################################################################################
ggsave2("Expected versus GeomAverage Amount of Habitat Across Disturbance Regimes Global Dispersal Only.jpg", 
        height=4.5, width=10, units="in", dpi=800)
        #height=9, width=20, units="in", dpi=800)

#RE-DEFINING ALL THE EXPERIMENT PARAMETERS THAT WERE RUN WITHIN INDIVIDUAL SIMULATIONS: 
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
#SINGLE SIMULATION FIGURES:
################################################################################
setwd("C:/Users/Administrator/Desktop/DEC 2022")
source("PStar_func.r")
source("lambda_M_func.r")
setwd("I:/600 reps 2023")
plot.data<-read_csv("plot_data_600reps2.csv")
head(summary.data)
#if want to pick a random iteration number (replicate)
#rx<-sample(summary.data$it.no, 1)
rx<-1 #cluster that generated the data
pcombos<-c(1:8,10:13,15:28,30:41,43:45) #parameter combos for 9,14,29 had no change in habitat over 10000 years!
pcombos<-1:45
#42 seems to also have something funny going on
pal<-c("#EC4E20", "#655A7C","#A4392A", "#AC8B41", "#E58507")
pal2<-c("#4B371C", "#586523", "#A59E2E", "#648E07", "#7C985C")
pcombos<-c(7,8,9,40,41,42,37,38,39,28,29,30) #for the example sims used in our global dispersal figures
#no change occured over 10,000 years in rx=1, w=1, use rx=2, w=1 instead for figure for illustrative purposes
#no change occured over 10,000 years in rx=1, w=45-7, use rx=2, w=45-7 instead for figure for illustrative purposes
#no change occured over 10,000 years in rx=1, w=45-43, use rx=2, w=45-43 instead for figure for illustrative purposes
#rx=2, w=45-1 is more interesting plot than with rx=1, might make a better example
for (w in pcombos){
  w=6 #to pick just one sim
  #pulling out the parameter values from the sim's data
  rx.summary.data<-summary.data%>%
    filter(it.no==rx)
  r.on<-rx.summary.data$r.on[w]
  r.off<-rx.summary.data$r.off[w]
  alpha<-rx.summary.data$alpha[w]
  delta<-rx.summary.data$delta[w]
  Initial.Lm<-rx.summary.data$Initial.Lm[w]
  #exp.QED.size<-rx.summary.data$exp.QED.size[w] #THIS NEEDS TO BE UPDATED to below:
  exp.QED.size<-unique(plot.data$exp.QED.size[plot.data$r.on==r.on & plot.data$r.off==r.off & plot.data$alpha==alpha]) #*****************************************************
  #setting the colour to plot the sim to match the colour of the boxplot it appears in
  simcolor.id=which(r.offs==r.off)
  #reading in the data for the sim
  configs.data<-read.csv(paste0("configdata Combined CTMC and SRLM simulation sample lmQED",Initial.Lm,
                                " ron",r.on," roff",r.off," a",alpha,"d",delta,"_",rx,".csv"))
  simu.data<-read.csv(paste0("Simudata Combined CTMC and SRLM simulation sample lmQED",Initial.Lm,
                             " ron",r.on," roff",r.off," a",alpha,"d",delta,"_",rx,".csv"))
  #if the simulation ended early in the last configuration not because it was the absorbing state but because the change in metapop
  #size was below 1x10-4
  if(simu.data$time[nrow(simu.data)]<10000 & sum(configs.data[nrow(configs.data),3:(n.patches+2)])>1){
    simu.data<-rbind(simu.data, simu.data[nrow(simu.data),]) #then let the metapop size at the end of 10000 be the same
    simu.data$time[nrow(simu.data)]<-10000
  }
  ##commented out code here is to be able to plot the cases where no habitat changes occured in 10000 years
  #if (length(t(as.matrix(configs.data==(2*n.patches+2))))){ 
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
  
  #calculates p*'s within each simulated config:
  e.rate<-delta
  c.rate<-1
  config.pstar<-rep(NA, nrow(configs.data))
  config.lm<-rep(NA, nrow(configs.data))
  for (i in 1:nrow(configs.data)){
    landscape.config<-landscape[just.configs.data[i,]!=0,]
    config.pstar[i]<-sum(pstar.function(landscape.config, alpha, e.rate, c.rate, self.rec, gamma, epsilon, iterations=1000))
    config.lm[i]<-get.lambda.M(landscape.config, alpha, gamma, epsilon, self.rec, e.rate, c.rate)/delta
  }
  configs.data$config.pstar<-config.pstar
  configs.data$config.lm<-config.lm
  configs.data$taus<-c(0,-(configs.data$V1[-length(configs.data$V1)]-configs.data$V1[-1]))
  configs.data
  
  
  #HISTOGRAMS: bin widths have to be manually adjusted to be most appropriate for the data
  ###########################################################################################
  #configs.data
  print(ggplot(configs.data) + theme_classic()
        + geom_histogram(aes(config.pstar), binwidth=1, fill=pal[simcolor.id], alpha=0.25)
        + geom_vline(aes(xintercept=mean(config.pstar)), color=pal[simcolor.id])#color="#B50A2AFF")
        + labs(x = "P*", y = "Frequency")
        + theme(text = element_text(size=15))
  )
  ggsave2(paste0("Pstar distribution lmQED20", 
                 " ron", r.on, " roff", r.off, " a", alpha,  " d", delta, ".jpeg"), 
          height=5, width=5, units="in", dpi=800)
  print(ggplot(configs.data) + theme_classic()
        + geom_histogram(aes(config.lm), binwidth=10, fill=pal[simcolor.id], alpha=0.25)
        + geom_vline(aes(xintercept=mean(config.lm)), color=pal[simcolor.id])#color="#B50A2AFF")
        + labs(x = "Persistence Capacity", y = "Frequency")
        + theme(text = element_text(size=15))
  )
  ggsave2(paste0("Lm distribution lmQED20", 
                 " ron", r.on, " roff", r.off, " a", alpha,  " d", delta, ".jpeg"), 
          height=5, width=5, units="in", dpi=800)
  #remove the initial condition as a time at which a habitat change occured since this is just the starting condition not a change
  configs.data2<-configs.data[-1,]
  print(ggplot(configs.data2) + theme_classic()
        + geom_histogram(aes(taus), binwidth=1, fill=pal2[simcolor.id], alpha=0.25)
        + geom_vline(aes(xintercept=mean(taus)), color=pal2[simcolor.id])#, color="#B50A2AFF")
        + labs(x = "Time Between Habitat Changes", y = "Frequency")
        + theme(text = element_text(size=15))
  )
  ggsave2(paste0("Tau distribution lmQED20", 
                 " ron", r.on, " roff", r.off, " a", alpha,  " d", delta, ".jpeg"), 
          height=5, width=5, units="in", dpi=800)
  ###########################################################################################
  
  #Simulation Figure:
  ########################################################################################################
  ##weight the occupancy values by area and plot metapop size over time
  #area.weights<-landscape$areas/sum(landscape$areas)
  #metapop.size<-sim.data[,-1]*area.weights *what Austin has, which gives the size of each patch 
  #but we don't wanna plot this for all of them, I'm just going to plot the sum
  metapop.size<-rowSums(simu.data[,2:(n.patches+1)]) #*area.weights) <-if we wanna add the area weights
  time<-simu.data$time
  n.on<-rowSums(configs.data[,2:(n.patches+1)])
  ts<-configs.data[,1]
  if (ts[length(ts)]<time[length(time)]){
    #so the bar of the last habitat configuration actually shows up if it ended from all habitat turning off
    #n.on<-c(n.on, n.on[length(n.on)]) 
    ts<-c(ts, time[length(time)])
  }
  
  #average simulated metaop size over time calculated using trapezoidal method (As used by Austin)
  avg.size=trapz(time,metapop.size)/time[length(time)]
  avg.size
  lead.time<-simu.data$time[-1]
  lag.time<-simu.data$time[-length(simu.data$time)]
  geom.avg.size=exp(sum(log(metapop.size[-length(metapop.size)])*(lead.time-lag.time))/lead.time[length(lead.time)])
  avg.size.reim<-sum(metapop.size[-length(metapop.size)]*(lead.time-lag.time))/lead.time[length(lead.time)]
    ((lead.time-lag.time)/lead.time[length(lead.time)]*(metapop.size[-length(metapop.size)]-metapop.size[-1])/2)+
    sum(metapop.size[-1]+(metapop.size[-length(metapop.size)]-metapop.size[-1])/lead.time[length(lead.time)])
  
  geom.avg.size
  
  #to calculate the avg. number of habitat patches on over time not excluding the end state:
  if (length(configs.data)==n.patches+1){ #if no configuration change happened in the time limit
    #n.on<-sum(configs.data[-1]) #the total number of patches is simply the sum of patches on in that config
    avg.n.on<-configs.data$n.on #so the average n.on is just n.on
  }else{
    n.on<-rowSums(configs.data[,2:(n.patches+1)])
    configs.data$time<-configs.data$V1
    lag.time<-configs.data$time
    lead.time<-c(configs.data$time[-1], time[length(time)])
    #reimann sum formula (rectangular area)
    avg.n.on<-sum(n.on*(lead.time-lag.time))/lead.time[length(lead.time)]
    geomavg.n.on<-exp(sum(log(n.on)*(lead.time-lag.time))/lead.time[length(lead.time)])
  }
  exp.n.on<-unique(plot.data$exp.n.on[(plot.data$r.on)==r.on&(plot.data$r.off)==r.off])
  
  #CALCULATING VARIANCE IN p* versus p
  var(configs.data$config.pstar)
  var(metapop.size)
  
  #figure
  plot(time, metapop.size)
  metapop.data<-data.frame(time, metapop.size)
  habitat.data<-data.frame(lead.time, lag.time,n.on)
  last.row<-c(time[length(time)], time[length(time)], n.on[length(n.on)]) #to include a point for how much habitat was present at the end of the simulation
  habitat.data<-rbind(habitat.data, last.row)
  print(ggplot(habitat.data) + theme_classic()
        #+ geom_rect(aes(xmin = tausmin, xmax = tausmax, 
        #                ymin = 0, ymax = n.on), 
        #            fill = pal2[simcolor.id], alpha=0.5)
        + geom_step(aes(x = lag.time, y = n.on),size=1, color=pal2[simcolor.id], alpha=0.5)
        + geom_hline(yintercept=geomavg.n.on, linetype="solid", size=1, color = pal2[simcolor.id])
        + geom_hline(yintercept=exp.n.on, linetype="solid", size=1, color = "Black")
        #+ geom_step(data = metapop.data, aes(x = time, y = metapop.size),size=1, color=pal[simcolor.id], alpha=0.7) 
        + labs(x = "time", y = "Number of Habitat Patches", 
               title = "Habitat Simulation Dynamics")
        + theme(text = element_text(size=15))
        #+ scale_y_continuous(sec.axis = sec_axis(~.*1, name = "Number of Habitable Patches"))
        #+ geom_hline(yintercept=exp.QED.size, linetype="solid", size=1, color = "Black")
        #+ geom_hline(yintercept=avg.size, linetype="solid", size=1, color = "#B50A2AFF")
        #+ xlim(0,taus[length(taus)])
        #+ xlim(0,10000)
        + ylim(0,n.patches)
  )
  ggsave2(paste0("Habitat Simulation Dynamics geomsample lmQED20", 
                 " ron", r.on, " roff", r.off, " a", alpha,  " d", delta, ".jpeg"), 
          height=5, width=10, units="in", dpi=800)
  if(is.na(geom.avg.size)){geom.avg.size<-0}
  print(ggplot() + theme_classic()
        + geom_step(data = metapop.data, aes(x = time, y = metapop.size),size=1, color=pal[simcolor.id], alpha=0.5) 
        + labs(x = "time", y = "Metapopulation Size", 
               title = "Metapopulation Simulation Dynamics")
        + theme(text = element_text(size=15))
        #+ scale_y_continuous(sec.axis = sec_axis(~.*1, name = "Number of Habitable Patches"))
        + geom_hline(yintercept=avg.size, linetype="solid", size=1, color = pal[simcolor.id])
        + geom_hline(yintercept=exp.QED.size, linetype="solid", size=1, color = "Black")
        #+ xlim(0,taus[length(taus)])
        #+ xlim(0,10000)
        + ylim(0,n.patches)
  )
  ggsave2(paste0("Metapopulation Simulation Dynamics geomsample lmQED20", 
                 " ron", r.on, " roff", r.off, " a", alpha,  " d", delta, ".jpeg"), 
          height=5, width=10, units="in", dpi=800)
  
}

  
  
  
  
  
  
  





# impact of local extinctions versus transients on metapop size plots
########################################################################################################
setwd("C:/Users/Administrator/Desktop/DEC 2022")
source("PStar_func.r")
source("lambda_M_func.r")
setwd("H:/1000 reps")
head(summary.data)

#pick a random iteration number (replicate)
#rx<-sample(summary.data$it.no, 1)
#rx
rx<-1 #cluster that generated the data
pcombos<-c(1:8,10:13,15:28,30:41,43:45) #parameter combos for 9,14,29 had no change in habitat over 10000 years!
#42 seems to also have something funny going on
rx<-1 #cluster that generated the data
#pcombos<-c(9,14,29)
plot.data<-matrix(rep(NA, length(pcombos)*11), nrow=length(pcombos), ncol=11)
z<-1
for (w in pcombos){
  #w=1
  #print(w)
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

  e.rate<-delta
  c.rate<-1
  config.pstar<-rep(NA, nrow(configs.data))
  config.lm<-rep(NA, nrow(configs.data))
  for (i in 1:nrow(configs.data)){
    landscape.config<-landscape[just.configs.data[i,]!=0,]
    config.pstar[i]<-sum(pstar.function(landscape.config, alpha, e.rate, c.rate, self.rec, gamma, epsilon, iterations=1000))
    config.lm[i]<-get.lambda.M(landscape.config, alpha, gamma, epsilon, self.rec, e.rate, c.rate)
    }
  configs.data$config.pstar<-config.pstar
  configs.data$config.lm<-config.lm

#configs.data

#now calculate average change in lm per config change per avg. time spent in a config for each sim
prior.lms<-configs.data$config.lm[-length(configs.data$config.lm)]
next.lms<-configs.data$config.lm[-1]
change.lm<-next.lms-prior.lms
avg.change.lm<-mean(change.lm)
avg.tau<-mean(configs.data$V1)
avg.change.lm.per.tau<-avg.change.lm/avg.tau
#and calculate average decline in Pstar per patch extinction for each sim
prior.pstars<-configs.data$config.pstar[-length(configs.data$config.pstar)]
next.pstars<-configs.data$config.pstar[-1]
avg.change.pstar<-mean(next.pstars-prior.pstars)
local.exts<-which(next.pstars<prior.pstars) #want only changes in p.star corresponding to local extinctions
ext.change.pstars<-next.pstars[local.exts]-prior.pstars[local.exts]
avg.ext.change.pstar<-mean(ext.change.pstars)

#now want the average metapopulation size in each sim
##weight the occupancy values by area and plot metapop size over time
#area.weights<-landscape$areas/sum(landscape$areas)
#metapop.size<-sim.data[,-1]*area.weights *what Austin has, which gives the size of each patch 
#but we don't wanna plot this for all of them, I'm just going to plot the sum
metapop.size<-rowSums(simu.data[,2:(n.patches+1)]) #*area.weights) <-if we wanna add the area weights
time<-simu.data$time
#average simulated metaop size over time calculated using trapezoidal method (As used by Austin)
avg.size=trapz(time,metapop.size)/time[length(time)]
avg.size

plot.data[z,]<-c(Initial.Lm, r.on, r.off, alpha, delta, avg.size, avg.change.lm.per.tau, avg.ext.change.pstar, avg.change.lm, avg.change.pstar, avg.tau)
z<-z+1
}

plot.data<-data.frame(plot.data)
colnames(plot.data)<-c("Initial.Lm", "r.on", "r.off", "alpha", "delta", "avg.size", "avg.change.lm.per.tau", "avg.ext.change.pstar", "avg.change.lm", "avg.change.pstar", "avg.tau")
print(ggplot(plot.data) + theme_classic()
      + geom_point(aes(x = avg.ext.change.pstar, y = avg.size), 
                   color=ghibli_palettes$MononokeMedium[3], alpha=0.5)
      + labs(x = "Avg. Loss of Occupancy Due to Local Habitat Loss", y = "Metapopulation Size", 
             title = "Occupancy Loss Due to Habitat Loss as a Predictor of Metapopulation Size")
      + theme(text = element_text(size=15))
)
ggsave2(paste0("Occupancy Loss Due to Habitat Loss as a Predictor of Metapopulation Size lmQED20", 
               " ron", r.on, " roff", r.off, " a", alpha,  " d", delta, ".jpeg"), 
        height=5, width=10, units="in", dpi=800)

print(ggplot(plot.data) + theme_classic()
      + geom_point(aes(x = avg.change.lm.per.tau, y = avg.size), 
                   color=ghibli_palettes$MononokeMedium[3], alpha=0.5)
      + labs(x = "Avg. Change in Persistence Capacity Relative per Avg. Time in a Configuration", y = "Metapopulation Size", 
             title = "Relative Duration of Transients as a Predictor of Metapopulation Size")
      + theme(text = element_text(size=15))
)
ggsave2(paste0("Relative Duration of Transients as a Predictor of Metapopulation Size lmQED20", 
               " ron", r.on, " roff", r.off, " a", alpha,  " d", delta, ".jpeg"), 
        height=5, width=10, units="in", dpi=800)






#######################################################################################################




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
