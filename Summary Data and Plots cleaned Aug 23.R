rm(list=ls()) #clear the workspace

#loading useful libraries:
################################################################################
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
library("Bolstad2") #for integration using simpson's rule
library("tidyverse") #for ggplot2 and data manip functions
library("data.table") #gives the efficient rbindlist function
library("pracma") #gives trapz function
library("cowplot") #to get ggsave2 function
setwd("I:/600 reps 2023") #current figure sims norm(x1)<1e-4
#setwd("H:/2023_07_03_reps")  #new sims norm(x1)<1e-10 -no difference in results
##################################################################################

# load in all the summary data files and put them into one data.frame
#takes about 30 min - 1 hr for ~10000 files and only needs to be run once to create summary data file
#######################################################################################################
#file.names <- list.files(pattern="summary data lmQED20")
#sim.summaries <- lapply(file.names, read_csv)
#summary.data<-rbindlist(sim.summaries)
#write.csv(summary.data, "summary_data_600reps.csv")
#######################################################################################################

#read in summary data and check number of rep sims run to acertain when if any parameter combinations
#sims failed 
#####################################################################################################
#setwd("H:/2023_07_03_reps")  #new sims norm(x1)<1e-10 - no diff in results, but can more clearly see 
# that it is only in the r.on=0.01, r.off=0.1, and alpha=10 (1/10th stepping stone scenario E) that 
#sims may fail -therefore best to exclude these from results since data could be biased
#summary.data<-read_csv("summary_data_600reps.csv")
#q<-group_by(summary.data, r.on, r.off, alpha) %>%
#  tally(is.na(avg.size))
#q #NA's indicating a problem with the ODE solver only happen for r.on=0.01, r.off=0.1, and alpha=10
##Exclude this parameter combination from my results
##summary.data<-na.omit(summary.data) # removes these NA's from from the data frame to make 
##plotting/data manip easier, but will need to remember to remove the r.on=0.01, r.off=0.1, alpha=10
##parameter combination from my plots
####################################################################################################

#read in the summary data and check the number of replicates run for each parameter combo
####################################################################################################
summary.data<-read_csv("summary_data_600reps.csv")
q<-group_by(summary.data, r.on, r.off, alpha)
tally(q)
####################################################################################################

#since I stopped the running of my simulations after approximately 200 reps were run (out of the 600 I 
#originally intended to run just because of time constraints), I take only the first 200 sims of each 
#parameter combination run just to ensure a common sample size
####################################################################################################
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
plot.data<-summary.data2 #this is the data from simulations to be used in my figures
###################################################################################################

# since the simulations were run we spotted a few errors in some of the values calculated -since the simulations
# don't depend on these values we were able to just recalculate them appropiately here and update the recorded data
#######################################################################################################
# A) loading the corrected required functions:
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
#setwd("H:/2023_07_03_reps")  # if want new sims norm(x1)<1e-10 - no diff in results
# B) redefining the parameter values used in the simulations for these calculations:
Initial.Lm<-20
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
# C) first recalculating the  equilibrium occupancy expected using metric provided by dewoody et al. 2005
# since we noticed delta wasn't cancelled out properly within it's calculation
plot.data$dewoody.p.star<-rep(NA, nrow(plot.data))
for (i in 1:nrow(plot.data)){
  plot.data$dewoody.p.star[i]<-sum(dewoody.pstar(landscape=landscape, alpha=plot.data$alpha[i], e.rate=plot.data$delta[i], c.rate=1, 
                self.rec=self.rec, gamma=gamma, epsilon=epsilon, iterations=iterations, 
                r.on=plot.data$r.on[i], r.off=plot.data$r.off[i]))
}
# D) calculating some additional values we wanted to check or add to our figures and tables for the paper or use in later calculations
# this takes a little while but only needs to be done once to add this to a new data file
plot.data$exp.n.on<-rep(NA, nrow(plot.data))
plot.data$T.absorp<-rep(NA, nrow(plot.data))
plot.data$rate.comparison<-rep(NA, nrow(plot.data))
r.offs<-c(1/100000,1/10000,1/1000,1/100,1/10)
r.ons<-c(1/100,1/10,1)
alphas<-c(0.01,0.1,10)
for (j in 1:length(r.offs)){
  for (k in 1:length(r.ons)){
    for (l in 1:length(alphas)){
      #extracting delta values from data:
      delta<-plot.data$delta[(plot.data$r.off==r.offs[j]&plot.data$r.on==r.ons[k]&plot.data$alpha==alphas[l])]
      delta<-delta[1]
      #recalculating QED metrics
      x<-QED.metrics(landscape=landscape, e.rate=delta, c.rate=1, gamma=gamma, epsilon=epsilon, self.rec=self.rec, alpha=alphas[l],r.on=r.ons[k], r.off=r.offs[j])
      # obtaining the number of patches expected to be on at QED in sims:
      s<-nrow(plot.data[(plot.data$r.off==r.offs[j]&plot.data$r.on==r.ons[k]&plot.data$alpha==alphas[l]),]) 
      # adding it to our data
      plot.data$exp.QED.size[(plot.data$r.off==r.offs[j]&plot.data$r.on==r.ons[k]&plot.data$alpha==alphas[l])]<-rep(x[[3]],s)
      new.configs.n.pstars<-cbind(x[[7]],x[[8]])
      write.csv(new.configs.n.pstars, paste0("new.configs.n.pstars",Initial.Lm,
                                             " ron",r.ons[k]," roff",r.offs[j]," a",alphas[l],".csv"))
    }
      # obtaining the QED, times to absorption, etc. of the CTMC and adding these to the plot data
      y<-quasi.eq(n.patches=10, r.on=r.ons[k], r.off=r.offs[j])
      r<-nrow(plot.data[(plot.data$r.off==r.offs[j]&plot.data$r.on==r.ons[k]),])
      plot.data$exp.n.on[(plot.data$r.off==r.offs[j]&plot.data$r.on==r.ons[k])]<-rep(x[[5]],r)
      plot.data$T.absorp[(plot.data$r.off==r.offs[j]&plot.data$r.on==r.ons[k])]<-rep(y[[2]],r)
      plot.data$rate.comparison[(plot.data$r.off==r.offs[j]&plot.data$r.on==r.ons[k])]<-rep(y[[3]],r)
      plot.data$med.exp.n.on[(plot.data$r.off==r.offs[j]&plot.data$r.on==r.ons[k])]<-rep(x[[9]],r)
      plot.data$QED.T.absorp[(plot.data$r.off==r.offs[j]&plot.data$r.on==r.ons[k])]<-rep(x[[10]],r)
  }}
#write.csv(plot.data, "plot_data_600reps2.csv")
######################################################################################################

# reading in this data with these new values added
######################################################################################################
plot.data<-read_csv("plot_data_600reps2.csv")
#####################################################################################################

# Generating a table of the the times to absorption
######################################################################################################
p<-plot.data%>%
  select(r.on, r.off, QED.T.absorp)%>%
  group_by(r.on, r.off)%>%
  reframe(unique(QED.T.absorp))
p # you'll notice there are some NA's the table. We know from when we tested the time to absorption 
# calculations for the CTMC that this can happen as values in our matrix become so small they are 
# r rounds them to 0, causing these calculations to fail. We can see this happens when the time to 
# absorption becomes very long >1x10^20, so we can conservatively say the time to absorption exceeds
# 1x10^20 in the table for our paper
write.csv(p, "Time to absorption table.csv")
#######################################################################################################

# converting some of the variables in our data to factors for the sake of plotting
########################################################################################################
plot.data$r.off<-as.factor(log10(plot.data$r.off))
plot.data$r.on<-as.factor(log10(plot.data$r.on))
plot.data$alpha<-factor(plot.data$alpha, 
                        levels=c("0.01", "0.1", "10"),
                        labels=c("Global Dispersal", "Stepping Stone", "1/10th Stepping Stone"))
plot.data$percent.occupied<-plot.data$avg.size/10*100
########################################################################################################

# recalculating the average number of patches on throughout time within sims to ensure this calculation does not include
# time spent in the no habitat state
# calculating the geometric mean number of patches on at QED to compare to the expected geometric mean number of patches
# on at QED for our figures and recalculating p*_QED since we noticed delta wasn't cancelled out properly within it's
# calculation also, then also calculating the variance in p versus p* within simulations
# again this takes a little while but only needs to be done once and saved in a new file for plotting:
########################################################################################################
r.offs<-c(1/100000, 1/10000,1/1000,1/100,1/10)
r.ons<-c(1/100,1/10,1)
alphas<-c(10,1/10,1/100)
Initial.Lm<-20
# new stuff to calculate:
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
z<-1
#reloading the data for individual sims to do these calculations and saving them in a new data frame:
for (j in 1:length(r.offs)){
  for (k in 1:length(r.ons)){
    for(l in 1:length(alphas)){
      #j<-1
      #k<-1
      #l<-1
      r.on<-r.ons[k] 
      r.off<-r.offs[j]
      alpha<-alphas[l]
      config.files <- list.files(pattern=paste0("configdata Combined CTMC and SRLM simulation sample lmQED", Initial.Lm, " ron", 
                                                r.on, " roff", r.off, " a", alpha,"*"))
      simu.files<-list.files(pattern=paste0("Simudata Combined CTMC and SRLM simulation sample lmQED", Initial.Lm,
                                            " ron",r.on," roff",r.off," a",alpha,"*"))
      configs.n.pstars.files<-list.files(pattern=paste0("new.configs.n.pstars",Initial.Lm,
                                                        " ron",r.on," roff",r.off," a",alpha,".csv"))
      #n<-length(config.files) # to do this for all data files
      n<-200 #to use just the first 200 files (the replicates we use in our figures)
      for (i in 1:n) {
        #i<-1
        configs.data<-read.csv(config.files[i])
        # *1 annoyingly in the cases where no change in habitat occured R output the config data as a vertical vector rather than
        #horizontal data frame, so these must be reformatted to match the others where multiple configs occured like so:
        ################################ *1
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
        ################################ *1 end
        simu.data<-read.csv(simu.files[i]) # reading in simulation data
        configs.n.pstars<-read.csv(paste0("new.configs.n.pstars",Initial.Lm,
                                          " ron",r.on," roff",r.off," a",alpha,".csv")) # reading in the new p* in all possible configs data
        # just to get all the possible configs to recalculate p* in 
        for (N in 1:n.patches){ #giving the configs data and pstar data common names to merge by
          names(configs.n.pstars)[names(configs.n.pstars) == paste0("V",N)] <- paste0("n",N)
          names(configs.data)[names(configs.data) == paste0("V",N+1)] <- paste0("n",N)
          names(configs.n.pstars)[names(configs.n.pstars) == paste0("X",N,".1")] <- paste0("pstar",N)
        }
        pstars.data<-merge(configs.data[-1], configs.n.pstars, by=paste0("n",1:n.patches))
        pstar.sizes<-rowSums(pstars.data[,(n.patches+2):(n.patches*2+1)])
        if(length(pstar.sizes)==1){ # if there was only on config and therefore only one value for pstar,
          p.star.var<-0 #there was zero variance in pstar
        }else{
          p.star.var<-var(pstar.sizes) #calculate variance in pstars within simulated configs throughout simulations
        }
        # X) CALCULATING AVERAGES OVER TIME IN SIMS: 
        ############################################# *X
        
        # X.1) AVERAGE OCCUPANCY THROUGH TIME:
        metapop.size<-rowSums(simu.data[,3:n.patches+2]) 
        lead.time<-simu.data$time[-1]
        lag.time<-simu.data$time[-length(simu.data$time)]
        p.var<-var(metapop.size) #calculate variance in p throught simulations
        n.on<-rowSums(configs.data[,2:(n.patches+1)])
        # X.1a) arithmetic avg. size using the Riemann sum formula (rectangular area)
        avg.size[z]<-sum(metapop.size[-length(metapop.size)]*(lead.time-lag.time))/lead.time[length(lead.time)]
        # X.1b) arithmetic mean occupancy using Trapezoidal Method:
        trapz.avg.size[z]<-trapz(simu.data$time, metapop.size)/lead.time[length(lead.time)]
        # X.1c) arithemtic mean occupancy through time using Simpson's Method:
        tiny.vec<-seq(0:0.0000000000000000000000000000000000000000000000001,length.out=length(metapop.size)) #adding the tiniest time difference
        # to every time point to ensure these are all unique (since sometimes this time difference was so small R rounded it to zero); this is 
        # necessary to use the sintegral function
        unique.time<-simu.data$time+tiny.vec
        simp.avg.size[z]<-sintegral(unique.time, metapop.size, n.pts=length(metapop.size)*10)$int/(simu.data$time[length(simu.data$time)])
        # X.1d) geometric mean occupancy through time using the Riemann sum formula
        if(any(metapop.size==0)){geom.avg.size[z]<-0
        }else{geom.avg.size[z]<-exp(sum(log(metapop.size[-length(metapop.size)])*(lead.time-lag.time))/lead.time[length(lead.time)])}
        
        # X.2) AVERAGE AMOUNT OF HABITAT THROUGH TIME::
        ts<-configs.data[,1]
        if (length(ts)==1){ #if no habitat change occured
          #average amount of habitat = the amount of habitat
          avg.n.on<-n.on[1]
          geomavg.n.on[z]<-n.on[1]
        }else{
          lag.ts<-ts
          lead.ts<-c(ts[-1],simu.data$time[length(simu.data$time)])
          # X.2a) arithemtic mean number of patches on using Reimann sum formula (rectangular area)
          avg.n.on<-sum(n.on*(lead.ts-lag.ts))/lead.ts[length(lead.ts)] 
          # X.2b) geometric mean number of patches on using Riemann sum formula (rectangular area)
          geomavg.n.on[z]<-exp(sum(log(n.on)*(lead.ts-lag.ts))/lead.ts[length(lead.ts)])
        }
        
        #recording these for new dataframe:
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
write.csv(new.avgs, "newavgs2.csv")
########################################################################################################

# reading in the dataframe with all these new averages calculated and reformatting variables as factors 
# for plotting
########################################################################################################
new.avgs<-read.csv("newavgs.csv")
new.avgs$r.off<-as.factor(log10(new.avgs$new.roffs))
new.avgs$r.on<-as.factor(log10(new.avgs$new.rons))
new.avgs$alpha<-factor(new.avgs$new.alphas, 
                       levels=c("0.01", "0.1", "10"),
                       labels=c("Global Dispersal", "Stepping Stone", "1/10th Stepping Stone"))
########################################################################################################

# data for table 2 in the paper
########################################################################################################
table2<-new.avgs%>%
  group_by(r.off, r.on, alpha)%>%
  summarise(mean.pstar.vars=mean(p.star.vars), mean.p.vars=mean(p.vars), mean.size=mean(geom.avg.size))
table2
write.csv(table2, "table2.csv")
########################################################################################################

# defining the custom colour pallete for the simulated occupancy figures (figs. 3-5)
########################################################################################################
pal<-c("#EC4E20", "#655A7C","#A4392A", "#AC8B41", "#E58507") #custom pallete "#462521","#D74E09", "#353531"
#"#016FB9"
########################################################################################################

# Plotting Figure 3 panel A):
########################################################################################################
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
ggsave2("Expected versus Average Simulated Metapopulation Sizes Across Disturbance Regimes Global Dispersal only.jpg", 
        height=4.5, width=10, units="in", dpi=800)
########################################################################################################
# Plotting Figure 4 panel A):
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
ggsave2("Expected versus Average Simulated Metapopulation Sizes Across Disturbance Regimes Stepping Stone only.jpg", 
        height=4.5, width=10, units="in", dpi=800)
########################################################################################################
# Plotting Figure 5 panel A): 
########################################################################################################
#*note boxplot for scenario E must be removed since results could be biased for
# since simulations with that parameter combination can fail due to errors in R and could bias results for 
# those parameter combos - I do this for the figure presented in the paper
p.data<-filter(plot.data, alpha=="1/10th Stepping Stone")
p2.data<-new.avgs[new.avgs$new.alphas==10,]
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
ggsave2("Expected versus Average Simulated Metapopulation Sizes Across Disturbance Regimes 1 10th Stepping Stone only.jpg", 
        height=4.5, width=10, units="in", dpi=800)
################################################################################################
# To make these plots without p*_DeWoody
################################################################################################
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
################################################################################################


# defining the custom colour pallete for the simulated number of habitat patches figure (fig. 2)
########################################################################################################
#pal<-c("#436969", "#444F30", "#858C6A", "#4F772D", "#90A955")#
#pal<-c("#4B371C", "#586523", "#A59E2E", "#648E07", "#7C985C")
pal<-c("#436969", "#4B371C", "#444F30", "#A59E2E", "#648E07")
#pal<-c("#70622D", "#2E4B6A","#585123", "#588157", "#A59E2E") #custom pallete "#042A2B","#3B252C", "#49496A", "#576CA8"
########################################################################################################

# Plotting Figure 2 panel A):
########################################################################################################
p2.data<-new.avgs[new.avgs$new.alphas==1/100,] #habitat data plotted from same sims as displayed in Figure 3. for comparison purposes
ggplot() + 
  geom_boxplot(data=p2.data, aes(x=r.on, y=geomavg.n.on, fill=r.off, col=r.off), alpha=0.25) +  
  geom_boxplot(data=p.data, aes(x=r.on, y=exp.n.on, fill=r.off, col="X"), alpha=0.25) +
  geom_line(data=p.data, aes(x=r.on,y=exp.n.on, linetype="Expected Avg. Number of Patches at QED"), alpha=0)+
  theme_classic() + scale_fill_manual(values=pal) + scale_color_manual(breaks=c("-5","-4","-3", "-2","-1","X"), values=c(pal,"black")) +
  labs(y = "Avg. Number of Habitatable Patches", x = "Log Rate of Recovery") + 
  guides(fill=guide_legend("Log Rate of Disturbance", override.aes = list(fill=pal, alpha=0.25, colour=pal)), 
         linetype=guide_legend("Estimate", override.aes = list(linetype=1, alpha=1, colour=c("black"))),
         colour="none")
ggsave2("Expected versus GeomAverage Amount of Habitat Across Disturbance Regimes Global Dispersal Only.jpg", 
        height=4.5, width=10, units="in", dpi=800)
########################################################################################################

# PLOTTING INFORMATION FROM INDIVIDUAL SIMULATIONS TO ADD TO THESE FIGURES 2-5 (panels B-D):
########################################################################################################
# A) redefining the parameters run during simulations: 
time.limit<-10000
r.offs<-c(1/100000, 1/10000,1/1000,1/100,1/10)
r.ons<-c(1/100,1/10,1)
alphas<-c(10,1/10,1/100)
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

# B) loading some necessary functions and the data used for the other figures
setwd("C:/Users/Administrator/Desktop/DEC 2022")
source("PStar_func.r")
source("lambda_M_func.r")
setwd("I:/600 reps 2023")
plot.data<-read_csv("plot_data_600reps2.csv")
head(summary.data)
#rx<-sample(summary.data$it.no, 1)  #if want to pick a random iteration number (replicate)
#pcombos<-1:45 #index of all parameter combinations

# C) defining which simulations to use in example plots
# the only time a rep other than rep 1 was chosen for these was just to choose simulations where enough habitat configuration
# changes occurred as to give a better visual sense of what simulation dynamics look like through time 
rxs<-c(rep(1,6), 14, 2, rep(1, 4)) #replicate numbers for chosen example simulation plots
pcombos<-c(43,44,45,22,23,24,7,8,9,4,5,6) #index of parameter combinations for chosen example simulation plots
examples<-data.frame(rx,pcombos)
# D) providing the appropiate plotting colours to match boxes in panel A) of figures 2-5:
pal<-c("#EC4E20", "#655A7C","#A4392A", "#AC8B41", "#E58507") #OCCUPANCY COLOUR PALLETE
pal2<-c("#4B371C", "#586523", "#A59E2E", "#648E07", "#7C985C") #HABITAT COLOUR PALLETE

for (w in 1:length(pcombos)){ #to create example plots for each parameter combination chosen above
  #w<-43 #to choose just a single set of parameters to create example plots for
  rx<-rxs[w]
  rx
  pcombo<-pcombos[w]
  pcombo
  
  # E) pulling out parameter values from the simulation data:
  rx.summary.data<-summary.data%>%
    filter(it.no==rx)
  r.on<-rx.summary.data$r.on[pcombo]
  r.off<-rx.summary.data$r.off[pcombo]
  alpha<-rx.summary.data$alpha[pcombo]
  delta<-rx.summary.data$delta[pcombo]
  Initial.Lm<-rx.summary.data$Initial.Lm[pcombo]
  exp.QED.size<-unique(plot.data$exp.QED.size[plot.data$r.on==r.on & plot.data$r.off==r.off & plot.data$alpha==alpha])
  
  # F) setting the colour to plot the sim to match the colour of the boxplot it appears in
  simcolor.id=which(r.offs==r.off)
  
  # G) reading in the data for the sim
  configs.data<-read.csv(paste0("configdata Combined CTMC and SRLM simulation sample lmQED",Initial.Lm,
                                " ron",r.on," roff",r.off," a",alpha,"d",delta,"_",rx,".csv"))
  simu.data<-read.csv(paste0("Simudata Combined CTMC and SRLM simulation sample lmQED",Initial.Lm,
                             " ron",r.on," roff",r.off," a",alpha,"d",delta,"_",rx,".csv"))
  
  # H) if the simulation ended early in the last configuration not because it was the absorbing state but because the change in metapop
  #size was below 1x10-4 then record the metapop size at time=10000 as it's last recorded occupancy: 
  if(simu.data$time[nrow(simu.data)]<10000 & sum(configs.data[nrow(configs.data),3:(n.patches+2)])>1){
    simu.data<-rbind(simu.data, simu.data[nrow(simu.data),]) #then let the metapop size at the end of 10000 be the same
    simu.data$time[nrow(simu.data)]<-10000
  }
  ## commented out code here is to be able to plot the cases where no habitat changes occured in 10000 years*
  ## *not necessary for the examples plotted for the paper of course since these were all chosen such that at least a few habitat changes occur
  ## to provide at least some visual representation of the dynamics seen across simulations
  #if (length(t(as.matrix(configs.data==(2*n.patches+2))))){ 
  #  configs.data<-t(as.matrix(configs.data))
  #  configs.data<-configs.data[-1,]
  #  configs.data<-data.frame(rbind(configs.data, configs.data))
  #  configs.data[2,1]<-10000
  #  just.configs.data<-configs.data[,-1]
  #  simu.data<-simu.data[,-1]
  #}else{
  just.configs.data<-configs.data[,-(1:2)] #extracting columns with just the config (on/off patch) data to use in calculations
  configs.data<-configs.data[,-1] #extracting columns with just the configs and their occurence times to use
  simu.data<-simu.data[,-1] #extracting just the simulated occupancy and time data to use
  #}
  #configs.data<-read.csv("configdata Combined CTMC and SRLM simulation sample lmQED20 ron1 roff0.001 a0.01d0.370311750341306_665.csv")
  #simu.data<-read.csv("Simudata Combined CTMC and SRLM simulation sample lmQED20 ron1 roff0.001 a0.01d0.370311750341306_665.csv")
  
  # I) recalculating p*'s within each simulated config, to correct for the minor p* calculation error mentioned above
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
  
  # J) PLOTTING THE HISTOGRAMS: bin widths have to be manually adjusted to be most appropriate for the data
  ############################################################################################
  # P* HISTOGRAM:
  ############################################
  #print(ggplot(configs.data) + theme_classic()
  #      + geom_histogram(aes(config.pstar), binwidth=1, fill=pal[simcolor.id], alpha=0.25)
  #      + geom_vline(aes(xintercept=mean(config.pstar)), color=pal[simcolor.id])#color="#B50A2AFF")
  #      + labs(x = "P*", y = "Frequency")
  #      + theme(text = element_text(size=15))
  #)
  #ggsave2(paste0("Pstar distribution lmQED20", 
  #               " ron", r.on, " roff", r.off, " a", alpha,  " d", delta, ".jpeg"), 
  #        height=5, width=5, units="in", dpi=800)
  ############################################
  # Lambda_M HISTOGRAM:
  ############################################
  #print(ggplot(configs.data) + theme_classic()
  #      + geom_histogram(aes(config.lm), binwidth=10, fill=pal[simcolor.id], alpha=0.25)
  #      + geom_vline(aes(xintercept=mean(config.lm)), color=pal[simcolor.id])#color="#B50A2AFF")
  #      + labs(x = "Persistence Capacity", y = "Frequency")
  #      + theme(text = element_text(size=15))
  #)
  #ggsave2(paste0("Lm distribution lmQED20", 
  #               " ron", r.on, " roff", r.off, " a", alpha,  " d", delta, ".jpeg"), 
  #        height=5, width=5, units="in", dpi=800)
  ############################################
  # TAUS HISTOGRAM:
  ############################################
  #remove the initial condition as a time at which a habitat change occured since this is just the starting condition not a change
  #configs.data2<-configs.data[-1,]
  #print(ggplot(configs.data2) + theme_classic()
  #      + geom_histogram(aes(taus), binwidth=1, fill=pal2[simcolor.id], alpha=0.25)
  #      + geom_vline(aes(xintercept=mean(taus)), color=pal2[simcolor.id])#, color="#B50A2AFF")
  #      + labs(x = "Time Between Habitat Changes", y = "Frequency")
  #      + theme(text = element_text(size=15))
  #)
  #ggsave2(paste0("Tau distribution lmQED20", 
  #               " ron", r.on, " roff", r.off, " a", alpha,  " d", delta, ".jpeg"), 
  #        height=5, width=5, units="in", dpi=800)
  #############################################
  
  # K) Simulation Figures:
  ########################################################################################################
  # K.1) Calculating average occupancy within the simulation to add to the figure:
  ##commented out code is to weight occupancies by patch areas/carrying capacities but isn't necessary for our simulations 
  ##since we keep these constant at 1 for all patches
  ##weight the occupancy values by area and plot metapop size over time
  #area.weights<-landscape$areas/sum(landscape$areas)
  #metapop.size<-sim.data[,-1]*area.weights
  metapop.size<-rowSums(simu.data[,2:(n.patches+1)]) #area.weights) <-if we wanna add the area weights
  time<-simu.data$time
  lead.time<-simu.data$time[-1]
  lag.time<-simu.data$time[-length(simu.data$time)]
  # K.1 a) arithmentic mean simulated metaop size over time calculated using trapezoidal method (as used by Austin)
  avg.size<-trapz(time,metapop.size)/time[length(time)]
  # K.1 b) geometric avg occupancy through time using Rieman sum method:
  geom.avg.size=exp(sum(log(metapop.size[-length(metapop.size)])*(lead.time-lag.time))/lead.time[length(lead.time)])
  # K.1 c) arithmentic mean simulated metaop size over time calculated using Riemann sum method
  avg.size.reim<-sum(metapop.size[-length(metapop.size)]*(lead.time-lag.time))/lead.time[length(lead.time)]
  # K.2) Calculating the avg. number of habitat patches on over time with the end state:
  n.on<-rowSums(configs.data[,2:(n.patches+1)])
  ts<-configs.data[,1]
  if (ts[length(ts)]<time[length(time)]){
    #so the time spent in the last habitat configuration actually shows up if it ended from all habitat turning off
    #n.on<-c(n.on, n.on[length(n.on)]) 
    ts<-c(ts, time[length(time)])
  }
  if (length(configs.data)==n.patches+1){ #if no configuration change happened in the time limit
    #n.on<-sum(configs.data[-1]) #the total number of patches is simply the sum of patches on in that config
    avg.n.on<-configs.data$n.on #so the average n.on is just n.on
  }else{
    n.on<-rowSums(configs.data[,2:(n.patches+1)])
    configs.data$time<-configs.data$V1
    lag.time<-configs.data$time
    lead.time<-c(configs.data$time[-1], time[length(time)])
    # K.2 a) arithemtic mean number of patches using the Riemann sum formula 
    avg.n.on<-sum(n.on*(lead.time-lag.time))/lead.time[length(lead.time)]
    # K.2 b) geometric mean number of patches using the Riemann sum formula 
    geomavg.n.on<-exp(sum(log(n.on)*(lead.time-lag.time))/lead.time[length(lead.time)])
  }
  # getting the expected amount of habitat at QED from the plot data used earlier
  exp.n.on<-unique(plot.data$exp.n.on[(plot.data$r.on)==r.on&(plot.data$r.off)==r.off])
  
  #CALCULATING VARIANCE IN p* versus p within the example simulations
  var(configs.data$config.pstar)
  var(metapop.size)
  
  #finally actually plotting the figure
  ############################################
  plot(time, metapop.size)
  metapop.data<-data.frame(time, metapop.size)
  habitat.data<-data.frame(lead.time, lag.time,n.on)
  last.row<-c(time[length(time)], time[length(time)], n.on[length(n.on)])
  #to include a point for how much habitat was present at the end of the simulation
  habitat.data<-rbind(habitat.data, last.row)
  ############################################
  # TO PLOT HABITAT DYNAMICS  for example simulation: (commented out lines allow the metapopulation dynamics to be added to this)
  ############################################
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
  ############################################
  # TO PLOT METAPOPULATION DYNAMICS for example simulation:
  ############################################
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
  ggsave2(paste0("Metapopulation Simulation Dynamics sample lmQED20", 
                 " ron", r.on, " roff", r.off, " a", alpha,  " d", delta, ".jpeg"), 
          height=5, width=10, units="in", dpi=800)
  ############################################
}

# even though it looks like the metapop size drops to 0 and comes back up sometimes, 
# it's just extremely small, never quite extinct in these instances, e.g. using the data of the last example sim
simu.data[rowSums(simu.data)==0,]