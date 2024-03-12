rm(list=ls()) #clear the workspace

library("pracma") #to get trapezoidal function
library("stringr") #to split strings (filenames)
library("TSA")
library("locfit")
#LOADING CODE SPECIFIC PACKAGES:
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
library("Matrix")
library("reticulate")
#LOADING REQUIRED FUNCTIONS
setwd("C:/Users/Administrator/Desktop/JUNE 2022")
source("quasi_eq_sparse_func.r")
source("CTMC_func.r")
source("SRLMODE_func.r")
source("at_equilibrium_rootfunc.r")
source("dispersal_kernel_func.r")
source("CTMC_and_SRLM_func.r")
source("lambda_M_func.r")
source("QED_metrics_func.r")
source("create_landscape_func.r")
source("Pstar_func.r")

#library("SparkR") #for lead and lag functions <- not avaialble in this version of r
setwd("C:/Users/Administrator/Desktop/HVM Parallel Data Batch 1") 

time.limit<-10000
#SETTING EXPERIMENT PARAMETERS:
r.offs<-c(1/100000, 1/10000,1/1000,1/100,1/10)
r.ons<-c(1/100,1/10,1)
alphas<-c(10,1/10,1/100) #avg. disp = 1/10 interpatch distance, the interpatch distance and global disp
n.patches<-10
landscape.type<-"linear"
landscape.limit<-100
patch.distribution<-"uniform"
areas.distribution<-"uniform"
areas.limit<-1
clustering.iters<-0
landscape<-create.landscape(n.patches=n.patches, landscape.type=landscape.type, landscape.limit=landscape.limit, 
                            patch.distribution=patch.distribution, areas.distribution=areas.distribution, areas.limit=areas.limit, 
                            clustering.iters=clustering.iters)
#1. create a dataframe with all possible configurations
#ordered according to the sequence they occur as rows in our transition matrix
#NOTE: updated to store these as sparse matrices to save memory :)
n<-n.patches
all.configs <- Matrix(0, nrow=2^n, ncol=n, sparse=TRUE) #sparse matrix version
#all.configs<-matrix(rep(NA, (2^n)*n), nrow=2^n, ncol=n) #no sparse matrix version
for (i in 1:n){
  #if (i==1){BB<-c(0,1) #no sparse vector version
  if (i==1){BB<-as(c(0,1), "sparseVector") #sparse vector version
  }else{BB<-c(BB[1:(2^(i-2))],BB,BB[((2^(i-2))+1):(2^(i-1))])}
  all.configs[,(n-i+1)]<-t(rep(BB,2^(n-i)))}
n.on.in.each.config<-rowSums(all.configs)

for (j in 1:length(r.offs)){
  for (k in 1:length(r.ons)){
    for (l in 1:length(alphas)){
      
      #j=1 #for testing
      #k=1 #for testing
      #l=1 #for testing
      #print(data.frame(j,k,l))
      
      r.on<-r.ons[k]
      r.off<-r.offs[j]
      alpha<-alphas[l]
      
      rep.sims <- list.files(pattern=paste0("Simudata Combined CTMC and SRLM simulation sample lmQED20 ron", 
                                            r.on, " roff", r.off, " a", alpha,"d*"))
      rep.configs <- list.files(pattern=paste0("configdata Combined CTMC and SRLM simulation sample lmQED20ron", 
                                               r.on, " roff", r.off, " a", alpha,"d*"))
      n.reps<-length(rep.sims)
      print(n.reps)
      
      #extract delta values from the file names
      x<-str_split(rep.sims, pattern=paste0("Simudata Combined CTMC and SRLM simulation sample lmQED20 ron", 
                                            r.on, " roff", r.off, " a", alpha,"d"))
      loc<-str_locate(rep.sims, pattern=paste0("Simudata Combined CTMC and SRLM simulation sample lmQED20 ron", 
                                               r.on, " roff", r.off, " a", alpha,"d"))
      start.d<-loc[,2]
      loc <- str_locate(rep.sims, pattern=paste0("_"))
      end.d<-loc[,1]
      delta<-substring(rep.sims,start.d+1,end.d-1)
      delta<-as.numeric(delta[1])
      
      #calculating the amount of habitat expected to be occupied at QED
      e.rate<-1/delta
      c.rate<-1
      self.rec<-1
      gamma<-0
      epsilon<-n.patches
      metrics<-QED.metrics(landscape=landscape, e.rate=e.rate, c.rate=c.rate, 
                                 gamma=gamma, epsilon=epsilon, self.rec=self.rec, 
                                 alpha=alpha, r.on=r.on, r.off=r.off)
      exp.QED.size<-metrics[[1]]
      exp.lm<-metrics[[2]]*delta
      
      #calculating expected amount of habitat at QED
      QED.output<-quasi.eq(n.patches, r.on, r.off)
      QED<-QED.output[[1]]#QED distribution
      #geometric mean habitat expected across configurations and time spent in them:
      exp.n.on<-exp(dot(QED,log(n.on.in.each.config[2:2^n.patches])))
      
      #Reading in the data and calculating summary statistics to plot:
      avg.sizes<-rep(NA, n.reps)
      avg.n.ons<-rep(NA, n.reps)
      for (i in 1:n.reps) {
        
        #i=1 #for testing
        #print(i)
        
        ##use for whole data frame:
        #data<-rbind(read.csv(rep.sims[i]))
        
        simu.data<-read.csv(rep.sims[i])
        configs.data<-read.csv(rep.configs[i])
        
        #getting the avg. metapop size over time.
        metapop.size<-rowSums(simu.data[,3:(n.patches+2)]) #*area.weights) <-if we wanna add area weights
        #but gotta get that data from the landscape
        time<-simu.data$time
        #average simulated metapop size over time calculated using trapezoidal method (As used by Austin)
        avg.size=trapz(time,metapop.size)/time[length(time)]
        avg.sizes[i]<-avg.size
        print(avg.sizes)
        
        #getting the avg. number of habitat patches on over time
        if(ncol(configs.data)==2){
          n.on<-10
          taus<-c(0,0)
        } else {
          n.on<-rowSums(configs.data[,3:(n.patches+2)])
          taus<-configs.data[,2]
        }
        ## if want average number of patches on over the 100000 years including once in the absorbing state
        #if (taus[length(taus)] < time.limit){
        #  lead.taus<-c(taus[-1], time.limit)
        #} else {
        #  taus[length(taus)]<-time.limit
        #  lead.taus<-c(taus[-1], time.limit)
        #}
        #to get average number of patches before hitting absorbing state
        lag.taus<-taus[-length(taus)]
        lead.taus<-taus[-1]
        #reimann sum formula (rectangular area)
        avg.n.on<-sum(n.on*(lead.taus-lag.taus))/lead.taus[length(lead.taus)]
        avg.n.ons[i]<-avg.n.on
        print(avg.n.ons)
      }
      
      
      #periodogram
      #p<-periodogram(simu.data$X1)
      
      
      for (i in 1:n.reps) {
      data<-data.frame(avg.sizes[i], avg.n.ons[i], r.on, r.off, alpha, delta, exp.QED.size, exp.lm, exp.n.on, i)
      if (j==1 & k==1 & l==1 & i==1){plot.data<-data}
      else {plot.data<-rbind(plot.data, data)}
      }
      
      ##use for whole data frame:
      #r.on
      #data<-cbind(data, r.on, r.off, alpha)
      #if (j==1 & k==1 & l==1){p.data<-data}
      #else{p.data<-rbind(p.data, data)}
      
    }}}
write.csv(plot.data, "plot_data_test1.csv")
