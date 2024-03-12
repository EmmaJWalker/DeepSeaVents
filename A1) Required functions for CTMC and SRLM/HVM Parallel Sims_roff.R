rm(list=ls()) #clear the workspace

#require(devtools)
#install_version("doParallel", version = "1.0.12", repos = "http://cran.us.r-project.org")
#packageVersion('doParallel')

# process in parallel
library("doParallel") 
cl <- makeCluster(detectCores(), type='PSOCK')
registerDoParallel(cl)

####################################################
single.batch<-function(i){
  
  reps.of.45<-22
  for (r in 0:(reps.of.45-1)){
    #i=1
    #r=1
  it.no=i+45*r
  print(it.no)
  
  #LOADING CODE SPECIFIC PACKAGES:
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
  #LOADING REQUIRED FUNCTIONS
  setwd("C:/Users/Administrator/Desktop/DEC 2022")
  source("create_landscape_func.r") #create.landscape()
  source("CTMC_and_SRLM_func_from_QED2_roff.r") #CTMC.SRLM()
  source("dispersal_kernel_func.r") #disp.kernel()
  source("lambda_M_func_roff.r") #get.lambda.M()
  source("QED_metrics_func_roff.r") #QED.metrics()
  source("dewoody_func.r") #dewoody()
  source("Pstar_func_roff.r") #pstar.function()
  source("quasi_eq_sparse_func.r") #quasi.eq()
  source("SRLMODE_func.r") #SRLM.ODE()
  source("at_equilibrium_rootfunc.r") #at.equilibrium()
  #source("single_parallel_batch_func.r") #single.batch()
  getwd()
  setwd("I:/45 reps 2023 roff")  
  

  #SETTING EXPERIMENT PARAMETERS:
  r.offs<-rev(c(1/100000, 1/10000,1/1000,1/100,1/10))
  r.ons<-rev(c(1/100,1/10,1))
  alphas<-c(10,1/10,1/100) #avg. disp = 1/10 interpatch distance, the interpatch distance and global disp
  n.patches<-10 #20
  landscape.type<-"linear"
  landscape.limit<-100
  patch.distribution<-"uniform"
  areas.distribution<-"uniform"
  areas.limit<-1
  clustering.iters<-0
  
  #CREATING EXPERIMENTAL LANDSCAPES
  landscape<-create.landscape(n.patches=n.patches, landscape.type=landscape.type, landscape.limit=landscape.limit, 
                              patch.distribution=patch.distribution, areas.distribution=areas.distribution, areas.limit=areas.limit, 
                              clustering.iters=clustering.iters)
  #print(ggplot(landscape, aes(x.coord,y.coord)) + theme_classic() + geom_point(aes(size = areas))) #plotting
  n.patches<-length(landscape$patch.ID) #for future use
  
  #RUNNING A BATCH OF SIMULATIONS OVER ALL THE PARAMETER COMBINATIONS WE PROVIDED
  for (j in 1:length(r.offs)){
    for (k in 1:length(r.ons)){
      for (l in 1:length(alphas)){
        
        #j=2
        #k=3
        #l=1
        
        run.check<-data.frame(j,k,l,it.no)
        write.csv(run.check, paste0("runcheck ", j, " ", k, " ", l, " ", it.no, ".csv"))
        r.off<-r.offs[j]
        r.on<-r.ons[k]
        
        e.rate<-1#0.1 #can pick just any number below 1 to start since we will just rescale this
        c.rate<-1  #can pick just any number below 1 to start since we will just rescale this
        alpha<-alphas[l] #global dispersal 1/landscape.limit
        self.rec<-1
        gamma<-0
        epsilon<-n.patches
        time.limit<-10000 #DOESN'T MAKE ANYTHING MUCH SLOWER OR FASTER
        
        #first let's scale our e.rate and c.rate to ensure a high enough persistence capacity the
        #metapopulation has the opportunity to experience some growth
        Initial.Lm<-20
        
        
        metrics<-QED.metrics(landscape=landscape, e.rate=e.rate, c.rate=c.rate, gamma=0, epsilon=n.patches, self.rec=1, alpha=alpha, r.on=r.on, r.off=r.off)
        delta<-metrics[[2]]/Initial.Lm
        e.rate<-delta-r.off
        c.rate<-1
        metrics<-QED.metrics(landscape=landscape, e.rate=e.rate, c.rate=c.rate, gamma=gamma, epsilon=epsilon, self.rec=self.rec, alpha=alpha, r.on=r.on, r.off=r.off)
        exp.lm<-metrics[[1]]
        exp.QED.size<-metrics[[3]]
        exp.n.on<-metrics[[5]]
        QED<-metrics[[6]]
        all.configs<-as.matrix(as(metrics[[7]], "sparseMatrix")) #IT'S FORMAT IS MESSED UP OTHERWISE
        all.configs<-all.configs[-1,] #all configs excluding the absorbing state (which we won't start in)
        pstar.configs<-metrics[[8]]
        pstar.configs<-pstar.configs[-1,] #all configs excluding the absorbing state (which we won't start in)
        #initial.config<-rbinom(n.patches, 1, 0.5) #to start in a random configuration *OR*
        #initial.config<-rep(1, n.patches) #to start with all patches on *OR*
        x<-sample(x=(1:nrow(all.configs)), size=1, prob=QED) #to choose a config to start with the probability of being in it at QED
        initial.config<-all.configs[x,]
        initial.p<-pstar.configs[x,]
        dewoody<-dewoody(e.rate=e.rate, c.rate=c.rate, r.on=r.on, r.off=r.off)
        dewoody.delta<-dewoody[[1]]
        dewoody.s<-dewoody[[2]]
        #exp.dewoody.size<-rowSums(pstar.configs/dewoody.s) #double check because they pull area out of the calculation of P* in each config and then later reweight by it... i think this should cancel out like so
        
        output<-CTMC.SRLM(total.t=time.limit, landscape=landscape, e.rate=e.rate, c.rate=c.rate, alpha=alpha,
                          gamma=gamma, epsilon=epsilon, self.rec=self.rec, r.on=r.on, r.off=r.off, initial.p=initial.p,
                          initial.config=initial.config) 
        
        simu.data<-output[[1]]
        configs.data<-output[[2]]
        #hist(configs.data[-1,1]-configs.data[-nrow(configs.data),1])
        
        if(any(is.na(simu.data))){ #if the ODE solver gave an NA
          summary.data<-data.frame(it.no, Initial.Lm, exp.lm, dewoody.delta, r.on, r.off, alpha, delta, 
                                   exp.n.on, avg.n.on=NA, exp.QED.size, avg.size=NA, dewoody.s) #output the parameter values and NA's for the calculated metrics
          write.csv(summary.data, paste0("summary data lmQED", Initial.Lm, " ron",                                
                                         r.on, " roff", r.off, " a", alpha,"d", delta, "_", it.no, ".csv"))
        }else{
        #getting the avg. metapop size over time.
        #simulated metapop size in each config:
        metapop.size<-rowSums(simu.data[,2:(n.patches+1)]) #*area.weights) <-if we wanna add area weights
        #expected metapop size in each config:
        
        #but gotta get that data from the landscape
        time<-simu.data$time
        #average simulated metapop size over time calculated using trapezoidal method (As used by Austin)
        avg.size=trapz(time,metapop.size)/time[length(time)]
        
        #getting the avg. number of habitat patches on over time
        if (length(configs.data)==n.patches+1){ #if no configuration change happened in the time limit
          n.on<-sum(configs.data[-1]) #the total number of patches is simply the sum of patches on in that config
          avg.n.on<-n.on #so the average n.on is just n.on
        }else{
          n.on<-rowSums(configs.data[,2:(n.patches+1)])
          taus<-configs.data[,1]
          if (taus[length(taus)] < time.limit){
            lead.taus<-c(taus[-1], time.limit)
          } else {
            taus[length(taus)]<-time.limit
            lead.taus<-c(taus[-1], time.limit)
          }
          #reimann sum formula (rectangular area)
          avg.n.on<-sum(n.on*(lead.taus-taus))/lead.taus[length(lead.taus)]
        }
          
        
        #periodogram
        #p<-periodogram(metapop.size/metapop.pstar) #need to get metapop.pstar by getting and 
        #using the indices configs sampled throughout the simulation and grabbing the P*'s in them
        #might be just as easy to just out put this necessary data and do this after running the sims
        
        summary.data<-data.frame(it.no, Initial.Lm, exp.lm, dewoody.delta, r.on, r.off, alpha, delta, exp.n.on, avg.n.on, exp.QED.size, avg.size, dewoody.s)
        write.csv(summary.data, paste0("summary data lmQED", Initial.Lm, " ron",                                
                                       r.on, " roff", r.off, " a", alpha,"d", delta, "_", it.no, ".csv"))
      }
        write.csv(simu.data, paste0("Simudata Combined CTMC and SRLM simulation sample lmQED", Initial.Lm, " ron", 
                                    r.on, " roff", r.off, " a", alpha,"d", delta, "_", it.no, ".csv"))
        write.csv(configs.data, paste0("configdata Combined CTMC and SRLM simulation sample lmQED", Initial.Lm, " ron", 
                                       r.on, " roff", r.off, " a", alpha,"d", delta, "_", it.no, ".csv"))
        
        configs.n.pstars<-data.frame(as.matrix(all.configs), pstar.configs)
        write.csv(configs.n.pstars, paste0("configs.n.pstars", Initial.Lm, " ron",                                
                                           r.on, " roff", r.off, " a", alpha,"d", delta, "_", it.no, ".csv"))
        }
      }
  }
  }
}
###############


#for 45 sims spread over 48 cores can have each core doing a rep of each parameter 
#combo and get 45 reps of each combo at once, should take about 2ish hours for each it.of.45
#can easily do 360 in a day, ~1000 in a weekend
#its.of.45=22 #22 = 990 reps total in ~44 Hrs ~2 days
#for (s in 0:its.of.45){
  foreach(i=1:45) %dopar% single.batch(i=(i))
#}

 

# turn parallel processing off and run sequentially again:
registerDoSEQ()
#stopCluster()
#stopImplicitCluster()
