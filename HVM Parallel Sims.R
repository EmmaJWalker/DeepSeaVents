rm(list=ls()) #clear the workspace

# process in parallel
library("doParallel") 
cl <- makeCluster(detectCores(), type='PSOCK')
registerDoParallel(cl)

#for 45 sims spread over 48 cores can have each core doing a rep of each parameter 
#combo and get 45 reps of each combo at once, should take about 2ish hours for each it.of.45
#can easily do 360 in a day, ~1000 in a weekend
its.of.45=1
  foreach(i=1:4) %dopar% 
    
    library("ggplot2")
  library("pracma")
  library("deSolve")
  library("RColorBrewer")
  library("colorspace")
  library("wesanderson") #SUPER NICE COLOUR SCHEMES!
  library("ggthemes")
  library("R.utils") #to get intToBin function
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
  source("sim_data_n_plot_function.r")
  getwd()
  setwd("C:/Users/abuga/OneDrive/Desktop/DeepSeaVents-master/Parallel test data") 
  getwd() #make sure output will go in the correct directory
  
  r.offs<-c(1/100000,1/10000,1/1000,1/100,1/10)
  r.ons<-c(1/100,1/10,1)
  alphas<-c(10,1/10,1/100) #avg. disp = 1/10 interpatch distance, the interpatch distance and global disp
  
  landscape<-create.landscape(n.patches=10, landscape.type="linear", landscape.limit=100, 
                              patch.distribution="uniform", areas.distribution="uniform", areas.limit=1, 
                              clustering.iters=0)
  print(ggplot(landscape, aes(x.coord,y.coord)) + theme_classic() + geom_point(aes(size = areas)))
  n.patches<-length(landscape$patch.ID) #for future use

for (j in 1:length(r.offs)){
  for (k in 1:length(r.ons)){
    for (l in 1:length(alphas)){
      
      sim.data.n.plot(landscape=landscape, lmQED.initial=20, gamma=0, 
                      epsilon=n.patches, self.rec=1, alpha=alphas[l], 
                      r.on=r.ons[k], r.off=r.offs[j], it.no=i)
      
    }
  }
}

# turn parallel processing off and run sequentially again:
registerDoSEQ()

stopCluster(cl)

stopImplicitCluster(cl)