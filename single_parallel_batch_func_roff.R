
single.batch<-function(i){
  
  it.no<-i
  
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
        
        run.check<-data.frame(j,k,l,it.no)
        write.csv(run.check, paste0("runcheck ", j, " ", k, " ", l, " ", it.no, ".csv"))
        r.off<-r.offs[j]
        r.on<-r.ons[k]
        
        e.rate<-1#0.1 #can pick just any number to start since we will just rescale this
        c.rate<-1#4  #can pick just any number to start since we will just rescale this
        alpha<-alphas[l] #global dispersal 1/landscape.limit
        self.rec<-1
        gamma<-0
        epsilon<-n.patches
        
        #first let's scale our e.rate and c.rate to ensure a high enough persistence capacity the
        #metapopulation has the opportunity to experience some growth
        Initial.Lm<-20
        
        
        metrics<-QED.metrics(landscape=landscape, e.rate=e.rate, c.rate=c.rate, gamma=0, epsilon=n.patches, self.rec=1, alpha=alpha, r.on=r.on, r.off=r.off)
        exp.QED.size<-metrics[[2]]
        delta<-metrics[[1]]/Initial.Lm
              e.rate<-delta-r.off
              c.rate<-1
              metrics<-QED.metrics(landscape=landscape, e.rate=e.rate, c.rate=c.rate, gamma=gamma, epsilon=epsilon, self.rec=self.rec, alpha=alpha, r.on=r.on, r.off=r.off)
              exp.lm<-metrics[[1]]
              QED<-metrics[[5]]
              all.configs<-as(metrics[[6]], "sparseMatrix") #it won't stay as a sparse matrix when I pull it out otherwise
              all.configs<-all.configs[-1,] #all configs excluding the absorbing state (which we won't start in)
              pstar.configs<-metrics[[7]]
              pstar.configs<-pstar.configs[-1,] #all configs excluding the absorbing state (which we won't start in)
              #initial.config<-rbinom(n.patches, 1, 0.5) #to start in a random configuration *OR*
              #initial.config<-rep(1, n.patches) #to start with all patches on *OR*
              x<-sample(x=(1:nrow(all.configs)), size=1, prob=QED) #to choose a config to start with the probability of being in it at QED
              initial.config<-all.configs[x,]
              initial.p<-pstar.configs[x,]

              
              output<-CTMC.SRLM(total.t=10000, landscape=landscape, e.rate=e.rate, c.rate=c.rate, alpha=alpha,
                                gamma=0, epsilon=epsilon, self.rec=self.rec, r.on=r.on, r.off=r.off, initial.p=initial.p,
                                initial.config=initial.config) 
  
              simu.data<-output[[1]]
              configs.data<-output[[2]]
              write.csv(simu.data, paste0("Simudata Combined CTMC and SRLM simulation sample lmQED", Initial.Lm, " ron", 
                                          r.on, " roff", r.off, " a", alpha,"d", delta, "_", it.no, ".csv"))
              write.csv(configs.data, paste0("configdata Combined CTMC and SRLM simulation sample lmQED", Initial.Lm, "ron", 
                                             r.on, " roff", r.off, " a", alpha,"d", delta, "_", it.no, ".csv"))
      }
    }
  }
}



