
#This function scales the expected persistence capacity at QED by delta 
#(the species specific extinction to colonization rate ratio) and then simulates the
#combined CTMC and SRLM 
#it outputs a graph of the simulated metapopulation size for each timestep and the total amount of 
#habitat after each transition with the value of lmQED, delta, r.on, r.off and alpha in the filename
#it also outputs the simulated data and configuration data in a file named with the same information
#note: the first row of data in the simu data file gives the expected metapop size at QED

sim.data.n.plot<-function(landscape, lmQED.initial, gamma, epsilon, self.rec, alpha, r.on, r.off, it.no){
  it.no=it.no
  landscape=landscape
  e.rate=1#these are rescaled so this is arbitrary
  c.rate=1#these are rescaled so this is arbitrary
  alpha=alpha
  self.rec=self.rec
  gamma=gamma
  epsilon=epsilon
  Initial.Lm=lmQED.initial
  
  #first let's scale our e.rate and c.rate to ensure a high enough persistence capacity the
  #metapopulation has the opportunity to experience some growth
  
  lm.n.pstar<-lm.n.pstar.QED(landscape=landscape, e.rate=e.rate, c.rate=c.rate, gamma=gamma, epsilon=epsilon, self.rec=self.rec, alpha=alpha, r1=r.on, r2=r.off)
  exp.QED.size<-lm.n.pstar[[1]]
  delta<-Initial.Lm/lm.n.pstar[[2]]
  delta
  e.rate<-1/delta
  c.rate<-1
  lm.n.pstar<-lm.n.pstar.QED(landscape=landscape, e.rate=e.rate, c.rate=c.rate, gamma=gamma, epsilon=epsilon, self.rec=self.rec, alpha=alpha, r1=r.on, r2=r.off)
  exp.QED.size<-lm.n.pstar[[1]]
  exp.QED.size
  exp.lm<-lm.n.pstar[[2]]*delta
  exp.lm
  
  output<-CTMC.SRLM(total.trans=2000, landscape=landscape, e.rate<-e.rate, c.rate<-c.rate, alpha<-alpha,
                    gamma<-0, epsilon=n.patches, self.rec<-1, r1<-r.on, r2<-r.off) 
  simu.data<-output[[1]]
  simu.data<-rbind(rep(exp.QED.size, length(simu.data)), simu.data) #the first line of the simu.data file 
  #contains the expected metapop size at QED
  configs.data<-output[[2]]
  write.csv(simu.data, paste0("Simudata Combined CTMC and SRLM simulation sample lmQED100 d", delta, 
                              " ron", r.on, " roff", r.off, " a", alpha,  it.no, ".csv"))
  write.csv(configs.data, paste0("configdata Combined CTMC and SRLM simulation sample lmQED100 d", delta, 
                                " ron", r.on, " roff", r.off, " a", alpha,  it.no, ".csv"))
  
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
               title = "Combined CTMC and SRLM Simulation")
        + theme(text = element_text(size=15))
        + scale_y_continuous(sec.axis = sec_axis(~.*1, name = "Number of Habitable Patches"))
        + geom_hline(yintercept=exp.QED.size, linetype="dashed", color = "red")
        + geom_hline(yintercept=avg.size, linetype="dashed", color = "black"))
  #setwd("C:/Users/abuga/Documents/HVM rscripts/Emma HVM cleaned/Integrated code for clustered landscapes/preliminary plots")
  dev.copy(png, paste0("Combined CTMC and SRLM simulation sample lmQED100 d", delta, 
                       " ron", r.on, " roff", r.off, " a", alpha,  it.no, ".png"))
  dev.off()
}
#check:
#landscape<-create.landscape(n.patches=10, landscape.type="linear", landscape.limit=100, 
#                            patch.distribution="uniform", areas.distribution="uniform", areas.limit=1, 
#                            clustering.iters=0)
#print(ggplot(landscape, aes(x.coord,y.coord)) + theme_classic() + geom_point(aes(size = areas)))
#n.patches<-length(landscape$patch.ID) #for future use

#test<-sim.data.n.plot(landscape=landscape, lmQED.initial=100, gamma=0, 
#                      epsilon=n.patches, self.rec=1, alpha=10, r.on=1/100000, r.off=1/100)
  
