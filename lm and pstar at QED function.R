#Finds the equilibrium metapopulation size for each configuration
###############################################################################
lm.n.pstar.QED<-function(landscape, e.rate, c.rate, gamma, epsilon, self.rec, alpha, r.on, r.off){
  
  #calculate the extinction rate for each patch
  #extinction.rates<-e.rate/areas #if want this to be patch dependent
  #calculate delta
  delta<-e.rate/c.rate
  
  lambda.M<-rep(0,2^n.patches)
  metapop.size<-rep(0,2^n.patches)
  
  #equilibrium metapopulation size: for each possible configuration
  #####################################################################################################
  for (i in 0:((2^n.patches)-1)){
    bin.val<-intToBin(i) #obtain the binary representation of i 
    bin.val<-strsplit(bin.val,"")
    config<-rep(0,n.patches)
    for (k in 1:length(bin.val[[1]])){ #for each binary digit in i
      config[length(bin.val[[1]])-k+1]<-as.numeric(as.character(bin.val[[1]][k]))}

    if (i==0){
      #when i=0, it's the config where everything is turned off 
      #so pulling out the matrix elements doesn't work
      lambda.M[i+1]<-0
      metapop.size[i+1]<-0
    } else { #if i is not zero
      ####################################################
      landscape.config<-landscape[config!=0,]
      lambda.M[i+1]<-get.lambda.M(landscape.config, alpha, gamma, epsilon, self.rec, e.rate, c.rate) 
      
      #calculate the persistence capacity within that configuration
      if(lambda.M[i+1]<delta){metapop.size[i+1]<-0 #if the persistence capacity is
      #below the extinction threshold, then the metapopultion is extinct, set size = 0
      } else {
        metapop.size[i+1]<-sum(pstar.function(landscape.config, alpha, delta, iterations=1000))
        #pstar gives expected occupancy of each patch at equilibrium, the total occupancy should be the sum
      }
    }
  }
  #get the time-averaged metapopulation sizes
  QED.output<-quasi.eq(n.patches, r.on, r.off)
  QED<-QED.output[[1]]#QED distribution
  T.absorp<-QED.output[[2]]#Time to absorption
  a.weighted.metapop.size<-dot(QED,metapop.size[2:2^n.patches]) #this is the expected metapop size in a given configuration x the porportion of time spent in it at QED and summed (therefore the arithmentic mean)
  #should actually be the geometric mean
  weighted.metapop.size<-exp(dot(QED,log(metapop.size[2:2^n.patches])))
  weighted.lm<-exp(dot(QED,log(lambda.M[2:2^n.patches])))
  return(list(weighted.metapop.size, weighted.lm, a.weighted.metapop.size))
}
#Check:
#landscape<-create.landscape(n.patches=10, landscape.type="linear", landscape.limit=100, 
#                            patch.distribution="uniform", areas.distribution="uniform", areas.limit=1, 
#                            clustering.iters=0)
#print(ggplot(landscape, aes(x.coord,y.coord)) + theme_classic() + geom_point(aes(size = areas)))
#n.patches<-length(landscape$patch.ID) #for future use
#
#test<-lm.n.pstar.QED(landscape=landscape, e.rate=0.1, c.rate=0.2, gamma=0, epsilon=n.patches,
#                     self.rec=1, alpha=10, r.on=1/100000, r.off=1/100)
#dewoody<-function(e.rate, c.rate, r.on, r.off){
#  dw.lm<-(e.rate+r.off)/(c.rate*(r.on/(r.on+r.off))) #Pr(patch can colonize another?)
#  s<-r.on/(r.on+r.off) #pr(patch is on)
#  r.o<-c.rate*dw.lm*s/(e.rate+r.off)
#  return(list(dw.lm, r.o))
#}
#test2<-dewoody(e.rate=0.1,c.rate=0.2,r.on=1/100000,r.off=1/100)

#exp.QED.size<-test[[1]]
#exp.QED.size
#exp.lm<-test[[2]]*e.rate/c.rate
#exp.lm

#test2
