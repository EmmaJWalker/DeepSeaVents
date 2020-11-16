SRLM.ODE<-function(t,p,parameters){
  with(c(p, parameters), {
    ###############################################################################################
    #parameters index just for reference
    ###############################
    #n.patches<-parameters[1]
    #c.rate<-parameters[2]
    #self.rec<-parameters[3]
    #E<-parameters[4:(3+n.patches)]
    #x.coord<-parameters[(4+n.patches):(3+n.patches*2)]
    #y.coord<-parameters[(4+n.patches*2):(3+n.patches*3)]
    #A<-parameters[(4+n.patches*3):(3+n.patches*4)]
    #config<-parameters[(4+n.patches*4):(3+n.patches*5)]
    #f.vals<-parameters[(4+n.patches*5):(3+n.patches*5+(n.patches*n.patches))]
    ###############################################################################################
    #Set up 1D and 2D matrices of parameters for input into the ODEs
    ###############################################################################################
    P=matrix(p, n.patches, 1)
    AREA=matrix(parameters[[1]][(4+n.patches*3):(3+n.patches*4)], n.patches, 1)
    CONFIG=matrix(parameters[[1]][(4+n.patches*4):(3+n.patches*5)], n.patches, 1)
    EXRATE=matrix(parameters[[1]][4:(3+n.patches)], n.patches, 1)
    DISPKERNEL=matrix(parameters[[1]][(4+n.patches*5):(3+n.patches*5+(n.patches*n.patches))], n.patches, n.patches)
    ###############################################################################################
    #create a 1D matrix to hold the right hand side of the ODE equations (dp)
    dP=matrix(rep(NA, n.patches), n.patches, 1) 
    interactions<-0 #we don't have any interactions in this model at the moment so ignore this term
    ###############################################################################################
    #Loops to construct the system of ODEs
    for (i in 1:n.patches){ #for each patch
      for (j in 1:n.patches){ #to each other patch (including itself)
        if (i!=j & CONFIG[i,1]!=0){ #if it's not itself and it's not off
          dP[i,1]<-interactions + c.rate*AREA[j,1]*DISPKERNEL[i,j]*P[j,1]*(1-P[i,1]) - EXRATE[i,1]*P[i,1]} #then the pop dynamics occur
        if (i==j & CONFIG[i,1]!=0){ #if the patch is itself and it's not off
          if (self.rec==1){ #if self recruitment is on 
            dP[i,1]<-interactions + c.rate*AREA[j,1]*DISPKERNEL[i,j]*P[j,1]*(1-P[i,1]) - EXRATE[i,1]*P[i,1]} #pop dynamics occur
        }
        else if ((i==j & self.rec==0) | CONFIG[i,]==0){ #otherwise, if self recruitment is off
          #or the patch is just off
          dP[i,1]<-0}#pop dynamics don't occur and the equation is just = 0
        interactions<-dP[i,1]
      }
    }
    ###############################################################################################
    dP.names<-c(rep(paste0("dP",1:n.patches)))
    dP.values<-dP
    RHS<-list(c(setNames(dP.values, dP.names)))
    return(RHS)
  }) # end with(as.list ...
}

