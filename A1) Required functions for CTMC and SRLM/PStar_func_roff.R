#P* FUNCTION
#################################################################################################################
#calculates the probability individual patches are occupied at equilibrium

#landscape = a dataframe created by the create landscape function which specifies patch locations, areas, 
#and and an interpatch distance matrix
#a = 1/(average dispersal distance) for a species
#delta = the ratio of the extinction to colonization rate of a species
#iterations = the number of iterations for which the iterative function f is iterated, 
#greater iterations = greater accuracy, less iterations = lower accuracy

#NOTE: delta can be factored out of this calculation by letting delta = 1
#################################################################################################################
pstar.function<-function(r.off, landscape, alpha, e.rate, c.rate, self.rec, gamma, epsilon, iterations){
  n.patches<-length(landscape$patch.ID)
  x.coord<-landscape$x.coord
  y.coord<-landscape$y.coord
  AREA<-landscape$areas
  delta<-(e.rate+r.off)/c.rate
  
  # construct the distance matrix
  dist.mat<-matrix(rep(0,n.patches*n.patches),n.patches,n.patches)
  for (i in 1:n.patches){
    for (j in 1:n.patches){
      dist.mat[i,j]<-sqrt((x.coord[i]-x.coord[j])^2+(y.coord[i]-y.coord[j])^2)
    }
  }
  #Getting the f(dij) values to plug into the rhs
  DISPKERNEL<-disp.kernel(dist.mat,alpha,gamma,epsilon,self.rec,n.patches)
  
  #STEP 1: SET UP OF MATRIX M
  M<-matrix(rep(NA, n.patches*n.patches), n.patches,n.patches)
  for (i in 1:n.patches){ #for each patch
    for (j in 1:n.patches){ #to each other patch (including itself)
      M[i,j]<-AREA[j]*DISPKERNEL[i,j]
    }
  }
  
  #ITERATING FUNCTION TO FIND P*
  p<-rep(NA, n.patches*iterations);dim(p)<-c(iterations,n.patches)
  p[1,1:n.patches]<-rep(0.1,n.patches)
  for (t in 1:(iterations-1)){ #iterate the following for a given number of iterations
    P<-p[t,] #use p of previous iteration
    g<-(M%*%P) #SRLM g function from Ovaskainen and Hanski, 2001
    f<-(g/(g+delta)) #SRLM f function from Ovaskainen and Hanski, 2001, 
    #which when iterated many times will give p*
    p[(t+1),]<-t(f)} 
  #set p for the next iteration = to the output value of p for this iteration
  p.star<-p[iterations,] #p* = the final iterations value of p after many iterations
  return(p.star)
}
#END P* FUNCTION
#################################################################################################################