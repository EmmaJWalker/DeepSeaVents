#stopping function for when the change in patch occupancy is less than 1e-4
#and therefore the SRLM is essentially at equilibrium

at.equilibrium <- function(t, p, parameters) {
  dP <- unlist(SRLM.ODE(t, p, parameters))
  Norm(dP) - 1e-4
  #if (p <= 1e-4){
  #  return(p - p)
  #}
  #if (Norm(dP) <= 1e-4){
  #  return(Norm(dP) - Norm(dP))
  #}
}

#too.small<-function(t, p, parameters) {
#  dP <- unlist(SRLM.ODE(t, p, parameters))
#  return(p - 1e-4)
#}