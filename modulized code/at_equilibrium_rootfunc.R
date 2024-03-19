# Provides a stopping condition to solving the SRLM when the change in patch 
# occupancy is less than 1e-4 and therefore the SRLM is essentially at 
# equilibrium

# INPUTS:
# t: a time sequence for the equations to be solved over
# p: occupancy of patches as described by the equations set up by this function
# parameters: a vector of parameters used to construct the equations (must be a 
# named list with unique names for every parameter, including all elements of 
# any vectors or matrices)

at.equilibrium <- function(t, p, parameters) {
  dP <- unlist(SRLM.ODE(t, p, parameters))
  Norm(dP) - 1e-4
}