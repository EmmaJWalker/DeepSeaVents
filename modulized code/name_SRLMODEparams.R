

name_SRLMODEparams<-function(parameter.values, f.names){
  n.patches<-parameter.values[1]
  parameter.names<-c("n.patches", "c.rate", "self.rec",
                     rep(paste0("extinction.rates",1:n.patches)),
                     rep(paste0("x.coord",1:n.patches)),
                     rep(paste0("y.coord",1:n.patches)),
                     rep(paste0("areas",1:n.patches)),
                     rep(paste0("config",1:n.patches))
                     ,f.names)
  parameters<-list(c(setNames(parameter.values, parameter.names)))

  return(parameters)
}