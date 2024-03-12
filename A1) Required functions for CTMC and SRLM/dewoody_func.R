dewoody<-function(e.rate, c.rate, r.on, r.off){
  dw.delta<-(e.rate+r.off)/(c.rate*(r.on/(r.on+r.off))) #lambda M threshold (as defined by delta) in a dynamic landscape as defined by dewoody et al. 2005 
  s<-r.on/(r.on+r.off) #pr(patch is on)
  r.o<-c.rate*dw.delta*s/(e.rate+r.off)
  return(list(dw.delta, s, r.o))
}
#test2<-dewoody(e.rate=0.1,c.rate=0.2,r.on=1/100000,r.off=1/100)
#test2