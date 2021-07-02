# setwd("C:/Users/Elisa/Dropbox/Docencia adscripción biomédicas/Temas Selectos - Lic en biología - Bio de sistemas/Prácticas computacionales/4 Redes booleanas estocásticas/esta es la buena")

library(BoolNet)
net=loadNetwork("RedEjemploJuanito.txt")
atr=getAttractors(net) # error pues no se puede calcular los atractores en un 
netD=loadNetwork("RedEjemploJuanitoDET.txt")
atrD=getAttractors(netD)
plotAttractors(atrD)

stateTransition(network=net, state=c(0,0,1), type="probabilistic") # partimos de un atractor


counter=0
for (iteration in 1:100){
  iteration
  XX=  stateTransition(network=net, state=c(0,0,1), type="probabilistic") # partimos de un atractor
 if (isTRUE(XX[[1]]==0 && XX[[2]]==0 && XX[[3]]==1)){
  counter=counter+1}
}

Prob_divergence=(100- counter)/100

## now we want to calculate the probability of movemente between attractors 
## division of labour calculated by the students 