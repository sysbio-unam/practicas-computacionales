# Título: Práctica 3: Redes booleanas estocásticas
# Nombre: Biología de sistemas
# Fecha: 10 de Abril del 2021
#######################################################

# Cargar librerías
library(Rlab)
library(BoolNet)
library(BoolFilter)
library(devtools)
install_github("mar-esther23/boolnet-perturb")
library(dplyr)
library(BoolNetPerturb)

# cargar red
data("p53net_DNAdsb1")

# obtener atractores de la red
at = getAttractors(p53net_DNAdsb1)

# mostrar red
plotAttractors(at)

# primera parte #########################################################
# 1. observar que el sistema determinista tiene un atractor cíclico 
# 2. Hacer simulaciones dinámicas con simulateNetwork sin ruido (q = 0) y con ruido
# variar q y p y observar resultados 

# simula error de observación 
# q --> qué tan grande es el ruido
# p --> probabilidad de que se aplique el ruido
# n.data --> pasos en el tiempo

data = simulateNetwork(p53net_DNAdsb1, n.data=10,
                       p=0.01, obsModel = list(type='Bernoulli', q=0.0))


plotTrajectory(data$X, labels = p53net_DNAdsb1$genes,  
               dataset2 = data$Y, 
               compare = TRUE)


# segunda parte ####################################################

## ahora con la red de arabidopsis
net = loadNetwork("ATH_flower_cell_fate_determination.net.txt")

# obtener atractores
attr = getAttractors(net)
attr

# mostrar atractores
par(mfrow = c(1,1))
plotAttractors(attr)

# ver como depende la diagonal del ruido (p = 0.0, sin ruido)
EL = epigeneticLandscape(net, p = 0.01)

# ¿cuál es la probabilidad de pasar de un atractor 5
# al 3 con un valor de p = 0.01?
EL[5,3]

# tercera parte  #######################################################

# tamaños de las cuencas de atracción
basinsize = seq(1, 10, by=1)

for (ii in 1:length(basinsize)){
  basinsize[ii] = attr$attractors[[ii]]$basinSize
}

n <- 13 # no.variables en el modelo
basinsize_norm <- basinsize/(2^n) # normalizar tamaños de las cuencas

# anadir una diagonal (que se vea que no hay relación lineal)
par(mfrow = c(1,1), din = c(1,1))
plot(diag(as.matrix(EL)), basinsize_norm, pch = 20, cex = 1, 
     xlab = "prob de permanecer en el atractor", 
     ylab = "tamaño de cuenca normalizado",cex.lab = 0.8,las=0.8)
lines(x = c(min(diag(as.matrix(EL))), max(diag(as.matrix(EL)))),
      y = c(min(basinsize_norm), max(basinsize_norm)))

# gráfica entre tamaños de las cuencas y prob de transición de atractor 
# j a atractor 1, j = 1, 2, ..., 10
# efecto de las transiciones no se atribuye directamente al tamaño de las cuencas
plot(EL[,1], basinsize_norm, pch = 20, cex = 1, 
     xlab = "probabilidad de transición",
     ylab = "tamaño de cuenca normalizado", cex.lab = 0.8)



# parte cuatro #########################################################

# si empiezas en atractor 1 después de 100 pasos cuál es la probabilidad de
# permanecer en 1

## simular una cadena de markov, 100 pasos
P = as.matrix(EL)
Nsteps <- 100 # number of steps
pi0 = c(1,0,0,0,0,0,0,0,0,0) # distribución de probabilidad inicial: x(0) =
v = vector("numeric", Nsteps) # creat un vector vacío de tamaño (prealocar)
r = length(pi0) # tamaño para la muestra de la distribución inicial 

v[1] = 1 #sample(r, 1, prob=pi0) 
# la primer entrada del vector es una muestra de 1,2 r
# la probabilidd de obtener cada uno de estos elementos está dada por pi0

# Una sola realización de la cadena
for (i in 2:Nsteps){
  v[i] = sample(r, 1, prob=P[v[i-1],]) # muestrea el nuevo valor: 
  # selecciona el renglón en la matriz de probabilidades que da el vector 
  # de probabilidad de acuerdo con el estado actual
}

# mostrar cadena
matplot(v, type="l", lwd=2.5, col=3,  xlab="t", ylab="Attractor")

## ahora, vamos a iterar varias veces, para sacar el promedio de 
# las veces que cada atractor fue visitado en cada paso de tiempo,
# es decir, la probabilidad del atractor i al tiempo N, para i = 1,2,3
iterations=500

# matriz para guardar las cadenas 
V = matrix(nrow = iterations, ncol = Nsteps)
V[,1] = rep(1, iterations)

for (jj in 1:iterations){
  for (i in 2:Nsteps){
    V[jj, i] = sample(r, 1, prob=P[v[i-1],]) 
  }
}

# probabilidad de estar en atractor 1 al tiempo 75
mean(V[, 75] == 1)


plot(seq(1,100), colMeans(V[,] == 3), type = "l")

