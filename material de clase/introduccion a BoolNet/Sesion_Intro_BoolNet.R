# Titulo: Introducción a BoolNet 

# Nombre: Biología de sistemas

# Fecha: Febrero 2021
####################################################################################

# entrar a R (Rstudio)
# así ponemos comentarios - esto no lo lo ve R

# revisa en qué directorio estás:
getwd() 

# ir a donde queremos estar (ojo con las diagonales / )
setwd("C:/Users/Elisa/Dropbox/Docencia adscripción biomédicas/Bio Mates PCBIOL 2021 1/Prácticas computacionales/1 Intro a BoolNet")

# utiliza la función help() para pedir ayuda
# ejemplo, si quieres pedir ayuda sobre cómo usar la función print() escribe:
help("print") 

# ojo - necesitamos internet para cargar paquetes
install.packages("BoolNet") 

# cargar (invocar) la librería
library(BoolNet) 


# cargar la red
net <- loadNetwork("red_ejemplo.txt")
net  

# mostrar la red
plotNetworkWiring(net)


# obtener atractores
attr <- getAttractors(net)
attr  

# mostrar atractores
plotAttractors(attr)

# perturbar la red
# asumir que el gen c siempre está apagado
mut = getAttractors(net, genesOFF=c(0,0,1))
mut

# mostrar atractores con la mutación
plotAttractors(mut,drawLegend = F)

# empezamos en (1,0,0) y preguntamos a qué atractor vamos a llegar
getPathToAttractor(net, c(1,0,0)) 


plotStateGraph(attr, piecewise=TRUE)

# nota: actualización asíncrona; fuente de ruido. repetir varias veces y ver si resultados
# coinciden con síncrono, para descartar que resultados sean artefactos del método computacional

att_asymchron = getAttractors(net, type="asynchronous") # asynchronament

plotAttractors(att_asymchron)

# %ginsim.org/
# % www.colomoto.org/


## Ahora veremos qué pasa con un atractor cíclico
cell_cyle <- loadNetwork("cellcycle.txt")
cell_cyle

# mostrar red
plotNetworkWiring(cell_cyle)


# obtener atractores
attrCC <- getAttractors(cell_cyle)
attrCC 

# mostrar atractores
plotAttractors(attrCC)
