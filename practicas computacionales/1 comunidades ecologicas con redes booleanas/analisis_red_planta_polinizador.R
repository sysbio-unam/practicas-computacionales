# Título: Comunidades ecologicas con redes booleanas

# Nombre: Biología de sistemas

# Fecha: Febreo 2021
#######################################################################################

# cargar librerías
library(BoolNet)

# establecer directorio de trabajo
# setwd("C:/Users/Elisa/Dropbox/Docencia adscripción biomédicas/EXAMENES Y RESPUESTAS A PRÁCTICAS PCBIOL 2021 1/ensamblaje de comunidades")

# cargar red
net <- loadNetwork("red_planta_polinizador.txt")
net

# mostrar la red
plotNetworkWiring(net)

# obtener atractores
attr <- getAttractors(net)
attr

# mostrar atractores
plotAttractors(attr)
plotStateGraph(attr) 

# empezamos en (0,1,1,1,1) y preguntamos a qué atractor vamos a llegar
path <- getPathToAttractor(net, c(0,1,1,1,1))
path

plotSequence(sequence=path)

