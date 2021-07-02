# Título: Adaptación perfecta

# Nombre: Biología de sistemas

# Fecha: 9 de diciembre del 2020
#######################################################################################

setwd("C:/Users/Elisa/Dropbox/Docencia adscripción biomédicas/Bio Mates PCBIOL 2021 1/Prácticas computacionales/Prácticas_Computacionales_ODEs_en_R")

library(deSolve)
library(sensitivity)
library(checkmate)

library(ODEnetwork)
library(ODEsensitivity)

##### Lotka-Volterra  #####
# modelo:
LVmod <- function(Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    Ingestion    <- rIng  * Prey * Predator
    GrowthPrey   <- rGrow * Prey * (1 - Prey/K)
    MortPredator <- rMort * Predator
    
    dPrey        <- GrowthPrey - Ingestion
    dPredator    <- Ingestion * assEff - MortPredator
    
    return(list(c(dPrey, dPredator)))
  })
}

# Los k = 5 parámetros rG, rI, rM, kAE y K se consideran variables de entrada para el análisis 
# de sensibilidad. Por lo tanto, analizaremos la sensibilidad de la población de presas y depredadores 
# frente a cambios en estos 5 parámetros.


# Los parámetros que se incluirán en el análisis de sensibilidad y sus límites superior e inferior:
LVpars  <- c("rIng", "rGrow", "rMort", "assEff", "K")
# normalmente alrededor del valor nominal (output de la optimizacion)
LVbinf <- c(0.05, 0.05, 0.05, 0.05, 1)
LVbsup <- c(1.00, 3.00, 0.95, 0.95, 20)

# El valor inicial de las variables de estado
LVinit  <- c(Prey = 1, Predator = 2)
# Los tiempos de interes
LVtimes <- c(0.01, seq(1, 5, by = 0.1))
set.seed(59281)
# Análisis de sensibilidad de Sobol '(aquí solo con n = 500, pero se recomienda n = 1000):
# Advertencia: ¡El siguiente código puede tardar mucho!

LVres_sobol <- ODEsobol(mod = LVmod,
                        pars = LVpars,
                        state_init = LVinit,
                        times = LVtimes,
                        n = 500,
                        rfuncs = "runif",
                        rargs = paste0("min = ", LVbinf,
                                       ", max = ", LVbsup),
                        sobol_method = "Martinez",
                        ode_method = "lsoda",
                        parallel_eval = TRUE,
                        parallel_eval_ncores = 2)

str(LVres_sobol, vec.len = 3, give.attr = FALSE)

# Es una lista de clase "ODEmorris" con un elemento para cada variable de estado (aquí, Prey y Predator).
# Esos elementos son matrices de 3 * longitud (LVpars) + 1 filas y columnas de longitud (LVtimes).

# La primera fila contiene una copia de todos los puntos de tiempo.
# Las otras filas contienen los 3 índices de sensibilidad de Sobol 
# para todos los parámetros en los 51 puntos de tiempo.

plot(LVres_sobol,  state_plot = "Predator")
plot(LVres_sobol,  state_plot = "Prey")


# ahora tracemos la distribución de índices por tiempo
library(vioplot)

x1 <- LVres_sobol$Prey$S[2, ]
x2 <- LVres_sobol$Prey$S[3, ]
x3 <- LVres_sobol$Prey$S[4, ]
x4 <- LVres_sobol$Prey$S[5, ]
x5 <- LVres_sobol$Prey$S[6, ]

par(pty="s") 
vioplot(x1, x2, x3,x4, x5, names=LVpars,  col="black")
title("Prey sensitivity")
LVres_sobol$Prey$S

y1 <- LVres_sobol$Predator$S[2, ]
y2 <- LVres_sobol$Predator$S[3, ]
y3 <- LVres_sobol$Predator$S[4, ]
y4 <- LVres_sobol$Predator$S[5, ]
y5 <- LVres_sobol$Predator$S[6, ]

par(pty="s") 
vioplot(y1, y2, y3,y4, y5, names=LVpars,  col="red")
title("Predator sensitivity")

