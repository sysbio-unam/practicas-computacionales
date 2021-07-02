# Título: Optimización paramétrica. Adaptación perfecta

# Nombre: Biología de sistemas

# Fecha: 9 de diciembre del 2020
#######################################################################################
## parameter estimation in R
# http://tbb.bio.uu.nl/rdb/grindR/
#rm(list=ls())


# (1) cargar Grind.R
source("Grind.R")

# (2) definir función (esto es, el modelo)
model <- function(time, state, parms) {
        with(as.list(c(state, parms)), {
                
                dR <- k1*S - k2*R*x
                dx <- k3*S - k4*x 
                
                return(list(c(dR, dx)))
        })
}

# (3) Declarar primer solución para el valor de los parámetros
p <- c(k1 = 5, k2 = 13, k3 = 1, k4 = 0.3, S = 2.5)

# Establecer condiciones iniciales (las cuales, efectivamente, son parámetros)
s <- c(R = 0, x = 0)

# (4) Cargar los datos
data <- data.frame(time = c(0, 0.5, 1, 2, 4, 6), 
                   R = c(0.0126,1.6059,0.8196,0.7323,0.5339,0.5613))
# El nombre de las columnas en el data.frame con los datos debe ser igual al nombre
# de las variables en el modelo

# (5) Correr la optimización
w <- c("k1", "k2", "k3", "k4") # nombres de los parámetros a optimizar
f <- fit(legend = FALSE, free = w, tstep = 0.0001, method = "BFGS")

# Podemos revisar los rangos de confianza, valor p, etc.
summary(f)

# Para obtener solo los parámetros óptimos podemos usar f$par
f$par

# Para obtener la sumar de los errores al cuadrado usamos f$ssr
f$ssr  
