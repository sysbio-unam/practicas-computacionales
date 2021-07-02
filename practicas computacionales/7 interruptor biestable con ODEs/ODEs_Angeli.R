# Título: Análisis de un sistma de ecuaciones diferenciales acopladas y no lineales

# Nombre: Biología de sistemas

# Fecha: Febrero 2021
#########################################################################################
##. Angeli, D., Ferrell, J. E. & Sontag, E. D. Detection of multistability, bifurcations, and hysteresis in a large class of biological positive-feedback systems. PNAS 101, 1822-7 (2004).

# Borramos todo
#rm(list=ls())

# Nos ubicamos donde queremos
# setwd(...)

# cargar la paquetería que necesitamos
library(deSolve)
library(phaseR)

# Declarar los valores de parámetros que permanecen constantes
alpha1 = 1
alpha2 = 1
beta1 = 200
beta2 = 10
gamma1 = 4
gamma2 = 4
K1 = 30
K2 = 1
v = 1

# establecer modelo
Angeli2004 <- function(t, y, parms) {  
    #         1      2       3      4       5     6      7   8   9
    # parms=(alpha1, alpha2, beta1, beta2, gamma1, gamma2, K1, K2, v)
         dX <- parms[1]*(1-y[1])-parms[3]*y[1]*(parms[9]*y[2])^parms[5]/(parms[7]+(parms[9]*y[2])^parms[5]);  
         dY <- parms[2]*(1-y[2])-parms[4]*y[2]*y[1]^parms[6]/(parms[8]+y[1]^parms[6]);
         list(c(dX,dY))
}

# condiciones iniciales
ini_1 <- c(0,0)
ini_2 <- c(0,0.9)

# tiempo de integración
tspan <- seq(from = 0, to = 10, by = 0.01)

############## PREGUNTA 1: DINÁMICA DEL SISTEMA ######################

# parámetro de bifurcación 
for (v in c(0.75, 1, 1.9)){
        
        # condiciones iniciales
        ini_1 <- c(0,0)
        ini_2 <- c(0,0.9)
 
        parms = c(alpha1, alpha2, beta1, beta2, gamma1, gamma2, K1, K2, v)

        # resolver EDOs
        out1 <- ode(y = ini_1, times = tspan, func = Angeli2004, parms = parms)
        
        # mostrar simulación
        plot(out1[,1], out1[,2],type = "l", ylim = c(0,1),
             col = "red", xlab = "Time", ylab = "X(t)", 
             main = paste("v = ", toString(v), sep=" "), lwd = 2)
        
        # ahora con la segunda condición inicial
        out2 <- ode(y = ini_2, times = tspan, func = Angeli2004, parms = parms)
        
        # mostrar sobre resultados sobre la gráfica anterior
        lines(out2[,1], out2[,2], type = "l", col = "blue", lwd = 2)
        
        legend("right", legend = c(paste0("IC 1: ","(",toString(ini_1),")"), paste0("IC 2: ","(",toString(ini_2),")")), 
               col = c(1,2), lty = 1, bty = "n", cex = 0.8, lwd = 2)
        
        # Añadamos a este diagrama de espacio fase un campo vectorial
        Angeli2004.flowField <- flowField(Angeli2004, xlim = c(0, 1), ylim = c(0, 1), parameters = parms, 
                                          points = 10, add = FALSE, main = paste0("v = ", v))
        Angeli2004.trajectory <- trajectory(Angeli2004, y0 = ini_1, tlim = c(0,10), parameters = parms, col = "blue")
        Angeli2004.trajectory <- trajectory(Angeli2004, y0 = ini_2, tlim = c(0,10), parameters = parms, col = "red")
}


############ PREGUNTA 2: Cuencas de atracción ###############

# parámetro de bifurcación 
for (v in c(1, 1.6)){
        
    parms["v"] <- v
    
    # Añadamos a este diagrama de espacio fase un campo vectorial
    flowField(Angeli2004, xlim = c(0, 1), ylim = c(0, 1), 
              parameters = parms, points = 10, add = FALSE,
              main = paste0("v = ", v))
    
    # Genera n condiciones iniciales al azar, pero sobre el márgen ([x=0,1; y=rand] y vice versa)
    for (ii in seq(1,20,1) ){
            
      # generar tres números aleatorios
      r1 = runif(1)
      r2=runif(1)
      r3=runif(1)
      
      if (r1 < 0.5) {
        ini = c(as.numeric(r2 < 0.5), r3)
        } else {
        ini = c(r3, as.numeric(r2 < 0.5))
        }
      
        trajectory(Angeli2004, y0 = ini, tlim = c(0,10), parameters = parms, col = "blue")
    }
}



############ PREGUNTA 3: Señales de alerta temprana ###############

LineWidth = 1
LineColor = 2

# parámetro de bifurcación - 
for (v in seq(0.2,1,0.1)){
        
        parms[9] <- v
        
        # ¡A integrar!
        out <- ode(y = ini_2, times = tspan, func = Angeli2004, parms = parms)
        
        if (v == 0.2){
                
                plot(out[,1], out[,3],type = "l", ylim=c(0,1),
                     col = 1, xlab = "Time", ylab = "X(t)", 
                     lwd = 1, main="Alentamiento crítico")
                
        } 
        
        else {
                
                lines(out[,1], out[,3],type = "l", col= LineColor, lwd = LineWidth)
        }
        
        LineWidth = LineWidth + 0.5
        LineColor = LineColor + 1
}

legend("right", legend = paste0("v = ",seq(0.2,1,0.1)), col = seq(1,11), cex = 0.8, 
       bty = "n", lty  = 1, lwd = 3)

############## PREGUNTA 4: Diagrama de bifurcación ######################

# cargar el código anexo:
source('Grind.r') 

model <- function(t, state, parms){  
  with(as.list(c(state,parms)), {
    dx = alpha1*(1-x)-beta1*x*(v*y)^gamma1/(K1+(v*y)^gamma1)
    dy = alpha2*(1-y)-beta2*y*x^gamma2/(K2+x^gamma2)
    return(list(c(dx, dy)))
  })
}
  
# declarar los valores de parámetros
p <-  c(alpha1=1, alpha2=1, beta1=200, beta2=10, gamma1=4, gamma2=4, K1=30, K2=1, v=1)

# condición inicial
s <- c(x=0,y=0)

# graficar plano de fase
plane(xmax=4)

mid <- newton(s,plot=T)
low <- newton(c(x=1,y=0),plot=T)
hig <- newton(c(x=0,y=1),plot=T)

continue(state=hig, parms=p, odes=model, x="v", step=0.001, xmin=0, xmax=2,y="y", ymin=0, ymax=1.1) # log="", time=0, positive=TRUE, add=TRUE)
continue(state=low, parms=p, odes=model, x="v", step=0.001, xmin=0, xmax=2,y="y", ymin=0, ymax=1.1, log="", time=0, positive=TRUE, add=TRUE)
continue(state=mid, parms=p, odes=model, x="v", step=0.001, xmin=0, xmax=2,y="y", ymin=0, ymax=1.1, log="", time=0, positive=TRUE, add=TRUE)
