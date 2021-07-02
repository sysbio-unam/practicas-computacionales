# Título: Simulación de ruido extrínseco e intrínseco con modelo de adaptación perfecta 

# Nombre: Biología de sistemas

# Fecha: 
#######################################################################################

# cargar Grind.R
#source("Grind.R")

# cargar librerías 
library(ggplot2)
library(GillespieSSA)


# 1. Calcular estado estacionario #####################################################
# calcular estado estacionario 
steady_state <- function(p) {
        with(as.list(p), {
                
                xss <- k3*S/k4
                
                Rss <- k1*k4/(k2*k3)
                
                return(list(xss = xss, Rss = Rss))
        })
}

# establecer parámetros 
p <- c(k1 = 5, k2 = 13, k3 = 1, k4 = 0.3, S = 2.5)

# calcular estado estacionario 
ss <- steady_state(p)
ss

# escribier modelo 
model <- function(time, state, parms) {
        with(as.list(c(state, parms)), {
                
                dR <- k1*S - k2*R*x
                dx <- k3*S - k4*x 
                
                return(list(c(dR, dx)))
        })
}

# condiciones iniciales 
s <- c(R = 0, x = 0)

# tiempo de integración 
times <- seq(0,20,0.1)

# resolver ED
out <- run(timeplot = F,table = T, times = times)
head(out)

# mostrar fase transitoria 
ggplot(out, aes(x = time, y = R)) + 
        geom_line(lwd = 2, col = "darkblue") +
        geom_line(data = data.frame(time = times, R = rep(ss$Rss, length(times))), 
                  col = "red", lty = 2, lwd = 2) +
        xlab("Tiempo") +
        ylab("R") + 
        theme_bw() + 
        theme(legend.position = "none",
              plot.title = element_text(size=25,face = "italic"),
              axis.text=element_text(size=20),
              axis.title=element_text(size=20),
              legend.title = element_text(size=20), 
              legend.text = element_text(size=20))

# mostrar fase transitoria 
ggplot(out, aes(x = time, y = x)) + 
        geom_line(lwd = 2, col = "darkblue") +
        geom_line(data = data.frame(time = times, x = rep(ss$xss, length(times))), 
                  col = "red", lty = 2, lwd = 2) +
        xlab("Tiempo") +
        ylab("x") + 
        theme_bw() + 
        theme(legend.position = "none",
              plot.title = element_text(size=25,face = "italic"),
              axis.text=element_text(size=20),
              axis.title=element_text(size=20),
              legend.title = element_text(size=20), 
              legend.text = element_text(size=20))

# 2. Seleccionar estado estacionario como condición inicial  #####################
s <- c(R = ss$Rss, x = ss$xss)

# resolver ED
out <- run(timeplot = F,table = T, times = times)
head(out)

# mostrar fase estacionaria 
ggplot(out, aes(x = time, y = R)) + 
        geom_line(lwd = 2, col = "darkblue") +
        geom_line(data = data.frame(time = times, R = rep(ss$Rss, length(times))), 
                  col = "red", lty = 2, lwd = 2) +
        xlab("Tiempo") +
        ylab("R") + 
        theme_bw() + 
        theme(legend.position = "none",
              plot.title = element_text(size=25,face = "italic"),
              axis.text=element_text(size=20),
              axis.title=element_text(size=20),
              legend.title = element_text(size=20), 
              legend.text = element_text(size=20))

# mostrar fase estacionaria 
ggplot(out, aes(x = time, y = x)) + 
        geom_line(lwd = 2, col = "darkblue") +
        geom_line(data = data.frame(time = times, x = rep(ss$xss, length(times))), 
                  col = "red", lty = 2, lwd = 2) +
        xlab("Tiempo") +
        ylab("x") + 
        theme_bw() + 
        theme(legend.position = "none",
              plot.title = element_text(size=25,face = "italic"),
              axis.text=element_text(size=20),
              axis.title=element_text(size=20),
              legend.title = element_text(size=20), 
              legend.text = element_text(size=20))

# 3. Simular ruido extr{inseco empleando el algoritmo de Euler-Mayurama ###########################
# simular utilizando el método de Euler-Maruyama
out <- run(state = s,after = "state<-state+rnorm(2,mean=0,sd=0.01)", timeplot = F, table = T)
head(out)

# mostrar fase estacionaria 
ggplot(out, aes(x = time, y = R)) + 
        geom_line(lwd = 2, col = "darkblue") +
        geom_line(data = data.frame(time = times, R = rep(ss$Rss, length(times))), 
                  col = "red", lty = 2, lwd = 2) +
        xlab("Tiempo") +
        ylab("R") + 
        ylim((c(0, 0.3)))+
        theme_bw() + 
        theme(legend.position = "none",
              plot.title = element_text(size=25,face = "italic"),
              axis.text=element_text(size=20),
              axis.title=element_text(size=20),
              legend.title = element_text(size=20), 
              legend.text = element_text(size=20))

# mostrar fase estacionaria 
ggplot(out, aes(x = time, y = x)) + 
        geom_line(lwd = 2, col = "darkblue") +
        geom_line(data = data.frame(time = times, x = rep(ss$xss, length(times))), 
                  col = "red", lty = 2, lwd = 2) +
        xlab("Tiempo") +
        ylab("x") +
        ylim(c(8,9))+
        theme_bw() + 
        theme(legend.position = "none",
              plot.title = element_text(size=25,face = "italic"),
              axis.text=element_text(size=20),
              axis.title=element_text(size=20),
              legend.title = element_text(size=20), 
              legend.text = element_text(size=20))

# 4. Realizar varias simulaciones para obtener la distribución de los valores finales de R y x

# número de iteraciones
iter <- 100

# matriz para guardar los valores finales de R y x en cada simulación
S_end <- matrix(nrow = iter, ncol = 2)

for (i in 1:nrow(S_end)) {
        # simular utilizando el método de Euler-Maruyama
        out <- run(state = s,after = "state<-state+rnorm(2,mean=0,sd=0.01)", timeplot = F, table = T)
        
        print(
                
                # mostrar fase estacionaria de R
                ggplot(out, aes(x = time, y = R)) + 
                        geom_line(lwd = 2, col = "darkblue") +
                        geom_line(data = data.frame(time = times, R = rep(ss$Rss, length(times))), 
                                  col = "red", lty = 2, lwd = 2) +
                        xlab("Tiempo") +
                        ylab("R") + 
                        ylim((c(0, 0.3)))+
                        theme_bw() + 
                        theme(legend.position = "none",
                              plot.title = element_text(size=25,face = "italic"),
                              axis.text=element_text(size=20),
                              axis.title=element_text(size=20),
                              legend.title = element_text(size=20), 
                              legend.text = element_text(size=20))
        )
        
        print(
                
                # mostrar fase estacionaria de x
                ggplot(out, aes(x = time, y = x)) + 
                        geom_line(lwd = 2, col = "darkblue") +
                        geom_line(data = data.frame(time = times, x = rep(ss$xss, length(times))), 
                                  col = "red", lty = 2, lwd = 2) +
                        xlab("Tiempo") +
                        ylab("x") +
                        ylim(c(8,9))+
                        theme_bw() + 
                        theme(legend.position = "none",
                              plot.title = element_text(size=25,face = "italic"),
                              axis.text=element_text(size=20),
                              axis.title=element_text(size=20),
                              legend.title = element_text(size=20), 
                              legend.text = element_text(size=20))
        )
        

        
        # guardar valor final de GFP
        S_end[i,] <- as.matrix(out[nrow(out), c(2,3)])
        
}

# convertir en data frame
S_end <- as.data.frame(S_end)


# mostrar distribución de valores finales de R
ggplot(S_end, aes(x = V1)) +
        geom_histogram(fill = "skyblue", col = "darkblue", alpha = 0.5,bins = 30) +
        geom_vline(xintercept = ss$Rss, col = "red", lty = 2, lwd = 2) +
        xlab(expression(R[final])) +
        ylab("frecuencia") + 
        theme_bw() + 
        theme(legend.position = "none",
              plot.title = element_text(size=25,face = "italic"),
              axis.text=element_text(size=20),
              axis.title=element_text(size=20),
              legend.title = element_text(size=20), 
              legend.text = element_text(size=20))

# mostrar distribución de valores finales de x
ggplot(S_end, aes(x = V2)) +
        geom_histogram(fill = "skyblue", col = "darkblue", alpha = 0.5,bins = 30) +
        geom_vline(xintercept = ss$xss, col = "red", lty = 2, lwd = 2) +
        xlab(expression(x[final])) +
        ylab("frecuencia") + 
        theme_bw() + 
        theme(legend.position = "none",
              plot.title = element_text(size=25,face = "italic"),
              axis.text=element_text(size=20),
              axis.title=element_text(size=20),
              legend.title = element_text(size=20), 
              legend.text = element_text(size=20))

####### Ahora con ruido intrínseco ###########################################
## Lotka predator-prey model

# dRdt = k1*S - k2*R*x
# dxdt = k3*S - k4*x

p <- c(k1 = 5, k2 = 13, k3 = 1, k4 = 0.3, S = 2.5)

nu <- matrix(c(+1,-1,0,0,0,0,+1,-1), nrow = 2, byrow=TRUE)

x0 <- c(R = 10, x = 10)

a <- c("k1*S","k2*R*x","k3*S","k4*x")

tf <- 10


out <- ssa(
        x0 = x0,
        a = a,
        nu = nu,
        parms = p,
        tf = tf,
        method = ssa.d(),
        simName = "dimero",
        verbose = FALSE,
        consoleInterval = 1
) 

ssa.plot(out, show.title = TRUE, show.legend = FALSE)

