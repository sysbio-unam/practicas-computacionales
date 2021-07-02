# Título: Señales de alerta temprana en sistemas biológicos bifurcantes

# Nombre: Biología de sistemas

# Fecha: Febrero 2021
#########################################################################################
##. Angeli, D., Ferrell, J. E. & Sontag, E. D. Detection of multistability, bifurcations, and hysteresis in a large class of biological positive-feedback systems. PNAS 101, 1822-7 (2004).
library(ggplot2)
#source("Grind.R")


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

# Diagrama de bifurcación para Y
continue(state=hig, parms=p, odes=model, x="v", step=0.001, xmin=0, xmax=2,y="y", ymin=0, ymax=1.1) # log="", time=0, positive=TRUE, add=TRUE)
continue(state=low, parms=p, odes=model, x="v", step=0.001, xmin=0, xmax=2,y="y", ymin=0, ymax=1.1, log="", time=0, positive=TRUE, add=TRUE)
continue(state=mid, parms=p, odes=model, x="v", step=0.001, xmin=0, xmax=2,y="y", ymin=0, ymax=1.1, log="", time=0, positive=TRUE, add=TRUE)

# Diagrama de bifurcación para Y
continue(state=hig, parms=p, odes=model, x="v", step=0.001, xmin=0, xmax=2,y="x", ymin=0, ymax=1.1) # log="", time=0, positive=TRUE, add=TRUE)
continue(state=low, parms=p, odes=model, x="v", step=0.001, xmin=0, xmax=2,y="x", ymin=0, ymax=1.1, log="", time=0, positive=TRUE, add=TRUE)
continue(state=mid, parms=p, odes=model, x="v", step=0.001, xmin=0, xmax=2,y="x", ymin=0, ymax=1.1, log="", time=0, positive=TRUE, add=TRUE)


# 2. Seleccionar estado estacionario como condición inicial 
s <- hig

# resolver ED
times <- seq(0,100,1)
out <- run(timeplot = F,table = T, times = times)
head(out)

# mostrar fase estacionaria 
ggplot(out, aes(x = time, y = x)) + 
  geom_line(lwd = 2, col = "darkblue") +
  geom_line(data = data.frame(time = times, x = rep(hig[1], length(times))), 
            col = "red", lty = 2, lwd = 2) +
  xlab("Tiempo") +
  ylab("x") + 
  ylim(0,0.2)+
  theme_bw() + 
  theme(legend.position = "none",
        plot.title = element_text(size=25,face = "italic"),
        axis.text=element_text(size=10),
        axis.title=element_text(size=10),
        legend.title = element_text(size=20), 
        legend.text = element_text(size=20))

# mostrar fase estacionaria 
ggplot(out, aes(x = time, y = y)) + 
  geom_line(lwd = 2, col = "darkblue") +
  geom_line(data = data.frame(time = times, y = rep(hig[2], length(times))), 
            col = "red", lty = 2, lwd = 2) +
  xlab("Tiempo") +
  ylab("y") + 
  ylim(0,1.1) +
  theme_bw() + 
  theme(legend.position = "none",
        plot.title = element_text(size=25,face = "italic"),
        axis.text=element_text(size=10),
        axis.title=element_text(size=10),
        legend.title = element_text(size=20), 
        legend.text = element_text(size=20))


# 3. Simular empleando el algoritmo de Euler-Mayurama 
# simular utilizando el método de Euler-Maruyama
out <- run(state = s,after = "state<-state+rnorm(2,mean=0,sd=0.01)", timeplot = F, table = T)
head(out)

# mostrar fase estacionaria 
ggplot(out, aes(x = time, y = x)) + 
  geom_line(lwd = 2, col = "darkblue") +
  geom_line(data = data.frame(time = times, x = rep(hig[1], length(times))), 
            col = "red", lty = 2, lwd = 2) +
  xlab("Tiempo") +
  ylab("x") + 
  theme_bw() + 
  theme(legend.position = "none",
        plot.title = element_text(size=25,face = "italic"),
        axis.text=element_text(size=10),
        axis.title=element_text(size=10),
        legend.title = element_text(size=20), 
        legend.text = element_text(size=20))

# mostrar fase estacionaria 
ggplot(out, aes(x = time, y = y)) + 
  geom_line(lwd = 2, col = "darkblue") +
  geom_line(data = data.frame(time = times, y = rep(hig[2], length(times))), 
            col = "red", lty = 2, lwd = 2) +
  xlab("Tiempo") +
  ylab("y") +
  theme_bw() + 
  theme(legend.position = "none",
        plot.title = element_text(size=25,face = "italic"),
        axis.text=element_text(size=10),
        axis.title=element_text(size=10),
        legend.title = element_text(size=20), 
        legend.text = element_text(size=20))

# 4. Realizar una simulación con un valor de tend muy grande para obtener la 
# distribución de los valores finales de x y y

#función para obtener la distribución de "x" y "y" en el equilibrio
get_distribution = function(t_end, delta_t, v, noise = 0.01) {
  
  p["v"] = v
  
  time = seq(0,t_end)
  
  after = paste0("state<-state+rnorm(2,mean=0,sd=",noise,")")
  
  out <- run(state = s, times = time, after = after, 
             timeplot = F, table = T)
  
  index = out[,1]%%delta_t == 0
  
  data_dist = out[index,]
  
  return(data_dist)
  
}
# 0 -> 10000
# delta_t = 100

data_dist = get_distribution(100,1,p["v"],1)
hist(data_dist[,2], xlab = "x", ylab = "fecuencia", main = "Distribución de xss")
var(data_dist[,2])
abline(v = hig[1], col = "red")

hist(data_dist[,3], xlab = "y", ylab = "fecuencia", main = "Distribución de yss")
var(data_dist[,3])
abline(v = hig[2], col = "blue")

# función para obtener la varianza de "x" y "y"
get_var = function(t_end, delta_t, v, noise) {
  data_dist = get_distribution(t_end, delta_t, v,noise)
  x = data_dist[,2]
  y = data_dist[,3]
  varience = c(var_x = var(x), var_y = var(y))
  return(varience)
}

get_var(1e2,1,p["v"],1)


# 5. Obtener la varianza de "xss" y "yss" para diferentes valores de v

v = seq(0,2,0.01)

VAR_XY = matrix(nrow = length(v), ncol = 3)
VAR_XY[,1] = v

for (i in 1:nrow(VAR_XY)) {
  
  VAR_XY[i,c(2,3)] = get_var(t_end = 100,delta_t = 1,v = v[i],noise = 0.1)
  
}

var_x = VAR_XY[,2]
var_y = VAR_XY[,3]
plot(v, var_x, xlab = "parámetro de bifurcación (v)", ylab = "varianza de xss")
abline(v = 1.796859, col = "red", lwd = 2)
abline(v = 0.8315781, col = "red", lwd = 2)
plot(v, var_y, xlab = "parámetro de bifurcación (v)", ylab = "varianza de yss")
abline(v = 1.796859, col = "red", lwd = 2)
abline(v = 0.8315781, col = "red", lwd = 2)
