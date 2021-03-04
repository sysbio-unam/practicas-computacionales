#learn how to infer parameters with Stan + and ODE model

#https://mc-stan.org/users/documentation/case-studies/lotka-volterra-predator-prey.html


getwd()

# load the data. I copy-pasted it in a Note pad and saved it in a Comma sepparated file from
#https://github.com/stan-dev/example-models/blob/master/knitr/lotka-volterra/hudson-bay-lynx-hare.csv
#it is a very simple file; time and values for X and Y variable

lynx_hare_df <-
  read.csv("hudson-bay-lynx-hare.csv",
           comment.char="#")

lynx_hare_df


# first just look at the data
library(ggplot2)
library(gridExtra)
library(knitr)
knitr::opts_chunk$set(cache = TRUE)
knitr::opts_chunk$set(tidy = FALSE, cache.extra = packageVersion('tufte'))
knitr::opts_chunk$set(comment = "")
library(reshape)
library(tufte)


lynx_hare_melted_df <- melt(as.matrix(lynx_hare_df[, 2:3]))

colnames(lynx_hare_melted_df) <- c("year", "species", "pelts")

lynx_hare_melted_df$year <-
  lynx_hare_melted_df$year +
  rep(1899, length(lynx_hare_melted_df$year))

ggplot(data = lynx_hare_melted_df,
       aes(x = year, y = pelts, color = species)) +
  geom_vline(xintercept = 0, color = "grey") +
  geom_hline(yintercept = 0, color = "grey") +
  geom_line(size = 0.75) +
  geom_point(size = 1.5) +
  ylab("pelts (thousands)") 

  
ggplot(data = lynx_hare_df,
       aes(x = Lynx, y = Hare, color = Year)) +
  geom_vline(xintercept = 0, color = "grey") +
  geom_hline(yintercept = 0, color = "grey") +
  geom_path(arrow = arrow(angle = 15, length = unit(0.15, "inches"))) +
  geom_point(size = 1.5) +
  xlab("lynx pelts (thousands)") +
  ylab("hare pelts (thousands)") +
 theme(legend.position="none")

#now we load rstam

library(rstan)
library(gridExtra)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores ())
set.seed(3) # for reproductibility


#contruct and stroe the model in lotka_volterra_model_stan.stan
# file -> new -> stan file 

## Statistical Model: Prior Knowledge and Unexplained Variation

#Solving the inverse problem
#For a given legal value of the model parameters and initial state, the Lotka-Volterra model 
#predicts population dynamics into the future (and into the past). 
#But given noisy data about population dynamics, how do we solve the inverse problem, 
#that of inferring the values of model parameters consistent with the data? 
#The general approach in Bayesian statistics is somewhat counterintuitive, 
#as it involves formulating the forward model then using general principles 
#to solve the inverse problem.

#####Specifically, a Bayesian model requires:
# 1) a mathematical model of what we know about the parameters (i.e., a prior) 
# 2) and a model of what we know about the data generating process given the parameters 
#(i.e., a sampling distribution).

#Mathematically, a prior density p(??) over the sequence of parameters ?? 
#encapsulates our knowledge of the parameters before seeing the data. 
#A sampling distribution (A sampling distribution p(??|y) is called a likelihood
#function when considered as a function L(??) of the parameters ?? for fixed data y
#which may have a continuous, discrete or mixed probability function, p(y|??)
#characterizes the distribution of observable data y given parameters ??.

#Bayes's rule gives us a general solution to the inverse problem, expressing the posterior
#p(??|y) in terms of the prior p(??) and likelihood p(y|??).

#Stan provides a form of Markov chain Monte Carlo (MCMC) sampling that draws a 
#sample ??(1),.,??(M) from the posterior to use for computational inference.

#Measurement Error and Unexplained Variation
#The Lotka-Volterra model is deterministic. 
#Given the system parameters and the initial conditions, the solutions to the 
#equations are fully determined. Population measurements are not so well behaved 
#that they exactly follow the model; there is residual, 
#unexplained variation, as well as measurement error.

#There are factors that impact predator and prey population size other than 
#the current population size. There are variable environmental effects, such as weather, 
#which will vary from season to season and year to year and affect population sizes. 
#Infectious diseases occasionally spread through a population, reducing its size 
#(Hewitt 1921). There are also more long-term environmental factors such as carrying capacity.

##Lotka-Volterra error model
#Solutions to the Lotka-Volterra equations replace the linear predictor xn?? (in a linear regression model),
#but we maintain the error term to compensate for measurement error and unexplained 
#variation in the data. In the case of population dynamics, 
#the data yn consists of measurements of the prey yn,1 and predator yn,2 populations at times tn.1818 This model makes the assumption that the underlying population sizes zn,k and measurements of it yn,k are continuous. This is a very tight approximation to counts when the numbers are in the thousands.

#The true population sizes at time t=0 are unknown-we only have measurements yinit1 and yinit2 of them. The true initial population sizes at time t=0 will be represented by a parameter zinit, so that

#zinit1zinit2==u(t=0)v(t=0).

#Next, let z1,.,zN be the solutions to the Lotka-Volterra differential equations at times t1,.,tN given initial conditions z(t=0)=zinit and parameters ??=(??,??,??,??). Each zn is a pair of prey and predator population sizes at the specified times,

#zn,1zn,2==u(tn)v(tn).

#The zn are random variables, but they are deterministic functions of the random variables for the initial state zinit and system parameters ??,??,??,??.

#The observed data is in the form of measurements yinit of the initial population of prey and predators, and subsequent measurements yn of the populations at times tn, where yinit and the yn consist of a pair of measured population sizes, for the prey and predator species.

#In summary, the measurements, yinit and yn, are drawn indepently 
#from a normal distribution centered at the underlying population sizes, 
#zinit and zn, with noise scales ??. Each quantity here, including the noise scale, 
#is a pair consisting of values for prey and predator.



###################### ok finally lets start

#Fitting the Hudson's Bay Company data
#First, the data is munged into a form suitable for Stan.
# basically as a list

N <- length(lynx_hare_df$Year) - 1
ts <- 1:N
y_init <- c(lynx_hare_df$Hare[1], lynx_hare_df$Lynx[1])
y <- as.matrix(lynx_hare_df[2:(N + 1), 2:3])
y <- cbind(y[ , 2], y[ , 1]); # hare, lynx order
lynx_hare_data <- list(N = N, ts = ts, y_init = y_init, y = y) 
lynx_hare_data

#Next, the model is translated to C++ and compiled. #again, this probably takes some time

model <- stan_model("lotka_volterra_model_stan.stan")

model


#Finally, the compiled model and data are used for sampling.

fit <- sampling(model, data = lynx_hare_data, seed = 123)


#The output can be displayed in tabular form, here limited to the median
#(0.5 quantile) and 80% interval (0.1 and 0.9 quantiles), 
#and restricted to the parameters of interest.

print(fit, pars=c("theta", "sigma", "z_init"),
      probs=c(0.1, 0.5, 0.9), digits = 3)


#The R^ values are all near 1, which is consistent with convergence. The effective sample size estimates for each parameter are sufficient for inference.2828 With effective sample sizes of roughly one thousand, standard errors are roughly one thirtieth the size of posterior standard deviations, being in an inverse square root relation. Thus we have reason to trust that Stan has produced an adequate approximation of the posterior.


#We can furthermore plot the marginal posterior densities and confirm the Markov chains are in agreement with one another.

stan_dens(fit, pars=c("theta", "sigma", "z_init"), separate_chains = TRUE)

#Posterior predictive checks
#We use posterior predictive checks to evaluate how well our model fits the data from which it was estimated.3333 This is "testing on the training data" in machine learning parlance, and while we would not trust it for final evaluation, it is an easy way to spot inconsistencies in the implementation of misspecification in the model.

#The basic idea is to take the posterior for the fitted model and use it to predict what the data should've looked like. That is, we will be replicating new y values that parallel the actual observations y. Becuase they are replicated values, we write them as as yrep. The distribution of these replicated values is given by the posterior predictive distribution,

#p(yrep|y) = ???p(yrep|??) p(??|y) d??,

#where ??=(??,??,??,??,zinit,??) is the vector of parameters for the model. Our two forms of uncertainty are represented in the two terms in the integral. The first is the sampling distribution for the replications, p(yrep|??), which is the distribution of observations yrep given parameters ??. This term encapsulates the unexplained variance and measurement error. The second term is the posterior p(??|y), which encapsulates our uncertainty in our parameter estimates ?? given the observations y. Here, the integral takes a weighted average of the sampling distribution, with weights given by the posterior. In statistical terms, we are calculating an expectation of a function of the parameters, f(??)=p(yrep|??), over the posterior p(??|y), which can be written concisely as a conditional expectation,

#Posterior predictive checks, including posterior means and 50% intervals along with the measured data.  If the model is well calibrated, as this one appears to be, 50% of the points are expected to fall in their 50% intervals.Posterior predictive checks, including posterior means and 50% intervals along with the measured data. If the model is well calibrated, as this one appears to be, 50% of the points are expected to fall in their 50% intervals.
#p(yrep|y) = E[p(yrep|??) ?????? y].


#We sample predictions, Ypred, from p(Ypred???Y) and use these samples to construct a fitted curve for students in bed, together with the uncertainty (90% interval, meaning observed data is expected to fall outside of this interval one in ten times). This posterior predictive check allows us to verify if the model captures the structure of the data. Here we see that the model gives a satisfying fit to the data, and that the model uncertainty is able to capture the variation of the data.


## Plot Posterior predictive checks, including posterior means and 50% intervals along with the measured data. If the model is well calibrated, as this one appears to be, 50% of the points are expected to fall in their 50% intervals.

z_init_draws <- extract(fit)$z_init
z_draws <- extract(fit)$z
y_init_rep_draws <- extract(fit)$y_init_rep
y_rep_draws <- extract(fit)$y_rep
predicted_pelts <- matrix(NA, 21, 2)
min_pelts <- matrix(NA, 21, 2)
max_pelts <- matrix(NA, 21, 2)
for (k in 1:2) {
  predicted_pelts[1, k] <- mean(y_init_rep_draws[ , k])
  min_pelts[1, k] <- quantile(y_init_rep_draws[ , k], 0.25)
  max_pelts[1, k] <- quantile(y_init_rep_draws[ , k], 0.75)
  for (n in 2:21) {
    predicted_pelts[n, k] <- mean(y_rep_draws[ , n - 1, k])
    min_pelts[n, k] <- quantile(y_rep_draws[ , n - 1, k], 0.25)
    max_pelts[n, k] <- quantile(y_rep_draws[ , n - 1, k], 0.75)
  }
}
lynx_hare_melted_df <- melt(as.matrix(lynx_hare_df[, 2:3]))
colnames(lynx_hare_melted_df) <- c("year", "species", "pelts")
lynx_hare_melted_df$year <-
  lynx_hare_melted_df$year +
  rep(1899, length(lynx_hare_melted_df$year))
Nmelt <- dim(lynx_hare_melted_df)[1]
lynx_hare_observe_df <- lynx_hare_melted_df
lynx_hare_observe_df$source <- rep("measurement", Nmelt)
lynx_hare_predict_df <-
  data.frame(year = rep(1900:1920, 2),
             species = c(rep("Lynx", 21), rep("Hare", 21)),
             pelts = c(predicted_pelts[, 2],
                       predicted_pelts[, 1]),
             min_pelts = c(min_pelts[, 2], min_pelts[, 1]),
             max_pelts = c(max_pelts[, 2], max_pelts[, 1]),
             source = rep("prediction", 42))
lynx_hare_observe_df$min_pelts = lynx_hare_predict_df$min_pelts
lynx_hare_observe_df$max_pelts = lynx_hare_predict_df$max_pelts
lynx_hare_observe_predict_df <-
  rbind(lynx_hare_observe_df, lynx_hare_predict_df)
population_plot2 <-
  ggplot(data = lynx_hare_observe_predict_df,
         aes(x = year, y = pelts, color = source)) +
  geom_vline(xintercept = 1900, color = "grey") +
  geom_hline(yintercept = 0, color = "grey") +
  facet_wrap( ~ species, ncol = 1) +
  geom_ribbon(aes(ymin = min_pelts, ymax = max_pelts),
              colour = NA, fill = "black", alpha = 0.1) +
  geom_line() +
  geom_point() +
  ylab("pelts (thousands)") +
  theme(legend.position="bottom")
population_plot2


#############

#Plot of expected population orbit for one hundred draws from the posterior. Each draw represents a different orbit determined by the differential equation system parameters. Together they provide a sketch of posterior uncertainty for the expected population dynamics. If the ODE solutions were extracted per month rather than per year, the resulting plots would appear fairly smooth.

ss <- extract(fit)
df <- data.frame(list(lynx = ss$z[1, 1:12 , 1], hare = ss$z[1, 1:12, 2], draw = 1))
for (m in 2:100) {
  df <- rbind(df, data.frame(list(lynx = ss$z[m, 1:12 , 1], hare = ss$z[m, 1:12, 2], draw = m)))
}
plot <- ggplot(df) +
  geom_vline(xintercept = 0, color = "grey") +
  geom_hline(yintercept = 0, color = "grey") +
  geom_path(aes(x = lynx, y = hare, colour = draw), size = 0.75, alpha = 0.2) +
  #  geom_point(size = 1.25) +
  xlab("lynx pelts (thousands)") +
  ylab("hare pelts (thousands)") +
  theme(legend.position="none")
plot


# Plot of expected pelt collection orbits for one hundred draws of system parameters from the posterior. Even if plotted at more fine-grained time intervals, error would remove any apparent smoothness. Extreme draws as seen here are typical when large values have high error on the multiplicative scale.

ss <- extract(fit)
df <- data.frame(list(lynx = ss$y_rep[1, 1:12 , 1], hare = ss$y_rep[1, 1:12, 2], draw = 1))
for (m in 2:100) {
  df <- rbind(df, data.frame(list(lynx = ss$y_rep[m, 1:12 , 1], hare = ss$y_rep[m, 1:12, 2], draw = m)))
}
plot <- ggplot(df) +
  geom_vline(xintercept = 0, color = "grey") +
  geom_hline(yintercept = 0, color = "grey") +
  geom_path(aes(x = lynx, y = hare, colour = draw), size = 0.75, alpha = 0.2) +
  xlab("lynx pelts (thousands)") +
  ylab("hare pelts (thousands)") +
  theme(legend.position="none")
# The uncertainty due to parameter estimation is rolled into the values of `z_init`, `z`, and `sigma`.  The uncertainty due to unexplained variation and measurement error is captured through the use of the lognormal pseudorandom number generator, `lognormal_rng`.  The additional noise in the measurements `y` over that of the underlying population predictions `z` is visualized in the plots.
plot

