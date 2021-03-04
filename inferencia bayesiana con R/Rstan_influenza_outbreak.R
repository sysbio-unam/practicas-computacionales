###### Learning to use stan to fit a ODE model


# example taken from 
# https://mc-stan.org/users/documentation/case-studies/boarding_school_case_study.html
# SIR model fitted to epidemiological data,


# the first time we do this we have to 
# install package 
#remove.packages("rstan")
#if (file.exists(".RData")) file.remove(".RData")


#pkgbuild::has_build_tools(debug = TRUE)

# make sure we are where we want to be
getwd()
setwd("C:/Users/Elisa/Dropbox/Modelo Fibrilación Auricular Rebeca y  Raquél Pacheco/Código/Stan archivos de prueba")

#n 
#install.packages("rstan", repos = "https://cloud.r-project.org/", dependencies = TRUE)


# load the packages, first simply to get the data
library(outbreaks)
library(tidyverse)

#and have a look of what is there

head(influenza_england_1978_school)


ggplot(data = influenza_england_1978_school) + 
  geom_point(mapping = aes(x = date, y = in_bed)) + 
  labs(y = "Number of students in bed")

library(rstan)
library(gridExtra)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores ())
set.seed(3) # for reproductibility


# create the model in a separate file, stored in the same directory
#File -> New File -> Stan File . Otherwise, open your favorite text editor. Either way, paste in the following and save your work to a file called schools.stan in R's working directory (which can be seen by executing getwd())
#setwd("C:/Users/josed/Downloads")
# sir_negbin.stan


## no we look again at the data,
# time series of cases
cases <- influenza_england_1978_school$in_bed  # Number of students in bed

# total count
N <- 763; ## total number of students

# times
n_days <- length(cases) 
t <- seq(0, n_days, by = 1)
t0 = 0 
t <- t[-1]

#initial conditions
i0 <- 1
s0 <- N - i0
r0 <- 0
y0 = c(S = s0, I = i0, R = r0)

# data for Stan: should be in a list format
data_sir <- list(n_days = n_days, y0 = y0, t0 = t0, ts = t, N = N, cases = cases)


# number of MCMC steps
niter <- 2000

#Next we compile the model, saved in the file sir_negbin.stan,

model <- stan_model("sir_negbin.stan") # this takes quite a while. 
#note, it might be necessary to restart r once after installing.
model


#and run MCMC. For this problem, it suffices to use Stan's defaults. Note that, as is standard practice, we run 4 Markov chains.

fit_sir_negbin <- sampling(model,
                           data = data_sir,
                           iter = niter,
                           chains = 4)

fit_sir_negbin


#Let's specify the parameters of interest.

pars=c('beta', 'gamma', "R0", "recovery_time")

#We start with a summary table of the results, which displays the posterior mean, standard error, quantiles, and some useful diagnostics.

print(fit_sir_negbin, pars = pars)


#Stan gives us a host of information to evaluate whether the inference is reliable. During sampling, warnings can tell us if something is wrong (here we have no warnings6). In the summary table, several quantities are available to check inference. Here we note that R^ is close to 1 (< 1.01), indicating the 4 Markov chains are in close agreement with one another. Furthermore the effective samples size, neff, is large (> 1007), which means the Markov chains were able to cohesively explore the parameter space. Conversely, large R^ and low neff would indicate that the Markov chains are poorly mixing. Apart from fixing coding errors, improving the mixing of the Markov chains almost always requires tweaking the model specification, for example with a reparameterization or stronger priors.

#We can furthermore plot the marginal posterior densities and confirm the Markov chains are in agreement with one another.

stan_dens(fit_sir_negbin, pars = pars, separate_chains = TRUE)

##Checking the model
#Now that we trust our inference, let us check the utility of our model. Utility is problem specific and can include the precise estimation of a quantity or predicting future behaviors. In general, it is good to check if our model, once fitted, produces simulations that are consistent with the observed data. This is the idea behind posterior predictive checks.

#We sample predictions, Ypred, from p(Ypred???Y) and use these samples to construct a fitted curve for students in bed, together with the uncertainty (90% interval, meaning observed data is expected to fall outside of this interval one in ten times). This posterior predictive check allows us to verify if the model captures the structure of the data. Here we see that the model gives a satisfying fit to the data, and that the model uncertainty is able to capture the variation of the data.

smr_pred <- cbind(as.data.frame(summary(
  fit_sir_negbin, pars = "pred_cases", probs = c(0.05, 0.5, 0.95))$summary), t, cases)

colnames(smr_pred) <- make.names(colnames(smr_pred)) # to remove % in the col names

ggplot(smr_pred, mapping = aes(x = t)) +
  geom_ribbon(aes(ymin = X5., ymax = X95.), fill = "orange", alpha = 0.6) +
  geom_line(mapping = aes(x = t, y = X50.)) + 
  geom_point(mapping = aes(y = cases)) +
  labs(x = "Day", y = "Number of students in bed")

#Maybe we also want to access the true number of infected people at each time, and not just the number of students in bed. This is a latent variable for which we have an estimation.

params <- lapply(t, function(i){sprintf("y[%s,2]", i)}) #number of infected for each day
smr_y <- as.data.frame(summary(fit_sir_negbin, 
                               pars = params, probs = c(0.05, 0.5, 0.95))$summary)
colnames(smr_y) <- make.names(colnames(smr_y)) # to remove % in the col names

ggplot(smr_y, mapping = aes(x = t)) +
  geom_ribbon(aes(ymin = X5., ymax = X95.), fill = "orange", alpha = 0.6) +
  geom_line(mapping = aes(x = t, y = X50.)) + 
  labs(x = "Day", y = "Number of infected students")



##Checking our priors
#We can check if our priors are sound by computing the a priori probability of various epidemiological parameters of interest. For instance for influenza, we know from domain knowledge that R0 is typically between 1 and 2, and that the recovery time is approximately 1 week. We want priors that allow for every reasonable configurations of the data but exclude pattently absurd scenarios, per our domain expertise. To check if our priors fulfill this role, we can do a prior predictive check.
#To conduct a prior predictive check, we take the same model as before, 
#put the parameters of interest in the generated_quantities code block, 
#and remove the sampling distribution term from the model. 
#Without the sampling distribution, 
#the parameters are not fitted to the data and are thus 
#sampled from their prior distribution. 
#The Stan code is thus the same as the final Stan code, 
#without the cases ~ neg_binomial_2(col(to_matrix(y), 2), phi); line. 
#A useful trick to make prior predictive check easy is to add a switch 
#compute_likelihood to the data. Then in the model code block :
  
# if (compute_likelihood == 1)
 #   cases ~ neg_binomial_2(col(to_matrix(y), 2), phi);
#This allows to do prior predictive check and inference with the same Stan file.

#We compile the model without the likelihood term


model <- stan_model("sir_negbin_check_priors.stan")

#and sample from it.

fit_sir_prior <- sampling(model,
                          data = data_sir, seed = 0, chains = 4)

#This gives us samples from the a priori distribution of parameters, which we can visualize. Here we show the distribution of the log of the recovery time, with the red bars showing loose bounds on the recovery time (1/2 day and 30 days). We observe that most of the probality mass is between the red bars but we still allow more extreme values, meaning our posterior can concentrate outside the bars, if the data warrants it.

s_prior <- rstan::extract(fit_sir_prior)
df_test <- tibble(r = s_prior$recovery_time)
ggplot(data = df_test) + 
  geom_histogram(mapping = aes(x = r), bins = 30) + 
  labs(x = "log(recovery time)") + 
  geom_vline(xintercept = 0.5, color = "red") + 
  geom_vline(xintercept = 30, color = "red") +
  scale_x_log10()


#We can do the same thing for R0 (again, on the log-scale), the loose bounds being 0.3 and 30.

df_test <- tibble(r = s_prior$R0)
ggplot(data = df_test) + 
  geom_histogram(mapping = aes(x = r), bins = 30) + 
  labs(x = "log(R0)") + 
  geom_vline(xintercept = 0.3, color = "red") + 
  geom_vline(xintercept = 30, color = "red") + 
  scale_x_log10()

#We thus see that these distributions are coherent with domain knowledge. 
# See:
#https://github.com/stan-dev/stan/wiki/Prior-Choice-Recommendations

# for more recommendations on prior choice.

#We can also plot trajectories of infection according to the prior, that is the number of infected people at each time accoring to prior distributions of parameters.

n_draws <- 1000
draws <- as_tibble(t(s_prior$y[,,2][1:n_draws,])) %>% add_column(t=t)
draws <-  pivot_longer(draws, c(1:1000) , names_to = "draw")
draws %>% 
  ggplot() + 
  geom_line(mapping = aes(x = t, y=value, group = draw), alpha = 0.6, size=0.1) +
  geom_hline(yintercept=763, color="red")  +
  geom_text(x=1.8, y=747, label="Population size", color="red") +
  labs(x = "Day", y="Number of infected students")

#And the median (black line) and 90% interval of the a priori number of student in bed (i.e the observed number of infected students).

smr_pred <- cbind(as.data.frame(summary(fit_sir_prior, pars="pred_cases", 
                                        probs=c(0.05, 0.5, 0.95))$summary), t)
colnames(smr_pred) <- make.names(colnames(smr_pred)) # to remove % in the col names

ggplot(smr_pred, mapping=aes(x=t)) +
  geom_ribbon(aes(ymin = X5., ymax = X95.), fill = "orange", alpha = 0.6) +
  geom_line(mapping=aes(x=t, y=X50.)) + 
  geom_hline(yintercept=763, color="red" ) +
  geom_text(x=1.8, y=747, label="Population size", color="red") +
  labs(x = "Day", y="Number of students in bed")


#It seems that most trajectories are reasonable and quite diverse. Still, some of the curves look a little bit funky and suggest we could refine our priors and make them more informative, although it may not be needed here.

#Typically, we can get away with priors that do not capture all our a priori knowledge, provided the data is informative enough. However when dealing with complicated models and relatively sparse data, we usually need well constructed priors to regularize our estimates and avoid non-identifiability.


##Can our inference algorithm recover the right parameters?
#  While there exist many theoretical guarantees for MCMC algorithms, modelers should realize that these rely on a set of assumptions which are not always easy to verify and that many of these guarantees are asymptotic. This means they are proven in the limit where we have an infinite number of samples from the posterior distribution. A very nice, if advanced, review on the subject can be found in Roberts and Rosenthal (2004). As practitioners, we must contend with finite computational resources and assumptions which may or may not hold. The diagnostics we reviewed earlier, e.g. R^, effective sample sizes, provide necessary conditions for the MCMC sampler to work but not sufficient ones. Nevertheless they are potent tools for diagnosing shortcomings in our inference. This section provides further such tools, from both a rigorous and a pragmatic perspective.

#Fitting the model to simulated data is, if done properly, an effective way to test whether our inference algorithm is reliable. If we cannot successfully fit the model in a controlled setting, it is unlikely we can do so with real data. This of course raises the question of what is meant by "successfully fitting" the model. In a Bayesian setting, this means our inference procedure recovers the correct posterior distribution. Unfortunately, even in a controlled setting, the posterior distribution is, but in the simplest cases, not tractable.

#A powerful method to check the accuracy of our Bayesian inference is simulation-based calibration (SBC) (Talts et al. 2018). SBC exploits a very nice consistency result. The intuition is the following: if we draw our parameters from our prior distribution
#??1,...,??2???p(??)
#and for each ??i simulate a data set Yi, we can by fitting the model multiple times recover the prior distribution from the estimated posteriors. This technique is a bit beyond the scope of this tutorial, though we vividly encourage the reader to consult the original paper, or to see here how this method fits in a principled Bayesian workflow.

##### For the time being, we focus on a simpler heuristic: fit the model to one simulated data set and check if we recover the correct parameter value. There are serious limitations with this approach: when do we consider that the estimated posterior distribution covers the correct value? How do we know if the variance of the posterior is properly estimated? etc. But the test is useful: in this controlled setting, do the chains converge? Is the computation of the log density numerically stable (e.g. are we using the right ODE integrator)? Do my priors prevent the chains from wondering into absurd regions of the parameter space? These are all questions this simple test can help us tackle.

#We take one arbitrary draw from the prior distribution

# one arbitrary draw from the prior distribution
draw <- 12 
# the number of predicted cases sampled from the prior distribution, which we will use as data
cases_simu <- s_prior$pred_cases[draw,] 

#And use it as data which we fit to our model.

data_simu <-  list (n_days  = n_days, y0 = y0, t0 = t0, ts = t, N=N, cases=cases_simu)

fit_simu <- sampling(model, data=data_simu, chains=4)

#####################################

#We can then examine the estimated posterior distribution.
params = c("beta", "gamma", "phi")
paste("true beta :", toString(s_prior$beta[draw]), 
      ", true gamma :", toString(s_prior$gamma[draw]), ", true phi :", toString(s_prior$phi[draw]))
## [1] "true beta : 2.53115910937212 , true gamma : 0.505381459021197 , true phi : 28.874355077517"
print(fit_simu, pars = params)

#We plot the posterior density (in red) to check if it matches the true value of the parameter (black line). The density is compatible with the true parameters, although not always centered on it. The latter is not alarming, especially if the model parameter ?? we sampled, lies on the tail of the prior distribution. We could repeat this process a few times to get a better sense of the performance of our inference algorithm.

plot_beta <- stan_dens(fit_simu, pars="beta") + geom_vline(xintercept =s_prior$beta[draw])
plot_gamma <- stan_dens(fit_simu, pars="gamma") + geom_vline(xintercept = s_prior$gamma[draw])
plot_phi <- stan_dens(fit_simu, pars="phi_inv") + geom_vline(xintercept = s_prior$phi_inv[draw])
grid.arrange(plot_beta, plot_gamma, plot_phi, nrow=1)