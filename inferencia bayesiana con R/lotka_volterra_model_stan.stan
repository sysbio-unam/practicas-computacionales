
// This is a comment
// this file reprsents the Stan model of the Lotka Volterra ODE together with the 
//
// the ODE model, stored as a function
functions {
  real[] dz_dt(real t,       // time
               real[] z,     // system state {prey, predator}
               real[] theta, // parameters
               real[] x_r,   // unused data
               int[] x_i) {
    real u = z[1];
    real v = z[2];

    real alpha = theta[1];
    real beta = theta[2];
    real gamma = theta[3];
    real delta = theta[4];

    real du_dt = (alpha - beta * v) * u;
    real dv_dt = (-gamma + delta * u) * v;

    return { du_dt, dv_dt };
  }
}
// the data
data {
  int<lower = 0> N;          // number of measurement times
  real ts[N];                // measurement times > 0
  real y_init[2];            // initial measured populations
  real<lower = 0> y[N, 2];   // measured populations
}
// parameters
parameters {
  real<lower = 0> theta[4];   // { alpha, beta, gamma, delta }
  real<lower = 0> z_init[2];  // initial population
  real<lower = 0> sigma[2];   // measurement errors
}
// transformed parameters - the ODE must be integrated
// uses runge kutta
transformed parameters {
  real z[N, 2]
    = integrate_ode_rk45(dz_dt, z_init, 0, ts, theta,
                         rep_array(0.0, 0), rep_array(0, 0),
                         1e-5, 1e-3, 5e2); 
}
// the actual model (bayesian model)
model {
  theta[{1, 3}] ~ normal(1, 0.5); // weakly informative prior distributions of first and 3rd parameter
  theta[{2, 4}] ~ normal(0.05, 0.05); // distr of second and 4th parameter
  sigma ~ lognormal(-1, 1); // prior for the measurement error
  z_init ~ lognormal(log(10), 1); // priori for the initial conditions for both state variables
  // for both state variables, the OBSERVED initial conditions and the dynamical traces are sampled from a 
  // log normal distribution with mean z (which corresponds to the "true" variable generated by the deterministic ode
  // and error sigma
  for (k in 1:2) {
    y_init[k] ~ lognormal(log(z_init[k]), sigma[k]); 
    y[ , k] ~ lognormal(log(z[, k]), sigma[k]);
  }
}
// this piece of code is for posterior predictive checks
// # Stan code for posterior predictive checks
// Stan defines predictive quantities in the generated quantities block, which is executed once per iteration.^[The log density and its gradient are typically evaluated many times per iteration to follow the Hamiltonian trajectory of the parameters given some initial momenta.]  The code declares variables at the top of the block, then defines them in a loop over the species, then over the times.
generated quantities {
  real y_init_rep[2];
  real y_rep[N, 2];
  for (k in 1:2) {
    y_init_rep[k] = lognormal_rng(log(z_init[k]), sigma[k]);
    for (n in 1:N)
      y_rep[n, k] = lognormal_rng(log(z[n, k]), sigma[k]);
  }
}
