data {
   int<lower=0> N; // Number of objects measured
   int<lower=0> J; // Number of judges
   int<lower=1> R; // Repetitions of measurements performed by each judge.
   array [N,J,R] real x;
   real l_mu,u_mu;
   real l_sigma_T,u_sigma_T;
   //real l_sigma_E,u_sigma_E;
   real l_sigma_J, u_sigma_J;
   real l_sigma_I, u_sigma_I;
 }
 
 parameters {
   real<lower=l_mu,upper=u_mu> mu;
   real<lower=l_sigma_T,upper=u_sigma_T> sigma_T; // Inherent object variance
   real<lower=l_sigma_J,upper=u_sigma_J> sigma_J; // Inter-observer variance
   //real<lower=l_sigma_E,upper=u_sigma_J> sigma_E;
   array [J] real<lower=l_sigma_I,upper=u_sigma_I> sigma_I; // Intra-observer variances
   array [J] real a; // Bias of the judges
   array [N] real b; // Movement of the object
 }
 
 model {
   mu ~ uniform(l_mu, u_mu);
   sigma_T ~ uniform(l_sigma_T, u_sigma_T);
   sigma_J ~ uniform(l_sigma_J, u_sigma_J);
   
   for (j in 1:J) {
     sigma_I[j] ~ uniform(l_sigma_T, u_sigma_T);
     a[j] ~ normal(0, sigma_J);
   }
   for (i in 1:N) {
     b[i] ~ normal(0, sigma_T);
   }
   for (i in 1:N) {
     for (j in 1:J) {
       for (k in 1:R) {
         x[i,j,k] ~ normal(mu + a[j] + b[i], sigma_I[j]);
         //x[i,j,k] ~ normal(mu + a[j] + b[i], sigma_E);
       }
     }
   }
 }