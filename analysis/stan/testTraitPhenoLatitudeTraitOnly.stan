data {
  int<lower = 1> n_spec; // number of random effect levels (species) 
  int<lower = 1> n_tran;
  // Traits
  int<lower = 1> N; // Sample size for trait data 
  int<lower = 1, upper = n_spec> trait_species[N]; // id of random effect (species)
  int<lower = 1, upper = n_tran> trait_transect[N];
  
  vector[N] lati;
  vector[N] yTraiti; // Observed trait

}

parameters{
  // Traits
  real mu_grand; // grand mean for trait value 
  vector[n_spec] muSp; // species offsets
  
  real mu_tran; // pop offsets
  vector[n_tran] b_tran;
  
  real<lower = 0> sigma_traity; // sd general
  real<lower = 0> sigma_sp; // sd species
  real<lower = 0> sigma_tran; // sd pop

}

transformed parameters{
  // Traits
  vector[N] y_hat; 
  vector[n_spec] mu_grand_sp;
  
  // Traits
  for(i in 1:n_spec){
    mu_grand_sp[i] = mu_grand + muSp[i];
  }
  for (i in 1:N){
    y_hat[i] = mu_grand + b_tran[trait_transect[i]] * lati[i] + muSp[trait_species[i]] ;
  }
  
}

model{
  // Traits
  //// likelihood
  yTraiti ~ normal(y_hat, sigma_traity);
  muSp ~ normal(0, sigma_sp);
  mu_tran ~ normal(0, 35);
  
  //// priors
  mu_grand ~ normal(10,10);
  sigma_sp ~ normal(4,10);
  sigma_tran ~ normal(0,5);
  sigma_traity ~ normal(3, 5);
  
  b_tran ~ normal(mu_tran, sigma_tran);
  

}


generated quantities {
} 
