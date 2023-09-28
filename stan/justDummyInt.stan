data {
  int<lower = 1> n_spec; // number of random effect levels (species) 
 
  // Traits
  int<lower = 1> N; // Sample size for trait data 
  int<lower = 1, upper = n_spec> trait_species[N]; // id of random effect (species)
  
  vector[N] lati;
  vector[N] tranE; //dummy variable for transect - 0/1
  vector[N] yTraiti; // Observed trait
  
}

transformed data{
  vector[N] interTranLat;

  interTranLat = tranE .* lati;
}

parameters{
  // Traits
 // real mu_grand; // grand mean for trait value 
  vector[n_spec] b_muSp; // species offsets
  real muSp;
  
  real b_tranE;
  
  //real mu_tranlat;
  real b_tranlat;
  
  //real<lower = 0> sigma_tranlat; // sd pop
  real<lower = 0> sigma_traity; // sd general
  real<lower = 0> sigma_sp; // sd species
  
}

transformed parameters{
  // Traits
  vector[N] y_hat; 
  
  for (i in 1:N){
    y_hat[i] =   b_muSp[trait_species[i]] + b_tranE * tranE[i]  + b_tranlat * interTranLat[i] 
    ;
  }
  
}

model{
  // Traits
  //// likelihood
  b_muSp ~ normal(muSp, sigma_sp);
  b_tranlat ~ normal(0,10);
  muSp ~ normal(10,5);
  //mu_grand ~ normal(10,5);
  b_tranE ~ normal(0,5);
  
 // mu_tranlat ~ normal(0,10);
 // sigma_tranlat ~ normal(0,10);
  
  //// priors
 // mu_grand ~ normal(10,10);
  sigma_sp ~ normal(4,5);
  sigma_traity ~ normal(3, 5);
  
  yTraiti ~ normal(y_hat, sigma_traity);
 
}


generated quantities {
} 
