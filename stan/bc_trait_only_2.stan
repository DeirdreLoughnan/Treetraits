//

data {
  int<lower = 1> n_spec; // number of random effect levels (species) 
  // Traits
  int<lower = 1> N; // Sample size for trait data 
  int<lower = 1, upper = n_spec> species[N]; // id of random effect (species)
  //int<lower = 1> n_study; // number of random effect levels (study) 
  //int<lower = 1, upper = n_study> study[N]; // id of random effect (study)
  vector[N] yTraiti; // Observed trait
  //Priors
  real prior_mu_grand_mu; 
  real prior_mu_grand_sigma;
  real prior_sigma_sp_mu;
  real prior_sigma_sp_sigma;
  // real prior_sigma_study_mu;
  // real prior_sigma_study_sigma;
  real prior_sigma_traity_mu;
  real prior_sigma_traity_sigma;
}

parameters{
  // Traits
  real mu_grand; // grand mean for trait value 
  vector[n_spec] muSp; // species offsets
  // vector[n_study] muStudy; // study offsets
  real<lower = 0> sigma_traity; // sd general
  real<lower = 0> sigma_sp; // sd species
  // real<lower = 0> sigma_study; // sd study
 
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
    y_hat[i] = mu_grand + muSp[species[i]] ;
    //+ muStudy[study[i]];
  }
  
}

model{
  // Traits
  //// likelihood
  yTraiti ~ normal(y_hat, sigma_traity);
  muSp ~ normal(0, sigma_sp);
  // muStudy ~ normal(0, sigma_study);
  //// priors
  mu_grand ~ normal(prior_mu_grand_mu, prior_mu_grand_sigma);
  sigma_sp ~ normal(prior_sigma_sp_mu, prior_sigma_sp_sigma);
  // sigma_study ~ normal(prior_sigma_study_mu, prior_sigma_study_sigma);
  sigma_traity ~ normal(prior_sigma_traity_mu, prior_sigma_traity_sigma);

  
}


generated quantities {
} 
