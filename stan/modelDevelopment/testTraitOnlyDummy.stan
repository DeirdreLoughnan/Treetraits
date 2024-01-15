data {
  int<lower = 1> n_spec; // number of random effect levels (species) 
  // Traits
  int<lower = 1> N; // Sample size for trait data 
  int<lower = 1, upper = n_spec> trait_species[N]; // id of random effect (species)
  vector [N] pop2;
  vector [N] pop3;
  vector [N] pop4;
  vector[N] yTraiti; // Observed trait
  
  // Phenology
  // int<lower = 1> Nph; // Sample size for phenology data
  // int<lower = 1, upper = n_spec> phenology_species[Nph]; // id of random effect (species)
  // vector[Nph] yPhenoi; // Outcome phenology
  // vector[Nph] forcei; // predictor forcing 
  // vector[Nph] chilli; // predictor chilling
  // vector[Nph] photoi; // predictor photoperiod
}

parameters{
  // Traits
  real mu_grand; // grand mean for trait value 
  vector[n_spec] muSp; // species offsets
  //vector[n_pop] mupop; // pop offsets
  real<lower = 0> sigma_traity; // sd general
  real<lower = 0> sigma_sp; // sd species
  //real<lower = 0> sigma_pop; // sd pop
  
  real b_pop2;
  real b_pop3;
  real b_pop4;
  // Phenology
  // real alphaForceSp[n_spec]; 
  // real muForceSp; 
  // real<lower = 0> sigmaForceSp;
  // real alphaChillSp[n_spec]; 
  // real muChillSp; 
  // real<lower = 0> sigmaChillSp;
  // real alphaPhotoSp[n_spec];
  // real muPhotoSp;
  // real<lower = 0> sigmaPhotoSp;
  // real alphaPhenoSp[n_spec];
  // real muPhenoSp;
  // real<lower = 0> sigmaPhenoSp; 
  // real betaTraitxForce; 
  // real betaTraitxChill;
  // real betaTraitxPhoto;
  // real<lower = 0> sigmapheno_y;
}

transformed parameters{
  // Traits
  vector[N] y_hat; 
  vector[n_spec] mu_grand_sp;
  // Phenology
  // real betaForceSp[n_spec];     //species level beta forcing 
  // real betaPhotoSp[n_spec];     //species level beta photoperiod
  // real betaChillSp[n_spec];     //species level beta chilling
  // Traits
  for(i in 1:n_spec){
    mu_grand_sp[i] = mu_grand + muSp[i];
  }
  for (i in 1:N){
    y_hat[i] = mu_grand + b_pop2 * pop2[i] + b_pop3 * pop3[i] + b_pop4 * pop4[i] + muSp[trait_species[i]] ;
  }
  // Phenology
  // for (isp in 1:n_spec){
  //   betaForceSp[isp] = alphaForceSp[isp] + betaTraitxForce * (mu_grand_sp[isp]);
  // }
  // for (isp in 1:n_spec){
  //   betaPhotoSp[isp] = alphaPhotoSp[isp] + betaTraitxPhoto * (mu_grand_sp[isp]);
  // }
  // for (isp in 1:n_spec){
  //   betaChillSp[isp] = alphaChillSp[isp] + betaTraitxChill * (mu_grand_sp[isp]);
  // }
}

model{
  // Traits
  //// likelihood
  yTraiti ~ normal(y_hat, sigma_traity);
  muSp ~ normal(0, sigma_sp);
  //mupop ~ normal(0, sigma_pop);
  b_pop2 ~ normal(0,10);
  b_pop3 ~ normal(0,10);
  b_pop4 ~ normal(0,10);
  //// priors
  mu_grand ~ normal(10,10);
  sigma_sp ~ normal(4,10);
  //sigma_pop ~ normal(2,5);
  sigma_traity ~ normal(3, 5);
  // Phenology
  //// likelihood
  // for (i in 1:Nph){
  //   yPhenoi[i] ~ normal(alphaPhenoSp[phenology_species[i]] + betaForceSp[phenology_species[i]] * forcei[i] + betaPhotoSp[phenology_species[i]] * photoi[i] + betaChillSp[phenology_species[i]] * chilli[i], sigmapheno_y);
  // }
  // alphaPhenoSp ~ normal(muPhenoSp, sigmaPhenoSp);
  // alphaForceSp ~ normal(muForceSp, sigmaForceSp); 
  // alphaChillSp ~ normal(muChillSp, sigmaChillSp);
  // alphaPhotoSp ~ normal(muPhotoSp, sigmaPhotoSp);
  // //// priors
  // muPhenoSp ~ normal(prior_muPhenoSp_mu, prior_muPhenoSp_sigma);
  // sigmaPhenoSp ~ normal(prior_sigmaPhenoSp_mu, prior_sigmaPhenoSp_sigma);
  // 
  // sigmapheno_y ~ normal(prior_sigmaphenoy_mu, prior_sigmaphenoy_sigma);

  // muForceSp ~ normal(prior_muForceSp_mu,  prior_muForceSp_sigma);
  // sigmaForceSp ~ normal(prior_sigmaForceSp_mu, prior_sigmaForceSp_sigma);
  // 
  // muChillSp ~ normal(prior_muChillSp_mu, prior_muChillSp_sigma);
  // sigmaChillSp ~ normal(prior_sigmaChillSp_mu, prior_sigmaChillSp_sigma);
  //  
  // muPhotoSp ~ normal(prior_muPhotoSp_mu, prior_muPhotoSp_sigma);
  // sigmaPhotoSp ~ normal(prior_sigmaPhotoSp_mu, prior_sigmaPhotoSp_sigma);
  //       
  // betaTraitxForce ~ normal(prior_betaTraitxForce_mu, prior_betaTraitxForce_sigma);
  // betaTraitxPhoto ~ normal(prior_betaTraitxPhoto_mu, prior_betaTraitxPhoto_sigma);
  // betaTraitxChill ~ normal(prior_betaTraitxChill_mu, prior_betaTraitxChill_sigma); 

}


generated quantities {
} 
