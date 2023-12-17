data {
  int<lower = 1> n_spec; // number of random effect levels (species) 
 
  // Traits
  int<lower = 1> N; // Sample size for trait data 
  int<lower = 1, upper = n_spec> trait_species[N]; // id of random effect (species)
  
  vector[N] lati;
  vector[N] tranE; //dummy variable for transect - 0/1
  vector[N] yTraiti; // Observed trait
  
   // Phenology
  int<lower = 1> Nph; // Sample size for phenology data
  int<lower = 1, upper = n_spec> phenology_species[Nph]; // id of random effect (species)
  vector[Nph] yPhenoi; // Outcome phenology
  // vector[Nph] forcei; // predictor forcing
  vector[Nph] chilli; // predictor chilling
  // vector[Nph] photoi; // predictor photoperiod
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
  
  real b_tranlat;
  
  real<lower = 0> sigma_traity; // sd general
  real<lower = 0> sigma_sp; // sd species
  
  //Phenology
 // real alphaForceSp[n_spec];
 // real muForceSp;
 // real<lower = 0> sigmaForceSp;
  real alphaChillSp[n_spec];
  real muChillSp;
  real<lower = 0> sigmaChillSp;
 // real alphaPhotoSp[n_spec];
 // real muPhotoSp;
 // real<lower = 0> sigmaPhotoSp;
  real alphaPhenoSp[n_spec];
  real muPhenoSp;
  real<lower = 0> sigmaPhenoSp;
  // real betaTraitxForce;
  real betaTraitxChill;
  // real betaTraitxPhoto;
  real<lower = 0> sigmapheno_y;
}

transformed parameters{
  // Traits
  vector[N] y_hat; 
  
  //Phenology
  real betaForceSp[n_spec];     //species level beta forcing
  real betaPhotoSp[n_spec];     //species level beta photoperiod
  real betaChillSp[n_spec];     //species level beta chilling

  for (i in 1:N){
    y_hat[i] =   b_muSp[trait_species[i]] + b_tranE * tranE[i]  + b_tranlat * interTranLat[i] 
    ;
  }
  
  // Phenology
  // for (isp in 1:n_spec){
  //   betaForceSp[isp] = alphaForceSp[isp] + betaTraitxForce * (b_muSp[isp]);
  // }
  // for (isp in 1:n_spec){
  //   betaPhotoSp[isp] = alphaPhotoSp[isp] + betaTraitxPhoto * (b_muSp[isp]);
  // }
  for (isp in 1:n_spec){
    betaChillSp[isp] = alphaChillSp[isp] + betaTraitxChill * (b_muSp[isp]);
  }

}

model{
  // Traits
  //// likelihood
  b_muSp ~ normal(muSp, sigma_sp);
  b_tranlat ~ normal(0,1);
  muSp ~ normal(0, 5);// beta(1,1);
  //mu_grand ~ normal(10,5);
  b_tranE ~ normal(0,5);

  //// priors
 // mu_grand ~ normal(10,10);
  sigma_sp ~ normal(0,5);
  sigma_traity ~ normal(0, 5);
  
  yTraiti ~ normal(y_hat, sigma_traity);
 
 // Phenology
  // likelihood
  for (i in 1:Nph){
    yPhenoi[i] ~ normal(alphaPhenoSp[phenology_species[i]] +
   // betaForceSp[phenology_species[i]] * forcei[i] + betaPhotoSp[phenology_species[i]] * photoi[i] 
    + betaChillSp[phenology_species[i]] * chilli[i], sigmapheno_y);
  }
  alphaPhenoSp ~ normal(muPhenoSp, sigmaPhenoSp);
  // alphaForceSp ~ normal(muForceSp, sigmaForceSp);
  alphaChillSp ~ normal(muChillSp, sigmaChillSp);
  // alphaPhotoSp ~ normal(muPhotoSp, sigmaPhotoSp);
  //// priors
  muPhenoSp ~ normal(40,20);
  sigmaPhenoSp ~ normal(5,5);

  sigmapheno_y ~ normal(10,5);
// 
  // muForceSp ~ normal(0,10);
  // sigmaForceSp ~ normal(5,5);

  muChillSp ~ normal(0,10);
  sigmaChillSp ~ normal(5,5);

  // muPhotoSp ~ normal(0,10);
  // sigmaPhotoSp ~ normal(5,5);

  // betaTraitxForce ~ normal(0,1);
  // betaTraitxPhoto ~ normal(0,1);
  betaTraitxChill ~ normal(0,1);


}


generated quantities {
} 
