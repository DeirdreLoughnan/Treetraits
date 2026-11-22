data {
  int<lower = 1> n_spec; // number of random effect levels (species) 
 
  // Traits
  int<lower = 1> N; // Sample size for trait data 
  array [N] int<lower = 1, upper = n_spec> trait_species; // id of random effect (species)
  
  vector[N] lati;
  vector[N] tranE; //dummy variable for transect - 0/1
  vector[N] treeShrub; //dummy variable for transect - 0/1
  vector[N] yTraiti; // Observed trait
  
   // Phenology
  int<lower = 1> Nph; // Sample size for phenology data
  array [Nph] int<lower = 1, upper = n_spec> phenology_species; // id of random effect (species)
  vector[Nph] yPhenoi; // Outcome phenology
  vector[Nph] forcei; // predictor forcing
  vector[Nph] chilli; // predictor chilling
  vector[Nph] photoi; // predictor photoperiod
}

transformed data{
  vector[N] interTranLat;

  interTranLat = tranE .* lati;
  //   
  //  vector[Nph] inter_wp;         // vector[N] inter_wp  = warm .* photo;
  //  vector[Nph] inter_wc;           
  //  vector[Nph] inter_pc;
  //  
  // inter_wp    = forcei .* photoi;
  // inter_wc    = forcei .* chilli;
  // inter_pc    = photoi .* chilli;
}

parameters{
  // Traits
  real mu_grand; // grand mean for trait value 
  vector[n_spec] muSp; // species offsets
  //real muSp;
  
  real b_tranE;  
  real b_treeShrub;  
  real b_lat;
  real b_tranlat;
  
  //real<lower = 0> sigma_tranlat; // sd pop
  real<lower = 0> sigma_traity; // sd general
  real<lower = 0> sigma_sp; // sd species
  
  //Phenology
  array [n_spec] real alphaForceSp;
  real muForceSp;
  real<lower = 0> sigmaForceSp;
  array [n_spec] real alphaChillSp;
  real muChillSp;
  real<lower = 0> sigmaChillSp;
  array [n_spec] real alphaPhotoSp;
  real muPhotoSp;
  real<lower = 0> sigmaPhotoSp;
  array [n_spec] real alphaPhenoSp;
  real muPhenoSp;
  real<lower = 0> sigmaPhenoSp;
  real betaTraitxForce;
  real betaTraitxChill;
  real betaTraitxPhoto;
  real<lower = 0> sigmapheno_y;
  
  // vector[n_spec] b_inter_wp;
  // vector[n_spec] b_inter_wc;
  // vector[n_spec] b_inter_pc;
  // 
  // real mu_b_inter_wp;
  // real mu_b_inter_wc;
  // real mu_b_inter_pc;
  // 
  // real<lower=0> sigma_b_inter_wp;
  // real<lower=0> sigma_b_inter_wc;
  // real<lower=0> sigma_b_inter_pc;
}

transformed parameters{
  // Traits
  vector[N] y_hat; 
  vector[n_spec] mu_grand_sp;
  
  //Phenology
  array [n_spec] real betaForceSp;     //species level beta forcing
  array [n_spec] real betaPhotoSp;     //species level beta photoperiod
  array [n_spec] real betaChillSp;     //species level beta chilling
  
  for(i in 1:n_spec){
    mu_grand_sp[i] = mu_grand + muSp[i];
  }
  
  for (i in 1:N){
    y_hat[i] =  mu_grand + muSp[trait_species[i]] + b_lat * lati[i] + b_tranE * tranE[i]  + b_tranlat * interTranLat[i] + b_treeShrub * treeShrub[i]
    ;
  }
  
  // Phenology
  for (isp in 1:n_spec){
    betaForceSp[isp] = alphaForceSp[isp] + betaTraitxForce * (mu_grand_sp[isp]);
  }
  for (isp in 1:n_spec){
    betaPhotoSp[isp] = alphaPhotoSp[isp] + betaTraitxPhoto * (mu_grand_sp[isp]);
  }
  for (isp in 1:n_spec){
    betaChillSp[isp] = alphaChillSp[isp] + betaTraitxChill * (mu_grand_sp[isp]);
  }
  
}

model{
  // Traits
  muSp ~ normal(0, sigma_sp);
  b_tranlat ~ normal(0,15);
  mu_grand ~ normal(0,30); 
  b_lat ~ normal(0,15);
  b_tranE ~ normal(0,15);
  b_treeShrub ~ normal(0,15);
 // mu_tranlat ~ normal(0,10);
 // sigma_tranlat ~ normal(0,10);
  
  //// priors
 // mu_grand ~ normal(10,10);
  sigma_sp ~ normal(0,15); 
  sigma_traity ~ normal(0,15);
  
  yTraiti ~ normal(y_hat, sigma_traity);
 

 // Phenology
  // likelihood
  for (i in 1:Nph){
    yPhenoi[i] ~ normal(alphaPhenoSp[phenology_species[i]] + 
    betaForceSp[phenology_species[i]] * forcei[i] + 
    betaPhotoSp[phenology_species[i]] * photoi[i] + 
    betaChillSp[phenology_species[i]] * chilli[i]
//    + b_inter_wp[phenology_species[i]] * inter_wp[i] +
// 		b_inter_wc[phenology_species[i]] * inter_wc[i] +
// 		b_inter_pc[phenology_species[i]] * inter_pc[i]
     , sigmapheno_y);
  }
		
  alphaPhenoSp ~ normal(muPhenoSp, sigmaPhenoSp);
  alphaForceSp ~ normal(muForceSp, sigmaForceSp);
  alphaChillSp ~ normal(muChillSp, sigmaChillSp);
  alphaPhotoSp ~ normal(muPhotoSp, sigmaPhotoSp);
  //// priors
  muPhenoSp ~ normal(40,30);
  sigmaPhenoSp ~ normal(5,15);

  sigmapheno_y ~ normal(10,15);

  muForceSp ~ normal(-15,30);
  sigmaForceSp ~ normal(5,15);

  muChillSp ~ normal(-15,30);
  sigmaChillSp ~ normal(5,15);

  muPhotoSp ~ normal(-15,30);
  sigmaPhotoSp ~ normal(5,15);

  betaTraitxForce ~ normal(0,20); //was 100
  betaTraitxPhoto ~ normal(0,20);
  betaTraitxChill ~ normal(0,20);

//   sigma_b_inter_wp ~ normal(0, 10);
// 	sigma_b_inter_wc ~ normal(0, 10);
// 	sigma_b_inter_pc ~ normal(0, 10);
// 	
//   b_inter_wp ~ normal(mu_b_inter_wp, sigma_b_inter_wp); // Delete because all in NCP now
//   b_inter_wc ~ normal(mu_b_inter_wc, sigma_b_inter_wc);		
//  	b_inter_pc ~ normal(mu_b_inter_pc, sigma_b_inter_pc);
//  	
//  	mu_b_inter_wp ~ normal(0, 35);
// 	mu_b_inter_wc ~ normal(0, 35);
// 	mu_b_inter_pc ~ normal(0, 35);
}


