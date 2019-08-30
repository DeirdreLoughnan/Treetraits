data{
    int<lower=1> N;
    int<lower=1> N_sp_in;
    real bday[N];
    real ht[N];
    real stom_d[N];
    real sla[N];
    real wood_den[N];
    real cn[N];
    int sp_in[N];
}
parameters{
    vector[N_sp_in] a;
    real mu_a;
    real bstom;
    real bsla;
    real bht;
    real bcn;
    real bwood;
    real<lower=0> sigma;
    real<lower=0> sigma_a;
}
model{
    vector[N] mu;
    sigma_a ~ normal( 0 , 10 );
    sigma ~ normal( 0 , 10 );
    bwood ~ normal( 0 , 50 );
    bcn ~ normal( 0 , 50 );
    bht ~ normal( 0 , 50 );
    bsla ~ normal( 0 , 50 );
    bstom ~ normal( 0 , 50 );
    mu_a ~ normal( 0 , 50 );
    a ~ normal( mu_a , sigma_a );
    for ( i in 1:N ) {
        mu[i] = a[sp_in[i]] + bstom * stom_d[i] + bsla * sla[i] + bht * ht[i] + bcn * cn[i] +      bwood * wood_den[i];
    }
    bday ~ normal( mu , sigma );
}
generated quantities{
    vector[N] mu;
    real dev;
    dev = 0;
    for ( i in 1:N ) {
        mu[i] = a[sp_in[i]] + bstom * stom_d[i] + bsla * sla[i] + bht * ht[i] + bcn * cn[i] +      bwood * wood_den[i];
    }
    dev = dev + (-2)*normal_lpdf( bday | mu , sigma );
}
