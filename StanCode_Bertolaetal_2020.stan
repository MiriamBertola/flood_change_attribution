//   Stan code of the flood change model presented in:

//   Bertola, M., Viglione, A., Vorogushyn, S., Lun, D., Merz, B., and Blöschl, G.: Do small and large floods have the same drivers of change? A regional attribution analysis in Europe, Hydrol. Earth Syst. Sci. Discuss., https://doi.org/10.5194/hess-2020-396, in review, 2020.


data {
int<lower=0> styr;         // n. of station years
int<lower=0> K;            // n. of stations
int<lower=0> k[styr];      // site identifier (1:K)
vector[styr] q;            // annual specific flood peaks (m3/s/km2)
vector[styr] lnS;          // log of catchment areas
vector[styr] lnX1;         // covariate 1
vector[styr] lnX2;         // covariate 2
vector[styr] lnX3;         // covariate 3
real<lower=0> lik_corr;    // likelihood correction factor
}
parameters {
  real lnalpha_20;
  real lnalpha_g0;
  real gamma_20;
  real gamma_g0;
  real<lower=0> alpha_21;
  real<lower=-alpha_21*2> alpha_g1;
  real<lower=0> alpha_22;
  real<lower=-alpha_22*2> alpha_g2;
  real<lower=0> alpha_23;
  real<lower=-alpha_23*2> alpha_g3;
  real lnsigma;
  vector[K] eps_lnq2;
}
transformed parameters {
real sigma;
sigma = exp(lnsigma);
}
model {
real a = -ln(-ln(0.5));
real b = -ln(-ln(0.99)) + ln(-ln(0.5));
vector[styr] q2;
vector[styr] x100p;
vector[styr] par2;
vector[styr] z;
// a-priori on parameters
lnalpha_20 ~ normal(0, 3);
lnalpha_g0 ~ normal(0, 3);
gamma_20 ~ normal(0, 2);
gamma_g0 ~ normal(0, 2);
alpha_21 ~ normal(0, 2);
alpha_g1 ~ normal(0, 2);
alpha_22 ~ normal(0, 2);
alpha_g2 ~ normal(0, 2);
alpha_23 ~ normal(0, 2);
alpha_g3 ~ normal(0, 2);
lnsigma ~ normal(0, 3);
eps_lnq2 ~ normal(0, exp(lnsigma));
// likelihood
q2 = exp(lnalpha_20 + gamma_20*lnS + alpha_21*lnX1 + alpha_22*lnX2 + alpha_23*lnX3 + eps_lnq2[k]);
x100p = exp(lnalpha_g0 + gamma_g0*lnS + alpha_g1*lnX1 + alpha_g2*lnX2 + alpha_g3*lnX3);
par2 = q2 .* x100p / b;
z = (q - (q2 .* (1 - x100p * a/b))) ./ par2;
target += sum(-log(par2) - z - exp(-z)) * lik_corr;
}
generated quantities {
real alpha_20;
real alpha_g0;
real a = -ln(-ln(0.5));
real b = -ln(-ln(0.99)) + ln(-ln(0.5));
vector[styr] q2;
vector[styr] x100p;
vector[styr] par2;
vector[styr] z; 
real log_lik[styr];
alpha_20 = exp(lnalpha_20);
alpha_g0 = exp(lnalpha_g0);
q2 = exp(lnalpha_20 + gamma_20*lnS + alpha_21*lnX1 + alpha_22*lnX2 + alpha_23*lnX3 + eps_lnq2[k]);
x100p = exp(lnalpha_g0 + gamma_g0*lnS + alpha_g1*lnX1 + alpha_g2*lnX2 + alpha_g3*lnX3);
par2 = q2 .* x100p / b;
z = (q - (q2 .* (1 - x100p * a/b))) ./ par2; 
for (i in 1:styr) {
log_lik[i] = (-log(par2[i]) - z[i] - exp(-z[i])) * lik_corr;
}
} 