// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include <RcppArmadillo.h>
#include <limits>
#include <omp.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]


#define minimum(x,y) (((x) < (y)) ? (x) : (y))

using namespace std;
using namespace Rcpp;
using namespace arma;



// ============================================================================-
// compute loglikelihoods ----
// ============================================================================-

// ----------------------------------------------------------------------------=
// Whale ----
// ----------------------------------------------------------------------------=

// compute true whale surface
// [[Rcpp::export]]
mat cppf_loglam_b0t_psi(List Xwhale, mat logZ, vec beta, vec psi, vec indtwhale){
  int n = indtwhale.size(), nG = logZ.n_rows, n_be = beta.size();
  mat res = zeros(nG, n), Xwhalei;
  
  for(int i = 0; i < n; i ++){
    Xwhalei = as<mat>(Xwhale[i]);
    res.col(i) = beta[i] + Xwhalei * beta.rows(n, n_be-2) + beta[n_be-1] * logZ.col(indtwhale[i]) + psi;
  }
  
  return res;
}


// [[Rcpp::export]]
mat cppf_loglam_cZ_b0t_psi(List Xwhale, mat logZ, vec beta, vec psi, vec indtwhale){
  int n = indtwhale.size(), nG = logZ.n_rows, n_be = beta.size();
  mat res = zeros(nG, n), Xwhalei;
  double meanlogZ = mean(mean(logZ));
  
  for(int i = 0; i < n; i ++){
    Xwhalei = as<mat>(Xwhale[i]);
    res.col(i) = beta[i] + Xwhalei * beta.rows(n, n_be-2) + beta[n_be-1] * (logZ.col(indtwhale[i]) - meanlogZ) + psi;
  }
  
  return res;
}


// [[Rcpp::export]]
mat cppf_loglam_b0t_psit(List Xwhale, mat logZ, vec beta, mat psi, vec indtwhale){
  int n = indtwhale.size(), nG = logZ.n_rows, n_be = beta.size();
  mat res = zeros(nG, n), Xwhalei;
  
  for(int i = 0; i < n; i ++){
    Xwhalei = as<mat>(Xwhale[i]);
    res.col(i) = beta[i] + Xwhalei * beta.rows(n, n_be-2) + beta[n_be-1] * logZ.col(indtwhale[i]) + psi.col(i);
  }
  
  return res;
}


// [[Rcpp::export]]
mat cppf_loglam_cZ_b0t_psit(List Xwhale, mat logZ, vec beta, mat psi, vec indtwhale){
  int n = indtwhale.size(), nG = logZ.n_rows, n_be = beta.size();
  mat res = zeros(nG, n), Xwhalei;
  double meanlogZ = mean(mean(logZ));
  
  for(int i = 0; i < n; i ++){
    Xwhalei = as<mat>(Xwhale[i]);
    res.col(i) = beta[i] + Xwhalei * beta.rows(n, n_be-2) + beta[n_be-1] * (logZ.col(indtwhale[i]) - meanlogZ) + psi.col(i);
  }
  
  return res;
}


// [[Rcpp::export]]
mat cppf_loglam_b0t_psit_woZ(List Xwhale, vec beta, mat psi){
  int n = psi.n_cols, nG = psi.n_rows, n_be = beta.size();
  mat res = zeros(nG, n), Xwhalei;
  
  for(int i = 0; i < n; i ++){
    Xwhalei = as<mat>(Xwhale[i]);
    res.col(i) = beta[i] + Xwhalei * beta.rows(n, n_be-1) + psi.col(i);
  }
  
  return res;
}


// [[Rcpp::export]]
mat cppf_loglam_b0t_psi_woZ(List Xwhale, vec beta, vec psi){
  int n = Xwhale.size(), nG = psi.size(), n_be = beta.size();
  mat res = zeros(nG, n), Xwhalei;
  
  for(int i = 0; i < n; i ++){
    Xwhalei = as<mat>(Xwhale[i]);
    res.col(i) = beta[i] + Xwhalei * beta.rows(n, n_be-1) + psi;
  }
  
  return res;
}


// distance sampling intensity
// [[Rcpp::export]]
double cppf_sumloglam_dist_obs(mat loglam, vec pii, vec f_dist_obs, mat inddist, vec indmdist){
  int n = inddist.n_rows;
  vec res = zeros(n);
  
  for(int i = 0; i < n; i ++){
    res[i] = log( pii[ indmdist[i] ] ) + log( f_dist_obs[i] ) + loglam(inddist(i,0), inddist(i,1));
  }
  
  return sum(res);
}




// integral of distance sampling intensity
// [[Rcpp::export]]
double cppf_sumintlam_dist(mat lam, vec pii, List f_dist, vec indmwhale, double cellarea){
  int n = indmwhale.size(), L;
  double res = 0;
  mat f_dist_t;
  
  for(int t = 0; t < n; t ++){
    f_dist_t = as<mat>(f_dist[t]);
    L = f_dist_t.n_cols;
    for(int l = 0; l < L; l ++){
      res = res + pii[ indmwhale[t] ] * sum( f_dist_t.col(l) % lam.col(t) ) * cellarea;
    }
  }
  return res;
}


// integral of distance sampling intensity
// [[Rcpp::export]]
vec cppf_intlam_dist(mat lam, vec pii, List f_dist, vec indmwhale, double cellarea){
  int n = indmwhale.size(), L;
  vec res = zeros(n);
  mat f_dist_t;
  
  for(int t = 0; t < n; t ++){
    f_dist_t = as<mat>(f_dist[t]);
    L = f_dist_t.n_cols;
    for(int l = 0; l < L; l ++){
      res[t] = res[t] + pii[ indmwhale[t] ] * sum( f_dist_t.col(l) % lam.col(t) ) * cellarea;
    }
  }
  return res;
}


// loglikelihood for distance sampling
// [[Rcpp::export]]
double cppf_loglik_dist(mat loglam, mat lam, vec pii, vec f_dist_obs, List f_dist, mat inddist, vec indmdist, vec indmwhale, double cellarea){
  double sumint = cppf_sumintlam_dist(lam, pii, f_dist, indmwhale, cellarea);
  double sumloglam = cppf_sumloglam_dist_obs(loglam, pii, f_dist_obs, inddist, indmdist);
  return - sumint + sumloglam;
}



// integral of distance sampling intensity
// [[Rcpp::export]]
vec cppf_sumintlam_dist_t(mat lam, vec pii, List f_dist, vec indmwhale, double cellarea){
  int n = indmwhale.size(), L;
  vec res = zeros(n);
  mat f_dist_t;
  
  for(int t = 0; t < n; t ++){
    f_dist_t = as<mat>(f_dist[t]);
    L = f_dist_t.n_cols;
    for(int l = 0; l < L; l ++){
      res[t] = pii[ indmwhale[t] ] * sum( f_dist_t.col(l) % lam.col(t) ) * cellarea;
    }
  }
  return res;
}


// distance sampling intensity
// [[Rcpp::export]]
vec cppf_sumloglam_dist_obs_t(mat loglam, vec pii, vec f_dist_obs, mat inddist, vec indmdist){
  int n = inddist.n_rows, nt = loglam.n_cols;
  vec resi = zeros(n), res = zeros(nt);
  uvec ind;
  
  for(int i = 0; i < n; i ++){
    resi[i] = log( pii[ indmdist[i] ] ) + log( f_dist_obs[i] ) + loglam(inddist(i,0), inddist(i,1));
  }
  
  for(int t = 0; t < nt; t ++){
    ind = find(indmdist == t);
    res[t] = sum(sum(resi.rows(ind)));
  }
  
  return res;
}


// loglikelihood for distance sampling for each day
// [[Rcpp::export]]
vec cppf_loglik_dist_t(mat loglam, mat lam, vec pii, vec f_dist_obs, List f_dist, mat inddist, vec indmdist, vec indmwhale, double cellarea){
  vec sumint = cppf_sumintlam_dist_t(lam, pii, f_dist, indmwhale, cellarea);
  vec sumloglam = cppf_sumloglam_dist_obs_t(loglam, pii, f_dist_obs, inddist, indmdist);
  return - sumint + sumloglam;
}




// integral of pam intensity
// [[Rcpp::export]]
List cppf_lam_pam(mat lam, vec dur_dist, double c, List p_pam, double cellarea){
  mat p_pam_i, dummy;
  int n = lam.n_cols, K;
  List res(n);
  
  for(int i = 0; i < n; i ++){
    p_pam_i = as<mat>(p_pam[i]);
    K = p_pam_i.n_cols;
    dummy = zeros(K);
    for(int k = 0; k < K; k ++){
      dummy[k] = dur_dist[i] * c * sum( p_pam_i.col(k) % lam.col(i) ) * cellarea;
    }
    res[i] = dummy;
  }
  return res;
}



// integral of pam intensity
// [[Rcpp::export]]
List cppf_lam_pam_ct(mat lam, vec dur_dist, vec c, List p_pam, double cellarea){
  mat p_pam_i, dummy;
  int n = lam.n_cols, K;
  List res(n);
  
  for(int i = 0; i < n; i ++){
    p_pam_i = as<mat>(p_pam[i]);
    K = p_pam_i.n_cols;
    dummy = zeros(K);
    for(int k = 0; k < K; k ++){
      dummy[k] = dur_dist[i] * c[i] * sum( p_pam_i.col(k) % lam.col(i) ) * cellarea;
    }
    res[i] = dummy;
  }
  return res;
}



// loglikelihood for pam
// [[Rcpp::export]]
double cppf_loglik_pam(List Wpam, mat lam, vec dur_dist, double c, List p_pam, double cellarea){
  int n = lam.n_cols;
  vec dummy;
  List lam_pam = cppf_lam_pam(lam, dur_dist, c, p_pam, cellarea);
  double res = 0;
  
  for(int i = 0; i < n; i ++){
    dummy = as<vec>(lam_pam[i]);
    res = res - sum( dummy ) + sum( as<vec>(Wpam[i]) % log(dummy) );
  }
  
  return res;
} 




// loglikelihood for pam
// [[Rcpp::export]]
double cppf_loglik_pam_ct(List Wpam, mat lam, vec dur_dist, vec c, List p_pam, double cellarea){
  int n = lam.n_cols;
  vec dummy;
  List lam_pam = cppf_lam_pam_ct(lam, dur_dist, c, p_pam, cellarea);
  double res = 0;
  
  for(int i = 0; i < n; i ++){
    dummy = as<vec>(lam_pam[i]);
    res = res - sum( dummy ) + sum( as<vec>(Wpam[i]) % log(dummy) );
  }
  
  return res;
} 


// loglikelihood for pam
// [[Rcpp::export]]
vec cppf_loglik_pam_ct_t(List Wpam, mat lam, vec dur_dist, vec c, List p_pam, double cellarea){
  int n = lam.n_cols;
  vec dummy;
  List lam_pam = cppf_lam_pam_ct(lam, dur_dist, c, p_pam, cellarea);
  vec res = zeros(n);
  
  for(int i = 0; i < n; i ++){
    dummy = as<vec>(lam_pam[i]);
    res[i] = - sum( dummy ) + sum( as<vec>(Wpam[i]) % log(dummy) );
  }
  
  return res;
} 



// ----------------------------------------------------------------------------=
// Zoop ----
// ----------------------------------------------------------------------------=

// compute true Zoop surface
// [[Rcpp::export]]
mat cppf_logZ_eta(List Xzoop, vec alpha, vec eta){
  int n = Xzoop.size(), nG = eta.n_rows;
  mat res = zeros(nG, n);
  mat Xzoopi;
  
  for(int i = 0; i < n; i ++){
    Xzoopi = as<mat>(Xzoop[i]);
    res.col(i) = Xzoopi * alpha + eta;
  }
  
  return res;
}


// [[Rcpp::export]]
mat cppf_logZ_a0t_etat(List Xzoop, vec alpha0, vec alpha, mat eta){
  int n = Xzoop.size(), nG = eta.n_rows;
  mat res = zeros(nG, n);
  mat Xzoopi;
  
  for(int i = 0; i < n; i ++){
    Xzoopi = as<mat>(Xzoop[i]);
    res.col(i) = alpha0[i] + Xzoopi * alpha + eta.col(i);
  }
  
  return res;
}

// ja: joint updating for alpha0 and alpha
// [[Rcpp::export]]
mat cppf_logZ_a0t_etat_ja(List Xzoop, vec alpha, mat eta){
  int n = Xzoop.size(), nG = eta.n_rows, n_al = alpha.size();
  mat res = zeros(nG, n);
  mat Xzoopi;
  
  for(int i = 0; i < n; i ++){
    Xzoopi = as<mat>(Xzoop[i]);
    res.col(i) = alpha[i] + Xzoopi * alpha.rows(n, n_al-1) + eta.col(i);
  }
  
  return res;
}


// [[Rcpp::export]]
mat cppf_logZ_a0t_eta_ja(List Xzoop, vec alpha, vec eta){
  int n = Xzoop.size(), nG = eta.size(), n_al = alpha.size();
  mat res = zeros(nG, n);
  mat Xzoopi;
  
  for(int i = 0; i < n; i ++){
    Xzoopi = as<mat>(Xzoop[i]);
    res.col(i) = alpha[i] + Xzoopi * alpha.rows(n, n_al-1) + eta;
  }
  
  return res;
}





// [[Rcpp::export]]
mat cppf_logZ_a0t_etat_pref(List Xzoop, vec alpha0, vec alpha, mat eta, vec delta){
  int n = Xzoop.size(), nG = eta.n_rows, n_al = alpha.size();
  mat res = zeros(nG, n);
  mat Xzoopi;
  
  for(int i = 0; i < n; i ++){
    Xzoopi = as<mat>(Xzoop[i]);
    res.col(i) = alpha0[i] + Xzoopi * alpha.rows(0, n_al-2) + eta.col(i) + alpha[n_al-1] * delta;
  }
  
  return res;
}


// [[Rcpp::export]]
mat cppf_logZ_a0t_etat_ja_pref(List Xzoop, vec alpha, mat eta, vec delta){
  int n = Xzoop.size(), nG = eta.n_rows, n_al = alpha.size();
  mat res = zeros(nG, n);
  mat Xzoopi;
  
  for(int i = 0; i < n; i ++){
    Xzoopi = as<mat>(Xzoop[i]);
    res.col(i) = alpha[i] + Xzoopi * alpha.rows(n, n_al-2) + eta.col(i) + alpha[n_al-1] * delta;
  }
  
  return res;
}


// true Zoop at observed oblique tow sites
// [[Rcpp::export]]
vec cppf_logZobl(mat logZ, mat indobl){
  int n = indobl.n_rows;
  vec res = zeros(n);
  
  for(int i = 0; i < n; i ++){
    res[i] = logZ(indobl(i,0), indobl(i,1));
  }
  return res; 
}


// loglikelihood for oblique tow
// [[Rcpp::export]]
double cppf_loglik_obl(vec logYobl, vec logZobl, double sig2_obl){
  return - 0.5 * sum( (logYobl - logZobl) % (logYobl - logZobl) ) / sig2_obl;
}


// loglikelihood for oblique tow for each day
// [[Rcpp::export]]
vec cppf_loglik_obl_t(vec logYobl, vec logZobl, double sig2_obl, vec indtobl, int nt){
  vec res = zeros(nt);
  uvec ind;
  
  for(int t = 0; t < nt; t ++){
    ind = find(indtobl == t);
    res[t] = - 0.5 * sum(sum( (logYobl.elem(ind) - logZobl.elem(ind)) % (logYobl.elem(ind) - logZobl.elem(ind)) )) / sig2_obl;
  }
  
  return res;
}


// true Zoop at observed surface tow sites
// [[Rcpp::export]]
vec cppf_logZsur(mat logZ, mat indsur){
  int n = indsur.n_rows;
  vec res = zeros(n);
  
  for(int i = 0; i < n; i ++){
    res[i] = logZ(indsur(i,0), indsur(i,1));
  }
  return res; 
}


// loglikelihood for surface tow
// [[Rcpp::export]]
double cppf_loglik_sur(vec logYsur, vec logZsur, double lam0, double lam1, double sig2_sur){
  vec dummy = logYsur - lam0  - lam1 * logZsur;
  return - 0.5 * sum( dummy % dummy ) / sig2_sur;
}



// loglikelihood for surface tow for each day
// [[Rcpp::export]]
vec cppf_loglik_sur_t(vec logYsur, vec logZsur, double lam0, double lam1, double sig2_sur, vec indtsur, int nt){
  vec res = zeros(nt), dummy = logYsur - lam0  - lam1 * logZsur;
  uvec ind;
  
  for(int t = 0; t < nt; t ++){
    ind = find(indtsur == t);
    res[t] = - 0.5 * sum(sum( dummy.elem(ind) % dummy.elem(ind) )) / sig2_sur;
  }
  
  return res;
}



// log intensity for preferential locations
// [[Rcpp::export]]
vec cppf_loglam_pref(double mu0, vec delta){
  return mu0 + delta;
}


// sum of log intensity at preferential locations
// [[Rcpp::export]]
double cppf_sumloglam_pref_obs(vec loglam_pref, vec indpref){
  int n = indpref.size();
  vec res = zeros(n);
  
  for(int i = 0; i < n; i ++){
    res[i] = loglam_pref[indpref[i]];
  }
  
  return sum(res);
}


// loglikelihood for preferential sampling
// [[Rcpp::export]]
double cppf_loglik_pref(vec loglam_pref, vec lam_pref, vec indpref, double cellarea){
  double integ = sum(lam_pref * cellarea);
  double sumloglam = cppf_sumloglam_pref_obs(loglam_pref, indpref);
  return - integ + sumloglam;
}




// ============================================================================-
// Gibbs update ----
// ============================================================================-

// Gibbs update for kappa_eta
// [[Rcpp::export]]
double cppf_gibbs_kappa_etat(mat eta, mat invReta, double a_kappa_eta, double b_kappa_eta){
  int n = eta.n_cols, nG = eta.n_rows;
  vec res = zeros(n);
  double shape = a_kappa_eta + 0.5 * nG * n;
  vec dummy_scale = zeros(n);
  
  for(int i = 0; i < n; i ++){
    dummy_scale[i] = sum( trans(eta.col(i)) * invReta * eta.col(i) );
  }
  double scale = b_kappa_eta + 0.5 * sum(dummy_scale);
  
  return 1/randg(distr_param(shape, 1/scale));
}


// [[Rcpp::export]]
double cppf_gibbs_kappa_eta(vec eta, mat invReta, double a_kappa_eta, double b_kappa_eta){
  double shape = a_kappa_eta + 0.5 * eta.n_rows;
  double scale = ( b_kappa_eta + 0.5 * trans(eta) * invReta * eta )[0];
  return 1/randg(distr_param(shape, 1/scale));
}

// Gibbs update for kappa_delta
// [[Rcpp::export]]
double cppf_gibbs_kappa_delta(vec delta, mat invRdelta, double a_kappa_delta, double b_kappa_delta){
  double shape = a_kappa_delta + 0.5 * delta.n_rows;
  double scale = ( b_kappa_delta + 0.5 * trans(delta) * invRdelta * delta )[0];
  return 1/randg(distr_param(shape, 1/scale));
}


// Gibbs update for sig2_obl
// [[Rcpp::export]]
double cppf_gibbs_sig2_obl(vec logYobl, vec logZobl, double a_sig2_obl, double b_sig2_obl){
  double shape = a_sig2_obl + 0.5 * logYobl.n_rows;
  double scale = b_sig2_obl + 0.5 * sum( (logYobl - logZobl) % (logYobl - logZobl) );
  return 1/randg(distr_param(shape, 1/scale));
}


// Gibbs update for alpha0tilde
// [[Rcpp::export]]
double cppf_gibbs_alpha0tilde(vec alpha0, double tau2){
  double var_alpha0tilde = 1 / ( alpha0.size() / tau2 + 1/100 );
  double mean_alpha0tilde = var_alpha0tilde / tau2 * sum(alpha0);
  return mean_alpha0tilde + sqrt(var_alpha0tilde) * randn();
}


// Gibbs update for lam0
// [[Rcpp::export]]
double cppf_gibbs_lam0(vec logYsur, vec logZsur, double lam1, double sig2_sur, double sig2_lam0){
  double var_lam0 = 1 / ( 1 / sig2_lam0 + logYsur.n_rows / sig2_sur );
  double mean_lam0 = sum( logYsur - lam1 * logZsur ) * var_lam0 / sig2_sur;
  return mean_lam0 + sqrt(var_lam0) * randn();
}


// Gibbs update for lam1
// [[Rcpp::export]]
double cppf_gibbs_lam1(vec logYsur, vec logZsur, double lam0, double sig2_sur, double sig2_lam1){
  double var_lam1 = 1 / ( 1 / sig2_lam1 + sum( logZsur % logZsur ) / sig2_sur );
  double mean_lam1 = sum( (logYsur - lam0) % logZsur ) * var_lam1 / sig2_sur;
  return mean_lam1 + sqrt(var_lam1) * randn();
}


// Gibbs update for sig2_sur
// [[Rcpp::export]]
double cppf_gibbs_sig2_sur(vec logYsur, vec logZsur, double lam0, double lam1, double a_sig2_sur, double b_sig2_sur){
  double shape = a_sig2_sur + 0.5 * logYsur.n_rows;
  vec dummy = logYsur - lam0 - lam1 * logZsur;
  double scale = b_sig2_sur + 0.5 * sum( dummy % dummy );
  return 1/randg(distr_param(shape, 1/scale));
}


// Gibbs update for kappa_psi
// [[Rcpp::export]]
double cppf_gibbs_kappa_psi(vec psi, mat invRpsi, double a_kappa_psi, double b_kappa_psi){
  int nG = psi.size();
  double shape = a_kappa_psi + 0.5 * nG;
  double scale = ( b_kappa_psi + 0.5 * trans(psi) * invRpsi * psi )[0];
  return 1/randg(distr_param(shape, 1/scale));
}


// [[Rcpp::export]]
double cppf_gibbs_kappa_psit(mat psi, mat invRpsi, double a_kappa_psi, double b_kappa_psi){
  int n = psi.n_cols, nG = psi.n_rows;
  vec res = zeros(n);
  double shape = a_kappa_psi + 0.5 * nG * n;
  vec dummy_scale = zeros(n);
  
  for(int i = 0; i < n; i ++){
    dummy_scale[i] = sum( trans(psi.col(i)) * invRpsi * psi.col(i) );
  }
  double scale = b_kappa_psi + 0.5 * sum(dummy_scale);
  
  return 1/randg(distr_param(shape, 1/scale));
}


// [[Rcpp::export]]
double cppf_gibbs_kappa_psis(mat psi, mat invRpsi, double a_kappa_psi, double b_kappa_psi){
  int n = psi.n_cols, nG = psi.n_rows;
  vec res = zeros(n);
  double shape = a_kappa_psi + 0.5 * nG * n;
  vec dummy_scale = zeros(n);
  
  for(int i = 0; i < n; i ++){
    dummy_scale[i] = sum( trans(psi.col(i)) * invRpsi * psi.col(i) );
  }
  double scale = b_kappa_psi + 0.5 * sum(dummy_scale);
  
  return 1/randg(distr_param(shape, 1/scale));
}



// constant call rate over time
// [[Rcpp::export]]
double cppf_gibbs(List Wpam, mat lam, vec dur_dist, List p_pam, double cellarea, double alpha_c, double beta_c){
  int nt = lam.n_cols, K, nind;
  mat p_pam_t;
  vec Wpam_t;
  uvec ind;
  double shape = alpha_c;
  double rate = 1 / beta_c;
  
  for(int t = 0; t < nt; t ++){
    
    Wpam_t = as<vec>(Wpam[t]);
    p_pam_t = as<mat>(p_pam[t]);
    K = p_pam_t.n_cols;
    
    for(int k = 0; k < K; k ++){
      shape = shape + Wpam_t[k];
      rate = rate + dur_dist[t] * sum( p_pam_t.col(k) % lam.col(t) ) * cellarea;
    }
  }
  
  return randg(distr_param(shape, 1/rate));
}


// daily call rates
// [[Rcpp::export]]
vec cppf_gibbs_ct(List Wpam, mat lam, vec dur_dist, List p_pam, double cellarea, double alpha_c, double beta_c){
  int nt = lam.n_cols, K, nind;
  mat p_pam_t;
  vec Wpam_t;
  uvec ind;
  vec res = zeros(nt);
  double shape, rate;
  
  for(int t = 0; t < nt; t ++){
    shape = alpha_c;
    rate = 1 / beta_c;
    
    Wpam_t = as<vec>(Wpam[t]);
    p_pam_t = as<mat>(p_pam[t]);
    K = p_pam_t.n_cols;
    
    for(int k = 0; k < K; k ++){
      shape = shape + Wpam_t[k];
      rate = rate + dur_dist[t] * sum( p_pam_t.col(k) % lam.col(t) ) * cellarea;
    }
    res[t] = randg(distr_param(shape, 1/rate));
  }
  
  return res;
}


// ============================================================================-
// post inference ----
// ============================================================================-

// dummy = cppf_logZ_a0t_etat_ja(Xzoop, postAlpha[iter,], sapply(1:length(postEta), function(t) postEta[[t]][iter,]))

// // ja: joint updating for alpha0 and alpha
// // [[Rcpp::export]]
// mat cppf_logZ_a0t_etat_ja(List Xzoop, vec alpha, mat eta){
//   int n = Xzoop.size(), nG = eta.n_rows, n_al = alpha.size();
//   mat res = zeros(nG, n);
//   mat Xzoopi;
//   
//   for(int i = 0; i < n; i ++){
//     Xzoopi = as<mat>(Xzoop[i]);
//     res.col(i) = alpha[i] + Xzoopi * alpha.rows(n, n_al-1) + eta.col(i);
//   }
//   
//   return res;
// }


// [[Rcpp::export]]
List cppf_logZ_a0t_etat_ja_list(List Xzoop, mat postAlpha, List postEta){
  int nzoop = Xzoop.size(), niters = postAlpha.n_rows, ngrid, n_al = postAlpha.n_cols;
  List res(nzoop);
  mat Xzoopt, postEtat, dummy;
  vec alpha;
  
  for(int t = 0; t < nzoop; t ++){
    Xzoopt = as<mat>(Xzoop[t]);
    postEtat = as<mat>(postEta[t]);
    
    ngrid = Xzoopt.n_rows;
    dummy = zeros(niters, ngrid);
    
    for(int iter = 0; iter < niters; iter ++){
      alpha = trans( postAlpha.row(iter) );
      dummy.row(iter) = trans( alpha[t] + Xzoopt * alpha.rows(nzoop, n_al-1) + trans(postEtat.row(iter)) );
    }
    res[t] = dummy;
  }
  
  return res;
}


// ============================================================================-
// Priors ----
// ============================================================================-

// Log of unnormalized MVN pdf
// [[Rcpp::export]]
double MVN_logh(vec y, vec mu, mat invSigma){
  vec result = - 0.5 * trans(y - mu) * invSigma * (y - mu);
  return result[0];
}


// Log of unnormalized normal
// [[Rcpp::export]]
double normal_logh(double y, double mu, double sig2){
  return - 0.5 / sig2 * (y - mu) * (y - mu);
}



// ============================================================================-
// MCMC ----
// ============================================================================-

// // a0t: daily zoop intercepts
// // etat: daily zoop GPs
// // b0t: daily whale intercepts
// // psit: daily whale GPs
// // ct: daily calling rates
// // ja: joint updating for alpha_0 and alpha
// // fixC: fix call rates to truth
// 
// 
// // [[Rcpp::export]]
// List cppf_fitJoint_cZ_a0t_etat_b0t_psi_ja(
//     int niter, mat distmat,
//     List Xzoop, vec logYobl, vec logYsur,
//     List Xwhale, mat Sdist, List Wpam, vec dur_dist,
//     vec f_dist_obs, List f_dist, List p_pam,
//     mat indobl, mat indsur, vec indtwhale, mat inddist, vec indmdist, vec indmwhale,
//     double alpha0tilde, vec alpha, mat eta, double lam0, double lam1,
//     double tau2, double kappa_eta, double phi_eta, double sig2_obl, double sig2_sur,
//     vec beta, vec psi, vec pii, double c, double shape_c, double scale_c,
//     double kappa_psi, double phi_psi,
//     double cellarea,
//     vec sigma2, mat COValpha, mat COVbeta,
//     bool updateCOV, int adaptInterval, double adaptFactorExponent, vec adapIter){
//   
//   double negativeInf = -std::numeric_limits<float>::infinity();;
//   
//   indobl = indobl-1, indsur = indsur-1, indtwhale = indtwhale-1;
//   inddist = inddist-1, indmdist = indmdist-1, indmwhale = indmwhale-1;
//   
//   int nG = distmat.n_rows, n_al = alpha.size(), n_be = beta.size(), ntzoop = eta.n_cols, n_sigma2 = sigma2.size();
//   
//   mat postAlpha0tilde = zeros(niter, 1), postAlpha = zeros(niter, n_al), postLam = zeros(niter, 2);
//   mat postSig2 = zeros(niter, 1), postKappa = zeros(niter, 2);
//   mat postBeta = zeros(niter, n_be), postC = zeros(niter, 1);
//   mat postPsi(niter, nG), postEtat = zeros(niter, nG);
//   List postEta(ntzoop);
//   for(int i = 0; i < ntzoop; i ++){
//     postEta[i] = zeros(niter, nG);
//   }
//   Rprintf("set up matrices for posterior samples\n");
//   
//   mat accprob = zeros(niter, n_sigma2);
//   vec rhat = zeros(n_sigma2), gamma1 = zeros(n_sigma2), gamma2 = zeros(n_sigma2), dummyaccprob;
//   double c0 = 1, c1 = adaptFactorExponent, ropt = 0.234, logprob;
//   mat cholCOValpha = trans( chol( sigma2[0] * ( COValpha + 0.0000000001 * diagmat(ones(n_al)) ) ) );
//   mat cholCOVbeta = trans( chol( sigma2[1] * ( COVbeta + 0.0000000001 * diagmat(ones(n_be)) ) ) );
//   Rprintf("set up adaptive proposals\n");
//   
//   vec newalpha, newbeta, newlogZobl, newlogZsur;
//   mat newlogZ, newloglam, newlam;
//   vec dummyvec;
//   
//   int count;
//   vec xvals;
//   double llprev, llthred, llnew, thetamin, thetamax, theta;
//   vec prioretat, priorpsi, newpsi;
//   mat neweta;
//   bool accept;
//   
//   mat Reta = exp(-distmat / phi_eta), Rpsi = exp(-distmat / phi_psi);
//   mat invReta = inv(Reta), invRpsi= inv(Rpsi), cholReta = chol(Reta), cholRpsi= chol(Rpsi);
//   rowvec oneinvReta = trans(ones(nG)) * invReta, oneinvRpsi = trans(ones(nG)) * invRpsi;
//   double oneinvRetaone = sum( oneinvReta * ones(nG) ), oneinvRpsione = sum( oneinvRpsi * ones(nG) );
//   Rprintf("set up GP covariance matrices\n");
//   
//   mat logZ = cppf_logZ_a0t_etat_ja(Xzoop, alpha, eta);
//   vec logZobl = cppf_logZobl(logZ, indobl);
//   vec logZsur = cppf_logZsur(logZ, indsur);
//   mat loglam = cppf_loglam_cZ_b0t_psi(Xwhale, logZ, beta, psi, indtwhale);
//   mat lam = exp(loglam);
//   Rprintf("set up initial values\n");
//   
//   
//   // Start MCMC
//   for(int s = 0; s < niter; s++) {
//     
//     // update alpha
//     if( updateCOV ){
//       if( (s+1 >= adaptInterval) && (s+1 - (adaptInterval * trunc((s+1) / adaptInterval)) == 0) ){
//         dummyaccprob = accprob.col(0);
//         rhat[0] = sum( dummyaccprob.rows(s+1-adaptInterval, s-1) ) / (adaptInterval-1);
//         gamma1[0] = 1 / pow(adapIter[0], c1);
//         gamma2[0] = c0 * gamma1[0];
//         sigma2[0] = exp( log(sigma2[0]) + gamma2[0] * (rhat[0] - ropt) );
//         
//         COValpha = COValpha + gamma1[0] * ( cov( postAlpha.rows(s+1-adaptInterval, s-1) ) - COValpha );
//         cholCOValpha = trans( chol( sigma2[0] * ( COValpha + 0.0000000001 * diagmat(ones(n_al)) ) ) );
//         
//         adapIter[0] = adapIter[0] + 1;
//       }
//     }
//     
//     newalpha = alpha + cholCOValpha * randn(n_al);
//     
//     newlogZ = cppf_logZ_a0t_etat_ja(Xzoop, newalpha, eta);
//     newlogZobl = cppf_logZobl(newlogZ, indobl);
//     newlogZsur = cppf_logZsur(newlogZ, indsur);
//     newloglam = cppf_loglam_cZ_b0t_psi(Xwhale, newlogZ, beta, psi, indtwhale);
//     newlam = exp(newloglam);
//     
//     logprob = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) -
//       cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
//       cppf_loglik_pam(Wpam, newlam, dur_dist, c, p_pam, cellarea) -
//       cppf_loglik_pam(Wpam, lam, dur_dist, c, p_pam, cellarea) +
//       cppf_loglik_obl(logYobl, newlogZobl, sig2_obl) -
//       cppf_loglik_obl(logYobl, logZobl, sig2_obl) +
//       cppf_loglik_sur(logYsur, newlogZsur, lam0, lam1, sig2_sur) -
//       cppf_loglik_sur(logYsur, logZsur, lam0, lam1, sig2_sur) +
//       MVN_logh(newalpha, zeros(n_al), diagmat(ones(n_al))/100) -
//       MVN_logh(alpha, zeros(n_al), diagmat(ones(n_al))/100);
//     
//     if (log(randu()) < logprob) {
//       
//       alpha = newalpha;
//       logZ = newlogZ;
//       logZobl = newlogZobl;
//       logZsur = newlogZsur;
//       loglam = newloglam;
//       lam = newlam;
//       
//       accprob(s, 0) = 1;
//     }
//     postAlpha.row(s) = trans(alpha);
//     if(s == 0){ Rprintf("sampled alpha\n"); }
//     
//     
//     // Gibbs update for alpha0tilde
//     alpha0tilde = cppf_gibbs_alpha0tilde(alpha.rows(0, ntzoop-1), tau2);
//     postAlpha0tilde[s] = alpha0tilde;
//     if(s == 0){ Rprintf("sampled alpha0tilde\n"); }
//     
//     
//     // update eta via Elliptical slice sampling
//     for(int t = 0; t < ntzoop; t ++){
//       neweta = eta;
//       
//       llprev = cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
//         cppf_loglik_pam(Wpam, lam, dur_dist, c, p_pam, cellarea) +
//         cppf_loglik_obl(logYobl, logZobl, sig2_obl) +
//         cppf_loglik_sur(logYsur, logZsur, lam0, lam1, sig2_sur);
//       prioretat = sqrt(kappa_eta) * trans(cholReta) * randn(nG);
//       
//       thetamin = 0;
//       thetamax = 2 * M_PI;
//       theta = thetamin + randu() * (thetamax - thetamin);
//       thetamin = theta - 2 * M_PI;
//       thetamax = theta;
//       
//       llthred = llprev + log(randu());
//       accept = false;
//       count = 0;
//       
//       while(accept == false){
//         count = count + 1;
//         
//         neweta.col(t) = eta.col(t) * cos(theta) + prioretat * sin(theta);
//         neweta.col(t) = neweta.col(t) - mean(neweta.col(t));
//         
//         newlogZ = cppf_logZ_a0t_etat_ja(Xzoop, alpha, neweta);
//         newlogZobl = cppf_logZobl(newlogZ, indobl);
//         newlogZsur = cppf_logZsur(newlogZ, indsur);
//         newloglam = cppf_loglam_cZ_b0t_psi(Xwhale, newlogZ, beta, psi, indtwhale);
//         newlam = exp(newloglam);
//         
//         llnew = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
//           cppf_loglik_pam(Wpam, newlam, dur_dist, c, p_pam, cellarea) +
//           cppf_loglik_obl(logYobl, newlogZobl, sig2_obl) +
//           cppf_loglik_sur(logYsur, newlogZsur, lam0, lam1, sig2_sur);
//         
//         if(llnew > llthred){
//           llprev = llnew;
//           accept = true;
//         } else {
//           if(theta < 0){ thetamin = theta; } else { thetamax = theta; }
//           theta = thetamin + randu() * (thetamax - thetamin);
//           if ( (count) % 1000 == 0 ) {
//             Rprintf("ESS for eta: %d iterations...\n", count);
//           }
//           if (count > 20000) stop("ESS is not converging...\n");
//         }
//       }
//       eta = neweta;
//       logZ = newlogZ;
//       logZobl = newlogZobl;
//       logZsur = newlogZsur;
//       loglam = newloglam;
//       lam = newlam;
//       
//       postEtat = as<mat>(postEta[t]);
//       postEtat.row(s) = trans(eta.col(t));
//       postEta[t] = postEtat;
//     }
//     if(s == 0){ Rprintf("sampled eta\n"); }
//     
//     
//     // update kappa_eta
//     kappa_eta = cppf_gibbs_kappa_etat(eta, invReta, 2, 2);
//     postKappa(s,0) = kappa_eta;
//     if(s == 0){ Rprintf("sampled kappa_eta\n"); }
//     
//     
//     // update lam0
//     lam0 = cppf_gibbs_lam0(logYsur, logZsur, lam1, sig2_sur, 100);
//     postLam(s, 0) = lam0;
//     if(s == 0){ Rprintf("sampled lam0\n"); }
//     
//     
//     // update lam1
//     lam1 = cppf_gibbs_lam1(logYsur, logZsur, lam0, sig2_sur, 100);
//     postLam(s, 1) = lam1;
//     if(s == 0){ Rprintf("sampled lam1\n"); }
//     
//     
//     // update sig2_sur
//     sig2_sur = cppf_gibbs_sig2_sur(logYsur, logZsur, lam0, lam1, 2, 2);
//     postSig2(s, 0) = sig2_sur;
//     if(s == 0){ Rprintf("sampled sig2_sur\n"); }
//     
//     
//     // update beta
//     if( updateCOV ){
//       if( (s+1 >= adaptInterval) && (s+1 - (adaptInterval * trunc((s+1) / adaptInterval)) == 0) ){
//         dummyaccprob = accprob.col(1);
//         rhat[1] = sum( dummyaccprob.rows(s+1-adaptInterval, s-1) ) / (adaptInterval-1);
//         gamma1[1] = 1 / pow(adapIter[1], c1);
//         gamma2[1] = c0 * gamma1[1];
//         sigma2[1] = exp( log(sigma2[1]) + gamma2[1] * (rhat[1] - ropt) );
//         
//         COVbeta = COVbeta + gamma1[1] * ( cov( postBeta.rows(s+1-adaptInterval, s-1)) - COVbeta );
//         cholCOVbeta = trans( chol( sigma2[1] * ( COVbeta + 0.0000000001 * diagmat(ones(n_be)) ) ) );
//         
//         adapIter[1] = adapIter[1] + 1;
//       }
//     }
//     
//     newbeta = beta + cholCOVbeta * randn(n_be);
//     
//     newloglam = cppf_loglam_cZ_b0t_psi(Xwhale, logZ, newbeta, psi, indtwhale);
//     newlam = exp(newloglam);
//     
//     logprob = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) -
//       cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
//       cppf_loglik_pam(Wpam, newlam, dur_dist, c, p_pam, cellarea) -
//       cppf_loglik_pam(Wpam, lam, dur_dist, c, p_pam, cellarea) +
//       MVN_logh(newbeta, zeros(n_be), diagmat(ones(n_be))/100) -
//       MVN_logh(beta, zeros(n_be), diagmat(ones(n_be))/100);
//     
//     if (log(randu()) < logprob) {
//       
//       beta = newbeta;
//       loglam = newloglam;
//       lam = newlam;
//       
//       accprob(s, 1) = 1;
//     }
//     
//     postBeta.row(s) = trans(beta);
//     if(s == 0){ Rprintf("sampled beta\n"); }
//     
//     
//     // update psi via elliptical slice sampling
//     llprev = cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
//       cppf_loglik_pam(Wpam, lam, dur_dist, c, p_pam, cellarea);
//     priorpsi = sqrt(kappa_psi) * trans(cholRpsi) * randn(nG);
//     
//     thetamin = 0;
//     thetamax = 2 * M_PI;
//     theta = thetamin + randu() * (thetamax - thetamin);
//     thetamin = theta - 2 * M_PI;
//     thetamax = theta;
//     
//     llthred = llprev + log(randu());
//     accept = false;
//     count = 0;
//     
//     while(accept == false){
//       count = count + 1;
//       
//       newpsi = psi * cos(theta) + priorpsi * sin(theta);
//       newpsi = newpsi - mean(newpsi);
//       
//       newloglam = cppf_loglam_cZ_b0t_psi(Xwhale, logZ, beta, newpsi, indtwhale);
//       newlam = exp(newloglam);
//       
//       llnew = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
//         cppf_loglik_pam(Wpam, newlam, dur_dist, c, p_pam, cellarea);
//       
//       if(llnew > llthred){
//         llprev = llnew;
//         accept = true;
//       } else {
//         if(theta < 0){ thetamin = theta; } else { thetamax = theta; }
//         theta = thetamin + randu() * (thetamax - thetamin);
//         if ( (count) % 1000 == 0 ) {
//           Rprintf("ESS for psi: %d iterations...\n", count);
//         }
//         if (count > 20000) stop("ESS is not converging...\n");
//       }
//     }
//     psi = newpsi;
//     loglam = newloglam;
//     lam = newlam;
//     
//     postPsi.row(s) = trans(psi);
//     if(s == 0){ Rprintf("sampled psi\n"); }
//     
//     
//     // update kappa_psi
//     kappa_psi = cppf_gibbs_kappa_psi(psi, invRpsi, 2, 2);
//     postKappa(s,1) = kappa_psi;
//     if(s == 0){ Rprintf("sampled kappa_psi\n"); }
//     
//     
//     // update c
//     c = cppf_gibbs(Wpam, lam, dur_dist, p_pam, cellarea, shape_c, scale_c);
//     postC(s,0) = c;
//     if(s == 0){ Rprintf("sampled c\n"); }
//     
//     
//     if ( (s+1) % 10 == 0 ) {
//       Rprintf("Generated %d samples...\n", s+1);
//     }
//   }
//   
//   
//   return Rcpp::List::create(Rcpp::Named("postAlpha0tilde") = postAlpha0tilde,
//                             Rcpp::Named("postAlpha") = postAlpha,
//                             Rcpp::Named("postEta") = postEta,
//                             Rcpp::Named("postLam") = postLam,
//                             Rcpp::Named("postKappa") = postKappa,
//                             Rcpp::Named("postSig2") = postSig2,
//                             Rcpp::Named("postBeta") = postBeta,
//                             Rcpp::Named("postC") = postC,
//                             Rcpp::Named("postPsi") = postPsi,
//                             Rcpp::Named("Accprob") = accprob,
//                             Rcpp::Named("sigma2") = sigma2,
//                             Rcpp::Named("adapIter") = adapIter,
//                             Rcpp::Named("COValpha") = COValpha,
//                             Rcpp::Named("COVbeta") = COVbeta);
//   
// }






// [[Rcpp::export]]
List cppf_fitJoint_etat_psit(
    int niter, mat distmat,
    List Xzoop, vec logYobl, vec logYsur,
    List Xwhale, mat Sdist, List Wpam, vec dur_dist,
    vec f_dist_obs, List f_dist, List p_pam,
    mat indobl, mat indsur, vec indtwhale, mat inddist, vec indmdist, vec indmwhale,
    double alpha0tilde, vec alpha, mat eta, double lam0, double lam1,
    double tau2, double kappa_eta, double phi_eta, double sig2_obl, double sig2_sur,
    vec beta, mat psi, vec pii, vec c, double shape_c, double scale_c,
    double kappa_psi, double phi_psi,
    double cellarea,
    vec sigma2, mat COValpha, mat COVbeta,
    bool updateCOV, int adaptInterval, double adaptFactorExponent, vec adapIter){

  double negativeInf = -std::numeric_limits<float>::infinity();;

  indobl = indobl-1, indsur = indsur-1, indtwhale = indtwhale-1;
  inddist = inddist-1, indmdist = indmdist-1, indmwhale = indmwhale-1;

  int nG = distmat.n_rows, n_al = alpha.size(), n_be = beta.size(), ntzoop = eta.n_cols, ntwhale = psi.n_cols, n_sigma2 = sigma2.size();

  mat postAlpha0tilde = zeros(niter, 1), postAlpha = zeros(niter, n_al), postLam = zeros(niter, 2);
  mat postSig2 = zeros(niter, 1), postKappa = zeros(niter, 2);
  mat postBeta = zeros(niter, n_be), postC = zeros(niter, c.size());
  mat postEtat = zeros(niter, nG), postPsit = zeros(niter, nG);
  List postEta(ntzoop), postPsi(ntwhale);
  for(int i = 0; i < ntzoop; i ++){ postEta[i] = zeros(niter, nG); }
  for(int i = 0; i < ntwhale; i ++){ postPsi[i] = zeros(niter, nG); }
  Rprintf("set up matrices for posterior samples\n");

  mat accprob = zeros(niter, n_sigma2);
  vec rhat = zeros(n_sigma2), gamma1 = zeros(n_sigma2), gamma2 = zeros(n_sigma2), dummyaccprob;
  double c0 = 1, c1 = adaptFactorExponent, ropt = 0.234, logprob;
  mat cholCOValpha = trans( chol( sigma2[0] * ( COValpha + 0.0000000001 * diagmat(ones(n_al)) ) ) );
  mat cholCOVbeta = trans( chol( sigma2[1] * ( COVbeta + 0.0000000001 * diagmat(ones(n_be)) ) ) );
  Rprintf("set up adaptive proposals\n");

  vec newalpha, newbeta, newlogZobl, newlogZsur;
  mat newlogZ, newloglam, newlam;
  vec dummyvec;

  int count;
  vec xvals;
  double llprev, llthred, llnew, thetamin, thetamax, theta;
  vec prioretat, priorpsit;
  mat neweta, newpsi;
  bool accept;

  mat Reta = exp(-distmat / phi_eta), Rpsi = exp(-distmat / phi_psi);
  mat invReta = inv(Reta), invRpsi= inv(Rpsi), cholReta = chol(Reta), cholRpsi= chol(Rpsi);
  rowvec oneinvReta = trans(ones(nG)) * invReta, oneinvRpsi = trans(ones(nG)) * invRpsi;
  double oneinvRetaone = sum( oneinvReta * ones(nG) ), oneinvRpsione = sum( oneinvRpsi * ones(nG) );
  Rprintf("set up GP covariance matrices\n");

  mat logZ = cppf_logZ_a0t_etat_ja(Xzoop, alpha, eta);
  vec logZobl = cppf_logZobl(logZ, indobl);
  vec logZsur = cppf_logZsur(logZ, indsur);
  mat loglam = cppf_loglam_cZ_b0t_psit(Xwhale, logZ, beta, psi, indtwhale);
  mat lam = exp(loglam);
  vec m_al = zeros(n_al), diag_al = ones(n_al)/100;
  Rprintf("set up initial values\n");


  // Start MCMC
  for(int s = 0; s < niter; s++) {

    // update alpha
    if( updateCOV ){
      if( (s+1 >= adaptInterval) && (s+1 - (adaptInterval * trunc((s+1) / adaptInterval)) == 0) ){
        dummyaccprob = accprob.col(0);
        rhat[0] = sum( dummyaccprob.rows(s+1-adaptInterval, s-1) ) / (adaptInterval-1);
        gamma1[0] = 1 / pow(adapIter[0], c1);
        gamma2[0] = c0 * gamma1[0];
        sigma2[0] = exp( log(sigma2[0]) + gamma2[0] * (rhat[0] - ropt) );

        COValpha = COValpha + gamma1[0] * ( cov( postAlpha.rows(s+1-adaptInterval, s-1) ) - COValpha );
        cholCOValpha = trans( chol( sigma2[0] * ( COValpha + 0.0000000001 * diagmat(ones(n_al)) ) ) );

        adapIter[0] = adapIter[0] + 1;
      }
    }

    newalpha = alpha + cholCOValpha * randn(n_al);

    newlogZ = cppf_logZ_a0t_etat_ja(Xzoop, newalpha, eta);
    newlogZobl = cppf_logZobl(newlogZ, indobl);
    newlogZsur = cppf_logZsur(newlogZ, indsur);
    newloglam = cppf_loglam_cZ_b0t_psit(Xwhale, newlogZ, beta, psi, indtwhale);
    newlam = exp(newloglam);

    for(int t = 0; t < ntzoop; t ++){
      m_al[t] = alpha0tilde;
      diag_al[t] = 1/tau2;
    }

    logprob = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) -
      cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
      cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) -
      cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
      cppf_loglik_obl(logYobl, newlogZobl, sig2_obl) -
      cppf_loglik_obl(logYobl, logZobl, sig2_obl) +
      cppf_loglik_sur(logYsur, newlogZsur, lam0, lam1, sig2_sur) -
      cppf_loglik_sur(logYsur, logZsur, lam0, lam1, sig2_sur) +
      MVN_logh(newalpha, m_al, diagmat(diag_al)) -
      MVN_logh(alpha, m_al, diagmat(diag_al));

    if (log(randu()) < logprob) {

      alpha = newalpha;
      logZ = newlogZ;
      logZobl = newlogZobl;
      logZsur = newlogZsur;
      loglam = newloglam;
      lam = newlam;

      accprob(s, 0) = 1;
    }
    postAlpha.row(s) = trans(alpha);
    if(s == 0){ Rprintf("sampled alpha\n"); }


    // Gibbs update for alpha0tilde
    alpha0tilde = cppf_gibbs_alpha0tilde(alpha.rows(0, ntzoop-1), tau2);
    postAlpha0tilde[s] = alpha0tilde;
    if(s == 0){ Rprintf("sampled alpha0tilde\n"); }



    // update eta via Elliptical slice sampling
    for(int t = 0; t < ntzoop; t ++){
      neweta = eta;

      llprev = cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
        cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
        cppf_loglik_obl(logYobl, logZobl, sig2_obl) +
        cppf_loglik_sur(logYsur, logZsur, lam0, lam1, sig2_sur);
      prioretat = sqrt(kappa_eta) * trans(cholReta) * randn(nG);

      thetamin = 0;
      thetamax = 2 * M_PI;
      theta = thetamin + randu() * (thetamax - thetamin);
      thetamin = theta - 2 * M_PI;
      thetamax = theta;

      llthred = llprev + log(randu());
      accept = false;
      count = 0;

      while(accept == false){
        count = count + 1;

        neweta.col(t) = eta.col(t) * cos(theta) + prioretat * sin(theta);
        neweta.col(t) = neweta.col(t) - mean(neweta.col(t));

        newlogZ = cppf_logZ_a0t_etat_ja(Xzoop, alpha, neweta);
        newlogZobl = cppf_logZobl(newlogZ, indobl);
        newlogZsur = cppf_logZsur(newlogZ, indsur);
        newloglam = cppf_loglam_cZ_b0t_psit(Xwhale, newlogZ, beta, psi, indtwhale);
        newlam = exp(newloglam);

        llnew = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
          cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) +
          cppf_loglik_obl(logYobl, newlogZobl, sig2_obl) +
          cppf_loglik_sur(logYsur, newlogZsur, lam0, lam1, sig2_sur);

        if(llnew > llthred){
          llprev = llnew;
          accept = true;
        } else {
          if(theta < 0){ thetamin = theta; } else { thetamax = theta; }
          theta = thetamin + randu() * (thetamax - thetamin);
          if ( (count) % 1000 == 0 ) {
            Rprintf("ESS for eta: %d iterations...\n", count);
          }
          if (count > 20000) stop("ESS is not converging...\n");
        }
      }
      eta = neweta;
      logZ = newlogZ;
      logZobl = newlogZobl;
      logZsur = newlogZsur;
      loglam = newloglam;
      lam = newlam;

      postEtat = as<mat>(postEta[t]);
      postEtat.row(s) = trans(eta.col(t));
      postEta[t] = postEtat;
    }
    if(s == 0){ Rprintf("sampled eta\n"); }


    // update kappa_eta
    kappa_eta = cppf_gibbs_kappa_etat(eta, invReta, 2, 2);
    postKappa(s,0) = kappa_eta;
    if(s == 0){ Rprintf("sampled kappa_eta\n"); }


    // update lam0
    lam0 = cppf_gibbs_lam0(logYsur, logZsur, lam1, sig2_sur, 100);
    postLam(s, 0) = lam0;
    if(s == 0){ Rprintf("sampled lam0\n"); }


    // update lam1
    lam1 = cppf_gibbs_lam1(logYsur, logZsur, lam0, sig2_sur, 100);
    postLam(s, 1) = lam1;
    if(s == 0){ Rprintf("sampled lam1\n"); }


    // update sig2_sur
    sig2_sur = cppf_gibbs_sig2_sur(logYsur, logZsur, lam0, lam1, 2, 2);
    postSig2(s, 0) = sig2_sur;
    if(s == 0){ Rprintf("sampled sig2_sur\n"); }


    // update beta
    if( updateCOV ){
      if( (s+1 >= adaptInterval) && (s+1 - (adaptInterval * trunc((s+1) / adaptInterval)) == 0) ){
        dummyaccprob = accprob.col(1);
        rhat[1] = sum( dummyaccprob.rows(s+1-adaptInterval, s-1) ) / (adaptInterval-1);
        gamma1[1] = 1 / pow(adapIter[1], c1);
        gamma2[1] = c0 * gamma1[1];
        sigma2[1] = exp( log(sigma2[1]) + gamma2[1] * (rhat[1] - ropt) );

        COVbeta = COVbeta + gamma1[1] * ( cov( postBeta.rows(s+1-adaptInterval, s-1)) - COVbeta );
        cholCOVbeta = trans( chol( sigma2[1] * ( COVbeta + 0.0000000001 * diagmat(ones(n_be)) ) ) );

        adapIter[1] = adapIter[1] + 1;
      }
    }

    newbeta = beta + cholCOVbeta * randn(n_be);

    newloglam = cppf_loglam_cZ_b0t_psit(Xwhale, logZ, newbeta, psi, indtwhale);
    newlam = exp(newloglam);

    logprob = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) -
      cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
      cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) -
      cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
      MVN_logh(newbeta, zeros(n_be), diagmat(ones(n_be))/100) -
      MVN_logh(beta, zeros(n_be), diagmat(ones(n_be))/100);

    if (log(randu()) < logprob) {

      beta = newbeta;
      loglam = newloglam;
      lam = newlam;

      accprob(s, 1) = 1;
    }

    postBeta.row(s) = trans(beta);
    if(s == 0){ Rprintf("sampled beta\n"); }



    // update psi via Elliptical slice sampling
    for(int t = 0; t < ntwhale; t ++){
      newpsi = psi;

      llprev = cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
        cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea);
      priorpsit = sqrt(kappa_psi) * trans(cholRpsi) * randn(nG);

      thetamin = 0;
      thetamax = 2 * M_PI;
      theta = thetamin + randu() * (thetamax - thetamin);
      thetamin = theta - 2 * M_PI;
      thetamax = theta;

      llthred = llprev + log(randu());
      accept = false;
      count = 0;

      while(accept == false){
        count = count + 1;

        newpsi.col(t) = psi.col(t) * cos(theta) + priorpsit * sin(theta);
        newpsi.col(t) = newpsi.col(t) - mean(newpsi.col(t));

        newloglam = cppf_loglam_cZ_b0t_psit(Xwhale, logZ, beta, newpsi, indtwhale);
        newlam = exp(newloglam);

        llnew = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
          cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea);

        if(llnew > llthred){
          llprev = llnew;
          accept = true;
        } else {
          if(theta < 0){ thetamin = theta; } else { thetamax = theta; }
          theta = thetamin + randu() * (thetamax - thetamin);
          if ( (count) % 1000 == 0 ) {
            Rprintf("ESS for psi: %d iterations...\n", count);
          }
          if (count > 20000) stop("ESS is not converging...\n");
        }
      }
      psi = newpsi;
      loglam = newloglam;
      lam = newlam;

      postPsit = as<mat>(postPsi[t]);
      postPsit.row(s) = trans(psi.col(t));
      postPsi[t] = postPsit;
    }
    if(s == 0){ Rprintf("sampled psi\n"); }


    // update kappa_psi
    kappa_psi = cppf_gibbs_kappa_psit(psi, invRpsi, 2, 2);
    postKappa(s,1) = kappa_psi;
    if(s == 0){ Rprintf("sampled kappa_psi\n"); }


    // update c
    c = cppf_gibbs_ct(Wpam, lam, dur_dist, p_pam, cellarea, shape_c, scale_c);
    postC.row(s) = trans(c);
    if(s == 0){ Rprintf("sampled c\n"); }


    if ( (s+1) % 10 == 0 ) {
      Rprintf("Generated %d samples...\n", s+1);
    }
  }


  return Rcpp::List::create(Rcpp::Named("postAlpha0tilde") = postAlpha0tilde,
                            Rcpp::Named("postAlpha") = postAlpha,
                            Rcpp::Named("postEta") = postEta,
                            Rcpp::Named("postLam") = postLam,
                            Rcpp::Named("postKappa") = postKappa,
                            Rcpp::Named("postSig2") = postSig2,
                            Rcpp::Named("postBeta") = postBeta,
                            Rcpp::Named("postC") = postC,
                            Rcpp::Named("postPsi") = postPsi,
                            Rcpp::Named("Accprob") = accprob,
                            Rcpp::Named("sigma2") = sigma2,
                            Rcpp::Named("adapIter") = adapIter,
                            Rcpp::Named("COValpha") = COValpha,
                            Rcpp::Named("COVbeta") = COVbeta);

}



// [[Rcpp::export]]
List cppf_fitJoint_etat_psi(
    int niter, mat distmat,
    List Xzoop, vec logYobl, vec logYsur,
    List Xwhale, mat Sdist, List Wpam, vec dur_dist,
    vec f_dist_obs, List f_dist, List p_pam,
    mat indobl, mat indsur, vec indtwhale, mat inddist, vec indmdist, vec indmwhale,
    double alpha0tilde, vec alpha, mat eta, double lam0, double lam1,
    double tau2, double kappa_eta, double phi_eta, double sig2_obl, double sig2_sur,
    vec beta, vec psi, vec pii, vec c, double shape_c, double scale_c,
    double kappa_psi, double phi_psi,
    double cellarea,
    vec sigma2, mat COValpha, mat COVbeta,
    bool updateCOV, int adaptInterval, double adaptFactorExponent, vec adapIter){
  
  double negativeInf = -std::numeric_limits<float>::infinity();;
  
  indobl = indobl-1, indsur = indsur-1, indtwhale = indtwhale-1;
  inddist = inddist-1, indmdist = indmdist-1, indmwhale = indmwhale-1;
  
  int nG = distmat.n_rows, n_al = alpha.size(), n_be = beta.size(), ntzoop = eta.n_cols, n_sigma2 = sigma2.size();
  
  mat postAlpha0tilde = zeros(niter, 1), postAlpha = zeros(niter, n_al), postLam = zeros(niter, 2);
  mat postSig2 = zeros(niter, 1), postKappa = zeros(niter, 2);
  mat postBeta = zeros(niter, n_be), postC = zeros(niter, c.size());
  mat postPsi(niter, nG), postEtat = zeros(niter, nG);
  List postEta(ntzoop);
  for(int i = 0; i < ntzoop; i ++){
    postEta[i] = zeros(niter, nG);
  }
  Rprintf("set up matrices for posterior samples\n");
  
  mat accprob = zeros(niter, n_sigma2);
  vec rhat = zeros(n_sigma2), gamma1 = zeros(n_sigma2), gamma2 = zeros(n_sigma2), dummyaccprob;
  double c0 = 1, c1 = adaptFactorExponent, ropt = 0.234, logprob;
  mat cholCOValpha = trans( chol( sigma2[0] * ( COValpha + 0.0000000001 * diagmat(ones(n_al)) ) ) );
  mat cholCOVbeta = trans( chol( sigma2[1] * ( COVbeta + 0.0000000001 * diagmat(ones(n_be)) ) ) );
  Rprintf("set up adaptive proposals\n");
  
  vec newalpha, newbeta, newlogZobl, newlogZsur;
  mat newlogZ, newloglam, newlam;
  vec dummyvec;
  
  int count;
  vec xvals;
  double llprev, llthred, llnew, thetamin, thetamax, theta;
  vec prioretat, priorpsi, newpsi;
  mat neweta;
  bool accept;
  
  mat Reta = exp(-distmat / phi_eta), Rpsi = exp(-distmat / phi_psi);
  mat invReta = inv(Reta), invRpsi= inv(Rpsi), cholReta = chol(Reta), cholRpsi= chol(Rpsi);
  rowvec oneinvReta = trans(ones(nG)) * invReta, oneinvRpsi = trans(ones(nG)) * invRpsi;
  double oneinvRetaone = sum( oneinvReta * ones(nG) ), oneinvRpsione = sum( oneinvRpsi * ones(nG) );
  Rprintf("set up GP covariance matrices\n");
  
  mat logZ = cppf_logZ_a0t_etat_ja(Xzoop, alpha, eta);
  vec logZobl = cppf_logZobl(logZ, indobl);
  vec logZsur = cppf_logZsur(logZ, indsur);
  mat loglam = cppf_loglam_cZ_b0t_psi(Xwhale, logZ, beta, psi, indtwhale);
  mat lam = exp(loglam);
  Rprintf("set up initial values\n");
  
  
  // Start MCMC
  for(int s = 0; s < niter; s++) {
    
    // update alpha
    if( updateCOV ){
      if( (s+1 >= adaptInterval) && (s+1 - (adaptInterval * trunc((s+1) / adaptInterval)) == 0) ){
        dummyaccprob = accprob.col(0);
        rhat[0] = sum( dummyaccprob.rows(s+1-adaptInterval, s-1) ) / (adaptInterval-1);
        gamma1[0] = 1 / pow(adapIter[0], c1);
        gamma2[0] = c0 * gamma1[0];
        sigma2[0] = exp( log(sigma2[0]) + gamma2[0] * (rhat[0] - ropt) );
        
        COValpha = COValpha + gamma1[0] * ( cov( postAlpha.rows(s+1-adaptInterval, s-1) ) - COValpha );
        cholCOValpha = trans( chol( sigma2[0] * ( COValpha + 0.0000000001 * diagmat(ones(n_al)) ) ) );
        
        adapIter[0] = adapIter[0] + 1;
      }
    }
    
    newalpha = alpha + cholCOValpha * randn(n_al);
    
    newlogZ = cppf_logZ_a0t_etat_ja(Xzoop, newalpha, eta);
    newlogZobl = cppf_logZobl(newlogZ, indobl);
    newlogZsur = cppf_logZsur(newlogZ, indsur);
    newloglam = cppf_loglam_cZ_b0t_psi(Xwhale, newlogZ, beta, psi, indtwhale);
    newlam = exp(newloglam);
    
    logprob = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) -
      cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
      cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) -
      cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
      cppf_loglik_obl(logYobl, newlogZobl, sig2_obl) -
      cppf_loglik_obl(logYobl, logZobl, sig2_obl) +
      cppf_loglik_sur(logYsur, newlogZsur, lam0, lam1, sig2_sur) -
      cppf_loglik_sur(logYsur, logZsur, lam0, lam1, sig2_sur) +
      MVN_logh(newalpha, zeros(n_al), diagmat(ones(n_al))/100) -
      MVN_logh(alpha, zeros(n_al), diagmat(ones(n_al))/100);
    
    if (log(randu()) < logprob) {
      
      alpha = newalpha;
      logZ = newlogZ;
      logZobl = newlogZobl;
      logZsur = newlogZsur;
      loglam = newloglam;
      lam = newlam;
      
      accprob(s, 0) = 1;
    }
    postAlpha.row(s) = trans(alpha);
    if(s == 0){ Rprintf("sampled alpha\n"); }
    
    
    // Gibbs update for alpha0tilde
    alpha0tilde = cppf_gibbs_alpha0tilde(alpha.rows(0, ntzoop-1), tau2);
    postAlpha0tilde[s] = alpha0tilde;
    if(s == 0){ Rprintf("sampled alpha0tilde\n"); }
    
    
    // update eta via Elliptical slice sampling
    for(int t = 0; t < ntzoop; t ++){
      neweta = eta;
      
      llprev = cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
        cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
        cppf_loglik_obl(logYobl, logZobl, sig2_obl) +
        cppf_loglik_sur(logYsur, logZsur, lam0, lam1, sig2_sur);
      prioretat = sqrt(kappa_eta) * trans(cholReta) * randn(nG);
      
      thetamin = 0;
      thetamax = 2 * M_PI;
      theta = thetamin + randu() * (thetamax - thetamin);
      thetamin = theta - 2 * M_PI;
      thetamax = theta;
      
      llthred = llprev + log(randu());
      accept = false;
      count = 0;
      
      while(accept == false){
        count = count + 1;
        
        neweta.col(t) = eta.col(t) * cos(theta) + prioretat * sin(theta);
        neweta.col(t) = neweta.col(t) - mean(neweta.col(t));
        
        newlogZ = cppf_logZ_a0t_etat_ja(Xzoop, alpha, neweta);
        newlogZobl = cppf_logZobl(newlogZ, indobl);
        newlogZsur = cppf_logZsur(newlogZ, indsur);
        newloglam = cppf_loglam_cZ_b0t_psi(Xwhale, newlogZ, beta, psi, indtwhale);
        newlam = exp(newloglam);
        
        llnew = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
          cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) +
          cppf_loglik_obl(logYobl, newlogZobl, sig2_obl) +
          cppf_loglik_sur(logYsur, newlogZsur, lam0, lam1, sig2_sur);
        
        if(llnew > llthred){
          llprev = llnew;
          accept = true;
        } else {
          if(theta < 0){ thetamin = theta; } else { thetamax = theta; }
          theta = thetamin + randu() * (thetamax - thetamin);
          if ( (count) % 1000 == 0 ) {
            Rprintf("ESS for eta: %d iterations...\n", count);
          }
          if (count > 20000) stop("ESS is not converging...\n");
        }
      }
      eta = neweta;
      logZ = newlogZ;
      logZobl = newlogZobl;
      logZsur = newlogZsur;
      loglam = newloglam;
      lam = newlam;
      
      postEtat = as<mat>(postEta[t]);
      postEtat.row(s) = trans(eta.col(t));
      postEta[t] = postEtat;
    }
    if(s == 0){ Rprintf("sampled eta\n"); }
    
    
    // update kappa_eta
    kappa_eta = cppf_gibbs_kappa_etat(eta, invReta, 2, 2);
    postKappa(s,0) = kappa_eta;
    if(s == 0){ Rprintf("sampled kappa_eta\n"); }
    
    
    // update lam0
    lam0 = cppf_gibbs_lam0(logYsur, logZsur, lam1, sig2_sur, 100);
    postLam(s, 0) = lam0;
    if(s == 0){ Rprintf("sampled lam0\n"); }
    
    
    // update lam1
    lam1 = cppf_gibbs_lam1(logYsur, logZsur, lam0, sig2_sur, 100);
    postLam(s, 1) = lam1;
    if(s == 0){ Rprintf("sampled lam1\n"); }
    
    
    // update sig2_sur
    sig2_sur = cppf_gibbs_sig2_sur(logYsur, logZsur, lam0, lam1, 2, 2);
    postSig2(s, 0) = sig2_sur;
    if(s == 0){ Rprintf("sampled sig2_sur\n"); }
    
    
    // update beta
    if( updateCOV ){
      if( (s+1 >= adaptInterval) && (s+1 - (adaptInterval * trunc((s+1) / adaptInterval)) == 0) ){
        dummyaccprob = accprob.col(1);
        rhat[1] = sum( dummyaccprob.rows(s+1-adaptInterval, s-1) ) / (adaptInterval-1);
        gamma1[1] = 1 / pow(adapIter[1], c1);
        gamma2[1] = c0 * gamma1[1];
        sigma2[1] = exp( log(sigma2[1]) + gamma2[1] * (rhat[1] - ropt) );
        
        COVbeta = COVbeta + gamma1[1] * ( cov( postBeta.rows(s+1-adaptInterval, s-1)) - COVbeta );
        cholCOVbeta = trans( chol( sigma2[1] * ( COVbeta + 0.0000000001 * diagmat(ones(n_be)) ) ) );
        
        adapIter[1] = adapIter[1] + 1;
      }
    }
    
    newbeta = beta + cholCOVbeta * randn(n_be);
    
    newloglam = cppf_loglam_cZ_b0t_psi(Xwhale, logZ, newbeta, psi, indtwhale);
    newlam = exp(newloglam);
    
    logprob = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) -
      cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
      cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) -
      cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
      MVN_logh(newbeta, zeros(n_be), diagmat(ones(n_be))/100) -
      MVN_logh(beta, zeros(n_be), diagmat(ones(n_be))/100);
    
    if (log(randu()) < logprob) {
      
      beta = newbeta;
      loglam = newloglam;
      lam = newlam;
      
      accprob(s, 1) = 1;
    }
    
    postBeta.row(s) = trans(beta);
    if(s == 0){ Rprintf("sampled beta\n"); }
    
    
    // update psi via elliptical slice sampling
    llprev = cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
      cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea);
    priorpsi = sqrt(kappa_psi) * trans(cholRpsi) * randn(nG);
    
    thetamin = 0;
    thetamax = 2 * M_PI;
    theta = thetamin + randu() * (thetamax - thetamin);
    thetamin = theta - 2 * M_PI;
    thetamax = theta;
    
    llthred = llprev + log(randu());
    accept = false;
    count = 0;
    
    while(accept == false){
      count = count + 1;
      
      newpsi = psi * cos(theta) + priorpsi * sin(theta);
      newpsi = newpsi - mean(newpsi);
      
      newloglam = cppf_loglam_cZ_b0t_psi(Xwhale, logZ, beta, newpsi, indtwhale);
      newlam = exp(newloglam);
      
      llnew = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
        cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea);
      
      if(llnew > llthred){
        llprev = llnew;
        accept = true;
      } else {
        if(theta < 0){ thetamin = theta; } else { thetamax = theta; }
        theta = thetamin + randu() * (thetamax - thetamin);
        if ( (count) % 1000 == 0 ) {
          Rprintf("ESS for psi: %d iterations...\n", count);
        }
        if (count > 20000) stop("ESS is not converging...\n");
      }
    }
    psi = newpsi;
    loglam = newloglam;
    lam = newlam;
    
    postPsi.row(s) = trans(psi);
    if(s == 0){ Rprintf("sampled psi\n"); }
    
    
    // update kappa_psi
    kappa_psi = cppf_gibbs_kappa_psi(psi, invRpsi, 2, 2);
    postKappa(s,1) = kappa_psi;
    if(s == 0){ Rprintf("sampled kappa_psi\n"); }
    
    
    // update c
    c = cppf_gibbs_ct(Wpam, lam, dur_dist, p_pam, cellarea, shape_c, scale_c);
    postC.row(s) = trans(c);
    if(s == 0){ Rprintf("sampled c\n"); }
    
    
    if ( (s+1) % 10 == 0 ) {
      Rprintf("Generated %d samples...\n", s+1);
    }
  }
  
  
  return Rcpp::List::create(Rcpp::Named("postAlpha0tilde") = postAlpha0tilde,
                            Rcpp::Named("postAlpha") = postAlpha,
                            Rcpp::Named("postEta") = postEta,
                            Rcpp::Named("postLam") = postLam,
                            Rcpp::Named("postKappa") = postKappa,
                            Rcpp::Named("postSig2") = postSig2,
                            Rcpp::Named("postBeta") = postBeta,
                            Rcpp::Named("postC") = postC,
                            Rcpp::Named("postPsi") = postPsi,
                            Rcpp::Named("Accprob") = accprob,
                            Rcpp::Named("sigma2") = sigma2,
                            Rcpp::Named("adapIter") = adapIter,
                            Rcpp::Named("COValpha") = COValpha,
                            Rcpp::Named("COVbeta") = COVbeta);
  
}





// [[Rcpp::export]]
List cppf_fitWhale_psit(
    int niter, mat distmat,
    List Xwhale, mat Sdist, List Wpam, vec dur_dist,
    vec f_dist_obs, List f_dist, List p_pam,
    mat inddist, vec indmdist, vec indmwhale,
    vec beta, mat psi, vec pii, vec c, double shape_c, double scale_c,
    double kappa_psi, double phi_psi,double cellarea,
    vec sigma2, mat COVbeta, bool updateCOV,
    int adaptInterval, double adaptFactorExponent, vec adapIter){

  double negativeInf = -std::numeric_limits<float>::infinity();;

  inddist = inddist-1, indmdist = indmdist-1, indmwhale = indmwhale-1;

  int nG = distmat.n_rows, n_be = beta.size(), ntwhale = psi.n_cols, n_sigma2 = sigma2.size();

  mat postKappa = zeros(niter, 1);
  mat postBeta = zeros(niter, n_be), postC = zeros(niter, c.size());
  mat postPsit = zeros(niter, nG);
  List postPsi(ntwhale);
  for(int i = 0; i < ntwhale; i ++){ postPsi[i] = zeros(niter, nG); }
  Rprintf("set up matrices for posterior samples\n");

  mat accprob = zeros(niter, n_sigma2);
  vec rhat = zeros(n_sigma2), gamma1 = zeros(n_sigma2), gamma2 = zeros(n_sigma2), dummyaccprob;
  double c0 = 1, c1 = adaptFactorExponent, ropt = 0.234, logprob;
  mat cholCOVbeta = trans( chol( sigma2[0] * ( COVbeta + 0.0000000001 * diagmat(ones(n_be)) ) ) );
  Rprintf("set up adaptive proposals\n");

  vec newbeta;
  mat newloglam, newlam;
  vec dummyvec;

  int count;
  vec xvals;
  double llprev, llthred, llnew, thetamin, thetamax, theta;
  vec priorpsit;
  mat newpsi;
  bool accept;

  mat Rpsi = exp(-distmat / phi_psi);
  mat invRpsi = inv(Rpsi), cholRpsi= chol(Rpsi);
  rowvec oneinvRpsi = trans(ones(nG)) * invRpsi;
  double oneinvRpsione = sum( oneinvRpsi * ones(nG) );
  Rprintf("set up GP covariance matrices\n");

  mat loglam = cppf_loglam_b0t_psit_woZ(Xwhale, beta, psi);
  mat lam = exp(loglam);
  Rprintf("set up initial values\n");


  // Start MCMC
  for(int s = 0; s < niter; s++) {

    // update beta
    if( updateCOV ){
      if( (s+1 >= adaptInterval) && (s+1 - (adaptInterval * trunc((s+1) / adaptInterval)) == 0) ){
        dummyaccprob = accprob.col(0);
        rhat[0] = sum( dummyaccprob.rows(s+1-adaptInterval, s-1) ) / (adaptInterval-1);
        gamma1[0] = 1 / pow(adapIter[0], c1);
        gamma2[0] = c0 * gamma1[0];
        sigma2[0] = exp( log(sigma2[0]) + gamma2[0] * (rhat[0] - ropt) );

        COVbeta = COVbeta + gamma1[0] * ( cov( postBeta.rows(s+1-adaptInterval, s-1)) - COVbeta );
        cholCOVbeta = trans( chol( sigma2[0] * ( COVbeta + 0.0000000001 * diagmat(ones(n_be)) ) ) );

        adapIter[0] = adapIter[0] + 1;
      }
    }

    newbeta = beta + cholCOVbeta * randn(n_be);

    newloglam = cppf_loglam_b0t_psit_woZ(Xwhale, newbeta, psi);
    newlam = exp(newloglam);

    logprob = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) -
      cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
      cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) -
      cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
      MVN_logh(newbeta, zeros(n_be), diagmat(ones(n_be))/100) -
      MVN_logh(beta, zeros(n_be), diagmat(ones(n_be))/100);

    if (log(randu()) < logprob) {

      beta = newbeta;
      loglam = newloglam;
      lam = newlam;

      accprob(s, 0) = 1;
    }

    postBeta.row(s) = trans(beta);
    if(s == 0){ Rprintf("sampled beta\n"); }



    // update psi via Elliptical slice sampling
    for(int t = 0; t < ntwhale; t ++){
      newpsi = psi;

      llprev = cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
        cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea);
      priorpsit = sqrt(kappa_psi) * trans(cholRpsi) * randn(nG);

      thetamin = 0;
      thetamax = 2 * M_PI;
      theta = thetamin + randu() * (thetamax - thetamin);
      thetamin = theta - 2 * M_PI;
      thetamax = theta;

      llthred = llprev + log(randu());
      accept = false;
      count = 0;

      while(accept == false){
        count = count + 1;

        newpsi.col(t) = psi.col(t) * cos(theta) + priorpsit * sin(theta);
        newpsi.col(t) = newpsi.col(t) - mean(newpsi.col(t));

        newloglam = cppf_loglam_b0t_psit_woZ(Xwhale, beta, newpsi);
        newlam = exp(newloglam);

        llnew = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
          cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea);

        if(llnew > llthred){
          llprev = llnew;
          accept = true;
        } else {
          if(theta < 0){ thetamin = theta; } else { thetamax = theta; }
          theta = thetamin + randu() * (thetamax - thetamin);
          if ( (count) % 1000 == 0 ) {
            Rprintf("ESS for psi: %d iterations...\n", count);
          }
          if (count > 20000) stop("ESS is not converging...\n");
        }
      }
      psi = newpsi;
      loglam = newloglam;
      lam = newlam;

      postPsit = as<mat>(postPsi[t]);
      postPsit.row(s) = trans(psi.col(t));
      postPsi[t] = postPsit;
    }
    if(s == 0){ Rprintf("sampled psi\n"); }


    // update kappa_psi
    kappa_psi = cppf_gibbs_kappa_psit(psi, invRpsi, 2, 2);
    postKappa(s,0) = kappa_psi;
    if(s == 0){ Rprintf("sampled kappa_psi\n"); }


    // update c
    c = cppf_gibbs_ct(Wpam, lam, dur_dist, p_pam, cellarea, shape_c, scale_c);
    postC.row(s) = trans(c);
    if(s == 0){ Rprintf("sampled c\n"); }


    if ( (s+1) % 10 == 0 ) {
      Rprintf("Generated %d samples...\n", s+1);
    }
  }


  return Rcpp::List::create(Rcpp::Named("postBeta") = postBeta,
                            Rcpp::Named("postKappa") = postKappa,
                            Rcpp::Named("postC") = postC,
                            Rcpp::Named("postPsi") = postPsi,
                            Rcpp::Named("Accprob") = accprob,
                            Rcpp::Named("sigma2") = sigma2,
                            Rcpp::Named("adapIter") = adapIter,
                            Rcpp::Named("COVbeta") = COVbeta);

}



// // [[Rcpp::export]]
// List cppf_fitWhale_b0t_psi_ct(
//     int niter, mat distmat,
//     List Xwhale, mat Sdist, List Wpam, vec dur_dist,
//     vec f_dist_obs, List f_dist, List p_pam,
//     mat inddist, vec indmdist, vec indmwhale,
//     vec beta, vec psi, vec pii, vec c, double shape_c, double scale_c,
//     double kappa_psi, double phi_psi,
//     double cellarea,
//     vec sigma2, mat COVbeta,
//     bool updateCOV, int adaptInterval, double adaptFactorExponent, vec adapIter){
//   
//   double negativeInf = -std::numeric_limits<float>::infinity();;
//   
//   inddist = inddist-1, indmdist = indmdist-1, indmwhale = indmwhale-1;
//   
//   int nG = distmat.n_rows, n_be = beta.size(), n_sigma2 = sigma2.size();
//   
//   mat postKappa = zeros(niter, 1);
//   mat postBeta = zeros(niter, n_be), postC = zeros(niter, c.size());
//   mat postPsi(niter, nG);
//   Rprintf("set up matrices for posterior samples\n");
//   
//   mat accprob = zeros(niter, n_sigma2);
//   vec rhat = zeros(n_sigma2), gamma1 = zeros(n_sigma2), gamma2 = zeros(n_sigma2), dummyaccprob;
//   double c0 = 1, c1 = adaptFactorExponent, ropt = 0.234, logprob;
//   mat cholCOVbeta = trans( chol( sigma2[0] * ( COVbeta + 0.0000000001 * diagmat(ones(n_be)) ) ) );
//   Rprintf("set up adaptive proposals\n");
//   
//   mat newloglam, newlam;
//   vec newbeta, dummyvec;
//   
//   int count;
//   vec xvals;
//   double llprev, llthred, llnew, thetamin, thetamax, theta;
//   vec priorpsi, newpsi;
//   bool accept;
//   
//   mat Rpsi = exp(-distmat / phi_psi);
//   mat invRpsi= inv(Rpsi), cholRpsi= chol(Rpsi);
//   rowvec oneinvRpsi = trans(ones(nG)) * invRpsi;
//   double oneinvRpsione = sum( oneinvRpsi * ones(nG) );
//   Rprintf("set up GP covariance matrices\n");
//   
//   mat loglam = cppf_loglam_b0t_psi_woZ(Xwhale, beta, psi);
//   mat lam = exp(loglam);
//   Rprintf("set up initial values\n");
//   
//   
//   // Start MCMC
//   for(int s = 0; s < niter; s++) {
//     
//     // update beta
//     if( updateCOV ){
//       if( (s+1 >= adaptInterval) && (s+1 - (adaptInterval * trunc((s+1) / adaptInterval)) == 0) ){
//         dummyaccprob = accprob.col(0);
//         rhat[0] = sum( dummyaccprob.rows(s+1-adaptInterval, s-1) ) / (adaptInterval-1);
//         gamma1[0] = 1 / pow(adapIter[0], c1);
//         gamma2[0] = c0 * gamma1[0];
//         sigma2[0] = exp( log(sigma2[0]) + gamma2[0] * (rhat[0] - ropt) );
//         
//         COVbeta = COVbeta + gamma1[0] * ( cov( postBeta.rows(s+1-adaptInterval, s-1)) - COVbeta );
//         cholCOVbeta = trans( chol( sigma2[0] * ( COVbeta + 0.0000000001 * diagmat(ones(n_be)) ) ) );
//         
//         adapIter[0] = adapIter[0] + 1;
//       }
//     }
//     
//     newbeta = beta + cholCOVbeta * randn(n_be);
//     
//     newloglam = cppf_loglam_b0t_psi_woZ(Xwhale, newbeta, psi);
//     newlam = exp(newloglam);
//     
//     logprob = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) -
//       cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
//       cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) -
//       cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
//       MVN_logh(newbeta, zeros(n_be), diagmat(ones(n_be))/100) -
//       MVN_logh(beta, zeros(n_be), diagmat(ones(n_be))/100);
//     
//     if (log(randu()) < logprob) {
//       
//       beta = newbeta;
//       loglam = newloglam;
//       lam = newlam;
//       
//       accprob(s, 0) = 1;
//     }
//     
//     postBeta.row(s) = trans(beta);
//     if(s == 0){ Rprintf("sampled beta\n"); }
//     
//     
//     // update psi via elliptical slice sampling
//     llprev = cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
//       cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea);
//     priorpsi = sqrt(kappa_psi) * trans(cholRpsi) * randn(nG);
//     
//     thetamin = 0;
//     thetamax = 2 * M_PI;
//     theta = thetamin + randu() * (thetamax - thetamin);
//     thetamin = theta - 2 * M_PI;
//     thetamax = theta;
//     
//     llthred = llprev + log(randu());
//     accept = false;
//     count = 0;
//     
//     while(accept == false){
//       count = count + 1;
//       
//       newpsi = psi * cos(theta) + priorpsi * sin(theta);
//       newpsi = newpsi - mean(newpsi);
//       
//       newloglam = cppf_loglam_b0t_psi_woZ(Xwhale, beta, newpsi);
//       newlam = exp(newloglam);
//       
//       llnew = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
//         cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea);
//       
//       if(llnew > llthred){
//         llprev = llnew;
//         accept = true;
//       } else {
//         if(theta < 0){ thetamin = theta; } else { thetamax = theta; }
//         theta = thetamin + randu() * (thetamax - thetamin);
//         if ( (count) % 1000 == 0 ) {
//           Rprintf("ESS for psi: %d iterations...\n", count);
//         }
//         if (count > 20000) stop("ESS is not converging...\n");
//       }
//     }
//     psi = newpsi;
//     loglam = newloglam;
//     lam = newlam;
//     
//     postPsi.row(s) = trans(psi);
//     if(s == 0){ Rprintf("sampled psi\n"); }
//     
//     
//     // update kappa_psi
//     kappa_psi = cppf_gibbs_kappa_psi(psi, invRpsi, 2, 2);
//     postKappa(s,0) = kappa_psi;
//     if(s == 0){ Rprintf("sampled kappa_psi\n"); }
//     
//     
//     // update c
//     c = cppf_gibbs_ct(Wpam, lam, dur_dist, p_pam, cellarea, shape_c, scale_c);
//     postC.row(s) = trans(c);
//     if(s == 0){ Rprintf("sampled c\n"); }
//     
//     
//     if ( (s+1) % 10 == 0 ) {
//       Rprintf("Generated %d samples...\n", s+1);
//     }
//   }
//   
//   
//   return Rcpp::List::create(Rcpp::Named("postKappa") = postKappa,
//                             Rcpp::Named("postBeta") = postBeta,
//                             Rcpp::Named("postC") = postC,
//                             Rcpp::Named("postPsi") = postPsi,
//                             Rcpp::Named("Accprob") = accprob,
//                             Rcpp::Named("sigma2") = sigma2,
//                             Rcpp::Named("adapIter") = adapIter,
//                             Rcpp::Named("COVbeta") = COVbeta);
//   
// }






// [[Rcpp::export]]
List cppf_fitZoop_etat(
    int niter, mat distmat,
    List Xzoop, vec logYobl, vec logYsur,
    mat indobl, mat indsur,
    double alpha0tilde, vec alpha, mat eta, double lam0, double lam1,
    double tau2, double kappa_eta, double phi_eta, double sig2_obl, double sig2_sur,
    double cellarea,
    vec sigma2, mat COValpha,
    bool updateCOV, int adaptInterval, double adaptFactorExponent, vec adapIter){

  double negativeInf = -std::numeric_limits<float>::infinity();;

  indobl = indobl-1, indsur = indsur-1;

  int nG = distmat.n_rows, n_al = alpha.size(), ntzoop = eta.n_cols, n_sigma2 = sigma2.size();

  mat postAlpha0tilde = zeros(niter, 1), postAlpha = zeros(niter, n_al), postLam = zeros(niter, 2);
  mat postSig2 = zeros(niter, 1), postKappa = zeros(niter, 1);
  mat postEtat = zeros(niter, nG);
  List postEta(ntzoop);
  for(int i = 0; i < ntzoop; i ++){
    postEta[i] = zeros(niter, nG);
  }
  Rprintf("set up matrices for posterior samples\n");

  mat accprob = zeros(niter, n_sigma2);
  vec rhat = zeros(n_sigma2), gamma1 = zeros(n_sigma2), gamma2 = zeros(n_sigma2), dummyaccprob;
  double c0 = 1, c1 = adaptFactorExponent, ropt = 0.234, logprob;
  mat cholCOValpha = trans( chol( sigma2[0] * ( COValpha + 0.0000000001 * diagmat(ones(n_al)) ) ) );
  Rprintf("set up adaptive proposals\n");

  vec newalpha, newlogZobl, newlogZsur;
  mat newlogZ;

  int count;
  vec xvals;
  double llprev, llthred, llnew, thetamin, thetamax, theta;
  vec prioretat;
  mat neweta;
  bool accept;

  mat Reta = exp(-distmat / phi_eta);
  mat invReta = inv(Reta), cholReta = chol(Reta);
  rowvec oneinvReta = trans(ones(nG)) * invReta;
  double oneinvRetaone = sum( oneinvReta * ones(nG) );
  Rprintf("set up GP covariance matrices\n");

  mat logZ = cppf_logZ_a0t_etat_ja(Xzoop, alpha, eta);
  vec logZobl = cppf_logZobl(logZ, indobl);
  vec logZsur = cppf_logZsur(logZ, indsur);
  Rprintf("set up initial values\n");


  // Start MCMC
  for(int s = 0; s < niter; s++) {

    // update alpha
    if( updateCOV ){
      if( (s+1 >= adaptInterval) && (s+1 - (adaptInterval * trunc((s+1) / adaptInterval)) == 0) ){
        dummyaccprob = accprob.col(0);
        rhat[0] = sum( dummyaccprob.rows(s+1-adaptInterval, s-1) ) / (adaptInterval-1);
        gamma1[0] = 1 / pow(adapIter[0], c1);
        gamma2[0] = c0 * gamma1[0];
        sigma2[0] = exp( log(sigma2[0]) + gamma2[0] * (rhat[0] - ropt) );

        COValpha = COValpha + gamma1[0] * ( cov( postAlpha.rows(s+1-adaptInterval, s-1) ) - COValpha );
        cholCOValpha = trans( chol( sigma2[0] * ( COValpha + 0.0000000001 * diagmat(ones(n_al)) ) ) );

        adapIter[0] = adapIter[0] + 1;
      }
    }

    newalpha = alpha + cholCOValpha * randn(n_al);

    newlogZ = cppf_logZ_a0t_etat_ja(Xzoop, newalpha, eta);
    newlogZobl = cppf_logZobl(newlogZ, indobl);
    newlogZsur = cppf_logZsur(newlogZ, indsur);

    logprob = cppf_loglik_obl(logYobl, newlogZobl, sig2_obl) -
      cppf_loglik_obl(logYobl, logZobl, sig2_obl) +
      cppf_loglik_sur(logYsur, newlogZsur, lam0, lam1, sig2_sur) -
      cppf_loglik_sur(logYsur, logZsur, lam0, lam1, sig2_sur) +
      MVN_logh(newalpha, zeros(n_al), diagmat(ones(n_al))/100) -
      MVN_logh(alpha, zeros(n_al), diagmat(ones(n_al))/100);

    if (log(randu()) < logprob) {

      alpha = newalpha;
      logZ = newlogZ;
      logZobl = newlogZobl;
      logZsur = newlogZsur;

      accprob(s, 0) = 1;
    }
    postAlpha.row(s) = trans(alpha);
    if(s == 0){ Rprintf("sampled alpha\n"); }


    // Gibbs update for alpha0tilde
    alpha0tilde = cppf_gibbs_alpha0tilde(alpha.rows(0, ntzoop-1), tau2);
    postAlpha0tilde[s] = alpha0tilde;
    if(s == 0){ Rprintf("sampled alpha0tilde\n"); }


    // update eta via Elliptical slice sampling
    for(int t = 0; t < ntzoop; t ++){
      neweta = eta;

      llprev = cppf_loglik_obl(logYobl, logZobl, sig2_obl) +
        cppf_loglik_sur(logYsur, logZsur, lam0, lam1, sig2_sur);
      prioretat = sqrt(kappa_eta) * trans(cholReta) * randn(nG);

      thetamin = 0;
      thetamax = 2 * M_PI;
      theta = thetamin + randu() * (thetamax - thetamin);
      thetamin = theta - 2 * M_PI;
      thetamax = theta;

      llthred = llprev + log(randu());
      accept = false;
      count = 0;

      while(accept == false){
        count = count + 1;

        neweta.col(t) = eta.col(t) * cos(theta) + prioretat * sin(theta);
        neweta.col(t) = neweta.col(t) - mean(neweta.col(t));

        newlogZ = cppf_logZ_a0t_etat_ja(Xzoop, alpha, neweta);
        newlogZobl = cppf_logZobl(newlogZ, indobl);
        newlogZsur = cppf_logZsur(newlogZ, indsur);

        llnew = cppf_loglik_obl(logYobl, newlogZobl, sig2_obl) +
          cppf_loglik_sur(logYsur, newlogZsur, lam0, lam1, sig2_sur);

        if(llnew > llthred){
          llprev = llnew;
          accept = true;
        } else {
          if(theta < 0){ thetamin = theta; } else { thetamax = theta; }
          theta = thetamin + randu() * (thetamax - thetamin);
          if ( (count) % 1000 == 0 ) {
            Rprintf("ESS for eta: %d iterations...\n", count);
          }
          if (count > 20000) stop("ESS is not converging...\n");
        }
      }
      eta = neweta;
      logZ = newlogZ;
      logZobl = newlogZobl;
      logZsur = newlogZsur;

      postEtat = as<mat>(postEta[t]);
      postEtat.row(s) = trans(eta.col(t));
      postEta[t] = postEtat;
    }
    if(s == 0){ Rprintf("sampled eta\n"); }


    // update kappa_eta
    kappa_eta = cppf_gibbs_kappa_etat(eta, invReta, 2, 2);
    postKappa(s,0) = kappa_eta;
    if(s == 0){ Rprintf("sampled kappa_eta\n"); }


    // update lam0
    lam0 = cppf_gibbs_lam0(logYsur, logZsur, lam1, sig2_sur, 100);
    postLam(s, 0) = lam0;
    if(s == 0){ Rprintf("sampled lam0\n"); }


    // update lam1
    lam1 = cppf_gibbs_lam1(logYsur, logZsur, lam0, sig2_sur, 100);
    postLam(s, 1) = lam1;
    if(s == 0){ Rprintf("sampled lam1\n"); }


    // update sig2_sur
    sig2_sur = cppf_gibbs_sig2_sur(logYsur, logZsur, lam0, lam1, 2, 2);
    postSig2(s, 0) = sig2_sur;
    if(s == 0){ Rprintf("sampled sig2_sur\n"); }


    if ( (s+1) % 10 == 0 ) {
      Rprintf("Generated %d samples...\n", s+1);
    }
  }


  return Rcpp::List::create(Rcpp::Named("postAlpha0tilde") = postAlpha0tilde,
                            Rcpp::Named("postAlpha") = postAlpha,
                            Rcpp::Named("postEta") = postEta,
                            Rcpp::Named("postLam") = postLam,
                            Rcpp::Named("postKappa") = postKappa,
                            Rcpp::Named("postSig2") = postSig2,
                            Rcpp::Named("Accprob") = accprob,
                            Rcpp::Named("sigma2") = sigma2,
                            Rcpp::Named("adapIter") = adapIter,
                            Rcpp::Named("COValpha") = COValpha);

}



// // ============================================================================-
// // uniform prior for c ----
// // ============================================================================-
// 
// // uc: uniform prior for c
// // [[Rcpp::export]]
// List cppf_fitJoint_cZ_a0t_etat_b0t_psi_ct_ja_uc(
//     int niter, mat distmat,
//     List Xzoop, vec logYobl, vec logYsur,
//     List Xwhale, mat Sdist, List Wpam, vec dur_dist,
//     vec f_dist_obs, List f_dist, List p_pam,
//     mat indobl, mat indsur, vec indtwhale, mat inddist, vec indmdist, vec indmwhale,
//     double alpha0tilde, vec alpha, mat eta, double lam0, double lam1,
//     double tau2, double kappa_eta, double phi_eta, double sig2_obl, double sig2_sur,
//     vec beta, vec psi, vec pii, vec c, double lb_c,
//     double kappa_psi, double phi_psi,
//     double cellarea,
//     vec sigma2, mat COValpha, mat COVbeta, mat COVc,
//     bool updateCOV, int adaptInterval, double adaptFactorExponent, vec adapIter){
//   
//   double negativeInf = -std::numeric_limits<float>::infinity();;
//   
//   indobl = indobl-1, indsur = indsur-1, indtwhale = indtwhale-1;
//   inddist = inddist-1, indmdist = indmdist-1, indmwhale = indmwhale-1;
//   
//   int nG = distmat.n_rows, n_al = alpha.size(), n_be = beta.size(), ntzoop = Xzoop.size(), ntwhale = Xwhale.size(), n_sigma2 = sigma2.size();
//   
//   mat postAlpha0tilde = zeros(niter, 1), postAlpha = zeros(niter, n_al), postLam = zeros(niter, 2);
//   mat postSig2 = zeros(niter, 1), postKappa = zeros(niter, 2);
//   mat postBeta = zeros(niter, n_be), postC = zeros(niter, c.size());
//   mat postPsi(niter, nG), postEtat = zeros(niter, nG);
//   List postEta(ntzoop);
//   for(int i = 0; i < ntzoop; i ++){
//     postEta[i] = zeros(niter, nG);
//   }
//   Rprintf("set up matrices for posterior samples\n");
//   
//   mat accprob = zeros(niter, n_sigma2);
//   vec rhat = zeros(n_sigma2), gamma1 = zeros(n_sigma2), gamma2 = zeros(n_sigma2), dummyaccprob;
//   double c0 = 1, c1 = adaptFactorExponent, ropt = 0.234, logprob;
//   mat cholCOValpha = trans( chol( sigma2[0] * ( COValpha + 0.0000000001 * diagmat(ones(n_al)) ) ) );
//   mat cholCOVbeta = trans( chol( sigma2[1] * ( COVbeta + 0.0000000001 * diagmat(ones(n_be)) ) ) );
//   mat cholCOVc = trans( chol( sigma2[2] * ( COVc + 0.0000000001 * diagmat(ones(ntwhale)) ) ) );
//   Rprintf("set up adaptive proposals\n");
//   
//   bool cWithinRange;
//   vec newalpha, newbeta, newlogZobl, newlogZsur, newc;
//   mat newlogZ, newloglam, newlam;
//   vec dummyvec;
//   
//   int count;
//   vec xvals;
//   double llprev, llthred, llnew, thetamin, thetamax, theta;
//   vec prioretat, priorpsi, newpsi;
//   mat neweta;
//   bool accept;
//   
//   mat Reta = exp(-distmat / phi_eta), Rpsi = exp(-distmat / phi_psi);
//   mat invReta = inv(Reta), invRpsi= inv(Rpsi), cholReta = chol(Reta), cholRpsi= chol(Rpsi);
//   rowvec oneinvReta = trans(ones(nG)) * invReta, oneinvRpsi = trans(ones(nG)) * invRpsi;
//   double oneinvRetaone = sum( oneinvReta * ones(nG) ), oneinvRpsione = sum( oneinvRpsi * ones(nG) );
//   Rprintf("set up GP covariance matrices\n");
//   
//   mat logZ = cppf_logZ_a0t_etat_ja(Xzoop, alpha, eta);
//   vec logZobl = cppf_logZobl(logZ, indobl);
//   vec logZsur = cppf_logZsur(logZ, indsur);
//   mat loglam = cppf_loglam_cZ_b0t_psi(Xwhale, logZ, beta, psi, indtwhale);
//   mat lam = exp(loglam);
//   Rprintf("set up initial values\n");
//   
//   
//   // Start MCMC
//   for(int s = 0; s < niter; s++) {
//     
//     // update alpha
//     if( updateCOV ){
//       if( (s+1 >= adaptInterval) && (s+1 - (adaptInterval * trunc((s+1) / adaptInterval)) == 0) ){
//         dummyaccprob = accprob.col(0);
//         rhat[0] = sum( dummyaccprob.rows(s+1-adaptInterval, s-1) ) / (adaptInterval-1);
//         gamma1[0] = 1 / pow(adapIter[0], c1);
//         gamma2[0] = c0 * gamma1[0];
//         sigma2[0] = exp( log(sigma2[0]) + gamma2[0] * (rhat[0] - ropt) );
//         
//         COValpha = COValpha + gamma1[0] * ( cov( postAlpha.rows(s+1-adaptInterval, s-1) ) - COValpha );
//         cholCOValpha = trans( chol( sigma2[0] * ( COValpha + 0.0000000001 * diagmat(ones(n_al)) ) ) );
//         
//         adapIter[0] = adapIter[0] + 1;
//       }
//     }
//     
//     newalpha = alpha + cholCOValpha * randn(n_al);
//     
//     newlogZ = cppf_logZ_a0t_etat_ja(Xzoop, newalpha, eta);
//     newlogZobl = cppf_logZobl(newlogZ, indobl);
//     newlogZsur = cppf_logZsur(newlogZ, indsur);
//     newloglam = cppf_loglam_cZ_b0t_psi(Xwhale, newlogZ, beta, psi, indtwhale);
//     newlam = exp(newloglam);
//     
//     logprob = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) -
//       cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
//       cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) -
//       cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
//       cppf_loglik_obl(logYobl, newlogZobl, sig2_obl) -
//       cppf_loglik_obl(logYobl, logZobl, sig2_obl) +
//       cppf_loglik_sur(logYsur, newlogZsur, lam0, lam1, sig2_sur) -
//       cppf_loglik_sur(logYsur, logZsur, lam0, lam1, sig2_sur) +
//       MVN_logh(newalpha, zeros(n_al), diagmat(ones(n_al))/100) -
//       MVN_logh(alpha, zeros(n_al), diagmat(ones(n_al))/100);
//     
//     if (log(randu()) < logprob) {
//       
//       alpha = newalpha;
//       logZ = newlogZ;
//       logZobl = newlogZobl;
//       logZsur = newlogZsur;
//       loglam = newloglam;
//       lam = newlam;
//       
//       accprob(s, 0) = 1;
//     }
//     postAlpha.row(s) = trans(alpha);
//     if(s == 0){ Rprintf("sampled alpha\n"); }
//     
//     
//     // Gibbs update for alpha0tilde
//     alpha0tilde = cppf_gibbs_alpha0tilde(alpha.rows(0, ntzoop-1), tau2);
//     postAlpha0tilde[s] = alpha0tilde;
//     if(s == 0){ Rprintf("sampled alpha0tilde\n"); }
//     
//     
//     // update eta via Elliptical slice sampling
//     for(int t = 0; t < ntzoop; t ++){
//       neweta = eta;
//       
//       llprev = cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
//         cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
//         cppf_loglik_obl(logYobl, logZobl, sig2_obl) +
//         cppf_loglik_sur(logYsur, logZsur, lam0, lam1, sig2_sur);
//       prioretat = sqrt(kappa_eta) * trans(cholReta) * randn(nG);
//       
//       thetamin = 0;
//       thetamax = 2 * M_PI;
//       theta = thetamin + randu() * (thetamax - thetamin);
//       thetamin = theta - 2 * M_PI;
//       thetamax = theta;
//       
//       llthred = llprev + log(randu());
//       accept = false;
//       count = 0;
//       
//       while(accept == false){
//         count = count + 1;
//         
//         neweta.col(t) = eta.col(t) * cos(theta) + prioretat * sin(theta);
//         neweta.col(t) = neweta.col(t) - mean(neweta.col(t));
//         
//         newlogZ = cppf_logZ_a0t_etat_ja(Xzoop, alpha, neweta);
//         newlogZobl = cppf_logZobl(newlogZ, indobl);
//         newlogZsur = cppf_logZsur(newlogZ, indsur);
//         newloglam = cppf_loglam_cZ_b0t_psi(Xwhale, newlogZ, beta, psi, indtwhale);
//         newlam = exp(newloglam);
//         
//         llnew = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
//           cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) +
//           cppf_loglik_obl(logYobl, newlogZobl, sig2_obl) +
//           cppf_loglik_sur(logYsur, newlogZsur, lam0, lam1, sig2_sur);
//         
//         if(llnew > llthred){
//           llprev = llnew;
//           accept = true;
//         } else {
//           if(theta < 0){ thetamin = theta; } else { thetamax = theta; }
//           theta = thetamin + randu() * (thetamax - thetamin);
//           if ( (count) % 1000 == 0 ) {
//             Rprintf("ESS for eta: %d iterations...\n", count);
//           }
//           if (count > 20000) stop("ESS is not converging...\n");
//         }
//       }
//       eta = neweta;
//       logZ = newlogZ;
//       logZobl = newlogZobl;
//       logZsur = newlogZsur;
//       loglam = newloglam;
//       lam = newlam;
//       
//       postEtat = as<mat>(postEta[t]);
//       postEtat.row(s) = trans(eta.col(t));
//       postEta[t] = postEtat;
//     }
//     if(s == 0){ Rprintf("sampled eta\n"); }
//     
//     
//     // update kappa_eta
//     kappa_eta = cppf_gibbs_kappa_etat(eta, invReta, 2, 2);
//     postKappa(s,0) = kappa_eta;
//     if(s == 0){ Rprintf("sampled kappa_eta\n"); }
//     
//     
//     // update lam0
//     lam0 = cppf_gibbs_lam0(logYsur, logZsur, lam1, sig2_sur, 100);
//     postLam(s, 0) = lam0;
//     if(s == 0){ Rprintf("sampled lam0\n"); }
//     
//     
//     // update lam1
//     lam1 = cppf_gibbs_lam1(logYsur, logZsur, lam0, sig2_sur, 100);
//     postLam(s, 1) = lam1;
//     if(s == 0){ Rprintf("sampled lam1\n"); }
//     
//     
//     // update sig2_sur
//     sig2_sur = cppf_gibbs_sig2_sur(logYsur, logZsur, lam0, lam1, 2, 2);
//     postSig2(s, 0) = sig2_sur;
//     if(s == 0){ Rprintf("sampled sig2_sur\n"); }
//     
//     
//     // update beta
//     if( updateCOV ){
//       if( (s+1 >= adaptInterval) && (s+1 - (adaptInterval * trunc((s+1) / adaptInterval)) == 0) ){
//         dummyaccprob = accprob.col(1);
//         rhat[1] = sum( dummyaccprob.rows(s+1-adaptInterval, s-1) ) / (adaptInterval-1);
//         gamma1[1] = 1 / pow(adapIter[1], c1);
//         gamma2[1] = c0 * gamma1[1];
//         sigma2[1] = exp( log(sigma2[1]) + gamma2[1] * (rhat[1] - ropt) );
//         
//         COVbeta = COVbeta + gamma1[1] * ( cov( postBeta.rows(s+1-adaptInterval, s-1)) - COVbeta );
//         cholCOVbeta = trans( chol( sigma2[1] * ( COVbeta + 0.0000000001 * diagmat(ones(n_be)) ) ) );
//         
//         adapIter[1] = adapIter[1] + 1;
//       }
//     }
//     
//     newbeta = beta + cholCOVbeta * randn(n_be);
//     
//     newloglam = cppf_loglam_cZ_b0t_psi(Xwhale, logZ, newbeta, psi, indtwhale);
//     newlam = exp(newloglam);
//     
//     logprob = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) -
//       cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
//       cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) -
//       cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
//       MVN_logh(newbeta, zeros(n_be), diagmat(ones(n_be))/100) -
//       MVN_logh(beta, zeros(n_be), diagmat(ones(n_be))/100);
//     
//     if (log(randu()) < logprob) {
//       
//       beta = newbeta;
//       loglam = newloglam;
//       lam = newlam;
//       
//       accprob(s, 1) = 1;
//     }
//     
//     postBeta.row(s) = trans(beta);
//     if(s == 0){ Rprintf("sampled beta\n"); }
//     
//     
//     // update psi via elliptical slice sampling
//     llprev = cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
//       cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea);
//     priorpsi = sqrt(kappa_psi) * trans(cholRpsi) * randn(nG);
//     
//     thetamin = 0;
//     thetamax = 2 * M_PI;
//     theta = thetamin + randu() * (thetamax - thetamin);
//     thetamin = theta - 2 * M_PI;
//     thetamax = theta;
//     
//     llthred = llprev + log(randu());
//     accept = false;
//     count = 0;
//     
//     while(accept == false){
//       count = count + 1;
//       
//       newpsi = psi * cos(theta) + priorpsi * sin(theta);
//       newpsi = newpsi - mean(newpsi);
//       
//       newloglam = cppf_loglam_cZ_b0t_psi(Xwhale, logZ, beta, newpsi, indtwhale);
//       newlam = exp(newloglam);
//       
//       llnew = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
//         cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea);
//       
//       if(llnew > llthred){
//         llprev = llnew;
//         accept = true;
//       } else {
//         if(theta < 0){ thetamin = theta; } else { thetamax = theta; }
//         theta = thetamin + randu() * (thetamax - thetamin);
//         if ( (count) % 1000 == 0 ) {
//           Rprintf("ESS for psi: %d iterations...\n", count);
//         }
//         if (count > 20000) stop("ESS is not converging...\n");
//       }
//     }
//     psi = newpsi;
//     loglam = newloglam;
//     lam = newlam;
//     
//     postPsi.row(s) = trans(psi);
//     if(s == 0){ Rprintf("sampled psi\n"); }
//     
//     
//     // update kappa_psi
//     kappa_psi = cppf_gibbs_kappa_psi(psi, invRpsi, 2, 2);
//     postKappa(s,1) = kappa_psi;
//     if(s == 0){ Rprintf("sampled kappa_psi\n"); }
//     
//     
//     // update c
//     if( updateCOV ){
//       if( (s+1 >= adaptInterval) && (s+1 - (adaptInterval * trunc((s+1) / adaptInterval)) == 0) ){
//         dummyaccprob = accprob.col(2);
//         rhat[2] = sum( dummyaccprob.rows(s+1-adaptInterval, s-1) ) / (adaptInterval-1);
//         gamma1[2] = 1 / pow(adapIter[2], c1);
//         gamma2[2] = c0 * gamma1[2];
//         sigma2[2] = exp( log(sigma2[2]) + gamma2[2] * (rhat[2] - ropt) );
//         
//         COVc = COVc + gamma1[2] * ( cov( postC.rows(s+1-adaptInterval, s-1)) - COVc );
//         cholCOVc = trans( chol( sigma2[2] * ( COVc + 0.0000000001 * diagmat(ones(ntwhale)) ) ) );
//         
//         adapIter[2] = adapIter[2] + 1;
//       }
//     }
//     
//     newc = c + cholCOVc * randn(ntwhale);
//     
//     cWithinRange = true;
//     for(int t = 0; t < ntwhale; t ++){
//       if( (newc[t] < lb_c) ){
//         cWithinRange = false;
//       }
//     }
//     
//     if(cWithinRange){
//       logprob = cppf_loglik_pam_ct(Wpam, lam, dur_dist, newc, p_pam, cellarea) -
//         cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea);
//     } else {
//       logprob = negativeInf;
//     }
//     
//     if (log(randu()) < logprob) {
//       c = newc;
//       accprob(s, 2) = 1;
//     }
//     postC.row(s) = trans(c);
//     if(s == 0){ Rprintf("sampled c\n"); }
//     
//     
//     
//     if ( (s+1) % 10 == 0 ) {
//       Rprintf("Generated %d samples...\n", s+1);
//     }
//   }
//   
//   
//   return Rcpp::List::create(Rcpp::Named("postAlpha0tilde") = postAlpha0tilde,
//                             Rcpp::Named("postAlpha") = postAlpha,
//                             Rcpp::Named("postEta") = postEta,
//                             Rcpp::Named("postLam") = postLam,
//                             Rcpp::Named("postKappa") = postKappa,
//                             Rcpp::Named("postSig2") = postSig2,
//                             Rcpp::Named("postBeta") = postBeta,
//                             Rcpp::Named("postC") = postC,
//                             Rcpp::Named("postPsi") = postPsi,
//                             Rcpp::Named("Accprob") = accprob,
//                             Rcpp::Named("sigma2") = sigma2,
//                             Rcpp::Named("adapIter") = adapIter,
//                             Rcpp::Named("COValpha") = COValpha,
//                             Rcpp::Named("COVbeta") = COVbeta,
//                             Rcpp::Named("COVc") = COVc);
//   
// }
// 
// 
// 
// // ============================================================================-
// // old ----
// // ============================================================================-
// // 
// // // [[Rcpp::export]]
// // List cppf_fitJoint_eta_b0t_psi_ct_fixC(
// //     int niter, mat distmat,
// //     List Xzoop, vec logYobl, vec logYsur,
// //     List Xwhale, mat Sdist, List Wpam, vec dur_dist,
// //     vec f_dist_obs, List f_dist, List p_pam,
// //     mat indobl, mat indsur, vec indtwhale, mat inddist, vec indmdist, vec indmwhale,
// //     vec alpha, vec eta, double lam0, double lam1,
// //     double kappa_eta, double phi_eta, double sig2_obl, double sig2_sur,
// //     vec beta, vec psi, vec pii, vec c,
// //     double kappa_psi, double phi_psi,
// //     double cellarea,
// //     vec sigma2, mat COValpha, mat COVbeta,
// //     bool updateCOV, int adaptInterval, double adaptFactorExponent, vec adapIter){
// //   
// //   double negativeInf = -std::numeric_limits<float>::infinity();;
// //   
// //   indobl = indobl-1, indsur = indsur-1, indtwhale = indtwhale-1;
// //   inddist = inddist-1, indmdist = indmdist-1, indmwhale = indmwhale-1;
// //   
// //   int nG = distmat.n_rows, n_al = alpha.size(), n_be = beta.size(), ntzoop = eta.n_cols, n_sigma2 = sigma2.size();
// //   
// //   mat postAlpha = zeros(niter, n_al), postLam = zeros(niter, 2);
// //   mat postSig2 = zeros(niter, 1), postKappa = zeros(niter, 2);
// //   mat postBeta = zeros(niter, n_be), postEta = zeros(niter, nG);
// //   mat postPsi = zeros(niter, nG);
// //   Rprintf("set up matrices for posterior samples\n");
// //   
// //   mat accprob = zeros(niter, n_sigma2);
// //   vec rhat = zeros(n_sigma2), gamma1 = zeros(n_sigma2), gamma2 = zeros(n_sigma2), dummyaccprob;
// //   double c0 = 1, c1 = adaptFactorExponent, ropt = 0.234, logprob;
// //   mat cholCOValpha = trans( chol( sigma2[0] * ( COValpha + 0.0000000001 * diagmat(ones(n_al)) ) ) );
// //   mat cholCOVbeta = trans( chol( sigma2[1] * ( COVbeta + 0.0000000001 * diagmat(ones(n_be)) ) ) );
// //   Rprintf("set up adaptive proposals\n");
// //   
// //   vec newalpha, newbeta, newlogZobl, newlogZsur;
// //   mat newlogZ, newloglam, newlam;
// //   vec dummyvec;
// //   
// //   int count;
// //   vec xvals;
// //   double llprev, llthred, llnew, thetamin, thetamax, theta;
// //   vec prioreta, neweta, priorpsi, newpsi;
// //   bool accept;
// //   
// //   mat Reta = exp(-distmat / phi_eta), Rpsi = exp(-distmat / phi_psi);
// //   mat invReta = inv(Reta), invRpsi= inv(Rpsi), cholReta = chol(Reta), cholRpsi= chol(Rpsi);
// //   rowvec oneinvReta = trans(ones(nG)) * invReta, oneinvRpsi = trans(ones(nG)) * invRpsi;
// //   double oneinvRetaone = sum( oneinvReta * ones(nG) ), oneinvRpsione = sum( oneinvRpsi * ones(nG) );
// //   Rprintf("set up GP covariance matrices\n");
// //   
// //   mat logZ = cppf_logZ_eta(Xzoop, alpha, eta);
// //   vec logZobl = cppf_logZobl(logZ, indobl);
// //   vec logZsur = cppf_logZsur(logZ, indsur);
// //   mat loglam = cppf_loglam_b0t_psi(Xwhale, logZ, beta, psi, indtwhale);
// //   mat lam = exp(loglam);
// //   vec m_al = zeros(n_al), diag_al = ones(n_al)/100;
// //   Rprintf("set up initial values\n");
// //   
// //   
// //   // Start MCMC
// //   for(int s = 0; s < niter; s++) {
// //     
// //     // update alpha
// //     if( updateCOV ){
// //       if( (s+1 >= adaptInterval) && (s+1 - (adaptInterval * trunc((s+1) / adaptInterval)) == 0) ){
// //         dummyaccprob = accprob.col(0);
// //         rhat[0] = sum( dummyaccprob.rows(s+1-adaptInterval, s-1) ) / (adaptInterval-1);
// //         gamma1[0] = 1 / pow(adapIter[0], c1);
// //         gamma2[0] = c0 * gamma1[0];
// //         sigma2[0] = exp( log(sigma2[0]) + gamma2[0] * (rhat[0] - ropt) );
// //         
// //         COValpha = COValpha + gamma1[0] * ( cov( postAlpha.rows(s+1-adaptInterval, s-1) ) - COValpha );
// //         cholCOValpha = trans( chol( sigma2[0] * ( COValpha + 0.0000000001 * diagmat(ones(n_al)) ) ) );
// //         
// //         adapIter[0] = adapIter[0] + 1;
// //       }
// //     }
// //     
// //     newalpha = alpha + cholCOValpha * randn(n_al);
// //     
// //     newlogZ = cppf_logZ_eta(Xzoop, newalpha, eta);
// //     newlogZobl = cppf_logZobl(newlogZ, indobl);
// //     newlogZsur = cppf_logZsur(newlogZ, indsur);
// //     newloglam = cppf_loglam_b0t_psi(Xwhale, newlogZ, beta, psi, indtwhale);
// //     newlam = exp(newloglam);
// //     
// //     logprob = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) -
// //       cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //       cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) -
// //       cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
// //       cppf_loglik_obl(logYobl, newlogZobl, sig2_obl) -
// //       cppf_loglik_obl(logYobl, logZobl, sig2_obl) +
// //       cppf_loglik_sur(logYsur, newlogZsur, lam0, lam1, sig2_sur) -
// //       cppf_loglik_sur(logYsur, logZsur, lam0, lam1, sig2_sur) +
// //       MVN_logh(newalpha, zeros(n_al), diagmat(ones(n_al))/100) -
// //       MVN_logh(alpha, zeros(n_al), diagmat(ones(n_al))/100);
// //     
// //     if (log(randu()) < logprob) {
// //       
// //       alpha = newalpha;
// //       logZ = newlogZ;
// //       logZobl = newlogZobl;
// //       logZsur = newlogZsur;
// //       loglam = newloglam;
// //       lam = newlam;
// //       
// //       accprob(s, 0) = 1;
// //     }
// //     postAlpha.row(s) = trans(alpha);
// //     if(s == 0){ Rprintf("sampled alpha\n"); }
// //     
// //     
// //     
// //     // update eta via Elliptical slice sampling
// //     llprev = cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //       cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
// //       cppf_loglik_obl(logYobl, logZobl, sig2_obl) +
// //       cppf_loglik_sur(logYsur, logZsur, lam0, lam1, sig2_sur);
// //     prioreta = sqrt(kappa_eta) * trans(cholReta) * randn(nG);
// //     
// //     thetamin = 0;
// //     thetamax = 2 * M_PI;
// //     theta = thetamin + randu() * (thetamax - thetamin);
// //     thetamin = theta - 2 * M_PI;
// //     thetamax = theta;
// //     
// //     llthred = llprev + log(randu());
// //     accept = false;
// //     count = 0;
// //     
// //     while(accept == false){
// //       count = count + 1;
// //       
// //       neweta = eta * cos(theta) + prioreta * sin(theta);
// //       neweta = neweta - mean(neweta);
// //       
// //       newlogZ = cppf_logZ_eta(Xzoop, alpha, neweta);
// //       newlogZobl = cppf_logZobl(newlogZ, indobl);
// //       newlogZsur = cppf_logZsur(newlogZ, indsur);
// //       newloglam = cppf_loglam_b0t_psi(Xwhale, newlogZ, beta, psi, indtwhale);
// //       newlam = exp(newloglam);
// //       
// //       llnew = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //         cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) +
// //         cppf_loglik_obl(logYobl, newlogZobl, sig2_obl) +
// //         cppf_loglik_sur(logYsur, newlogZsur, lam0, lam1, sig2_sur);
// //       
// //       if(llnew > llthred){
// //         llprev = llnew;
// //         accept = true;
// //       } else {
// //         if(theta < 0){ thetamin = theta; } else { thetamax = theta; }
// //         theta = thetamin + randu() * (thetamax - thetamin);
// //         if ( (count) % 1000 == 0 ) {
// //           Rprintf("ESS for eta: %d iterations...\n", count);
// //         }
// //         if (count > 20000) stop("ESS is not converging...\n");
// //       }
// //     }
// //     eta = neweta;
// //     logZ = newlogZ;
// //     logZobl = newlogZobl;
// //     logZsur = newlogZsur;
// //     loglam = newloglam;
// //     lam = newlam;
// //     
// //     postEta.row(s) = trans(eta);
// //     if(s == 0){ Rprintf("sampled eta\n"); }
// //     
// //     
// //     // update kappa_eta
// //     kappa_eta = cppf_gibbs_kappa_eta(eta, invReta, 2, 2);
// //     postKappa(s,0) = kappa_eta;
// //     if(s == 0){ Rprintf("sampled kappa_eta\n"); }
// //     
// //     
// //     // update lam0
// //     lam0 = cppf_gibbs_lam0(logYsur, logZsur, lam1, sig2_sur, 100);
// //     postLam(s, 0) = lam0;
// //     if(s == 0){ Rprintf("sampled lam0\n"); }
// //     
// //     
// //     // update lam1
// //     lam1 = cppf_gibbs_lam1(logYsur, logZsur, lam0, sig2_sur, 100);
// //     postLam(s, 1) = lam1;
// //     if(s == 0){ Rprintf("sampled lam1\n"); }
// //     
// //     
// //     // update sig2_sur
// //     sig2_sur = cppf_gibbs_sig2_sur(logYsur, logZsur, lam0, lam1, 2, 2);
// //     postSig2(s, 0) = sig2_sur;
// //     if(s == 0){ Rprintf("sampled sig2_sur\n"); }
// //     
// //     
// //     // update beta
// //     if( updateCOV ){
// //       if( (s+1 >= adaptInterval) && (s+1 - (adaptInterval * trunc((s+1) / adaptInterval)) == 0) ){
// //         dummyaccprob = accprob.col(1);
// //         rhat[1] = sum( dummyaccprob.rows(s+1-adaptInterval, s-1) ) / (adaptInterval-1);
// //         gamma1[1] = 1 / pow(adapIter[1], c1);
// //         gamma2[1] = c0 * gamma1[1];
// //         sigma2[1] = exp( log(sigma2[1]) + gamma2[1] * (rhat[1] - ropt) );
// //         
// //         COVbeta = COVbeta + gamma1[1] * ( cov( postBeta.rows(s+1-adaptInterval, s-1)) - COVbeta );
// //         cholCOVbeta = trans( chol( sigma2[1] * ( COVbeta + 0.0000000001 * diagmat(ones(n_be)) ) ) );
// //         
// //         adapIter[1] = adapIter[1] + 1;
// //       }
// //     }
// //     
// //     newbeta = beta + cholCOVbeta * randn(n_be);
// //     
// //     newloglam = cppf_loglam_b0t_psi(Xwhale, logZ, newbeta, psi, indtwhale);
// //     newlam = exp(newloglam);
// //     
// //     logprob = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) -
// //       cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //       cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) -
// //       cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
// //       MVN_logh(newbeta, zeros(n_be), diagmat(ones(n_be))/100) -
// //       MVN_logh(beta, zeros(n_be), diagmat(ones(n_be))/100);
// //     
// //     if (log(randu()) < logprob) {
// //       
// //       beta = newbeta;
// //       loglam = newloglam;
// //       lam = newlam;
// //       
// //       accprob(s, 1) = 1;
// //     }
// //     
// //     postBeta.row(s) = trans(beta);
// //     if(s == 0){ Rprintf("sampled beta\n"); }
// //     
// //     
// //     
// //     // update psi via elliptical slice sampling
// //     llprev = cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //       cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea);
// //     priorpsi = sqrt(kappa_psi) * trans(cholRpsi) * randn(nG);
// //     
// //     thetamin = 0;
// //     thetamax = 2 * M_PI;
// //     theta = thetamin + randu() * (thetamax - thetamin);
// //     thetamin = theta - 2 * M_PI;
// //     thetamax = theta;
// //     
// //     llthred = llprev + log(randu());
// //     accept = false;
// //     count = 0;
// //     
// //     while(accept == false){
// //       count = count + 1;
// //       
// //       newpsi = psi * cos(theta) + priorpsi * sin(theta);
// //       newpsi = newpsi - mean(newpsi);
// //       
// //       newloglam = cppf_loglam_b0t_psi(Xwhale, logZ, beta, newpsi, indtwhale);
// //       newlam = exp(newloglam);
// //       
// //       llnew = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //         cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea);
// //       
// //       if(llnew > llthred){
// //         llprev = llnew;
// //         accept = true;
// //       } else {
// //         if(theta < 0){ thetamin = theta; } else { thetamax = theta; }
// //         theta = thetamin + randu() * (thetamax - thetamin);
// //         if ( (count) % 1000 == 0 ) {
// //           Rprintf("ESS for psi: %d iterations...\n", count);
// //         }
// //         if (count > 20000) stop("ESS is not converging...\n");
// //       }
// //     }
// //     psi = newpsi;
// //     loglam = newloglam;
// //     lam = newlam;
// //     
// //     postPsi.row(s) = trans(psi);
// //     if(s == 0){ Rprintf("sampled psi\n"); }
// //     
// //     
// //     // update kappa_psi
// //     kappa_psi = cppf_gibbs_kappa_psi(psi, invRpsi, 2, 2);
// //     postKappa(s,1) = kappa_psi;
// //     if(s == 0){ Rprintf("sampled kappa_psi\n"); }
// //     
// //     
// //     if ( (s+1) % 10 == 0 ) {
// //       Rprintf("Generated %d samples...\n", s+1);
// //     }
// //   }
// //   
// //   
// //   return Rcpp::List::create(Rcpp::Named("postAlpha") = postAlpha,
// //                             Rcpp::Named("postEta") = postEta,
// //                             Rcpp::Named("postLam") = postLam,
// //                             Rcpp::Named("postKappa") = postKappa,
// //                             Rcpp::Named("postSig2") = postSig2,
// //                             Rcpp::Named("postBeta") = postBeta,
// //                             Rcpp::Named("postPsi") = postPsi,
// //                             Rcpp::Named("Accprob") = accprob,
// //                             Rcpp::Named("sigma2") = sigma2,
// //                             Rcpp::Named("adapIter") = adapIter,
// //                             Rcpp::Named("COValpha") = COValpha,
// //                             Rcpp::Named("COVbeta") = COVbeta);
// //   
// // }
// // 
// // 
// // 
// // // [[Rcpp::export]]
// // List cppf_fitJoint_eta_b0t_psit_ct_fixC(
// //     int niter, mat distmat,
// //     List Xzoop, vec logYobl, vec logYsur,
// //     List Xwhale, mat Sdist, List Wpam, vec dur_dist,
// //     vec f_dist_obs, List f_dist, List p_pam,
// //     mat indobl, mat indsur, vec indtwhale, mat inddist, vec indmdist, vec indmwhale,
// //     vec alpha, vec eta, double lam0, double lam1,
// //     double kappa_eta, double phi_eta, double sig2_obl, double sig2_sur,
// //     vec beta, mat psi, vec pii, vec c,
// //     double kappa_psi, double phi_psi,
// //     double cellarea,
// //     vec sigma2, mat COValpha, mat COVbeta,
// //     bool updateCOV, int adaptInterval, double adaptFactorExponent, vec adapIter){
// //   
// //   double negativeInf = -std::numeric_limits<float>::infinity();;
// //   
// //   indobl = indobl-1, indsur = indsur-1, indtwhale = indtwhale-1;
// //   inddist = inddist-1, indmdist = indmdist-1, indmwhale = indmwhale-1;
// //   
// //   int nG = distmat.n_rows, n_al = alpha.size(), n_be = beta.size(), ntzoop = eta.n_cols, ntwhale = psi.n_cols, n_sigma2 = sigma2.size();
// //   
// //   mat postAlpha = zeros(niter, n_al), postLam = zeros(niter, 2);
// //   mat postSig2 = zeros(niter, 1), postKappa = zeros(niter, 2);
// //   mat postBeta = zeros(niter, n_be), postEta = zeros(niter, nG);
// //   mat postPsit = zeros(niter, nG);
// //   List postPsi(ntwhale);
// //   for(int i = 0; i < ntwhale; i ++){ postPsi[i] = zeros(niter, nG); }
// //   Rprintf("set up matrices for posterior samples\n");
// //   
// //   mat accprob = zeros(niter, n_sigma2);
// //   vec rhat = zeros(n_sigma2), gamma1 = zeros(n_sigma2), gamma2 = zeros(n_sigma2), dummyaccprob;
// //   double c0 = 1, c1 = adaptFactorExponent, ropt = 0.234, logprob;
// //   mat cholCOValpha = trans( chol( sigma2[0] * ( COValpha + 0.0000000001 * diagmat(ones(n_al)) ) ) );
// //   mat cholCOVbeta = trans( chol( sigma2[1] * ( COVbeta + 0.0000000001 * diagmat(ones(n_be)) ) ) );
// //   Rprintf("set up adaptive proposals\n");
// //   
// //   vec newalpha, newbeta, newlogZobl, newlogZsur;
// //   mat newlogZ, newloglam, newlam;
// //   vec dummyvec;
// //   
// //   int count;
// //   vec xvals;
// //   double llprev, llthred, llnew, thetamin, thetamax, theta;
// //   vec prioreta, neweta, priorpsit;
// //   mat newpsi;
// //   bool accept;
// //   
// //   mat Reta = exp(-distmat / phi_eta), Rpsi = exp(-distmat / phi_psi);
// //   mat invReta = inv(Reta), invRpsi= inv(Rpsi), cholReta = chol(Reta), cholRpsi= chol(Rpsi);
// //   rowvec oneinvReta = trans(ones(nG)) * invReta, oneinvRpsi = trans(ones(nG)) * invRpsi;
// //   double oneinvRetaone = sum( oneinvReta * ones(nG) ), oneinvRpsione = sum( oneinvRpsi * ones(nG) );
// //   Rprintf("set up GP covariance matrices\n");
// //   
// //   mat logZ = cppf_logZ_eta(Xzoop, alpha, eta);
// //   vec logZobl = cppf_logZobl(logZ, indobl);
// //   vec logZsur = cppf_logZsur(logZ, indsur);
// //   mat loglam = cppf_loglam_b0t_psit(Xwhale, logZ, beta, psi, indtwhale);
// //   mat lam = exp(loglam);
// //   vec m_al = zeros(n_al), diag_al = ones(n_al)/100;
// //   Rprintf("set up initial values\n");
// //   
// //   
// //   // Start MCMC
// //   for(int s = 0; s < niter; s++) {
// //     
// //     // update alpha
// //     if( updateCOV ){
// //       if( (s+1 >= adaptInterval) && (s+1 - (adaptInterval * trunc((s+1) / adaptInterval)) == 0) ){
// //         dummyaccprob = accprob.col(0);
// //         rhat[0] = sum( dummyaccprob.rows(s+1-adaptInterval, s-1) ) / (adaptInterval-1);
// //         gamma1[0] = 1 / pow(adapIter[0], c1);
// //         gamma2[0] = c0 * gamma1[0];
// //         sigma2[0] = exp( log(sigma2[0]) + gamma2[0] * (rhat[0] - ropt) );
// //         
// //         COValpha = COValpha + gamma1[0] * ( cov( postAlpha.rows(s+1-adaptInterval, s-1) ) - COValpha );
// //         cholCOValpha = trans( chol( sigma2[0] * ( COValpha + 0.0000000001 * diagmat(ones(n_al)) ) ) );
// //         
// //         adapIter[0] = adapIter[0] + 1;
// //       }
// //     }
// //     
// //     newalpha = alpha + cholCOValpha * randn(n_al);
// //     
// //     newlogZ = cppf_logZ_eta(Xzoop, newalpha, eta);
// //     newlogZobl = cppf_logZobl(newlogZ, indobl);
// //     newlogZsur = cppf_logZsur(newlogZ, indsur);
// //     newloglam = cppf_loglam_b0t_psit(Xwhale, newlogZ, beta, psi, indtwhale);
// //     newlam = exp(newloglam);
// //     
// //     logprob = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) -
// //       cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //       cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) -
// //       cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
// //       cppf_loglik_obl(logYobl, newlogZobl, sig2_obl) -
// //       cppf_loglik_obl(logYobl, logZobl, sig2_obl) +
// //       cppf_loglik_sur(logYsur, newlogZsur, lam0, lam1, sig2_sur) -
// //       cppf_loglik_sur(logYsur, logZsur, lam0, lam1, sig2_sur) +
// //       MVN_logh(newalpha, zeros(n_al), diagmat(ones(n_al))/100) -
// //       MVN_logh(alpha, zeros(n_al), diagmat(ones(n_al))/100);
// //     
// //     if (log(randu()) < logprob) {
// //       
// //       alpha = newalpha;
// //       logZ = newlogZ;
// //       logZobl = newlogZobl;
// //       logZsur = newlogZsur;
// //       loglam = newloglam;
// //       lam = newlam;
// //       
// //       accprob(s, 0) = 1;
// //     }
// //     postAlpha.row(s) = trans(alpha);
// //     if(s == 0){ Rprintf("sampled alpha\n"); }
// //     
// //     
// //     
// //     // update eta via Elliptical slice sampling
// //     llprev = cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //       cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
// //       cppf_loglik_obl(logYobl, logZobl, sig2_obl) +
// //       cppf_loglik_sur(logYsur, logZsur, lam0, lam1, sig2_sur);
// //     prioreta = sqrt(kappa_eta) * trans(cholReta) * randn(nG);
// //     
// //     thetamin = 0;
// //     thetamax = 2 * M_PI;
// //     theta = thetamin + randu() * (thetamax - thetamin);
// //     thetamin = theta - 2 * M_PI;
// //     thetamax = theta;
// //     
// //     llthred = llprev + log(randu());
// //     accept = false;
// //     count = 0;
// //     
// //     while(accept == false){
// //       count = count + 1;
// //       
// //       neweta = eta * cos(theta) + prioreta * sin(theta);
// //       neweta = neweta - mean(neweta);
// //       
// //       newlogZ = cppf_logZ_eta(Xzoop, alpha, neweta);
// //       newlogZobl = cppf_logZobl(newlogZ, indobl);
// //       newlogZsur = cppf_logZsur(newlogZ, indsur);
// //       newloglam = cppf_loglam_b0t_psit(Xwhale, newlogZ, beta, psi, indtwhale);
// //       newlam = exp(newloglam);
// //       
// //       llnew = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //         cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) +
// //         cppf_loglik_obl(logYobl, newlogZobl, sig2_obl) +
// //         cppf_loglik_sur(logYsur, newlogZsur, lam0, lam1, sig2_sur);
// //       
// //       if(llnew > llthred){
// //         llprev = llnew;
// //         accept = true;
// //       } else {
// //         if(theta < 0){ thetamin = theta; } else { thetamax = theta; }
// //         theta = thetamin + randu() * (thetamax - thetamin);
// //         if ( (count) % 1000 == 0 ) {
// //           Rprintf("ESS for eta: %d iterations...\n", count);
// //         }
// //         if (count > 20000) stop("ESS is not converging...\n");
// //       }
// //     }
// //     eta = neweta;
// //     logZ = newlogZ;
// //     logZobl = newlogZobl;
// //     logZsur = newlogZsur;
// //     loglam = newloglam;
// //     lam = newlam;
// //     
// //     postEta.row(s) = trans(eta);
// //     if(s == 0){ Rprintf("sampled eta\n"); }
// //     
// //     
// //     // update kappa_eta
// //     kappa_eta = cppf_gibbs_kappa_eta(eta, invReta, 2, 2);
// //     postKappa(s,0) = kappa_eta;
// //     if(s == 0){ Rprintf("sampled kappa_eta\n"); }
// //     
// //     
// //     // update lam0
// //     lam0 = cppf_gibbs_lam0(logYsur, logZsur, lam1, sig2_sur, 100);
// //     postLam(s, 0) = lam0;
// //     if(s == 0){ Rprintf("sampled lam0\n"); }
// //     
// //     
// //     // update lam1
// //     lam1 = cppf_gibbs_lam1(logYsur, logZsur, lam0, sig2_sur, 100);
// //     postLam(s, 1) = lam1;
// //     if(s == 0){ Rprintf("sampled lam1\n"); }
// //     
// //     
// //     // update sig2_sur
// //     sig2_sur = cppf_gibbs_sig2_sur(logYsur, logZsur, lam0, lam1, 2, 2);
// //     postSig2(s, 0) = sig2_sur;
// //     if(s == 0){ Rprintf("sampled sig2_sur\n"); }
// //     
// //     
// //     // update beta
// //     if( updateCOV ){
// //       if( (s+1 >= adaptInterval) && (s+1 - (adaptInterval * trunc((s+1) / adaptInterval)) == 0) ){
// //         dummyaccprob = accprob.col(1);
// //         rhat[1] = sum( dummyaccprob.rows(s+1-adaptInterval, s-1) ) / (adaptInterval-1);
// //         gamma1[1] = 1 / pow(adapIter[1], c1);
// //         gamma2[1] = c0 * gamma1[1];
// //         sigma2[1] = exp( log(sigma2[1]) + gamma2[1] * (rhat[1] - ropt) );
// //         
// //         COVbeta = COVbeta + gamma1[1] * ( cov( postBeta.rows(s+1-adaptInterval, s-1)) - COVbeta );
// //         cholCOVbeta = trans( chol( sigma2[1] * ( COVbeta + 0.0000000001 * diagmat(ones(n_be)) ) ) );
// //         
// //         adapIter[1] = adapIter[1] + 1;
// //       }
// //     }
// //     
// //     newbeta = beta + cholCOVbeta * randn(n_be);
// //     
// //     newloglam = cppf_loglam_b0t_psit(Xwhale, logZ, newbeta, psi, indtwhale);
// //     newlam = exp(newloglam);
// //     
// //     logprob = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) -
// //       cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //       cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) -
// //       cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
// //       MVN_logh(newbeta, zeros(n_be), diagmat(ones(n_be))/100) -
// //       MVN_logh(beta, zeros(n_be), diagmat(ones(n_be))/100);
// //     
// //     if (log(randu()) < logprob) {
// //       
// //       beta = newbeta;
// //       loglam = newloglam;
// //       lam = newlam;
// //       
// //       accprob(s, 1) = 1;
// //     }
// //     
// //     postBeta.row(s) = trans(beta);
// //     if(s == 0){ Rprintf("sampled beta\n"); }
// //     
// //     
// //     
// //     // update psi via Elliptical slice sampling
// //     for(int t = 0; t < ntwhale; t ++){
// //       newpsi = psi;
// //       
// //       llprev = cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //         cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea);
// //       priorpsit = sqrt(kappa_psi) * trans(cholRpsi) * randn(nG);
// //       
// //       thetamin = 0;
// //       thetamax = 2 * M_PI;
// //       theta = thetamin + randu() * (thetamax - thetamin);
// //       thetamin = theta - 2 * M_PI;
// //       thetamax = theta;
// //       
// //       llthred = llprev + log(randu());
// //       accept = false;
// //       count = 0;
// //       
// //       while(accept == false){
// //         count = count + 1;
// //         
// //         newpsi.col(t) = psi.col(t) * cos(theta) + priorpsit * sin(theta);
// //         newpsi.col(t) = newpsi.col(t) - mean(newpsi.col(t));
// //         
// //         newloglam = cppf_loglam_b0t_psit(Xwhale, logZ, beta, newpsi, indtwhale);
// //         newlam = exp(newloglam);
// //         
// //         llnew = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //           cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea);
// //         
// //         if(llnew > llthred){
// //           llprev = llnew;
// //           accept = true;
// //         } else {
// //           if(theta < 0){ thetamin = theta; } else { thetamax = theta; }
// //           theta = thetamin + randu() * (thetamax - thetamin);
// //           if ( (count) % 1000 == 0 ) {
// //             Rprintf("ESS for psi: %d iterations...\n", count);
// //           }
// //           if (count > 20000) stop("ESS is not converging...\n");
// //         }
// //       }
// //       psi = newpsi;
// //       loglam = newloglam;
// //       lam = newlam;
// //       
// //       postPsit = as<mat>(postPsi[t]);
// //       postPsit.row(s) = trans(psi.col(t));
// //       postPsi[t] = postPsit;
// //     }
// //     if(s == 0){ Rprintf("sampled psi\n"); }
// //     
// //     
// //     // update kappa_psi
// //     kappa_psi = cppf_gibbs_kappa_psit(psi, invRpsi, 2, 2);
// //     postKappa(s,1) = kappa_psi;
// //     if(s == 0){ Rprintf("sampled kappa_psi\n"); }
// //     
// //     
// //     if ( (s+1) % 10 == 0 ) {
// //       Rprintf("Generated %d samples...\n", s+1);
// //     }
// //   }
// //   
// //   
// //   return Rcpp::List::create(Rcpp::Named("postAlpha") = postAlpha,
// //                             Rcpp::Named("postEta") = postEta,
// //                             Rcpp::Named("postLam") = postLam,
// //                             Rcpp::Named("postKappa") = postKappa,
// //                             Rcpp::Named("postSig2") = postSig2,
// //                             Rcpp::Named("postBeta") = postBeta,
// //                             Rcpp::Named("postPsi") = postPsi,
// //                             Rcpp::Named("Accprob") = accprob,
// //                             Rcpp::Named("sigma2") = sigma2,
// //                             Rcpp::Named("adapIter") = adapIter,
// //                             Rcpp::Named("COValpha") = COValpha,
// //                             Rcpp::Named("COVbeta") = COVbeta);
// //   
// // }
// // 
// // 
// // // [[Rcpp::export]]
// // List cppf_fitJoint_a0t_etat_b0t_psi_ct_ja_fixC(
// //     int niter, mat distmat,
// //     List Xzoop, vec logYobl, vec logYsur,
// //     List Xwhale, mat Sdist, List Wpam, vec dur_dist,
// //     vec f_dist_obs, List f_dist, List p_pam,
// //     mat indobl, mat indsur, vec indtwhale, mat inddist, vec indmdist, vec indmwhale,
// //     double alpha0tilde, vec alpha, mat eta, double lam0, double lam1,
// //     double tau2, double kappa_eta, double phi_eta, double sig2_obl, double sig2_sur,
// //     vec beta, vec psi, vec pii, vec c, 
// //     double kappa_psi, double phi_psi,
// //     double cellarea,
// //     vec sigma2, mat COValpha, mat COVbeta,
// //     bool updateCOV, int adaptInterval, double adaptFactorExponent, vec adapIter){
// //   
// //   double negativeInf = -std::numeric_limits<float>::infinity();;
// //   
// //   indobl = indobl-1, indsur = indsur-1, indtwhale = indtwhale-1;
// //   inddist = inddist-1, indmdist = indmdist-1, indmwhale = indmwhale-1;
// //   
// //   int nG = distmat.n_rows, n_al = alpha.size(), n_be = beta.size(), ntzoop = eta.n_cols, n_sigma2 = sigma2.size();
// //   
// //   mat postAlpha0tilde = zeros(niter, 1), postAlpha = zeros(niter, n_al), postLam = zeros(niter, 2);
// //   mat postSig2 = zeros(niter, 1), postKappa = zeros(niter, 2);
// //   mat postBeta = zeros(niter, n_be);
// //   mat postPsi(niter, nG), postEtat = zeros(niter, nG);
// //   List postEta(ntzoop);
// //   for(int i = 0; i < ntzoop; i ++){
// //     postEta[i] = zeros(niter, nG);
// //   }
// //   Rprintf("set up matrices for posterior samples\n");
// //   
// //   mat accprob = zeros(niter, n_sigma2);
// //   vec rhat = zeros(n_sigma2), gamma1 = zeros(n_sigma2), gamma2 = zeros(n_sigma2), dummyaccprob;
// //   double c0 = 1, c1 = adaptFactorExponent, ropt = 0.234, logprob;
// //   mat cholCOValpha = trans( chol( sigma2[0] * ( COValpha + 0.0000000001 * diagmat(ones(n_al)) ) ) );
// //   mat cholCOVbeta = trans( chol( sigma2[1] * ( COVbeta + 0.0000000001 * diagmat(ones(n_be)) ) ) );
// //   Rprintf("set up adaptive proposals\n");
// //   
// //   vec newalpha, newbeta, newlogZobl, newlogZsur;
// //   mat newlogZ, newloglam, newlam;
// //   vec dummyvec;
// //   
// //   int count;
// //   vec xvals;
// //   double llprev, llthred, llnew, thetamin, thetamax, theta;
// //   vec prioretat, priorpsi, newpsi;
// //   mat neweta;
// //   bool accept;
// //   
// //   mat Reta = exp(-distmat / phi_eta), Rpsi = exp(-distmat / phi_psi);
// //   mat invReta = inv(Reta), invRpsi= inv(Rpsi), cholReta = chol(Reta), cholRpsi= chol(Rpsi);
// //   rowvec oneinvReta = trans(ones(nG)) * invReta, oneinvRpsi = trans(ones(nG)) * invRpsi;
// //   double oneinvRetaone = sum( oneinvReta * ones(nG) ), oneinvRpsione = sum( oneinvRpsi * ones(nG) );
// //   Rprintf("set up GP covariance matrices\n");
// //   
// //   mat logZ = cppf_logZ_a0t_etat_ja(Xzoop, alpha, eta);
// //   vec logZobl = cppf_logZobl(logZ, indobl);
// //   vec logZsur = cppf_logZsur(logZ, indsur);
// //   mat loglam = cppf_loglam_b0t_psi(Xwhale, logZ, beta, psi, indtwhale);
// //   mat lam = exp(loglam);
// //   Rprintf("set up initial values\n");
// //   
// //   
// //   // Start MCMC
// //   for(int s = 0; s < niter; s++) {
// //     
// //     // update alpha
// //     if( updateCOV ){
// //       if( (s+1 >= adaptInterval) && (s+1 - (adaptInterval * trunc((s+1) / adaptInterval)) == 0) ){
// //         dummyaccprob = accprob.col(0);
// //         rhat[0] = sum( dummyaccprob.rows(s+1-adaptInterval, s-1) ) / (adaptInterval-1);
// //         gamma1[0] = 1 / pow(adapIter[0], c1);
// //         gamma2[0] = c0 * gamma1[0];
// //         sigma2[0] = exp( log(sigma2[0]) + gamma2[0] * (rhat[0] - ropt) );
// //         
// //         COValpha = COValpha + gamma1[0] * ( cov( postAlpha.rows(s+1-adaptInterval, s-1) ) - COValpha );
// //         cholCOValpha = trans( chol( sigma2[0] * ( COValpha + 0.0000000001 * diagmat(ones(n_al)) ) ) );
// //         
// //         adapIter[0] = adapIter[0] + 1;
// //       }
// //     }
// //     
// //     newalpha = alpha + cholCOValpha * randn(n_al);
// //     
// //     newlogZ = cppf_logZ_a0t_etat_ja(Xzoop, newalpha, eta);
// //     newlogZobl = cppf_logZobl(newlogZ, indobl);
// //     newlogZsur = cppf_logZsur(newlogZ, indsur);
// //     newloglam = cppf_loglam_b0t_psi(Xwhale, newlogZ, beta, psi, indtwhale);
// //     newlam = exp(newloglam);
// //     
// //     logprob = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) -
// //       cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //       cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) -
// //       cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
// //       cppf_loglik_obl(logYobl, newlogZobl, sig2_obl) -
// //       cppf_loglik_obl(logYobl, logZobl, sig2_obl) +
// //       cppf_loglik_sur(logYsur, newlogZsur, lam0, lam1, sig2_sur) -
// //       cppf_loglik_sur(logYsur, logZsur, lam0, lam1, sig2_sur) +
// //       MVN_logh(newalpha, zeros(n_al), diagmat(ones(n_al))/100) -
// //       MVN_logh(alpha, zeros(n_al), diagmat(ones(n_al))/100);
// //     
// //     if (log(randu()) < logprob) {
// //       
// //       alpha = newalpha;
// //       logZ = newlogZ;
// //       logZobl = newlogZobl;
// //       logZsur = newlogZsur;
// //       loglam = newloglam;
// //       lam = newlam;
// //       
// //       accprob(s, 0) = 1;
// //     }
// //     postAlpha.row(s) = trans(alpha);
// //     if(s == 0){ Rprintf("sampled alpha\n"); }
// //     
// //     
// //     // Gibbs update for alpha0tilde
// //     alpha0tilde = cppf_gibbs_alpha0tilde(alpha.rows(0, ntzoop-1), tau2);
// //     postAlpha0tilde[s] = alpha0tilde;
// //     if(s == 0){ Rprintf("sampled alpha0tilde\n"); }
// //     
// //     
// //     // update eta via Elliptical slice sampling
// //     for(int t = 0; t < ntzoop; t ++){
// //       neweta = eta;
// //       
// //       llprev = cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //         cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
// //         cppf_loglik_obl(logYobl, logZobl, sig2_obl) +
// //         cppf_loglik_sur(logYsur, logZsur, lam0, lam1, sig2_sur);
// //       prioretat = sqrt(kappa_eta) * trans(cholReta) * randn(nG);
// //       
// //       thetamin = 0;
// //       thetamax = 2 * M_PI;
// //       theta = thetamin + randu() * (thetamax - thetamin);
// //       thetamin = theta - 2 * M_PI;
// //       thetamax = theta;
// //       
// //       llthred = llprev + log(randu());
// //       accept = false;
// //       count = 0;
// //       
// //       while(accept == false){
// //         count = count + 1;
// //         
// //         neweta.col(t) = eta.col(t) * cos(theta) + prioretat * sin(theta);
// //         neweta.col(t) = neweta.col(t) - mean(neweta.col(t));
// //         
// //         newlogZ = cppf_logZ_a0t_etat_ja(Xzoop, alpha, neweta);
// //         newlogZobl = cppf_logZobl(newlogZ, indobl);
// //         newlogZsur = cppf_logZsur(newlogZ, indsur);
// //         newloglam = cppf_loglam_b0t_psi(Xwhale, newlogZ, beta, psi, indtwhale);
// //         newlam = exp(newloglam);
// //         
// //         llnew = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //           cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) +
// //           cppf_loglik_obl(logYobl, newlogZobl, sig2_obl) +
// //           cppf_loglik_sur(logYsur, newlogZsur, lam0, lam1, sig2_sur);
// //         
// //         if(llnew > llthred){
// //           llprev = llnew;
// //           accept = true;
// //         } else {
// //           if(theta < 0){ thetamin = theta; } else { thetamax = theta; }
// //           theta = thetamin + randu() * (thetamax - thetamin);
// //           if ( (count) % 1000 == 0 ) {
// //             Rprintf("ESS for eta: %d iterations...\n", count);
// //           }
// //           if (count > 20000) stop("ESS is not converging...\n");
// //         }
// //       }
// //       eta = neweta;
// //       logZ = newlogZ;
// //       logZobl = newlogZobl;
// //       logZsur = newlogZsur;
// //       loglam = newloglam;
// //       lam = newlam;
// //       
// //       postEtat = as<mat>(postEta[t]);
// //       postEtat.row(s) = trans(eta.col(t));
// //       postEta[t] = postEtat;
// //     }
// //     if(s == 0){ Rprintf("sampled eta\n"); }
// //     
// //     
// //     // update kappa_eta
// //     kappa_eta = cppf_gibbs_kappa_etat(eta, invReta, 2, 2);
// //     postKappa(s,0) = kappa_eta;
// //     if(s == 0){ Rprintf("sampled kappa_eta\n"); }
// //     
// //     
// //     // update lam0
// //     lam0 = cppf_gibbs_lam0(logYsur, logZsur, lam1, sig2_sur, 100);
// //     postLam(s, 0) = lam0;
// //     if(s == 0){ Rprintf("sampled lam0\n"); }
// //     
// //     
// //     // update lam1
// //     lam1 = cppf_gibbs_lam1(logYsur, logZsur, lam0, sig2_sur, 100);
// //     postLam(s, 1) = lam1;
// //     if(s == 0){ Rprintf("sampled lam1\n"); }
// //     
// //     
// //     // update sig2_sur
// //     sig2_sur = cppf_gibbs_sig2_sur(logYsur, logZsur, lam0, lam1, 2, 2);
// //     postSig2(s, 0) = sig2_sur;
// //     if(s == 0){ Rprintf("sampled sig2_sur\n"); }
// //     
// //     
// //     // update beta
// //     if( updateCOV ){
// //       if( (s+1 >= adaptInterval) && (s+1 - (adaptInterval * trunc((s+1) / adaptInterval)) == 0) ){
// //         dummyaccprob = accprob.col(1);
// //         rhat[1] = sum( dummyaccprob.rows(s+1-adaptInterval, s-1) ) / (adaptInterval-1);
// //         gamma1[1] = 1 / pow(adapIter[1], c1);
// //         gamma2[1] = c0 * gamma1[1];
// //         sigma2[1] = exp( log(sigma2[1]) + gamma2[1] * (rhat[1] - ropt) );
// //         
// //         COVbeta = COVbeta + gamma1[1] * ( cov( postBeta.rows(s+1-adaptInterval, s-1)) - COVbeta );
// //         cholCOVbeta = trans( chol( sigma2[1] * ( COVbeta + 0.0000000001 * diagmat(ones(n_be)) ) ) );
// //         
// //         adapIter[1] = adapIter[1] + 1;
// //       }
// //     }
// //     
// //     newbeta = beta + cholCOVbeta * randn(n_be);
// //     
// //     newloglam = cppf_loglam_b0t_psi(Xwhale, logZ, newbeta, psi, indtwhale);
// //     newlam = exp(newloglam);
// //     
// //     logprob = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) -
// //       cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //       cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) -
// //       cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
// //       MVN_logh(newbeta, zeros(n_be), diagmat(ones(n_be))/100) -
// //       MVN_logh(beta, zeros(n_be), diagmat(ones(n_be))/100);
// //     
// //     if (log(randu()) < logprob) {
// //       
// //       beta = newbeta;
// //       loglam = newloglam;
// //       lam = newlam;
// //       
// //       accprob(s, 1) = 1;
// //     }
// //     
// //     postBeta.row(s) = trans(beta);
// //     if(s == 0){ Rprintf("sampled beta\n"); }
// //     
// //     
// //     // update psi via elliptical slice sampling
// //     llprev = cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //       cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea);
// //     priorpsi = sqrt(kappa_psi) * trans(cholRpsi) * randn(nG);
// //     
// //     thetamin = 0;
// //     thetamax = 2 * M_PI;
// //     theta = thetamin + randu() * (thetamax - thetamin);
// //     thetamin = theta - 2 * M_PI;
// //     thetamax = theta;
// //     
// //     llthred = llprev + log(randu());
// //     accept = false;
// //     count = 0;
// //     
// //     while(accept == false){
// //       count = count + 1;
// //       
// //       newpsi = psi * cos(theta) + priorpsi * sin(theta);
// //       newpsi = newpsi - mean(newpsi);
// //       
// //       newloglam = cppf_loglam_b0t_psi(Xwhale, logZ, beta, newpsi, indtwhale);
// //       newlam = exp(newloglam);
// //       
// //       llnew = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //         cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea);
// //       
// //       if(llnew > llthred){
// //         llprev = llnew;
// //         accept = true;
// //       } else {
// //         if(theta < 0){ thetamin = theta; } else { thetamax = theta; }
// //         theta = thetamin + randu() * (thetamax - thetamin);
// //         if ( (count) % 1000 == 0 ) {
// //           Rprintf("ESS for psi: %d iterations...\n", count);
// //         }
// //         if (count > 20000) stop("ESS is not converging...\n");
// //       }
// //     }
// //     psi = newpsi;
// //     loglam = newloglam;
// //     lam = newlam;
// //     
// //     postPsi.row(s) = trans(psi);
// //     if(s == 0){ Rprintf("sampled psi\n"); }
// //     
// //     
// //     // update kappa_psi
// //     kappa_psi = cppf_gibbs_kappa_psi(psi, invRpsi, 2, 2);
// //     postKappa(s,1) = kappa_psi;
// //     if(s == 0){ Rprintf("sampled kappa_psi\n"); }
// //     
// //     
// //     if ( (s+1) % 10 == 0 ) {
// //       Rprintf("Generated %d samples...\n", s+1);
// //     }
// //   }
// //   
// //   
// //   return Rcpp::List::create(Rcpp::Named("postAlpha0tilde") = postAlpha0tilde,
// //                             Rcpp::Named("postAlpha") = postAlpha,
// //                             Rcpp::Named("postEta") = postEta,
// //                             Rcpp::Named("postLam") = postLam,
// //                             Rcpp::Named("postKappa") = postKappa,
// //                             Rcpp::Named("postSig2") = postSig2,
// //                             Rcpp::Named("postBeta") = postBeta,
// //                             Rcpp::Named("postPsi") = postPsi,
// //                             Rcpp::Named("Accprob") = accprob,
// //                             Rcpp::Named("sigma2") = sigma2,
// //                             Rcpp::Named("adapIter") = adapIter,
// //                             Rcpp::Named("COValpha") = COValpha,
// //                             Rcpp::Named("COVbeta") = COVbeta);
// //   
// // }
// // 
// // 
// // 
// // // [[Rcpp::export]]
// // List cppf_fitJoint_a0t_etat_b0t_psit_ct_ja_fixC(
// //     int niter, mat distmat,
// //     List Xzoop, vec logYobl, vec logYsur,
// //     List Xwhale, mat Sdist, List Wpam, vec dur_dist,
// //     vec f_dist_obs, List f_dist, List p_pam,
// //     mat indobl, mat indsur, vec indtwhale, mat inddist, vec indmdist, vec indmwhale,
// //     double alpha0tilde, vec alpha, mat eta, double lam0, double lam1,
// //     double tau2, double kappa_eta, double phi_eta, double sig2_obl, double sig2_sur,
// //     vec beta, mat psi, vec pii, vec c,
// //     double kappa_psi, double phi_psi,
// //     double cellarea,
// //     vec sigma2, mat COValpha, mat COVbeta,
// //     bool updateCOV, int adaptInterval, double adaptFactorExponent, vec adapIter){
// //   
// //   double negativeInf = -std::numeric_limits<float>::infinity();;
// //   
// //   indobl = indobl-1, indsur = indsur-1, indtwhale = indtwhale-1;
// //   inddist = inddist-1, indmdist = indmdist-1, indmwhale = indmwhale-1;
// //   
// //   int nG = distmat.n_rows, n_al = alpha.size(), n_be = beta.size(), ntzoop = eta.n_cols, ntwhale = psi.n_cols, n_sigma2 = sigma2.size();
// //   
// //   mat postAlpha0tilde = zeros(niter, 1), postAlpha = zeros(niter, n_al), postLam = zeros(niter, 2);
// //   mat postSig2 = zeros(niter, 1), postKappa = zeros(niter, 2);
// //   mat postBeta = zeros(niter, n_be);
// //   mat postEtat = zeros(niter, nG), postPsit = zeros(niter, nG);
// //   List postEta(ntzoop), postPsi(ntwhale);
// //   for(int i = 0; i < ntzoop; i ++){ postEta[i] = zeros(niter, nG); }
// //   for(int i = 0; i < ntwhale; i ++){ postPsi[i] = zeros(niter, nG); }
// //   Rprintf("set up matrices for posterior samples\n");
// //   
// //   mat accprob = zeros(niter, n_sigma2);
// //   vec rhat = zeros(n_sigma2), gamma1 = zeros(n_sigma2), gamma2 = zeros(n_sigma2), dummyaccprob;
// //   double c0 = 1, c1 = adaptFactorExponent, ropt = 0.234, logprob;
// //   mat cholCOValpha = trans( chol( sigma2[0] * ( COValpha + 0.0000000001 * diagmat(ones(n_al)) ) ) );
// //   mat cholCOVbeta = trans( chol( sigma2[1] * ( COVbeta + 0.0000000001 * diagmat(ones(n_be)) ) ) );
// //   Rprintf("set up adaptive proposals\n");
// //   
// //   vec newalpha, newbeta, newlogZobl, newlogZsur;
// //   mat newlogZ, newloglam, newlam;
// //   vec dummyvec;
// //   
// //   int count;
// //   vec xvals;
// //   double llprev, llthred, llnew, thetamin, thetamax, theta;
// //   vec prioretat, priorpsit;
// //   mat neweta, newpsi;
// //   bool accept;
// //   
// //   mat Reta = exp(-distmat / phi_eta), Rpsi = exp(-distmat / phi_psi);
// //   mat invReta = inv(Reta), invRpsi= inv(Rpsi), cholReta = chol(Reta), cholRpsi= chol(Rpsi);
// //   rowvec oneinvReta = trans(ones(nG)) * invReta, oneinvRpsi = trans(ones(nG)) * invRpsi;
// //   double oneinvRetaone = sum( oneinvReta * ones(nG) ), oneinvRpsione = sum( oneinvRpsi * ones(nG) );
// //   Rprintf("set up GP covariance matrices\n");
// //   
// //   mat logZ = cppf_logZ_a0t_etat_ja(Xzoop, alpha, eta);
// //   vec logZobl = cppf_logZobl(logZ, indobl);
// //   vec logZsur = cppf_logZsur(logZ, indsur);
// //   mat loglam = cppf_loglam_b0t_psit(Xwhale, logZ, beta, psi, indtwhale);
// //   mat lam = exp(loglam);
// //   vec m_al = zeros(n_al), diag_al = ones(n_al)/100;
// //   Rprintf("set up initial values\n");
// //   
// //   
// //   // Start MCMC
// //   for(int s = 0; s < niter; s++) {
// //     
// //     // update alpha
// //     if( updateCOV ){
// //       if( (s+1 >= adaptInterval) && (s+1 - (adaptInterval * trunc((s+1) / adaptInterval)) == 0) ){
// //         dummyaccprob = accprob.col(0);
// //         rhat[0] = sum( dummyaccprob.rows(s+1-adaptInterval, s-1) ) / (adaptInterval-1);
// //         gamma1[0] = 1 / pow(adapIter[0], c1);
// //         gamma2[0] = c0 * gamma1[0];
// //         sigma2[0] = exp( log(sigma2[0]) + gamma2[0] * (rhat[0] - ropt) );
// //         
// //         COValpha = COValpha + gamma1[0] * ( cov( postAlpha.rows(s+1-adaptInterval, s-1) ) - COValpha );
// //         cholCOValpha = trans( chol( sigma2[0] * ( COValpha + 0.0000000001 * diagmat(ones(n_al)) ) ) );
// //         
// //         adapIter[0] = adapIter[0] + 1;
// //       }
// //     }
// //     
// //     newalpha = alpha + cholCOValpha * randn(n_al);
// //     
// //     newlogZ = cppf_logZ_a0t_etat_ja(Xzoop, newalpha, eta);
// //     newlogZobl = cppf_logZobl(newlogZ, indobl);
// //     newlogZsur = cppf_logZsur(newlogZ, indsur);
// //     newloglam = cppf_loglam_b0t_psit(Xwhale, newlogZ, beta, psi, indtwhale);
// //     newlam = exp(newloglam);
// //     
// //     for(int t = 0; t < ntzoop; t ++){
// //       m_al[t] = alpha0tilde;
// //       diag_al[t] = 1/tau2;
// //     }
// //     
// //     logprob = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) -
// //       cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //       cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) -
// //       cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
// //       cppf_loglik_obl(logYobl, newlogZobl, sig2_obl) -
// //       cppf_loglik_obl(logYobl, logZobl, sig2_obl) +
// //       cppf_loglik_sur(logYsur, newlogZsur, lam0, lam1, sig2_sur) -
// //       cppf_loglik_sur(logYsur, logZsur, lam0, lam1, sig2_sur) +
// //       MVN_logh(newalpha, m_al, diagmat(diag_al)) -
// //       MVN_logh(alpha, m_al, diagmat(diag_al));
// //     
// //     if (log(randu()) < logprob) {
// //       
// //       alpha = newalpha;
// //       logZ = newlogZ;
// //       logZobl = newlogZobl;
// //       logZsur = newlogZsur;
// //       loglam = newloglam;
// //       lam = newlam;
// //       
// //       accprob(s, 0) = 1;
// //     }
// //     postAlpha.row(s) = trans(alpha);
// //     if(s == 0){ Rprintf("sampled alpha\n"); }
// //     
// //     
// //     // Gibbs update for alpha0tilde
// //     alpha0tilde = cppf_gibbs_alpha0tilde(alpha.rows(0, ntzoop-1), tau2);
// //     postAlpha0tilde[s] = alpha0tilde;
// //     if(s == 0){ Rprintf("sampled alpha0tilde\n"); }
// //     
// //     
// //     
// //     // update eta via Elliptical slice sampling
// //     for(int t = 0; t < ntzoop; t ++){
// //       neweta = eta;
// //       
// //       llprev = cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //         cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
// //         cppf_loglik_obl(logYobl, logZobl, sig2_obl) +
// //         cppf_loglik_sur(logYsur, logZsur, lam0, lam1, sig2_sur);
// //       prioretat = sqrt(kappa_eta) * trans(cholReta) * randn(nG);
// //       
// //       thetamin = 0;
// //       thetamax = 2 * M_PI;
// //       theta = thetamin + randu() * (thetamax - thetamin);
// //       thetamin = theta - 2 * M_PI;
// //       thetamax = theta;
// //       
// //       llthred = llprev + log(randu());
// //       accept = false;
// //       count = 0;
// //       
// //       while(accept == false){
// //         count = count + 1;
// //         
// //         neweta.col(t) = eta.col(t) * cos(theta) + prioretat * sin(theta);
// //         neweta.col(t) = neweta.col(t) - mean(neweta.col(t));
// //         
// //         newlogZ = cppf_logZ_a0t_etat_ja(Xzoop, alpha, neweta);
// //         newlogZobl = cppf_logZobl(newlogZ, indobl);
// //         newlogZsur = cppf_logZsur(newlogZ, indsur);
// //         newloglam = cppf_loglam_b0t_psit(Xwhale, newlogZ, beta, psi, indtwhale);
// //         newlam = exp(newloglam);
// //         
// //         llnew = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //           cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) +
// //           cppf_loglik_obl(logYobl, newlogZobl, sig2_obl) +
// //           cppf_loglik_sur(logYsur, newlogZsur, lam0, lam1, sig2_sur);
// //         
// //         if(llnew > llthred){
// //           llprev = llnew;
// //           accept = true;
// //         } else {
// //           if(theta < 0){ thetamin = theta; } else { thetamax = theta; }
// //           theta = thetamin + randu() * (thetamax - thetamin);
// //           if ( (count) % 1000 == 0 ) {
// //             Rprintf("ESS for eta: %d iterations...\n", count);
// //           }
// //           if (count > 20000) stop("ESS is not converging...\n");
// //         }
// //       }
// //       eta = neweta;
// //       logZ = newlogZ;
// //       logZobl = newlogZobl;
// //       logZsur = newlogZsur;
// //       loglam = newloglam;
// //       lam = newlam;
// //       
// //       postEtat = as<mat>(postEta[t]);
// //       postEtat.row(s) = trans(eta.col(t));
// //       postEta[t] = postEtat;
// //     }
// //     if(s == 0){ Rprintf("sampled eta\n"); }
// //     
// //     
// //     // update kappa_eta
// //     kappa_eta = cppf_gibbs_kappa_etat(eta, invReta, 2, 2);
// //     postKappa(s,0) = kappa_eta;
// //     if(s == 0){ Rprintf("sampled kappa_eta\n"); }
// //     
// //     
// //     // update lam0
// //     lam0 = cppf_gibbs_lam0(logYsur, logZsur, lam1, sig2_sur, 100);
// //     postLam(s, 0) = lam0;
// //     if(s == 0){ Rprintf("sampled lam0\n"); }
// //     
// //     
// //     // update lam1
// //     lam1 = cppf_gibbs_lam1(logYsur, logZsur, lam0, sig2_sur, 100);
// //     postLam(s, 1) = lam1;
// //     if(s == 0){ Rprintf("sampled lam1\n"); }
// //     
// //     
// //     // update sig2_sur
// //     sig2_sur = cppf_gibbs_sig2_sur(logYsur, logZsur, lam0, lam1, 2, 2);
// //     postSig2(s, 0) = sig2_sur;
// //     if(s == 0){ Rprintf("sampled sig2_sur\n"); }
// //     
// //     
// //     // update beta
// //     if( updateCOV ){
// //       if( (s+1 >= adaptInterval) && (s+1 - (adaptInterval * trunc((s+1) / adaptInterval)) == 0) ){
// //         dummyaccprob = accprob.col(1);
// //         rhat[1] = sum( dummyaccprob.rows(s+1-adaptInterval, s-1) ) / (adaptInterval-1);
// //         gamma1[1] = 1 / pow(adapIter[1], c1);
// //         gamma2[1] = c0 * gamma1[1];
// //         sigma2[1] = exp( log(sigma2[1]) + gamma2[1] * (rhat[1] - ropt) );
// //         
// //         COVbeta = COVbeta + gamma1[1] * ( cov( postBeta.rows(s+1-adaptInterval, s-1)) - COVbeta );
// //         cholCOVbeta = trans( chol( sigma2[1] * ( COVbeta + 0.0000000001 * diagmat(ones(n_be)) ) ) );
// //         
// //         adapIter[1] = adapIter[1] + 1;
// //       }
// //     }
// //     
// //     newbeta = beta + cholCOVbeta * randn(n_be);
// //     
// //     newloglam = cppf_loglam_b0t_psit(Xwhale, logZ, newbeta, psi, indtwhale);
// //     newlam = exp(newloglam);
// //     
// //     logprob = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) -
// //       cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //       cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) -
// //       cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
// //       MVN_logh(newbeta, zeros(n_be), diagmat(ones(n_be))/100) -
// //       MVN_logh(beta, zeros(n_be), diagmat(ones(n_be))/100);
// //     
// //     if (log(randu()) < logprob) {
// //       
// //       beta = newbeta;
// //       loglam = newloglam;
// //       lam = newlam;
// //       
// //       accprob(s, 1) = 1;
// //     }
// //     
// //     postBeta.row(s) = trans(beta);
// //     if(s == 0){ Rprintf("sampled beta\n"); }
// //     
// //     
// //     
// //     // update psi via Elliptical slice sampling
// //     for(int t = 0; t < ntwhale; t ++){
// //       newpsi = psi;
// //       
// //       llprev = cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //         cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea);
// //       priorpsit = sqrt(kappa_psi) * trans(cholRpsi) * randn(nG);
// //       
// //       thetamin = 0;
// //       thetamax = 2 * M_PI;
// //       theta = thetamin + randu() * (thetamax - thetamin);
// //       thetamin = theta - 2 * M_PI;
// //       thetamax = theta;
// //       
// //       llthred = llprev + log(randu());
// //       accept = false;
// //       count = 0;
// //       
// //       while(accept == false){
// //         count = count + 1;
// //         
// //         newpsi.col(t) = psi.col(t) * cos(theta) + priorpsit * sin(theta);
// //         newpsi.col(t) = newpsi.col(t) - mean(newpsi.col(t));
// //         
// //         newloglam = cppf_loglam_b0t_psit(Xwhale, logZ, beta, newpsi, indtwhale);
// //         newlam = exp(newloglam);
// //         
// //         llnew = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //           cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea);
// //         
// //         if(llnew > llthred){
// //           llprev = llnew;
// //           accept = true;
// //         } else {
// //           if(theta < 0){ thetamin = theta; } else { thetamax = theta; }
// //           theta = thetamin + randu() * (thetamax - thetamin);
// //           if ( (count) % 1000 == 0 ) {
// //             Rprintf("ESS for psi: %d iterations...\n", count);
// //           }
// //           if (count > 20000) stop("ESS is not converging...\n");
// //         }
// //       }
// //       psi = newpsi;
// //       loglam = newloglam;
// //       lam = newlam;
// //       
// //       postPsit = as<mat>(postPsi[t]);
// //       postPsit.row(s) = trans(psi.col(t));
// //       postPsi[t] = postPsit;
// //     }
// //     if(s == 0){ Rprintf("sampled psi\n"); }
// //     
// //     
// //     // update kappa_psi
// //     kappa_psi = cppf_gibbs_kappa_psit(psi, invRpsi, 2, 2);
// //     postKappa(s,1) = kappa_psi;
// //     if(s == 0){ Rprintf("sampled kappa_psi\n"); }
// //     
// //     
// //     if ( (s+1) % 10 == 0 ) {
// //       Rprintf("Generated %d samples...\n", s+1);
// //     }
// //   }
// //   
// //   
// //   return Rcpp::List::create(Rcpp::Named("postAlpha0tilde") = postAlpha0tilde,
// //                             Rcpp::Named("postAlpha") = postAlpha,
// //                             Rcpp::Named("postEta") = postEta,
// //                             Rcpp::Named("postLam") = postLam,
// //                             Rcpp::Named("postKappa") = postKappa,
// //                             Rcpp::Named("postSig2") = postSig2,
// //                             Rcpp::Named("postBeta") = postBeta,
// //                             Rcpp::Named("postPsi") = postPsi,
// //                             Rcpp::Named("Accprob") = accprob,
// //                             Rcpp::Named("sigma2") = sigma2,
// //                             Rcpp::Named("adapIter") = adapIter,
// //                             Rcpp::Named("COValpha") = COValpha,
// //                             Rcpp::Named("COVbeta") = COVbeta);
// //   
// // }
// // 
// // 
// // 
// // 
// // 
// // 
// // 
// // 
// // // ============================================================================-
// // // centering zoop covariate ----
// // // ============================================================================-
// // 
// // 
// // // [[Rcpp::export]]
// // List cppf_fitJoint_cZ_eta_b0t_psi_ct_fixC(
// //     int niter, mat distmat,
// //     List Xzoop, vec logYobl, vec logYsur,
// //     List Xwhale, mat Sdist, List Wpam, vec dur_dist,
// //     vec f_dist_obs, List f_dist, List p_pam,
// //     mat indobl, mat indsur, vec indtwhale, mat inddist, vec indmdist, vec indmwhale,
// //     vec alpha, vec eta, double lam0, double lam1,
// //     double kappa_eta, double phi_eta, double sig2_obl, double sig2_sur,
// //     vec beta, vec psi, vec pii, vec c,
// //     double kappa_psi, double phi_psi,
// //     double cellarea,
// //     vec sigma2, mat COValpha, mat COVbeta,
// //     bool updateCOV, int adaptInterval, double adaptFactorExponent, vec adapIter){
// //   
// //   double negativeInf = -std::numeric_limits<float>::infinity();;
// //   
// //   indobl = indobl-1, indsur = indsur-1, indtwhale = indtwhale-1;
// //   inddist = inddist-1, indmdist = indmdist-1, indmwhale = indmwhale-1;
// //   
// //   int nG = distmat.n_rows, n_al = alpha.size(), n_be = beta.size(), ntzoop = eta.n_cols, n_sigma2 = sigma2.size();
// //   
// //   mat postAlpha = zeros(niter, n_al), postLam = zeros(niter, 2);
// //   mat postSig2 = zeros(niter, 1), postKappa = zeros(niter, 2);
// //   mat postBeta = zeros(niter, n_be), postEta = zeros(niter, nG);
// //   mat postPsi = zeros(niter, nG);
// //   Rprintf("set up matrices for posterior samples\n");
// //   
// //   mat accprob = zeros(niter, n_sigma2);
// //   vec rhat = zeros(n_sigma2), gamma1 = zeros(n_sigma2), gamma2 = zeros(n_sigma2), dummyaccprob;
// //   double c0 = 1, c1 = adaptFactorExponent, ropt = 0.234, logprob;
// //   mat cholCOValpha = trans( chol( sigma2[0] * ( COValpha + 0.0000000001 * diagmat(ones(n_al)) ) ) );
// //   mat cholCOVbeta = trans( chol( sigma2[1] * ( COVbeta + 0.0000000001 * diagmat(ones(n_be)) ) ) );
// //   Rprintf("set up adaptive proposals\n");
// //   
// //   vec newalpha, newbeta, newlogZobl, newlogZsur;
// //   mat newlogZ, newloglam, newlam;
// //   vec dummyvec;
// //   
// //   int count;
// //   vec xvals;
// //   double llprev, llthred, llnew, thetamin, thetamax, theta;
// //   vec prioreta, neweta, priorpsi, newpsi;
// //   bool accept;
// //   
// //   mat Reta = exp(-distmat / phi_eta), Rpsi = exp(-distmat / phi_psi);
// //   mat invReta = inv(Reta), invRpsi= inv(Rpsi), cholReta = chol(Reta), cholRpsi= chol(Rpsi);
// //   rowvec oneinvReta = trans(ones(nG)) * invReta, oneinvRpsi = trans(ones(nG)) * invRpsi;
// //   double oneinvRetaone = sum( oneinvReta * ones(nG) ), oneinvRpsione = sum( oneinvRpsi * ones(nG) );
// //   Rprintf("set up GP covariance matrices\n");
// //   
// //   mat logZ = cppf_logZ_eta(Xzoop, alpha, eta);
// //   vec logZobl = cppf_logZobl(logZ, indobl);
// //   vec logZsur = cppf_logZsur(logZ, indsur);
// //   mat loglam = cppf_loglam_cZ_b0t_psi(Xwhale, logZ, beta, psi, indtwhale);
// //   mat lam = exp(loglam);
// //   vec m_al = zeros(n_al), diag_al = ones(n_al)/100;
// //   Rprintf("set up initial values\n");
// //   
// //   
// //   // Start MCMC
// //   for(int s = 0; s < niter; s++) {
// //     
// //     // update alpha
// //     if( updateCOV ){
// //       if( (s+1 >= adaptInterval) && (s+1 - (adaptInterval * trunc((s+1) / adaptInterval)) == 0) ){
// //         dummyaccprob = accprob.col(0);
// //         rhat[0] = sum( dummyaccprob.rows(s+1-adaptInterval, s-1) ) / (adaptInterval-1);
// //         gamma1[0] = 1 / pow(adapIter[0], c1);
// //         gamma2[0] = c0 * gamma1[0];
// //         sigma2[0] = exp( log(sigma2[0]) + gamma2[0] * (rhat[0] - ropt) );
// //         
// //         COValpha = COValpha + gamma1[0] * ( cov( postAlpha.rows(s+1-adaptInterval, s-1) ) - COValpha );
// //         cholCOValpha = trans( chol( sigma2[0] * ( COValpha + 0.0000000001 * diagmat(ones(n_al)) ) ) );
// //         
// //         adapIter[0] = adapIter[0] + 1;
// //       }
// //     }
// //     
// //     newalpha = alpha + cholCOValpha * randn(n_al);
// //     
// //     newlogZ = cppf_logZ_eta(Xzoop, newalpha, eta);
// //     newlogZobl = cppf_logZobl(newlogZ, indobl);
// //     newlogZsur = cppf_logZsur(newlogZ, indsur);
// //     newloglam = cppf_loglam_cZ_b0t_psi(Xwhale, newlogZ, beta, psi, indtwhale);
// //     newlam = exp(newloglam);
// //     
// //     logprob = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) -
// //       cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //       cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) -
// //       cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
// //       cppf_loglik_obl(logYobl, newlogZobl, sig2_obl) -
// //       cppf_loglik_obl(logYobl, logZobl, sig2_obl) +
// //       cppf_loglik_sur(logYsur, newlogZsur, lam0, lam1, sig2_sur) -
// //       cppf_loglik_sur(logYsur, logZsur, lam0, lam1, sig2_sur) +
// //       MVN_logh(newalpha, zeros(n_al), diagmat(ones(n_al))/100) -
// //       MVN_logh(alpha, zeros(n_al), diagmat(ones(n_al))/100);
// //     
// //     if (log(randu()) < logprob) {
// //       
// //       alpha = newalpha;
// //       logZ = newlogZ;
// //       logZobl = newlogZobl;
// //       logZsur = newlogZsur;
// //       loglam = newloglam;
// //       lam = newlam;
// //       
// //       accprob(s, 0) = 1;
// //     }
// //     postAlpha.row(s) = trans(alpha);
// //     if(s == 0){ Rprintf("sampled alpha\n"); }
// //     
// //     
// //     
// //     // update eta via Elliptical slice sampling
// //     llprev = cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //       cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
// //       cppf_loglik_obl(logYobl, logZobl, sig2_obl) +
// //       cppf_loglik_sur(logYsur, logZsur, lam0, lam1, sig2_sur);
// //     prioreta = sqrt(kappa_eta) * trans(cholReta) * randn(nG);
// //     
// //     thetamin = 0;
// //     thetamax = 2 * M_PI;
// //     theta = thetamin + randu() * (thetamax - thetamin);
// //     thetamin = theta - 2 * M_PI;
// //     thetamax = theta;
// //     
// //     llthred = llprev + log(randu());
// //     accept = false;
// //     count = 0;
// //     
// //     while(accept == false){
// //       count = count + 1;
// //       
// //       neweta = eta * cos(theta) + prioreta * sin(theta);
// //       neweta = neweta - mean(neweta);
// //       
// //       newlogZ = cppf_logZ_eta(Xzoop, alpha, neweta);
// //       newlogZobl = cppf_logZobl(newlogZ, indobl);
// //       newlogZsur = cppf_logZsur(newlogZ, indsur);
// //       newloglam = cppf_loglam_cZ_b0t_psi(Xwhale, newlogZ, beta, psi, indtwhale);
// //       newlam = exp(newloglam);
// //       
// //       llnew = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //         cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) +
// //         cppf_loglik_obl(logYobl, newlogZobl, sig2_obl) +
// //         cppf_loglik_sur(logYsur, newlogZsur, lam0, lam1, sig2_sur);
// //       
// //       if(llnew > llthred){
// //         llprev = llnew;
// //         accept = true;
// //       } else {
// //         if(theta < 0){ thetamin = theta; } else { thetamax = theta; }
// //         theta = thetamin + randu() * (thetamax - thetamin);
// //         if ( (count) % 1000 == 0 ) {
// //           Rprintf("ESS for eta: %d iterations...\n", count);
// //         }
// //         if (count > 20000) stop("ESS is not converging...\n");
// //       }
// //     }
// //     eta = neweta;
// //     logZ = newlogZ;
// //     logZobl = newlogZobl;
// //     logZsur = newlogZsur;
// //     loglam = newloglam;
// //     lam = newlam;
// //     
// //     postEta.row(s) = trans(eta);
// //     if(s == 0){ Rprintf("sampled eta\n"); }
// //     
// //     
// //     // update kappa_eta
// //     kappa_eta = cppf_gibbs_kappa_eta(eta, invReta, 2, 2);
// //     postKappa(s,0) = kappa_eta;
// //     if(s == 0){ Rprintf("sampled kappa_eta\n"); }
// //     
// //     
// //     // update lam0
// //     lam0 = cppf_gibbs_lam0(logYsur, logZsur, lam1, sig2_sur, 100);
// //     postLam(s, 0) = lam0;
// //     if(s == 0){ Rprintf("sampled lam0\n"); }
// //     
// //     
// //     // update lam1
// //     lam1 = cppf_gibbs_lam1(logYsur, logZsur, lam0, sig2_sur, 100);
// //     postLam(s, 1) = lam1;
// //     if(s == 0){ Rprintf("sampled lam1\n"); }
// //     
// //     
// //     // update sig2_sur
// //     sig2_sur = cppf_gibbs_sig2_sur(logYsur, logZsur, lam0, lam1, 2, 2);
// //     postSig2(s, 0) = sig2_sur;
// //     if(s == 0){ Rprintf("sampled sig2_sur\n"); }
// //     
// //     
// //     // update beta
// //     if( updateCOV ){
// //       if( (s+1 >= adaptInterval) && (s+1 - (adaptInterval * trunc((s+1) / adaptInterval)) == 0) ){
// //         dummyaccprob = accprob.col(1);
// //         rhat[1] = sum( dummyaccprob.rows(s+1-adaptInterval, s-1) ) / (adaptInterval-1);
// //         gamma1[1] = 1 / pow(adapIter[1], c1);
// //         gamma2[1] = c0 * gamma1[1];
// //         sigma2[1] = exp( log(sigma2[1]) + gamma2[1] * (rhat[1] - ropt) );
// //         
// //         COVbeta = COVbeta + gamma1[1] * ( cov( postBeta.rows(s+1-adaptInterval, s-1)) - COVbeta );
// //         cholCOVbeta = trans( chol( sigma2[1] * ( COVbeta + 0.0000000001 * diagmat(ones(n_be)) ) ) );
// //         
// //         adapIter[1] = adapIter[1] + 1;
// //       }
// //     }
// //     
// //     newbeta = beta + cholCOVbeta * randn(n_be);
// //     
// //     newloglam = cppf_loglam_cZ_b0t_psi(Xwhale, logZ, newbeta, psi, indtwhale);
// //     newlam = exp(newloglam);
// //     
// //     logprob = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) -
// //       cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //       cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) -
// //       cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
// //       MVN_logh(newbeta, zeros(n_be), diagmat(ones(n_be))/100) -
// //       MVN_logh(beta, zeros(n_be), diagmat(ones(n_be))/100);
// //     
// //     if (log(randu()) < logprob) {
// //       
// //       beta = newbeta;
// //       loglam = newloglam;
// //       lam = newlam;
// //       
// //       accprob(s, 1) = 1;
// //     }
// //     
// //     postBeta.row(s) = trans(beta);
// //     if(s == 0){ Rprintf("sampled beta\n"); }
// //     
// //     
// //     
// //     // update psi via elliptical slice sampling
// //     llprev = cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //       cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea);
// //     priorpsi = sqrt(kappa_psi) * trans(cholRpsi) * randn(nG);
// //     
// //     thetamin = 0;
// //     thetamax = 2 * M_PI;
// //     theta = thetamin + randu() * (thetamax - thetamin);
// //     thetamin = theta - 2 * M_PI;
// //     thetamax = theta;
// //     
// //     llthred = llprev + log(randu());
// //     accept = false;
// //     count = 0;
// //     
// //     while(accept == false){
// //       count = count + 1;
// //       
// //       newpsi = psi * cos(theta) + priorpsi * sin(theta);
// //       newpsi = newpsi - mean(newpsi);
// //       
// //       newloglam = cppf_loglam_cZ_b0t_psi(Xwhale, logZ, beta, newpsi, indtwhale);
// //       newlam = exp(newloglam);
// //       
// //       llnew = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //         cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea);
// //       
// //       if(llnew > llthred){
// //         llprev = llnew;
// //         accept = true;
// //       } else {
// //         if(theta < 0){ thetamin = theta; } else { thetamax = theta; }
// //         theta = thetamin + randu() * (thetamax - thetamin);
// //         if ( (count) % 1000 == 0 ) {
// //           Rprintf("ESS for psi: %d iterations...\n", count);
// //         }
// //         if (count > 20000) stop("ESS is not converging...\n");
// //       }
// //     }
// //     psi = newpsi;
// //     loglam = newloglam;
// //     lam = newlam;
// //     
// //     postPsi.row(s) = trans(psi);
// //     if(s == 0){ Rprintf("sampled psi\n"); }
// //     
// //     
// //     // update kappa_psi
// //     kappa_psi = cppf_gibbs_kappa_psi(psi, invRpsi, 2, 2);
// //     postKappa(s,1) = kappa_psi;
// //     if(s == 0){ Rprintf("sampled kappa_psi\n"); }
// //     
// //     
// //     if ( (s+1) % 10 == 0 ) {
// //       Rprintf("Generated %d samples...\n", s+1);
// //     }
// //   }
// //   
// //   
// //   return Rcpp::List::create(Rcpp::Named("postAlpha") = postAlpha,
// //                             Rcpp::Named("postEta") = postEta,
// //                             Rcpp::Named("postLam") = postLam,
// //                             Rcpp::Named("postKappa") = postKappa,
// //                             Rcpp::Named("postSig2") = postSig2,
// //                             Rcpp::Named("postBeta") = postBeta,
// //                             Rcpp::Named("postPsi") = postPsi,
// //                             Rcpp::Named("Accprob") = accprob,
// //                             Rcpp::Named("sigma2") = sigma2,
// //                             Rcpp::Named("adapIter") = adapIter,
// //                             Rcpp::Named("COValpha") = COValpha,
// //                             Rcpp::Named("COVbeta") = COVbeta);
// //   
// // }
// // 
// // 
// // 
// // // [[Rcpp::export]]
// // List cppf_fitJoint_cZ_eta_b0t_psit_ct_fixC(
// //     int niter, mat distmat,
// //     List Xzoop, vec logYobl, vec logYsur,
// //     List Xwhale, mat Sdist, List Wpam, vec dur_dist,
// //     vec f_dist_obs, List f_dist, List p_pam,
// //     mat indobl, mat indsur, vec indtwhale, mat inddist, vec indmdist, vec indmwhale,
// //     vec alpha, vec eta, double lam0, double lam1,
// //     double kappa_eta, double phi_eta, double sig2_obl, double sig2_sur,
// //     vec beta, mat psi, vec pii, vec c,
// //     double kappa_psi, double phi_psi,
// //     double cellarea,
// //     vec sigma2, mat COValpha, mat COVbeta,
// //     bool updateCOV, int adaptInterval, double adaptFactorExponent, vec adapIter){
// //   
// //   double negativeInf = -std::numeric_limits<float>::infinity();;
// //   
// //   indobl = indobl-1, indsur = indsur-1, indtwhale = indtwhale-1;
// //   inddist = inddist-1, indmdist = indmdist-1, indmwhale = indmwhale-1;
// //   
// //   int nG = distmat.n_rows, n_al = alpha.size(), n_be = beta.size(), ntzoop = eta.n_cols, ntwhale = psi.n_cols, n_sigma2 = sigma2.size();
// //   
// //   mat postAlpha = zeros(niter, n_al), postLam = zeros(niter, 2);
// //   mat postSig2 = zeros(niter, 1), postKappa = zeros(niter, 2);
// //   mat postBeta = zeros(niter, n_be), postEta = zeros(niter, nG);
// //   mat postPsit = zeros(niter, nG);
// //   List postPsi(ntwhale);
// //   for(int i = 0; i < ntwhale; i ++){ postPsi[i] = zeros(niter, nG); }
// //   Rprintf("set up matrices for posterior samples\n");
// //   
// //   mat accprob = zeros(niter, n_sigma2);
// //   vec rhat = zeros(n_sigma2), gamma1 = zeros(n_sigma2), gamma2 = zeros(n_sigma2), dummyaccprob;
// //   double c0 = 1, c1 = adaptFactorExponent, ropt = 0.234, logprob;
// //   mat cholCOValpha = trans( chol( sigma2[0] * ( COValpha + 0.0000000001 * diagmat(ones(n_al)) ) ) );
// //   mat cholCOVbeta = trans( chol( sigma2[1] * ( COVbeta + 0.0000000001 * diagmat(ones(n_be)) ) ) );
// //   Rprintf("set up adaptive proposals\n");
// //   
// //   vec newalpha, newbeta, newlogZobl, newlogZsur;
// //   mat newlogZ, newloglam, newlam;
// //   vec dummyvec;
// //   
// //   int count;
// //   vec xvals;
// //   double llprev, llthred, llnew, thetamin, thetamax, theta;
// //   vec prioreta, neweta, priorpsit;
// //   mat newpsi;
// //   bool accept;
// //   
// //   mat Reta = exp(-distmat / phi_eta), Rpsi = exp(-distmat / phi_psi);
// //   mat invReta = inv(Reta), invRpsi= inv(Rpsi), cholReta = chol(Reta), cholRpsi= chol(Rpsi);
// //   rowvec oneinvReta = trans(ones(nG)) * invReta, oneinvRpsi = trans(ones(nG)) * invRpsi;
// //   double oneinvRetaone = sum( oneinvReta * ones(nG) ), oneinvRpsione = sum( oneinvRpsi * ones(nG) );
// //   Rprintf("set up GP covariance matrices\n");
// //   
// //   mat logZ = cppf_logZ_eta(Xzoop, alpha, eta);
// //   vec logZobl = cppf_logZobl(logZ, indobl);
// //   vec logZsur = cppf_logZsur(logZ, indsur);
// //   mat loglam = cppf_loglam_cZ_b0t_psit(Xwhale, logZ, beta, psi, indtwhale);
// //   mat lam = exp(loglam);
// //   vec m_al = zeros(n_al), diag_al = ones(n_al)/100;
// //   Rprintf("set up initial values\n");
// //   
// //   
// //   // Start MCMC
// //   for(int s = 0; s < niter; s++) {
// //     
// //     // update alpha
// //     if( updateCOV ){
// //       if( (s+1 >= adaptInterval) && (s+1 - (adaptInterval * trunc((s+1) / adaptInterval)) == 0) ){
// //         dummyaccprob = accprob.col(0);
// //         rhat[0] = sum( dummyaccprob.rows(s+1-adaptInterval, s-1) ) / (adaptInterval-1);
// //         gamma1[0] = 1 / pow(adapIter[0], c1);
// //         gamma2[0] = c0 * gamma1[0];
// //         sigma2[0] = exp( log(sigma2[0]) + gamma2[0] * (rhat[0] - ropt) );
// //         
// //         COValpha = COValpha + gamma1[0] * ( cov( postAlpha.rows(s+1-adaptInterval, s-1) ) - COValpha );
// //         cholCOValpha = trans( chol( sigma2[0] * ( COValpha + 0.0000000001 * diagmat(ones(n_al)) ) ) );
// //         
// //         adapIter[0] = adapIter[0] + 1;
// //       }
// //     }
// //     
// //     newalpha = alpha + cholCOValpha * randn(n_al);
// //     
// //     newlogZ = cppf_logZ_eta(Xzoop, newalpha, eta);
// //     newlogZobl = cppf_logZobl(newlogZ, indobl);
// //     newlogZsur = cppf_logZsur(newlogZ, indsur);
// //     newloglam = cppf_loglam_cZ_b0t_psit(Xwhale, newlogZ, beta, psi, indtwhale);
// //     newlam = exp(newloglam);
// //     
// //     logprob = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) -
// //       cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //       cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) -
// //       cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
// //       cppf_loglik_obl(logYobl, newlogZobl, sig2_obl) -
// //       cppf_loglik_obl(logYobl, logZobl, sig2_obl) +
// //       cppf_loglik_sur(logYsur, newlogZsur, lam0, lam1, sig2_sur) -
// //       cppf_loglik_sur(logYsur, logZsur, lam0, lam1, sig2_sur) +
// //       MVN_logh(newalpha, zeros(n_al), diagmat(ones(n_al))/100) -
// //       MVN_logh(alpha, zeros(n_al), diagmat(ones(n_al))/100);
// //     
// //     if (log(randu()) < logprob) {
// //       
// //       alpha = newalpha;
// //       logZ = newlogZ;
// //       logZobl = newlogZobl;
// //       logZsur = newlogZsur;
// //       loglam = newloglam;
// //       lam = newlam;
// //       
// //       accprob(s, 0) = 1;
// //     }
// //     postAlpha.row(s) = trans(alpha);
// //     if(s == 0){ Rprintf("sampled alpha\n"); }
// //     
// //     
// //     
// //     // update eta via Elliptical slice sampling
// //     llprev = cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //       cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
// //       cppf_loglik_obl(logYobl, logZobl, sig2_obl) +
// //       cppf_loglik_sur(logYsur, logZsur, lam0, lam1, sig2_sur);
// //     prioreta = sqrt(kappa_eta) * trans(cholReta) * randn(nG);
// //     
// //     thetamin = 0;
// //     thetamax = 2 * M_PI;
// //     theta = thetamin + randu() * (thetamax - thetamin);
// //     thetamin = theta - 2 * M_PI;
// //     thetamax = theta;
// //     
// //     llthred = llprev + log(randu());
// //     accept = false;
// //     count = 0;
// //     
// //     while(accept == false){
// //       count = count + 1;
// //       
// //       neweta = eta * cos(theta) + prioreta * sin(theta);
// //       neweta = neweta - mean(neweta);
// //       
// //       newlogZ = cppf_logZ_eta(Xzoop, alpha, neweta);
// //       newlogZobl = cppf_logZobl(newlogZ, indobl);
// //       newlogZsur = cppf_logZsur(newlogZ, indsur);
// //       newloglam = cppf_loglam_cZ_b0t_psit(Xwhale, newlogZ, beta, psi, indtwhale);
// //       newlam = exp(newloglam);
// //       
// //       llnew = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //         cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) +
// //         cppf_loglik_obl(logYobl, newlogZobl, sig2_obl) +
// //         cppf_loglik_sur(logYsur, newlogZsur, lam0, lam1, sig2_sur);
// //       
// //       if(llnew > llthred){
// //         llprev = llnew;
// //         accept = true;
// //       } else {
// //         if(theta < 0){ thetamin = theta; } else { thetamax = theta; }
// //         theta = thetamin + randu() * (thetamax - thetamin);
// //         if ( (count) % 1000 == 0 ) {
// //           Rprintf("ESS for eta: %d iterations...\n", count);
// //         }
// //         if (count > 20000) stop("ESS is not converging...\n");
// //       }
// //     }
// //     eta = neweta;
// //     logZ = newlogZ;
// //     logZobl = newlogZobl;
// //     logZsur = newlogZsur;
// //     loglam = newloglam;
// //     lam = newlam;
// //     
// //     postEta.row(s) = trans(eta);
// //     if(s == 0){ Rprintf("sampled eta\n"); }
// //     
// //     
// //     // update kappa_eta
// //     kappa_eta = cppf_gibbs_kappa_eta(eta, invReta, 2, 2);
// //     postKappa(s,0) = kappa_eta;
// //     if(s == 0){ Rprintf("sampled kappa_eta\n"); }
// //     
// //     
// //     // update lam0
// //     lam0 = cppf_gibbs_lam0(logYsur, logZsur, lam1, sig2_sur, 100);
// //     postLam(s, 0) = lam0;
// //     if(s == 0){ Rprintf("sampled lam0\n"); }
// //     
// //     
// //     // update lam1
// //     lam1 = cppf_gibbs_lam1(logYsur, logZsur, lam0, sig2_sur, 100);
// //     postLam(s, 1) = lam1;
// //     if(s == 0){ Rprintf("sampled lam1\n"); }
// //     
// //     
// //     // update sig2_sur
// //     sig2_sur = cppf_gibbs_sig2_sur(logYsur, logZsur, lam0, lam1, 2, 2);
// //     postSig2(s, 0) = sig2_sur;
// //     if(s == 0){ Rprintf("sampled sig2_sur\n"); }
// //     
// //     
// //     // update beta
// //     if( updateCOV ){
// //       if( (s+1 >= adaptInterval) && (s+1 - (adaptInterval * trunc((s+1) / adaptInterval)) == 0) ){
// //         dummyaccprob = accprob.col(1);
// //         rhat[1] = sum( dummyaccprob.rows(s+1-adaptInterval, s-1) ) / (adaptInterval-1);
// //         gamma1[1] = 1 / pow(adapIter[1], c1);
// //         gamma2[1] = c0 * gamma1[1];
// //         sigma2[1] = exp( log(sigma2[1]) + gamma2[1] * (rhat[1] - ropt) );
// //         
// //         COVbeta = COVbeta + gamma1[1] * ( cov( postBeta.rows(s+1-adaptInterval, s-1)) - COVbeta );
// //         cholCOVbeta = trans( chol( sigma2[1] * ( COVbeta + 0.0000000001 * diagmat(ones(n_be)) ) ) );
// //         
// //         adapIter[1] = adapIter[1] + 1;
// //       }
// //     }
// //     
// //     newbeta = beta + cholCOVbeta * randn(n_be);
// //     
// //     newloglam = cppf_loglam_cZ_b0t_psit(Xwhale, logZ, newbeta, psi, indtwhale);
// //     newlam = exp(newloglam);
// //     
// //     logprob = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) -
// //       cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //       cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) -
// //       cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
// //       MVN_logh(newbeta, zeros(n_be), diagmat(ones(n_be))/100) -
// //       MVN_logh(beta, zeros(n_be), diagmat(ones(n_be))/100);
// //     
// //     if (log(randu()) < logprob) {
// //       
// //       beta = newbeta;
// //       loglam = newloglam;
// //       lam = newlam;
// //       
// //       accprob(s, 1) = 1;
// //     }
// //     
// //     postBeta.row(s) = trans(beta);
// //     if(s == 0){ Rprintf("sampled beta\n"); }
// //     
// //     
// //     
// //     // update psi via Elliptical slice sampling
// //     for(int t = 0; t < ntwhale; t ++){
// //       newpsi = psi;
// //       
// //       llprev = cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //         cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea);
// //       priorpsit = sqrt(kappa_psi) * trans(cholRpsi) * randn(nG);
// //       
// //       thetamin = 0;
// //       thetamax = 2 * M_PI;
// //       theta = thetamin + randu() * (thetamax - thetamin);
// //       thetamin = theta - 2 * M_PI;
// //       thetamax = theta;
// //       
// //       llthred = llprev + log(randu());
// //       accept = false;
// //       count = 0;
// //       
// //       while(accept == false){
// //         count = count + 1;
// //         
// //         newpsi.col(t) = psi.col(t) * cos(theta) + priorpsit * sin(theta);
// //         newpsi.col(t) = newpsi.col(t) - mean(newpsi.col(t));
// //         
// //         newloglam = cppf_loglam_cZ_b0t_psit(Xwhale, logZ, beta, newpsi, indtwhale);
// //         newlam = exp(newloglam);
// //         
// //         llnew = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //           cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea);
// //         
// //         if(llnew > llthred){
// //           llprev = llnew;
// //           accept = true;
// //         } else {
// //           if(theta < 0){ thetamin = theta; } else { thetamax = theta; }
// //           theta = thetamin + randu() * (thetamax - thetamin);
// //           if ( (count) % 1000 == 0 ) {
// //             Rprintf("ESS for psi: %d iterations...\n", count);
// //           }
// //           if (count > 20000) stop("ESS is not converging...\n");
// //         }
// //       }
// //       psi = newpsi;
// //       loglam = newloglam;
// //       lam = newlam;
// //       
// //       postPsit = as<mat>(postPsi[t]);
// //       postPsit.row(s) = trans(psi.col(t));
// //       postPsi[t] = postPsit;
// //     }
// //     if(s == 0){ Rprintf("sampled psi\n"); }
// //     
// //     
// //     // update kappa_psi
// //     kappa_psi = cppf_gibbs_kappa_psit(psi, invRpsi, 2, 2);
// //     postKappa(s,1) = kappa_psi;
// //     if(s == 0){ Rprintf("sampled kappa_psi\n"); }
// //     
// //     
// //     if ( (s+1) % 10 == 0 ) {
// //       Rprintf("Generated %d samples...\n", s+1);
// //     }
// //   }
// //   
// //   
// //   return Rcpp::List::create(Rcpp::Named("postAlpha") = postAlpha,
// //                             Rcpp::Named("postEta") = postEta,
// //                             Rcpp::Named("postLam") = postLam,
// //                             Rcpp::Named("postKappa") = postKappa,
// //                             Rcpp::Named("postSig2") = postSig2,
// //                             Rcpp::Named("postBeta") = postBeta,
// //                             Rcpp::Named("postPsi") = postPsi,
// //                             Rcpp::Named("Accprob") = accprob,
// //                             Rcpp::Named("sigma2") = sigma2,
// //                             Rcpp::Named("adapIter") = adapIter,
// //                             Rcpp::Named("COValpha") = COValpha,
// //                             Rcpp::Named("COVbeta") = COVbeta);
// //   
// // }
// // 
// // 
// // // [[Rcpp::export]]
// // List cppf_fitJoint_cZ_a0t_etat_b0t_psi_ct_ja_fixC(
// //     int niter, mat distmat,
// //     List Xzoop, vec logYobl, vec logYsur,
// //     List Xwhale, mat Sdist, List Wpam, vec dur_dist,
// //     vec f_dist_obs, List f_dist, List p_pam,
// //     mat indobl, mat indsur, vec indtwhale, mat inddist, vec indmdist, vec indmwhale,
// //     double alpha0tilde, vec alpha, mat eta, double lam0, double lam1,
// //     double tau2, double kappa_eta, double phi_eta, double sig2_obl, double sig2_sur,
// //     vec beta, vec psi, vec pii, vec c, 
// //     double kappa_psi, double phi_psi,
// //     double cellarea,
// //     vec sigma2, mat COValpha, mat COVbeta,
// //     bool updateCOV, int adaptInterval, double adaptFactorExponent, vec adapIter){
// //   
// //   double negativeInf = -std::numeric_limits<float>::infinity();;
// //   
// //   indobl = indobl-1, indsur = indsur-1, indtwhale = indtwhale-1;
// //   inddist = inddist-1, indmdist = indmdist-1, indmwhale = indmwhale-1;
// //   
// //   int nG = distmat.n_rows, n_al = alpha.size(), n_be = beta.size(), ntzoop = eta.n_cols, n_sigma2 = sigma2.size();
// //   
// //   mat postAlpha0tilde = zeros(niter, 1), postAlpha = zeros(niter, n_al), postLam = zeros(niter, 2);
// //   mat postSig2 = zeros(niter, 1), postKappa = zeros(niter, 2);
// //   mat postBeta = zeros(niter, n_be);
// //   mat postPsi(niter, nG), postEtat = zeros(niter, nG);
// //   List postEta(ntzoop);
// //   for(int i = 0; i < ntzoop; i ++){
// //     postEta[i] = zeros(niter, nG);
// //   }
// //   Rprintf("set up matrices for posterior samples\n");
// //   
// //   mat accprob = zeros(niter, n_sigma2);
// //   vec rhat = zeros(n_sigma2), gamma1 = zeros(n_sigma2), gamma2 = zeros(n_sigma2), dummyaccprob;
// //   double c0 = 1, c1 = adaptFactorExponent, ropt = 0.234, logprob;
// //   mat cholCOValpha = trans( chol( sigma2[0] * ( COValpha + 0.0000000001 * diagmat(ones(n_al)) ) ) );
// //   mat cholCOVbeta = trans( chol( sigma2[1] * ( COVbeta + 0.0000000001 * diagmat(ones(n_be)) ) ) );
// //   Rprintf("set up adaptive proposals\n");
// //   
// //   vec newalpha, newbeta, newlogZobl, newlogZsur;
// //   mat newlogZ, newloglam, newlam;
// //   vec dummyvec;
// //   
// //   int count;
// //   vec xvals;
// //   double llprev, llthred, llnew, thetamin, thetamax, theta;
// //   vec prioretat, priorpsi, newpsi;
// //   mat neweta;
// //   bool accept;
// //   
// //   mat Reta = exp(-distmat / phi_eta), Rpsi = exp(-distmat / phi_psi);
// //   mat invReta = inv(Reta), invRpsi= inv(Rpsi), cholReta = chol(Reta), cholRpsi= chol(Rpsi);
// //   rowvec oneinvReta = trans(ones(nG)) * invReta, oneinvRpsi = trans(ones(nG)) * invRpsi;
// //   double oneinvRetaone = sum( oneinvReta * ones(nG) ), oneinvRpsione = sum( oneinvRpsi * ones(nG) );
// //   Rprintf("set up GP covariance matrices\n");
// //   
// //   mat logZ = cppf_logZ_a0t_etat_ja(Xzoop, alpha, eta);
// //   vec logZobl = cppf_logZobl(logZ, indobl);
// //   vec logZsur = cppf_logZsur(logZ, indsur);
// //   mat loglam = cppf_loglam_cZ_b0t_psi(Xwhale, logZ, beta, psi, indtwhale);
// //   mat lam = exp(loglam);
// //   Rprintf("set up initial values\n");
// //   
// //   
// //   // Start MCMC
// //   for(int s = 0; s < niter; s++) {
// //     
// //     // update alpha
// //     if( updateCOV ){
// //       if( (s+1 >= adaptInterval) && (s+1 - (adaptInterval * trunc((s+1) / adaptInterval)) == 0) ){
// //         dummyaccprob = accprob.col(0);
// //         rhat[0] = sum( dummyaccprob.rows(s+1-adaptInterval, s-1) ) / (adaptInterval-1);
// //         gamma1[0] = 1 / pow(adapIter[0], c1);
// //         gamma2[0] = c0 * gamma1[0];
// //         sigma2[0] = exp( log(sigma2[0]) + gamma2[0] * (rhat[0] - ropt) );
// //         
// //         COValpha = COValpha + gamma1[0] * ( cov( postAlpha.rows(s+1-adaptInterval, s-1) ) - COValpha );
// //         cholCOValpha = trans( chol( sigma2[0] * ( COValpha + 0.0000000001 * diagmat(ones(n_al)) ) ) );
// //         
// //         adapIter[0] = adapIter[0] + 1;
// //       }
// //     }
// //     
// //     newalpha = alpha + cholCOValpha * randn(n_al);
// //     
// //     newlogZ = cppf_logZ_a0t_etat_ja(Xzoop, newalpha, eta);
// //     newlogZobl = cppf_logZobl(newlogZ, indobl);
// //     newlogZsur = cppf_logZsur(newlogZ, indsur);
// //     newloglam = cppf_loglam_cZ_b0t_psi(Xwhale, newlogZ, beta, psi, indtwhale);
// //     newlam = exp(newloglam);
// //     
// //     logprob = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) -
// //       cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //       cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) -
// //       cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
// //       cppf_loglik_obl(logYobl, newlogZobl, sig2_obl) -
// //       cppf_loglik_obl(logYobl, logZobl, sig2_obl) +
// //       cppf_loglik_sur(logYsur, newlogZsur, lam0, lam1, sig2_sur) -
// //       cppf_loglik_sur(logYsur, logZsur, lam0, lam1, sig2_sur) +
// //       MVN_logh(newalpha, zeros(n_al), diagmat(ones(n_al))/100) -
// //       MVN_logh(alpha, zeros(n_al), diagmat(ones(n_al))/100);
// //     
// //     if (log(randu()) < logprob) {
// //       
// //       alpha = newalpha;
// //       logZ = newlogZ;
// //       logZobl = newlogZobl;
// //       logZsur = newlogZsur;
// //       loglam = newloglam;
// //       lam = newlam;
// //       
// //       accprob(s, 0) = 1;
// //     }
// //     postAlpha.row(s) = trans(alpha);
// //     if(s == 0){ Rprintf("sampled alpha\n"); }
// //     
// //     
// //     // Gibbs update for alpha0tilde
// //     alpha0tilde = cppf_gibbs_alpha0tilde(alpha.rows(0, ntzoop-1), tau2);
// //     postAlpha0tilde[s] = alpha0tilde;
// //     if(s == 0){ Rprintf("sampled alpha0tilde\n"); }
// //     
// //     
// //     // update eta via Elliptical slice sampling
// //     for(int t = 0; t < ntzoop; t ++){
// //       neweta = eta;
// //       
// //       llprev = cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //         cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
// //         cppf_loglik_obl(logYobl, logZobl, sig2_obl) +
// //         cppf_loglik_sur(logYsur, logZsur, lam0, lam1, sig2_sur);
// //       prioretat = sqrt(kappa_eta) * trans(cholReta) * randn(nG);
// //       
// //       thetamin = 0;
// //       thetamax = 2 * M_PI;
// //       theta = thetamin + randu() * (thetamax - thetamin);
// //       thetamin = theta - 2 * M_PI;
// //       thetamax = theta;
// //       
// //       llthred = llprev + log(randu());
// //       accept = false;
// //       count = 0;
// //       
// //       while(accept == false){
// //         count = count + 1;
// //         
// //         neweta.col(t) = eta.col(t) * cos(theta) + prioretat * sin(theta);
// //         neweta.col(t) = neweta.col(t) - mean(neweta.col(t));
// //         
// //         newlogZ = cppf_logZ_a0t_etat_ja(Xzoop, alpha, neweta);
// //         newlogZobl = cppf_logZobl(newlogZ, indobl);
// //         newlogZsur = cppf_logZsur(newlogZ, indsur);
// //         newloglam = cppf_loglam_cZ_b0t_psi(Xwhale, newlogZ, beta, psi, indtwhale);
// //         newlam = exp(newloglam);
// //         
// //         llnew = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //           cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) +
// //           cppf_loglik_obl(logYobl, newlogZobl, sig2_obl) +
// //           cppf_loglik_sur(logYsur, newlogZsur, lam0, lam1, sig2_sur);
// //         
// //         if(llnew > llthred){
// //           llprev = llnew;
// //           accept = true;
// //         } else {
// //           if(theta < 0){ thetamin = theta; } else { thetamax = theta; }
// //           theta = thetamin + randu() * (thetamax - thetamin);
// //           if ( (count) % 1000 == 0 ) {
// //             Rprintf("ESS for eta: %d iterations...\n", count);
// //           }
// //           if (count > 20000) stop("ESS is not converging...\n");
// //         }
// //       }
// //       eta = neweta;
// //       logZ = newlogZ;
// //       logZobl = newlogZobl;
// //       logZsur = newlogZsur;
// //       loglam = newloglam;
// //       lam = newlam;
// //       
// //       postEtat = as<mat>(postEta[t]);
// //       postEtat.row(s) = trans(eta.col(t));
// //       postEta[t] = postEtat;
// //     }
// //     if(s == 0){ Rprintf("sampled eta\n"); }
// //     
// //     
// //     // update kappa_eta
// //     kappa_eta = cppf_gibbs_kappa_etat(eta, invReta, 2, 2);
// //     postKappa(s,0) = kappa_eta;
// //     if(s == 0){ Rprintf("sampled kappa_eta\n"); }
// //     
// //     
// //     // update lam0
// //     lam0 = cppf_gibbs_lam0(logYsur, logZsur, lam1, sig2_sur, 100);
// //     postLam(s, 0) = lam0;
// //     if(s == 0){ Rprintf("sampled lam0\n"); }
// //     
// //     
// //     // update lam1
// //     lam1 = cppf_gibbs_lam1(logYsur, logZsur, lam0, sig2_sur, 100);
// //     postLam(s, 1) = lam1;
// //     if(s == 0){ Rprintf("sampled lam1\n"); }
// //     
// //     
// //     // update sig2_sur
// //     sig2_sur = cppf_gibbs_sig2_sur(logYsur, logZsur, lam0, lam1, 2, 2);
// //     postSig2(s, 0) = sig2_sur;
// //     if(s == 0){ Rprintf("sampled sig2_sur\n"); }
// //     
// //     
// //     // update beta
// //     if( updateCOV ){
// //       if( (s+1 >= adaptInterval) && (s+1 - (adaptInterval * trunc((s+1) / adaptInterval)) == 0) ){
// //         dummyaccprob = accprob.col(1);
// //         rhat[1] = sum( dummyaccprob.rows(s+1-adaptInterval, s-1) ) / (adaptInterval-1);
// //         gamma1[1] = 1 / pow(adapIter[1], c1);
// //         gamma2[1] = c0 * gamma1[1];
// //         sigma2[1] = exp( log(sigma2[1]) + gamma2[1] * (rhat[1] - ropt) );
// //         
// //         COVbeta = COVbeta + gamma1[1] * ( cov( postBeta.rows(s+1-adaptInterval, s-1)) - COVbeta );
// //         cholCOVbeta = trans( chol( sigma2[1] * ( COVbeta + 0.0000000001 * diagmat(ones(n_be)) ) ) );
// //         
// //         adapIter[1] = adapIter[1] + 1;
// //       }
// //     }
// //     
// //     newbeta = beta + cholCOVbeta * randn(n_be);
// //     
// //     newloglam = cppf_loglam_cZ_b0t_psi(Xwhale, logZ, newbeta, psi, indtwhale);
// //     newlam = exp(newloglam);
// //     
// //     logprob = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) -
// //       cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //       cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) -
// //       cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
// //       MVN_logh(newbeta, zeros(n_be), diagmat(ones(n_be))/100) -
// //       MVN_logh(beta, zeros(n_be), diagmat(ones(n_be))/100);
// //     
// //     if (log(randu()) < logprob) {
// //       
// //       beta = newbeta;
// //       loglam = newloglam;
// //       lam = newlam;
// //       
// //       accprob(s, 1) = 1;
// //     }
// //     
// //     postBeta.row(s) = trans(beta);
// //     if(s == 0){ Rprintf("sampled beta\n"); }
// //     
// //     
// //     // update psi via elliptical slice sampling
// //     llprev = cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //       cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea);
// //     priorpsi = sqrt(kappa_psi) * trans(cholRpsi) * randn(nG);
// //     
// //     thetamin = 0;
// //     thetamax = 2 * M_PI;
// //     theta = thetamin + randu() * (thetamax - thetamin);
// //     thetamin = theta - 2 * M_PI;
// //     thetamax = theta;
// //     
// //     llthred = llprev + log(randu());
// //     accept = false;
// //     count = 0;
// //     
// //     while(accept == false){
// //       count = count + 1;
// //       
// //       newpsi = psi * cos(theta) + priorpsi * sin(theta);
// //       newpsi = newpsi - mean(newpsi);
// //       
// //       newloglam = cppf_loglam_cZ_b0t_psi(Xwhale, logZ, beta, newpsi, indtwhale);
// //       newlam = exp(newloglam);
// //       
// //       llnew = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //         cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea);
// //       
// //       if(llnew > llthred){
// //         llprev = llnew;
// //         accept = true;
// //       } else {
// //         if(theta < 0){ thetamin = theta; } else { thetamax = theta; }
// //         theta = thetamin + randu() * (thetamax - thetamin);
// //         if ( (count) % 1000 == 0 ) {
// //           Rprintf("ESS for psi: %d iterations...\n", count);
// //         }
// //         if (count > 20000) stop("ESS is not converging...\n");
// //       }
// //     }
// //     psi = newpsi;
// //     loglam = newloglam;
// //     lam = newlam;
// //     
// //     postPsi.row(s) = trans(psi);
// //     if(s == 0){ Rprintf("sampled psi\n"); }
// //     
// //     
// //     // update kappa_psi
// //     kappa_psi = cppf_gibbs_kappa_psi(psi, invRpsi, 2, 2);
// //     postKappa(s,1) = kappa_psi;
// //     if(s == 0){ Rprintf("sampled kappa_psi\n"); }
// //     
// //     
// //     if ( (s+1) % 10 == 0 ) {
// //       Rprintf("Generated %d samples...\n", s+1);
// //     }
// //   }
// //   
// //   
// //   return Rcpp::List::create(Rcpp::Named("postAlpha0tilde") = postAlpha0tilde,
// //                             Rcpp::Named("postAlpha") = postAlpha,
// //                             Rcpp::Named("postEta") = postEta,
// //                             Rcpp::Named("postLam") = postLam,
// //                             Rcpp::Named("postKappa") = postKappa,
// //                             Rcpp::Named("postSig2") = postSig2,
// //                             Rcpp::Named("postBeta") = postBeta,
// //                             Rcpp::Named("postPsi") = postPsi,
// //                             Rcpp::Named("Accprob") = accprob,
// //                             Rcpp::Named("sigma2") = sigma2,
// //                             Rcpp::Named("adapIter") = adapIter,
// //                             Rcpp::Named("COValpha") = COValpha,
// //                             Rcpp::Named("COVbeta") = COVbeta);
// //   
// // }
// // 
// // 
// // 
// // // ============================================================================-
// // // estimate call rates ----
// // // ============================================================================-
// // 
// // // [[Rcpp::export]]
// // List cppf_fitJoint_eta_b0t_psi_ct(
// //     int niter, mat distmat,
// //     List Xzoop, vec logYobl, vec logYsur,
// //     List Xwhale, mat Sdist, List Wpam, vec dur_dist,
// //     vec f_dist_obs, List f_dist, List p_pam,
// //     mat indobl, mat indsur, vec indtwhale, mat inddist, vec indmdist, vec indmwhale,
// //     vec alpha, vec eta, double lam0, double lam1,
// //     double kappa_eta, double phi_eta, double sig2_obl, double sig2_sur,
// //     vec beta, vec psi, vec pii, vec c, double shape_c, double scale_c,
// //     double kappa_psi, double phi_psi,
// //     double cellarea,
// //     vec sigma2, mat COValpha, mat COVbeta,
// //     bool updateCOV, int adaptInterval, double adaptFactorExponent, vec adapIter){
// //   
// //   double negativeInf = -std::numeric_limits<float>::infinity();;
// //   
// //   indobl = indobl-1, indsur = indsur-1, indtwhale = indtwhale-1;
// //   inddist = inddist-1, indmdist = indmdist-1, indmwhale = indmwhale-1;
// //   
// //   int nG = distmat.n_rows, n_al = alpha.size(), n_be = beta.size(), ntzoop = eta.n_cols, n_sigma2 = sigma2.size();
// //   
// //   mat postAlpha = zeros(niter, n_al), postLam = zeros(niter, 2);
// //   mat postSig2 = zeros(niter, 1), postKappa = zeros(niter, 2);
// //   mat postBeta = zeros(niter, n_be), postEta = zeros(niter, nG);
// //   mat postPsi = zeros(niter, nG), postC = zeros(niter, c.size());
// //   Rprintf("set up matrices for posterior samples\n");
// //   
// //   mat accprob = zeros(niter, n_sigma2);
// //   vec rhat = zeros(n_sigma2), gamma1 = zeros(n_sigma2), gamma2 = zeros(n_sigma2), dummyaccprob;
// //   double c0 = 1, c1 = adaptFactorExponent, ropt = 0.234, logprob;
// //   mat cholCOValpha = trans( chol( sigma2[0] * ( COValpha + 0.0000000001 * diagmat(ones(n_al)) ) ) );
// //   mat cholCOVbeta = trans( chol( sigma2[1] * ( COVbeta + 0.0000000001 * diagmat(ones(n_be)) ) ) );
// //   Rprintf("set up adaptive proposals\n");
// //   
// //   vec newalpha, newbeta, newlogZobl, newlogZsur;
// //   mat newlogZ, newloglam, newlam;
// //   vec dummyvec;
// //   
// //   int count;
// //   vec xvals;
// //   double llprev, llthred, llnew, thetamin, thetamax, theta;
// //   vec prioreta, neweta, priorpsi, newpsi;
// //   bool accept;
// //   
// //   mat Reta = exp(-distmat / phi_eta), Rpsi = exp(-distmat / phi_psi);
// //   mat invReta = inv(Reta), invRpsi= inv(Rpsi), cholReta = chol(Reta), cholRpsi= chol(Rpsi);
// //   rowvec oneinvReta = trans(ones(nG)) * invReta, oneinvRpsi = trans(ones(nG)) * invRpsi;
// //   double oneinvRetaone = sum( oneinvReta * ones(nG) ), oneinvRpsione = sum( oneinvRpsi * ones(nG) );
// //   Rprintf("set up GP covariance matrices\n");
// //   
// //   mat logZ = cppf_logZ_eta(Xzoop, alpha, eta);
// //   vec logZobl = cppf_logZobl(logZ, indobl);
// //   vec logZsur = cppf_logZsur(logZ, indsur);
// //   mat loglam = cppf_loglam_b0t_psi(Xwhale, logZ, beta, psi, indtwhale);
// //   mat lam = exp(loglam);
// //   vec m_al = zeros(n_al), diag_al = ones(n_al)/100;
// //   Rprintf("set up initial values\n");
// //   
// //   
// //   // Start MCMC
// //   for(int s = 0; s < niter; s++) {
// //     
// //     // update alpha
// //     if( updateCOV ){
// //       if( (s+1 >= adaptInterval) && (s+1 - (adaptInterval * trunc((s+1) / adaptInterval)) == 0) ){
// //         dummyaccprob = accprob.col(0);
// //         rhat[0] = sum( dummyaccprob.rows(s+1-adaptInterval, s-1) ) / (adaptInterval-1);
// //         gamma1[0] = 1 / pow(adapIter[0], c1);
// //         gamma2[0] = c0 * gamma1[0];
// //         sigma2[0] = exp( log(sigma2[0]) + gamma2[0] * (rhat[0] - ropt) );
// //         
// //         COValpha = COValpha + gamma1[0] * ( cov( postAlpha.rows(s+1-adaptInterval, s-1) ) - COValpha );
// //         cholCOValpha = trans( chol( sigma2[0] * ( COValpha + 0.0000000001 * diagmat(ones(n_al)) ) ) );
// //         
// //         adapIter[0] = adapIter[0] + 1;
// //       }
// //     }
// //     
// //     newalpha = alpha + cholCOValpha * randn(n_al);
// //     
// //     newlogZ = cppf_logZ_eta(Xzoop, newalpha, eta);
// //     newlogZobl = cppf_logZobl(newlogZ, indobl);
// //     newlogZsur = cppf_logZsur(newlogZ, indsur);
// //     newloglam = cppf_loglam_b0t_psi(Xwhale, newlogZ, beta, psi, indtwhale);
// //     newlam = exp(newloglam);
// //     
// //     logprob = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) -
// //       cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //       cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) -
// //       cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
// //       cppf_loglik_obl(logYobl, newlogZobl, sig2_obl) -
// //       cppf_loglik_obl(logYobl, logZobl, sig2_obl) +
// //       cppf_loglik_sur(logYsur, newlogZsur, lam0, lam1, sig2_sur) -
// //       cppf_loglik_sur(logYsur, logZsur, lam0, lam1, sig2_sur) +
// //       MVN_logh(newalpha, zeros(n_al), diagmat(ones(n_al))/100) -
// //       MVN_logh(alpha, zeros(n_al), diagmat(ones(n_al))/100);
// //     
// //     if (log(randu()) < logprob) {
// //       
// //       alpha = newalpha;
// //       logZ = newlogZ;
// //       logZobl = newlogZobl;
// //       logZsur = newlogZsur;
// //       loglam = newloglam;
// //       lam = newlam;
// //       
// //       accprob(s, 0) = 1;
// //     }
// //     postAlpha.row(s) = trans(alpha);
// //     if(s == 0){ Rprintf("sampled alpha\n"); }
// //     
// //     
// //     
// //     // update eta via Elliptical slice sampling
// //     llprev = cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //       cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
// //       cppf_loglik_obl(logYobl, logZobl, sig2_obl) +
// //       cppf_loglik_sur(logYsur, logZsur, lam0, lam1, sig2_sur);
// //     prioreta = sqrt(kappa_eta) * trans(cholReta) * randn(nG);
// //     
// //     thetamin = 0;
// //     thetamax = 2 * M_PI;
// //     theta = thetamin + randu() * (thetamax - thetamin);
// //     thetamin = theta - 2 * M_PI;
// //     thetamax = theta;
// //     
// //     llthred = llprev + log(randu());
// //     accept = false;
// //     count = 0;
// //     
// //     while(accept == false){
// //       count = count + 1;
// //       
// //       neweta = eta * cos(theta) + prioreta * sin(theta);
// //       neweta = neweta - mean(neweta);
// //       
// //       newlogZ = cppf_logZ_eta(Xzoop, alpha, neweta);
// //       newlogZobl = cppf_logZobl(newlogZ, indobl);
// //       newlogZsur = cppf_logZsur(newlogZ, indsur);
// //       newloglam = cppf_loglam_b0t_psi(Xwhale, newlogZ, beta, psi, indtwhale);
// //       newlam = exp(newloglam);
// //       
// //       llnew = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //         cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) +
// //         cppf_loglik_obl(logYobl, newlogZobl, sig2_obl) +
// //         cppf_loglik_sur(logYsur, newlogZsur, lam0, lam1, sig2_sur);
// //       
// //       if(llnew > llthred){
// //         llprev = llnew;
// //         accept = true;
// //       } else {
// //         if(theta < 0){ thetamin = theta; } else { thetamax = theta; }
// //         theta = thetamin + randu() * (thetamax - thetamin);
// //         if ( (count) % 1000 == 0 ) {
// //           Rprintf("ESS for eta: %d iterations...\n", count);
// //         }
// //         if (count > 20000) stop("ESS is not converging...\n");
// //       }
// //     }
// //     eta = neweta;
// //     logZ = newlogZ;
// //     logZobl = newlogZobl;
// //     logZsur = newlogZsur;
// //     loglam = newloglam;
// //     lam = newlam;
// //     
// //     postEta.row(s) = trans(eta);
// //     if(s == 0){ Rprintf("sampled eta\n"); }
// //     
// //     
// //     // update kappa_eta
// //     kappa_eta = cppf_gibbs_kappa_eta(eta, invReta, 2, 2);
// //     postKappa(s,0) = kappa_eta;
// //     if(s == 0){ Rprintf("sampled kappa_eta\n"); }
// //     
// //     
// //     // update lam0
// //     lam0 = cppf_gibbs_lam0(logYsur, logZsur, lam1, sig2_sur, 100);
// //     postLam(s, 0) = lam0;
// //     if(s == 0){ Rprintf("sampled lam0\n"); }
// //     
// //     
// //     // update lam1
// //     lam1 = cppf_gibbs_lam1(logYsur, logZsur, lam0, sig2_sur, 100);
// //     postLam(s, 1) = lam1;
// //     if(s == 0){ Rprintf("sampled lam1\n"); }
// //     
// //     
// //     // update sig2_sur
// //     sig2_sur = cppf_gibbs_sig2_sur(logYsur, logZsur, lam0, lam1, 2, 2);
// //     postSig2(s, 0) = sig2_sur;
// //     if(s == 0){ Rprintf("sampled sig2_sur\n"); }
// //     
// //     
// //     // update beta
// //     if( updateCOV ){
// //       if( (s+1 >= adaptInterval) && (s+1 - (adaptInterval * trunc((s+1) / adaptInterval)) == 0) ){
// //         dummyaccprob = accprob.col(1);
// //         rhat[1] = sum( dummyaccprob.rows(s+1-adaptInterval, s-1) ) / (adaptInterval-1);
// //         gamma1[1] = 1 / pow(adapIter[1], c1);
// //         gamma2[1] = c0 * gamma1[1];
// //         sigma2[1] = exp( log(sigma2[1]) + gamma2[1] * (rhat[1] - ropt) );
// //         
// //         COVbeta = COVbeta + gamma1[1] * ( cov( postBeta.rows(s+1-adaptInterval, s-1)) - COVbeta );
// //         cholCOVbeta = trans( chol( sigma2[1] * ( COVbeta + 0.0000000001 * diagmat(ones(n_be)) ) ) );
// //         
// //         adapIter[1] = adapIter[1] + 1;
// //       }
// //     }
// //     
// //     newbeta = beta + cholCOVbeta * randn(n_be);
// //     
// //     newloglam = cppf_loglam_b0t_psi(Xwhale, logZ, newbeta, psi, indtwhale);
// //     newlam = exp(newloglam);
// //     
// //     logprob = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) -
// //       cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //       cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) -
// //       cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
// //       MVN_logh(newbeta, zeros(n_be), diagmat(ones(n_be))/100) -
// //       MVN_logh(beta, zeros(n_be), diagmat(ones(n_be))/100);
// //     
// //     if (log(randu()) < logprob) {
// //       
// //       beta = newbeta;
// //       loglam = newloglam;
// //       lam = newlam;
// //       
// //       accprob(s, 1) = 1;
// //     }
// //     
// //     postBeta.row(s) = trans(beta);
// //     if(s == 0){ Rprintf("sampled beta\n"); }
// //     
// //     
// //     
// //     // update psi via elliptical slice sampling
// //     llprev = cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //       cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea);
// //     priorpsi = sqrt(kappa_psi) * trans(cholRpsi) * randn(nG);
// //     
// //     thetamin = 0;
// //     thetamax = 2 * M_PI;
// //     theta = thetamin + randu() * (thetamax - thetamin);
// //     thetamin = theta - 2 * M_PI;
// //     thetamax = theta;
// //     
// //     llthred = llprev + log(randu());
// //     accept = false;
// //     count = 0;
// //     
// //     while(accept == false){
// //       count = count + 1;
// //       
// //       newpsi = psi * cos(theta) + priorpsi * sin(theta);
// //       newpsi = newpsi - mean(newpsi);
// //       
// //       newloglam = cppf_loglam_b0t_psi(Xwhale, logZ, beta, newpsi, indtwhale);
// //       newlam = exp(newloglam);
// //       
// //       llnew = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //         cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea);
// //       
// //       if(llnew > llthred){
// //         llprev = llnew;
// //         accept = true;
// //       } else {
// //         if(theta < 0){ thetamin = theta; } else { thetamax = theta; }
// //         theta = thetamin + randu() * (thetamax - thetamin);
// //         if ( (count) % 1000 == 0 ) {
// //           Rprintf("ESS for psi: %d iterations...\n", count);
// //         }
// //         if (count > 20000) stop("ESS is not converging...\n");
// //       }
// //     }
// //     psi = newpsi;
// //     loglam = newloglam;
// //     lam = newlam;
// //     
// //     postPsi.row(s) = trans(psi);
// //     if(s == 0){ Rprintf("sampled psi\n"); }
// //     
// //     
// //     // update kappa_psi
// //     kappa_psi = cppf_gibbs_kappa_psi(psi, invRpsi, 2, 2);
// //     postKappa(s,1) = kappa_psi;
// //     if(s == 0){ Rprintf("sampled kappa_psi\n"); }
// //     
// //     
// //     // update c
// //     c = cppf_gibbs_ct(Wpam, lam, dur_dist, p_pam, cellarea, shape_c, scale_c);
// //     postC.row(s) = trans(c);
// //     if(s == 0){ Rprintf("sampled c\n"); }
// //     
// //     
// //     if ( (s+1) % 10 == 0 ) {
// //       Rprintf("Generated %d samples...\n", s+1);
// //     }
// //   }
// //   
// //   
// //   return Rcpp::List::create(Rcpp::Named("postAlpha") = postAlpha,
// //                             Rcpp::Named("postEta") = postEta,
// //                             Rcpp::Named("postLam") = postLam,
// //                             Rcpp::Named("postKappa") = postKappa,
// //                             Rcpp::Named("postSig2") = postSig2,
// //                             Rcpp::Named("postBeta") = postBeta,
// //                             Rcpp::Named("postC") = postC,
// //                             Rcpp::Named("postPsi") = postPsi,
// //                             Rcpp::Named("Accprob") = accprob,
// //                             Rcpp::Named("sigma2") = sigma2,
// //                             Rcpp::Named("adapIter") = adapIter,
// //                             Rcpp::Named("COValpha") = COValpha,
// //                             Rcpp::Named("COVbeta") = COVbeta);
// //   
// // }
// // 
// // 
// // 
// // // [[Rcpp::export]]
// // List cppf_fitJoint_eta_b0t_psit_ct(
// //     int niter, mat distmat,
// //     List Xzoop, vec logYobl, vec logYsur,
// //     List Xwhale, mat Sdist, List Wpam, vec dur_dist,
// //     vec f_dist_obs, List f_dist, List p_pam,
// //     mat indobl, mat indsur, vec indtwhale, mat inddist, vec indmdist, vec indmwhale,
// //     vec alpha, vec eta, double lam0, double lam1,
// //     double kappa_eta, double phi_eta, double sig2_obl, double sig2_sur,
// //     vec beta, mat psi, vec pii, vec c, double shape_c, double scale_c,
// //     double kappa_psi, double phi_psi,
// //     double cellarea,
// //     vec sigma2, mat COValpha, mat COVbeta,
// //     bool updateCOV, int adaptInterval, double adaptFactorExponent, vec adapIter){
// //   
// //   double negativeInf = -std::numeric_limits<float>::infinity();;
// //   
// //   indobl = indobl-1, indsur = indsur-1, indtwhale = indtwhale-1;
// //   inddist = inddist-1, indmdist = indmdist-1, indmwhale = indmwhale-1;
// //   
// //   int nG = distmat.n_rows, n_al = alpha.size(), n_be = beta.size(), ntzoop = eta.n_cols, ntwhale = psi.n_cols, n_sigma2 = sigma2.size();
// //   
// //   mat postAlpha = zeros(niter, n_al), postLam = zeros(niter, 2);
// //   mat postSig2 = zeros(niter, 1), postKappa = zeros(niter, 2);
// //   mat postBeta = zeros(niter, n_be), postEta = zeros(niter, nG);
// //   mat postPsit = zeros(niter, nG), postC = zeros(niter, c.size());
// //   List postPsi(ntwhale);
// //   for(int i = 0; i < ntwhale; i ++){ postPsi[i] = zeros(niter, nG); }
// //   Rprintf("set up matrices for posterior samples\n");
// //   
// //   mat accprob = zeros(niter, n_sigma2);
// //   vec rhat = zeros(n_sigma2), gamma1 = zeros(n_sigma2), gamma2 = zeros(n_sigma2), dummyaccprob;
// //   double c0 = 1, c1 = adaptFactorExponent, ropt = 0.234, logprob;
// //   mat cholCOValpha = trans( chol( sigma2[0] * ( COValpha + 0.0000000001 * diagmat(ones(n_al)) ) ) );
// //   mat cholCOVbeta = trans( chol( sigma2[1] * ( COVbeta + 0.0000000001 * diagmat(ones(n_be)) ) ) );
// //   Rprintf("set up adaptive proposals\n");
// //   
// //   vec newalpha, newbeta, newlogZobl, newlogZsur;
// //   mat newlogZ, newloglam, newlam;
// //   vec dummyvec;
// //   
// //   int count;
// //   vec xvals;
// //   double llprev, llthred, llnew, thetamin, thetamax, theta;
// //   vec prioreta, neweta, priorpsit;
// //   mat newpsi;
// //   bool accept;
// //   
// //   mat Reta = exp(-distmat / phi_eta), Rpsi = exp(-distmat / phi_psi);
// //   mat invReta = inv(Reta), invRpsi= inv(Rpsi), cholReta = chol(Reta), cholRpsi= chol(Rpsi);
// //   rowvec oneinvReta = trans(ones(nG)) * invReta, oneinvRpsi = trans(ones(nG)) * invRpsi;
// //   double oneinvRetaone = sum( oneinvReta * ones(nG) ), oneinvRpsione = sum( oneinvRpsi * ones(nG) );
// //   Rprintf("set up GP covariance matrices\n");
// //   
// //   mat logZ = cppf_logZ_eta(Xzoop, alpha, eta);
// //   vec logZobl = cppf_logZobl(logZ, indobl);
// //   vec logZsur = cppf_logZsur(logZ, indsur);
// //   mat loglam = cppf_loglam_b0t_psit(Xwhale, logZ, beta, psi, indtwhale);
// //   mat lam = exp(loglam);
// //   vec m_al = zeros(n_al), diag_al = ones(n_al)/100;
// //   Rprintf("set up initial values\n");
// //   
// //   
// //   // Start MCMC
// //   for(int s = 0; s < niter; s++) {
// //     
// //     // update alpha
// //     if( updateCOV ){
// //       if( (s+1 >= adaptInterval) && (s+1 - (adaptInterval * trunc((s+1) / adaptInterval)) == 0) ){
// //         dummyaccprob = accprob.col(0);
// //         rhat[0] = sum( dummyaccprob.rows(s+1-adaptInterval, s-1) ) / (adaptInterval-1);
// //         gamma1[0] = 1 / pow(adapIter[0], c1);
// //         gamma2[0] = c0 * gamma1[0];
// //         sigma2[0] = exp( log(sigma2[0]) + gamma2[0] * (rhat[0] - ropt) );
// //         
// //         COValpha = COValpha + gamma1[0] * ( cov( postAlpha.rows(s+1-adaptInterval, s-1) ) - COValpha );
// //         cholCOValpha = trans( chol( sigma2[0] * ( COValpha + 0.0000000001 * diagmat(ones(n_al)) ) ) );
// //         
// //         adapIter[0] = adapIter[0] + 1;
// //       }
// //     }
// //     
// //     newalpha = alpha + cholCOValpha * randn(n_al);
// //     
// //     newlogZ = cppf_logZ_eta(Xzoop, newalpha, eta);
// //     newlogZobl = cppf_logZobl(newlogZ, indobl);
// //     newlogZsur = cppf_logZsur(newlogZ, indsur);
// //     newloglam = cppf_loglam_b0t_psit(Xwhale, newlogZ, beta, psi, indtwhale);
// //     newlam = exp(newloglam);
// //     
// //     logprob = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) -
// //       cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //       cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) -
// //       cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
// //       cppf_loglik_obl(logYobl, newlogZobl, sig2_obl) -
// //       cppf_loglik_obl(logYobl, logZobl, sig2_obl) +
// //       cppf_loglik_sur(logYsur, newlogZsur, lam0, lam1, sig2_sur) -
// //       cppf_loglik_sur(logYsur, logZsur, lam0, lam1, sig2_sur) +
// //       MVN_logh(newalpha, zeros(n_al), diagmat(ones(n_al))/100) -
// //       MVN_logh(alpha, zeros(n_al), diagmat(ones(n_al))/100);
// //     
// //     if (log(randu()) < logprob) {
// //       
// //       alpha = newalpha;
// //       logZ = newlogZ;
// //       logZobl = newlogZobl;
// //       logZsur = newlogZsur;
// //       loglam = newloglam;
// //       lam = newlam;
// //       
// //       accprob(s, 0) = 1;
// //     }
// //     postAlpha.row(s) = trans(alpha);
// //     if(s == 0){ Rprintf("sampled alpha\n"); }
// //     
// //     
// //     
// //     // update eta via Elliptical slice sampling
// //     llprev = cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //       cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
// //       cppf_loglik_obl(logYobl, logZobl, sig2_obl) +
// //       cppf_loglik_sur(logYsur, logZsur, lam0, lam1, sig2_sur);
// //     prioreta = sqrt(kappa_eta) * trans(cholReta) * randn(nG);
// //     
// //     thetamin = 0;
// //     thetamax = 2 * M_PI;
// //     theta = thetamin + randu() * (thetamax - thetamin);
// //     thetamin = theta - 2 * M_PI;
// //     thetamax = theta;
// //     
// //     llthred = llprev + log(randu());
// //     accept = false;
// //     count = 0;
// //     
// //     while(accept == false){
// //       count = count + 1;
// //       
// //       neweta = eta * cos(theta) + prioreta * sin(theta);
// //       neweta = neweta - mean(neweta);
// //       
// //       newlogZ = cppf_logZ_eta(Xzoop, alpha, neweta);
// //       newlogZobl = cppf_logZobl(newlogZ, indobl);
// //       newlogZsur = cppf_logZsur(newlogZ, indsur);
// //       newloglam = cppf_loglam_b0t_psit(Xwhale, newlogZ, beta, psi, indtwhale);
// //       newlam = exp(newloglam);
// //       
// //       llnew = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //         cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) +
// //         cppf_loglik_obl(logYobl, newlogZobl, sig2_obl) +
// //         cppf_loglik_sur(logYsur, newlogZsur, lam0, lam1, sig2_sur);
// //       
// //       if(llnew > llthred){
// //         llprev = llnew;
// //         accept = true;
// //       } else {
// //         if(theta < 0){ thetamin = theta; } else { thetamax = theta; }
// //         theta = thetamin + randu() * (thetamax - thetamin);
// //         if ( (count) % 1000 == 0 ) {
// //           Rprintf("ESS for eta: %d iterations...\n", count);
// //         }
// //         if (count > 20000) stop("ESS is not converging...\n");
// //       }
// //     }
// //     eta = neweta;
// //     logZ = newlogZ;
// //     logZobl = newlogZobl;
// //     logZsur = newlogZsur;
// //     loglam = newloglam;
// //     lam = newlam;
// //     
// //     postEta.row(s) = trans(eta);
// //     if(s == 0){ Rprintf("sampled eta\n"); }
// //     
// //     
// //     // update kappa_eta
// //     kappa_eta = cppf_gibbs_kappa_eta(eta, invReta, 2, 2);
// //     postKappa(s,0) = kappa_eta;
// //     if(s == 0){ Rprintf("sampled kappa_eta\n"); }
// //     
// //     
// //     // update lam0
// //     lam0 = cppf_gibbs_lam0(logYsur, logZsur, lam1, sig2_sur, 100);
// //     postLam(s, 0) = lam0;
// //     if(s == 0){ Rprintf("sampled lam0\n"); }
// //     
// //     
// //     // update lam1
// //     lam1 = cppf_gibbs_lam1(logYsur, logZsur, lam0, sig2_sur, 100);
// //     postLam(s, 1) = lam1;
// //     if(s == 0){ Rprintf("sampled lam1\n"); }
// //     
// //     
// //     // update sig2_sur
// //     sig2_sur = cppf_gibbs_sig2_sur(logYsur, logZsur, lam0, lam1, 2, 2);
// //     postSig2(s, 0) = sig2_sur;
// //     if(s == 0){ Rprintf("sampled sig2_sur\n"); }
// //     
// //     
// //     // update beta
// //     if( updateCOV ){
// //       if( (s+1 >= adaptInterval) && (s+1 - (adaptInterval * trunc((s+1) / adaptInterval)) == 0) ){
// //         dummyaccprob = accprob.col(1);
// //         rhat[1] = sum( dummyaccprob.rows(s+1-adaptInterval, s-1) ) / (adaptInterval-1);
// //         gamma1[1] = 1 / pow(adapIter[1], c1);
// //         gamma2[1] = c0 * gamma1[1];
// //         sigma2[1] = exp( log(sigma2[1]) + gamma2[1] * (rhat[1] - ropt) );
// //         
// //         COVbeta = COVbeta + gamma1[1] * ( cov( postBeta.rows(s+1-adaptInterval, s-1)) - COVbeta );
// //         cholCOVbeta = trans( chol( sigma2[1] * ( COVbeta + 0.0000000001 * diagmat(ones(n_be)) ) ) );
// //         
// //         adapIter[1] = adapIter[1] + 1;
// //       }
// //     }
// //     
// //     newbeta = beta + cholCOVbeta * randn(n_be);
// //     
// //     newloglam = cppf_loglam_b0t_psit(Xwhale, logZ, newbeta, psi, indtwhale);
// //     newlam = exp(newloglam);
// //     
// //     logprob = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) -
// //       cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //       cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) -
// //       cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
// //       MVN_logh(newbeta, zeros(n_be), diagmat(ones(n_be))/100) -
// //       MVN_logh(beta, zeros(n_be), diagmat(ones(n_be))/100);
// //     
// //     if (log(randu()) < logprob) {
// //       
// //       beta = newbeta;
// //       loglam = newloglam;
// //       lam = newlam;
// //       
// //       accprob(s, 1) = 1;
// //     }
// //     
// //     postBeta.row(s) = trans(beta);
// //     if(s == 0){ Rprintf("sampled beta\n"); }
// //     
// //     
// //     
// //     // update psi via Elliptical slice sampling
// //     for(int t = 0; t < ntwhale; t ++){
// //       newpsi = psi;
// //       
// //       llprev = cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //         cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea);
// //       priorpsit = sqrt(kappa_psi) * trans(cholRpsi) * randn(nG);
// //       
// //       thetamin = 0;
// //       thetamax = 2 * M_PI;
// //       theta = thetamin + randu() * (thetamax - thetamin);
// //       thetamin = theta - 2 * M_PI;
// //       thetamax = theta;
// //       
// //       llthred = llprev + log(randu());
// //       accept = false;
// //       count = 0;
// //       
// //       while(accept == false){
// //         count = count + 1;
// //         
// //         newpsi.col(t) = psi.col(t) * cos(theta) + priorpsit * sin(theta);
// //         newpsi.col(t) = newpsi.col(t) - mean(newpsi.col(t));
// //         
// //         newloglam = cppf_loglam_b0t_psit(Xwhale, logZ, beta, newpsi, indtwhale);
// //         newlam = exp(newloglam);
// //         
// //         llnew = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //           cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea);
// //         
// //         if(llnew > llthred){
// //           llprev = llnew;
// //           accept = true;
// //         } else {
// //           if(theta < 0){ thetamin = theta; } else { thetamax = theta; }
// //           theta = thetamin + randu() * (thetamax - thetamin);
// //           if ( (count) % 1000 == 0 ) {
// //             Rprintf("ESS for psi: %d iterations...\n", count);
// //           }
// //           if (count > 20000) stop("ESS is not converging...\n");
// //         }
// //       }
// //       psi = newpsi;
// //       loglam = newloglam;
// //       lam = newlam;
// //       
// //       postPsit = as<mat>(postPsi[t]);
// //       postPsit.row(s) = trans(psi.col(t));
// //       postPsi[t] = postPsit;
// //     }
// //     if(s == 0){ Rprintf("sampled psi\n"); }
// //     
// //     
// //     // update kappa_psi
// //     kappa_psi = cppf_gibbs_kappa_psit(psi, invRpsi, 2, 2);
// //     postKappa(s,1) = kappa_psi;
// //     if(s == 0){ Rprintf("sampled kappa_psi\n"); }
// //     
// //     
// //     // update c
// //     c = cppf_gibbs_ct(Wpam, lam, dur_dist, p_pam, cellarea, shape_c, scale_c);
// //     postC.row(s) = trans(c);
// //     if(s == 0){ Rprintf("sampled c\n"); }
// //     
// //     
// //     
// //     if ( (s+1) % 10 == 0 ) {
// //       Rprintf("Generated %d samples...\n", s+1);
// //     }
// //   }
// //   
// //   
// //   return Rcpp::List::create(Rcpp::Named("postAlpha") = postAlpha,
// //                             Rcpp::Named("postEta") = postEta,
// //                             Rcpp::Named("postLam") = postLam,
// //                             Rcpp::Named("postKappa") = postKappa,
// //                             Rcpp::Named("postSig2") = postSig2,
// //                             Rcpp::Named("postBeta") = postBeta,
// //                             Rcpp::Named("postC") = postC,
// //                             Rcpp::Named("postPsi") = postPsi,
// //                             Rcpp::Named("Accprob") = accprob,
// //                             Rcpp::Named("sigma2") = sigma2,
// //                             Rcpp::Named("adapIter") = adapIter,
// //                             Rcpp::Named("COValpha") = COValpha,
// //                             Rcpp::Named("COVbeta") = COVbeta);
// //   
// // }
// // 
// // 
// // // [[Rcpp::export]]
// // List cppf_fitJoint_a0t_etat_b0t_psi_ct_ja(
// //     int niter, mat distmat,
// //     List Xzoop, vec logYobl, vec logYsur,
// //     List Xwhale, mat Sdist, List Wpam, vec dur_dist,
// //     vec f_dist_obs, List f_dist, List p_pam,
// //     mat indobl, mat indsur, vec indtwhale, mat inddist, vec indmdist, vec indmwhale,
// //     double alpha0tilde, vec alpha, mat eta, double lam0, double lam1,
// //     double tau2, double kappa_eta, double phi_eta, double sig2_obl, double sig2_sur,
// //     vec beta, vec psi, vec pii, vec c, double shape_c, double scale_c,
// //     double kappa_psi, double phi_psi,
// //     double cellarea,
// //     vec sigma2, mat COValpha, mat COVbeta,
// //     bool updateCOV, int adaptInterval, double adaptFactorExponent, vec adapIter){
// //   
// //   double negativeInf = -std::numeric_limits<float>::infinity();;
// //   
// //   indobl = indobl-1, indsur = indsur-1, indtwhale = indtwhale-1;
// //   inddist = inddist-1, indmdist = indmdist-1, indmwhale = indmwhale-1;
// //   
// //   int nG = distmat.n_rows, n_al = alpha.size(), n_be = beta.size(), ntzoop = eta.n_cols, n_sigma2 = sigma2.size();
// //   
// //   mat postAlpha0tilde = zeros(niter, 1), postAlpha = zeros(niter, n_al), postLam = zeros(niter, 2);
// //   mat postSig2 = zeros(niter, 1), postKappa = zeros(niter, 2);
// //   mat postBeta = zeros(niter, n_be), postC = zeros(niter, c.size());
// //   mat postPsi(niter, nG), postEtat = zeros(niter, nG);
// //   List postEta(ntzoop);
// //   for(int i = 0; i < ntzoop; i ++){
// //     postEta[i] = zeros(niter, nG);
// //   }
// //   Rprintf("set up matrices for posterior samples\n");
// //   
// //   mat accprob = zeros(niter, n_sigma2);
// //   vec rhat = zeros(n_sigma2), gamma1 = zeros(n_sigma2), gamma2 = zeros(n_sigma2), dummyaccprob;
// //   double c0 = 1, c1 = adaptFactorExponent, ropt = 0.234, logprob;
// //   mat cholCOValpha = trans( chol( sigma2[0] * ( COValpha + 0.0000000001 * diagmat(ones(n_al)) ) ) );
// //   mat cholCOVbeta = trans( chol( sigma2[1] * ( COVbeta + 0.0000000001 * diagmat(ones(n_be)) ) ) );
// //   Rprintf("set up adaptive proposals\n");
// //   
// //   vec newalpha, newbeta, newlogZobl, newlogZsur;
// //   mat newlogZ, newloglam, newlam;
// //   vec dummyvec;
// //   
// //   int count;
// //   vec xvals;
// //   double llprev, llthred, llnew, thetamin, thetamax, theta;
// //   vec prioretat, priorpsi, newpsi;
// //   mat neweta;
// //   bool accept;
// //   
// //   mat Reta = exp(-distmat / phi_eta), Rpsi = exp(-distmat / phi_psi);
// //   mat invReta = inv(Reta), invRpsi= inv(Rpsi), cholReta = chol(Reta), cholRpsi= chol(Rpsi);
// //   rowvec oneinvReta = trans(ones(nG)) * invReta, oneinvRpsi = trans(ones(nG)) * invRpsi;
// //   double oneinvRetaone = sum( oneinvReta * ones(nG) ), oneinvRpsione = sum( oneinvRpsi * ones(nG) );
// //   Rprintf("set up GP covariance matrices\n");
// //   
// //   mat logZ = cppf_logZ_a0t_etat_ja(Xzoop, alpha, eta);
// //   vec logZobl = cppf_logZobl(logZ, indobl);
// //   vec logZsur = cppf_logZsur(logZ, indsur);
// //   mat loglam = cppf_loglam_b0t_psi(Xwhale, logZ, beta, psi, indtwhale);
// //   mat lam = exp(loglam);
// //   Rprintf("set up initial values\n");
// //   
// //   
// //   // Start MCMC
// //   for(int s = 0; s < niter; s++) {
// //     
// //     // update alpha
// //     if( updateCOV ){
// //       if( (s+1 >= adaptInterval) && (s+1 - (adaptInterval * trunc((s+1) / adaptInterval)) == 0) ){
// //         dummyaccprob = accprob.col(0);
// //         rhat[0] = sum( dummyaccprob.rows(s+1-adaptInterval, s-1) ) / (adaptInterval-1);
// //         gamma1[0] = 1 / pow(adapIter[0], c1);
// //         gamma2[0] = c0 * gamma1[0];
// //         sigma2[0] = exp( log(sigma2[0]) + gamma2[0] * (rhat[0] - ropt) );
// //         
// //         COValpha = COValpha + gamma1[0] * ( cov( postAlpha.rows(s+1-adaptInterval, s-1) ) - COValpha );
// //         cholCOValpha = trans( chol( sigma2[0] * ( COValpha + 0.0000000001 * diagmat(ones(n_al)) ) ) );
// //         
// //         adapIter[0] = adapIter[0] + 1;
// //       }
// //     }
// //     
// //     newalpha = alpha + cholCOValpha * randn(n_al);
// //     
// //     newlogZ = cppf_logZ_a0t_etat_ja(Xzoop, newalpha, eta);
// //     newlogZobl = cppf_logZobl(newlogZ, indobl);
// //     newlogZsur = cppf_logZsur(newlogZ, indsur);
// //     newloglam = cppf_loglam_b0t_psi(Xwhale, newlogZ, beta, psi, indtwhale);
// //     newlam = exp(newloglam);
// //     
// //     logprob = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) -
// //       cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //       cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) -
// //       cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
// //       cppf_loglik_obl(logYobl, newlogZobl, sig2_obl) -
// //       cppf_loglik_obl(logYobl, logZobl, sig2_obl) +
// //       cppf_loglik_sur(logYsur, newlogZsur, lam0, lam1, sig2_sur) -
// //       cppf_loglik_sur(logYsur, logZsur, lam0, lam1, sig2_sur) +
// //       MVN_logh(newalpha, zeros(n_al), diagmat(ones(n_al))/100) -
// //       MVN_logh(alpha, zeros(n_al), diagmat(ones(n_al))/100);
// //     
// //     if (log(randu()) < logprob) {
// //       
// //       alpha = newalpha;
// //       logZ = newlogZ;
// //       logZobl = newlogZobl;
// //       logZsur = newlogZsur;
// //       loglam = newloglam;
// //       lam = newlam;
// //       
// //       accprob(s, 0) = 1;
// //     }
// //     postAlpha.row(s) = trans(alpha);
// //     if(s == 0){ Rprintf("sampled alpha\n"); }
// //     
// //     
// //     // Gibbs update for alpha0tilde
// //     alpha0tilde = cppf_gibbs_alpha0tilde(alpha.rows(0, ntzoop-1), tau2);
// //     postAlpha0tilde[s] = alpha0tilde;
// //     if(s == 0){ Rprintf("sampled alpha0tilde\n"); }
// //     
// //     
// //     // update eta via Elliptical slice sampling
// //     for(int t = 0; t < ntzoop; t ++){
// //       neweta = eta;
// //       
// //       llprev = cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //         cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
// //         cppf_loglik_obl(logYobl, logZobl, sig2_obl) +
// //         cppf_loglik_sur(logYsur, logZsur, lam0, lam1, sig2_sur);
// //       prioretat = sqrt(kappa_eta) * trans(cholReta) * randn(nG);
// //       
// //       thetamin = 0;
// //       thetamax = 2 * M_PI;
// //       theta = thetamin + randu() * (thetamax - thetamin);
// //       thetamin = theta - 2 * M_PI;
// //       thetamax = theta;
// //       
// //       llthred = llprev + log(randu());
// //       accept = false;
// //       count = 0;
// //       
// //       while(accept == false){
// //         count = count + 1;
// //         
// //         neweta.col(t) = eta.col(t) * cos(theta) + prioretat * sin(theta);
// //         neweta.col(t) = neweta.col(t) - mean(neweta.col(t));
// //         
// //         newlogZ = cppf_logZ_a0t_etat_ja(Xzoop, alpha, neweta);
// //         newlogZobl = cppf_logZobl(newlogZ, indobl);
// //         newlogZsur = cppf_logZsur(newlogZ, indsur);
// //         newloglam = cppf_loglam_b0t_psi(Xwhale, newlogZ, beta, psi, indtwhale);
// //         newlam = exp(newloglam);
// //         
// //         llnew = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //           cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) +
// //           cppf_loglik_obl(logYobl, newlogZobl, sig2_obl) +
// //           cppf_loglik_sur(logYsur, newlogZsur, lam0, lam1, sig2_sur);
// //         
// //         if(llnew > llthred){
// //           llprev = llnew;
// //           accept = true;
// //         } else {
// //           if(theta < 0){ thetamin = theta; } else { thetamax = theta; }
// //           theta = thetamin + randu() * (thetamax - thetamin);
// //           if ( (count) % 1000 == 0 ) {
// //             Rprintf("ESS for eta: %d iterations...\n", count);
// //           }
// //           if (count > 20000) stop("ESS is not converging...\n");
// //         }
// //       }
// //       eta = neweta;
// //       logZ = newlogZ;
// //       logZobl = newlogZobl;
// //       logZsur = newlogZsur;
// //       loglam = newloglam;
// //       lam = newlam;
// //       
// //       postEtat = as<mat>(postEta[t]);
// //       postEtat.row(s) = trans(eta.col(t));
// //       postEta[t] = postEtat;
// //     }
// //     if(s == 0){ Rprintf("sampled eta\n"); }
// //     
// //     
// //     // update kappa_eta
// //     kappa_eta = cppf_gibbs_kappa_etat(eta, invReta, 2, 2);
// //     postKappa(s,0) = kappa_eta;
// //     if(s == 0){ Rprintf("sampled kappa_eta\n"); }
// //     
// //     
// //     // update lam0
// //     lam0 = cppf_gibbs_lam0(logYsur, logZsur, lam1, sig2_sur, 100);
// //     postLam(s, 0) = lam0;
// //     if(s == 0){ Rprintf("sampled lam0\n"); }
// //     
// //     
// //     // update lam1
// //     lam1 = cppf_gibbs_lam1(logYsur, logZsur, lam0, sig2_sur, 100);
// //     postLam(s, 1) = lam1;
// //     if(s == 0){ Rprintf("sampled lam1\n"); }
// //     
// //     
// //     // update sig2_sur
// //     sig2_sur = cppf_gibbs_sig2_sur(logYsur, logZsur, lam0, lam1, 2, 2);
// //     postSig2(s, 0) = sig2_sur;
// //     if(s == 0){ Rprintf("sampled sig2_sur\n"); }
// //     
// //     
// //     // update beta
// //     if( updateCOV ){
// //       if( (s+1 >= adaptInterval) && (s+1 - (adaptInterval * trunc((s+1) / adaptInterval)) == 0) ){
// //         dummyaccprob = accprob.col(1);
// //         rhat[1] = sum( dummyaccprob.rows(s+1-adaptInterval, s-1) ) / (adaptInterval-1);
// //         gamma1[1] = 1 / pow(adapIter[1], c1);
// //         gamma2[1] = c0 * gamma1[1];
// //         sigma2[1] = exp( log(sigma2[1]) + gamma2[1] * (rhat[1] - ropt) );
// //         
// //         COVbeta = COVbeta + gamma1[1] * ( cov( postBeta.rows(s+1-adaptInterval, s-1)) - COVbeta );
// //         cholCOVbeta = trans( chol( sigma2[1] * ( COVbeta + 0.0000000001 * diagmat(ones(n_be)) ) ) );
// //         
// //         adapIter[1] = adapIter[1] + 1;
// //       }
// //     }
// //     
// //     newbeta = beta + cholCOVbeta * randn(n_be);
// //     
// //     newloglam = cppf_loglam_b0t_psi(Xwhale, logZ, newbeta, psi, indtwhale);
// //     newlam = exp(newloglam);
// //     
// //     logprob = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) -
// //       cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //       cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) -
// //       cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
// //       MVN_logh(newbeta, zeros(n_be), diagmat(ones(n_be))/100) -
// //       MVN_logh(beta, zeros(n_be), diagmat(ones(n_be))/100);
// //     
// //     if (log(randu()) < logprob) {
// //       
// //       beta = newbeta;
// //       loglam = newloglam;
// //       lam = newlam;
// //       
// //       accprob(s, 1) = 1;
// //     }
// //     
// //     postBeta.row(s) = trans(beta);
// //     if(s == 0){ Rprintf("sampled beta\n"); }
// //     
// //     
// //     // update psi via elliptical slice sampling
// //     llprev = cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //       cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea);
// //     priorpsi = sqrt(kappa_psi) * trans(cholRpsi) * randn(nG);
// //     
// //     thetamin = 0;
// //     thetamax = 2 * M_PI;
// //     theta = thetamin + randu() * (thetamax - thetamin);
// //     thetamin = theta - 2 * M_PI;
// //     thetamax = theta;
// //     
// //     llthred = llprev + log(randu());
// //     accept = false;
// //     count = 0;
// //     
// //     while(accept == false){
// //       count = count + 1;
// //       
// //       newpsi = psi * cos(theta) + priorpsi * sin(theta);
// //       newpsi = newpsi - mean(newpsi);
// //       
// //       newloglam = cppf_loglam_b0t_psi(Xwhale, logZ, beta, newpsi, indtwhale);
// //       newlam = exp(newloglam);
// //       
// //       llnew = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //         cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea);
// //       
// //       if(llnew > llthred){
// //         llprev = llnew;
// //         accept = true;
// //       } else {
// //         if(theta < 0){ thetamin = theta; } else { thetamax = theta; }
// //         theta = thetamin + randu() * (thetamax - thetamin);
// //         if ( (count) % 1000 == 0 ) {
// //           Rprintf("ESS for psi: %d iterations...\n", count);
// //         }
// //         if (count > 20000) stop("ESS is not converging...\n");
// //       }
// //     }
// //     psi = newpsi;
// //     loglam = newloglam;
// //     lam = newlam;
// //     
// //     postPsi.row(s) = trans(psi);
// //     if(s == 0){ Rprintf("sampled psi\n"); }
// //     
// //     
// //     // update kappa_psi
// //     kappa_psi = cppf_gibbs_kappa_psi(psi, invRpsi, 2, 2);
// //     postKappa(s,1) = kappa_psi;
// //     if(s == 0){ Rprintf("sampled kappa_psi\n"); }
// //     
// //     
// //     // update c
// //     c = cppf_gibbs_ct(Wpam, lam, dur_dist, p_pam, cellarea, shape_c, scale_c);
// //     postC.row(s) = trans(c);
// //     if(s == 0){ Rprintf("sampled c\n"); }
// //     
// //     
// //     if ( (s+1) % 10 == 0 ) {
// //       Rprintf("Generated %d samples...\n", s+1);
// //     }
// //   }
// //   
// //   
// //   return Rcpp::List::create(Rcpp::Named("postAlpha0tilde") = postAlpha0tilde,
// //                             Rcpp::Named("postAlpha") = postAlpha,
// //                             Rcpp::Named("postEta") = postEta,
// //                             Rcpp::Named("postLam") = postLam,
// //                             Rcpp::Named("postKappa") = postKappa,
// //                             Rcpp::Named("postSig2") = postSig2,
// //                             Rcpp::Named("postBeta") = postBeta,
// //                             Rcpp::Named("postC") = postC,
// //                             Rcpp::Named("postPsi") = postPsi,
// //                             Rcpp::Named("Accprob") = accprob,
// //                             Rcpp::Named("sigma2") = sigma2,
// //                             Rcpp::Named("adapIter") = adapIter,
// //                             Rcpp::Named("COValpha") = COValpha,
// //                             Rcpp::Named("COVbeta") = COVbeta);
// //   
// // }
// // 
// // 
// // 
// // // ============================================================================-
// // // estimate call rates and centering zoop covariate ----
// // // ============================================================================-
// // 
// // 
// // // [[Rcpp::export]]
// // List cppf_fitJoint_cZ_eta_b0t_psi_ct(
// //     int niter, mat distmat,
// //     List Xzoop, vec logYobl, vec logYsur,
// //     List Xwhale, mat Sdist, List Wpam, vec dur_dist,
// //     vec f_dist_obs, List f_dist, List p_pam,
// //     mat indobl, mat indsur, vec indtwhale, mat inddist, vec indmdist, vec indmwhale,
// //     vec alpha, vec eta, double lam0, double lam1,
// //     double kappa_eta, double phi_eta, double sig2_obl, double sig2_sur,
// //     vec beta, vec psi, vec pii, vec c, double shape_c, double scale_c,
// //     double kappa_psi, double phi_psi,
// //     double cellarea,
// //     vec sigma2, mat COValpha, mat COVbeta,
// //     bool updateCOV, int adaptInterval, double adaptFactorExponent, vec adapIter){
// //   
// //   double negativeInf = -std::numeric_limits<float>::infinity();;
// //   
// //   indobl = indobl-1, indsur = indsur-1, indtwhale = indtwhale-1;
// //   inddist = inddist-1, indmdist = indmdist-1, indmwhale = indmwhale-1;
// //   
// //   int nG = distmat.n_rows, n_al = alpha.size(), n_be = beta.size(), ntzoop = eta.n_cols, n_sigma2 = sigma2.size();
// //   
// //   mat postAlpha = zeros(niter, n_al), postLam = zeros(niter, 2);
// //   mat postSig2 = zeros(niter, 1), postKappa = zeros(niter, 2);
// //   mat postBeta = zeros(niter, n_be), postEta = zeros(niter, nG);
// //   mat postPsi = zeros(niter, nG), postC = zeros(niter, c.size());
// //   Rprintf("set up matrices for posterior samples\n");
// //   
// //   mat accprob = zeros(niter, n_sigma2);
// //   vec rhat = zeros(n_sigma2), gamma1 = zeros(n_sigma2), gamma2 = zeros(n_sigma2), dummyaccprob;
// //   double c0 = 1, c1 = adaptFactorExponent, ropt = 0.234, logprob;
// //   mat cholCOValpha = trans( chol( sigma2[0] * ( COValpha + 0.0000000001 * diagmat(ones(n_al)) ) ) );
// //   mat cholCOVbeta = trans( chol( sigma2[1] * ( COVbeta + 0.0000000001 * diagmat(ones(n_be)) ) ) );
// //   Rprintf("set up adaptive proposals\n");
// //   
// //   vec newalpha, newbeta, newlogZobl, newlogZsur;
// //   mat newlogZ, newloglam, newlam;
// //   vec dummyvec;
// //   
// //   int count;
// //   vec xvals;
// //   double llprev, llthred, llnew, thetamin, thetamax, theta;
// //   vec prioreta, neweta, priorpsi, newpsi;
// //   bool accept;
// //   
// //   mat Reta = exp(-distmat / phi_eta), Rpsi = exp(-distmat / phi_psi);
// //   mat invReta = inv(Reta), invRpsi= inv(Rpsi), cholReta = chol(Reta), cholRpsi= chol(Rpsi);
// //   rowvec oneinvReta = trans(ones(nG)) * invReta, oneinvRpsi = trans(ones(nG)) * invRpsi;
// //   double oneinvRetaone = sum( oneinvReta * ones(nG) ), oneinvRpsione = sum( oneinvRpsi * ones(nG) );
// //   Rprintf("set up GP covariance matrices\n");
// //   
// //   mat logZ = cppf_logZ_eta(Xzoop, alpha, eta);
// //   vec logZobl = cppf_logZobl(logZ, indobl);
// //   vec logZsur = cppf_logZsur(logZ, indsur);
// //   mat loglam = cppf_loglam_cZ_b0t_psi(Xwhale, logZ, beta, psi, indtwhale);
// //   mat lam = exp(loglam);
// //   vec m_al = zeros(n_al), diag_al = ones(n_al)/100;
// //   Rprintf("set up initial values\n");
// //   
// //   
// //   // Start MCMC
// //   for(int s = 0; s < niter; s++) {
// //     
// //     // update alpha
// //     if( updateCOV ){
// //       if( (s+1 >= adaptInterval) && (s+1 - (adaptInterval * trunc((s+1) / adaptInterval)) == 0) ){
// //         dummyaccprob = accprob.col(0);
// //         rhat[0] = sum( dummyaccprob.rows(s+1-adaptInterval, s-1) ) / (adaptInterval-1);
// //         gamma1[0] = 1 / pow(adapIter[0], c1);
// //         gamma2[0] = c0 * gamma1[0];
// //         sigma2[0] = exp( log(sigma2[0]) + gamma2[0] * (rhat[0] - ropt) );
// //         
// //         COValpha = COValpha + gamma1[0] * ( cov( postAlpha.rows(s+1-adaptInterval, s-1) ) - COValpha );
// //         cholCOValpha = trans( chol( sigma2[0] * ( COValpha + 0.0000000001 * diagmat(ones(n_al)) ) ) );
// //         
// //         adapIter[0] = adapIter[0] + 1;
// //       }
// //     }
// //     
// //     newalpha = alpha + cholCOValpha * randn(n_al);
// //     
// //     newlogZ = cppf_logZ_eta(Xzoop, newalpha, eta);
// //     newlogZobl = cppf_logZobl(newlogZ, indobl);
// //     newlogZsur = cppf_logZsur(newlogZ, indsur);
// //     newloglam = cppf_loglam_cZ_b0t_psi(Xwhale, newlogZ, beta, psi, indtwhale);
// //     newlam = exp(newloglam);
// //     
// //     logprob = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) -
// //       cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //       cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) -
// //       cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
// //       cppf_loglik_obl(logYobl, newlogZobl, sig2_obl) -
// //       cppf_loglik_obl(logYobl, logZobl, sig2_obl) +
// //       cppf_loglik_sur(logYsur, newlogZsur, lam0, lam1, sig2_sur) -
// //       cppf_loglik_sur(logYsur, logZsur, lam0, lam1, sig2_sur) +
// //       MVN_logh(newalpha, zeros(n_al), diagmat(ones(n_al))/100) -
// //       MVN_logh(alpha, zeros(n_al), diagmat(ones(n_al))/100);
// //     
// //     if (log(randu()) < logprob) {
// //       
// //       alpha = newalpha;
// //       logZ = newlogZ;
// //       logZobl = newlogZobl;
// //       logZsur = newlogZsur;
// //       loglam = newloglam;
// //       lam = newlam;
// //       
// //       accprob(s, 0) = 1;
// //     }
// //     postAlpha.row(s) = trans(alpha);
// //     if(s == 0){ Rprintf("sampled alpha\n"); }
// //     
// //     
// //     
// //     // update eta via Elliptical slice sampling
// //     llprev = cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //       cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
// //       cppf_loglik_obl(logYobl, logZobl, sig2_obl) +
// //       cppf_loglik_sur(logYsur, logZsur, lam0, lam1, sig2_sur);
// //     prioreta = sqrt(kappa_eta) * trans(cholReta) * randn(nG);
// //     
// //     thetamin = 0;
// //     thetamax = 2 * M_PI;
// //     theta = thetamin + randu() * (thetamax - thetamin);
// //     thetamin = theta - 2 * M_PI;
// //     thetamax = theta;
// //     
// //     llthred = llprev + log(randu());
// //     accept = false;
// //     count = 0;
// //     
// //     while(accept == false){
// //       count = count + 1;
// //       
// //       neweta = eta * cos(theta) + prioreta * sin(theta);
// //       neweta = neweta - mean(neweta);
// //       
// //       newlogZ = cppf_logZ_eta(Xzoop, alpha, neweta);
// //       newlogZobl = cppf_logZobl(newlogZ, indobl);
// //       newlogZsur = cppf_logZsur(newlogZ, indsur);
// //       newloglam = cppf_loglam_cZ_b0t_psi(Xwhale, newlogZ, beta, psi, indtwhale);
// //       newlam = exp(newloglam);
// //       
// //       llnew = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //         cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) +
// //         cppf_loglik_obl(logYobl, newlogZobl, sig2_obl) +
// //         cppf_loglik_sur(logYsur, newlogZsur, lam0, lam1, sig2_sur);
// //       
// //       if(llnew > llthred){
// //         llprev = llnew;
// //         accept = true;
// //       } else {
// //         if(theta < 0){ thetamin = theta; } else { thetamax = theta; }
// //         theta = thetamin + randu() * (thetamax - thetamin);
// //         if ( (count) % 1000 == 0 ) {
// //           Rprintf("ESS for eta: %d iterations...\n", count);
// //         }
// //         if (count > 20000) stop("ESS is not converging...\n");
// //       }
// //     }
// //     eta = neweta;
// //     logZ = newlogZ;
// //     logZobl = newlogZobl;
// //     logZsur = newlogZsur;
// //     loglam = newloglam;
// //     lam = newlam;
// //     
// //     postEta.row(s) = trans(eta);
// //     if(s == 0){ Rprintf("sampled eta\n"); }
// //     
// //     
// //     // update kappa_eta
// //     kappa_eta = cppf_gibbs_kappa_eta(eta, invReta, 2, 2);
// //     postKappa(s,0) = kappa_eta;
// //     if(s == 0){ Rprintf("sampled kappa_eta\n"); }
// //     
// //     
// //     // update lam0
// //     lam0 = cppf_gibbs_lam0(logYsur, logZsur, lam1, sig2_sur, 100);
// //     postLam(s, 0) = lam0;
// //     if(s == 0){ Rprintf("sampled lam0\n"); }
// //     
// //     
// //     // update lam1
// //     lam1 = cppf_gibbs_lam1(logYsur, logZsur, lam0, sig2_sur, 100);
// //     postLam(s, 1) = lam1;
// //     if(s == 0){ Rprintf("sampled lam1\n"); }
// //     
// //     
// //     // update sig2_sur
// //     sig2_sur = cppf_gibbs_sig2_sur(logYsur, logZsur, lam0, lam1, 2, 2);
// //     postSig2(s, 0) = sig2_sur;
// //     if(s == 0){ Rprintf("sampled sig2_sur\n"); }
// //     
// //     
// //     // update beta
// //     if( updateCOV ){
// //       if( (s+1 >= adaptInterval) && (s+1 - (adaptInterval * trunc((s+1) / adaptInterval)) == 0) ){
// //         dummyaccprob = accprob.col(1);
// //         rhat[1] = sum( dummyaccprob.rows(s+1-adaptInterval, s-1) ) / (adaptInterval-1);
// //         gamma1[1] = 1 / pow(adapIter[1], c1);
// //         gamma2[1] = c0 * gamma1[1];
// //         sigma2[1] = exp( log(sigma2[1]) + gamma2[1] * (rhat[1] - ropt) );
// //         
// //         COVbeta = COVbeta + gamma1[1] * ( cov( postBeta.rows(s+1-adaptInterval, s-1)) - COVbeta );
// //         cholCOVbeta = trans( chol( sigma2[1] * ( COVbeta + 0.0000000001 * diagmat(ones(n_be)) ) ) );
// //         
// //         adapIter[1] = adapIter[1] + 1;
// //       }
// //     }
// //     
// //     newbeta = beta + cholCOVbeta * randn(n_be);
// //     
// //     newloglam = cppf_loglam_cZ_b0t_psi(Xwhale, logZ, newbeta, psi, indtwhale);
// //     newlam = exp(newloglam);
// //     
// //     logprob = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) -
// //       cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //       cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) -
// //       cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
// //       MVN_logh(newbeta, zeros(n_be), diagmat(ones(n_be))/100) -
// //       MVN_logh(beta, zeros(n_be), diagmat(ones(n_be))/100);
// //     
// //     if (log(randu()) < logprob) {
// //       
// //       beta = newbeta;
// //       loglam = newloglam;
// //       lam = newlam;
// //       
// //       accprob(s, 1) = 1;
// //     }
// //     
// //     postBeta.row(s) = trans(beta);
// //     if(s == 0){ Rprintf("sampled beta\n"); }
// //     
// //     
// //     
// //     // update psi via elliptical slice sampling
// //     llprev = cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //       cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea);
// //     priorpsi = sqrt(kappa_psi) * trans(cholRpsi) * randn(nG);
// //     
// //     thetamin = 0;
// //     thetamax = 2 * M_PI;
// //     theta = thetamin + randu() * (thetamax - thetamin);
// //     thetamin = theta - 2 * M_PI;
// //     thetamax = theta;
// //     
// //     llthred = llprev + log(randu());
// //     accept = false;
// //     count = 0;
// //     
// //     while(accept == false){
// //       count = count + 1;
// //       
// //       newpsi = psi * cos(theta) + priorpsi * sin(theta);
// //       newpsi = newpsi - mean(newpsi);
// //       
// //       newloglam = cppf_loglam_cZ_b0t_psi(Xwhale, logZ, beta, newpsi, indtwhale);
// //       newlam = exp(newloglam);
// //       
// //       llnew = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //         cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea);
// //       
// //       if(llnew > llthred){
// //         llprev = llnew;
// //         accept = true;
// //       } else {
// //         if(theta < 0){ thetamin = theta; } else { thetamax = theta; }
// //         theta = thetamin + randu() * (thetamax - thetamin);
// //         if ( (count) % 1000 == 0 ) {
// //           Rprintf("ESS for psi: %d iterations...\n", count);
// //         }
// //         if (count > 20000) stop("ESS is not converging...\n");
// //       }
// //     }
// //     psi = newpsi;
// //     loglam = newloglam;
// //     lam = newlam;
// //     
// //     postPsi.row(s) = trans(psi);
// //     if(s == 0){ Rprintf("sampled psi\n"); }
// //     
// //     
// //     // update kappa_psi
// //     kappa_psi = cppf_gibbs_kappa_psi(psi, invRpsi, 2, 2);
// //     postKappa(s,1) = kappa_psi;
// //     if(s == 0){ Rprintf("sampled kappa_psi\n"); }
// //     
// //     
// //     // update c
// //     c = cppf_gibbs_ct(Wpam, lam, dur_dist, p_pam, cellarea, shape_c, scale_c);
// //     postC.row(s) = trans(c);
// //     if(s == 0){ Rprintf("sampled c\n"); }
// //     
// //     
// //     if ( (s+1) % 10 == 0 ) {
// //       Rprintf("Generated %d samples...\n", s+1);
// //     }
// //   }
// //   
// //   
// //   return Rcpp::List::create(Rcpp::Named("postAlpha") = postAlpha,
// //                             Rcpp::Named("postEta") = postEta,
// //                             Rcpp::Named("postLam") = postLam,
// //                             Rcpp::Named("postKappa") = postKappa,
// //                             Rcpp::Named("postSig2") = postSig2,
// //                             Rcpp::Named("postBeta") = postBeta,
// //                             Rcpp::Named("postC") = postC,
// //                             Rcpp::Named("postPsi") = postPsi,
// //                             Rcpp::Named("Accprob") = accprob,
// //                             Rcpp::Named("sigma2") = sigma2,
// //                             Rcpp::Named("adapIter") = adapIter,
// //                             Rcpp::Named("COValpha") = COValpha,
// //                             Rcpp::Named("COVbeta") = COVbeta);
// //   
// // }
// // 
// // 
// // 
// // // [[Rcpp::export]]
// // List cppf_fitJoint_cZ_eta_b0t_psit_ct(
// //     int niter, mat distmat,
// //     List Xzoop, vec logYobl, vec logYsur,
// //     List Xwhale, mat Sdist, List Wpam, vec dur_dist,
// //     vec f_dist_obs, List f_dist, List p_pam,
// //     mat indobl, mat indsur, vec indtwhale, mat inddist, vec indmdist, vec indmwhale,
// //     vec alpha, vec eta, double lam0, double lam1,
// //     double kappa_eta, double phi_eta, double sig2_obl, double sig2_sur,
// //     vec beta, mat psi, vec pii, vec c, double shape_c, double scale_c,
// //     double kappa_psi, double phi_psi,
// //     double cellarea,
// //     vec sigma2, mat COValpha, mat COVbeta,
// //     bool updateCOV, int adaptInterval, double adaptFactorExponent, vec adapIter){
// //   
// //   double negativeInf = -std::numeric_limits<float>::infinity();;
// //   
// //   indobl = indobl-1, indsur = indsur-1, indtwhale = indtwhale-1;
// //   inddist = inddist-1, indmdist = indmdist-1, indmwhale = indmwhale-1;
// //   
// //   int nG = distmat.n_rows, n_al = alpha.size(), n_be = beta.size(), ntzoop = eta.n_cols, ntwhale = psi.n_cols, n_sigma2 = sigma2.size();
// //   
// //   mat postAlpha = zeros(niter, n_al), postLam = zeros(niter, 2);
// //   mat postSig2 = zeros(niter, 1), postKappa = zeros(niter, 2);
// //   mat postBeta = zeros(niter, n_be), postEta = zeros(niter, nG);
// //   mat postPsit = zeros(niter, nG), postC = zeros(niter, c.size());
// //   List postPsi(ntwhale);
// //   for(int i = 0; i < ntwhale; i ++){ postPsi[i] = zeros(niter, nG); }
// //   Rprintf("set up matrices for posterior samples\n");
// //   
// //   mat accprob = zeros(niter, n_sigma2);
// //   vec rhat = zeros(n_sigma2), gamma1 = zeros(n_sigma2), gamma2 = zeros(n_sigma2), dummyaccprob;
// //   double c0 = 1, c1 = adaptFactorExponent, ropt = 0.234, logprob;
// //   mat cholCOValpha = trans( chol( sigma2[0] * ( COValpha + 0.0000000001 * diagmat(ones(n_al)) ) ) );
// //   mat cholCOVbeta = trans( chol( sigma2[1] * ( COVbeta + 0.0000000001 * diagmat(ones(n_be)) ) ) );
// //   Rprintf("set up adaptive proposals\n");
// //   
// //   vec newalpha, newbeta, newlogZobl, newlogZsur;
// //   mat newlogZ, newloglam, newlam;
// //   vec dummyvec;
// //   
// //   int count;
// //   vec xvals;
// //   double llprev, llthred, llnew, thetamin, thetamax, theta;
// //   vec prioreta, neweta, priorpsit;
// //   mat newpsi;
// //   bool accept;
// //   
// //   mat Reta = exp(-distmat / phi_eta), Rpsi = exp(-distmat / phi_psi);
// //   mat invReta = inv(Reta), invRpsi= inv(Rpsi), cholReta = chol(Reta), cholRpsi= chol(Rpsi);
// //   rowvec oneinvReta = trans(ones(nG)) * invReta, oneinvRpsi = trans(ones(nG)) * invRpsi;
// //   double oneinvRetaone = sum( oneinvReta * ones(nG) ), oneinvRpsione = sum( oneinvRpsi * ones(nG) );
// //   Rprintf("set up GP covariance matrices\n");
// //   
// //   mat logZ = cppf_logZ_eta(Xzoop, alpha, eta);
// //   vec logZobl = cppf_logZobl(logZ, indobl);
// //   vec logZsur = cppf_logZsur(logZ, indsur);
// //   mat loglam = cppf_loglam_cZ_b0t_psit(Xwhale, logZ, beta, psi, indtwhale);
// //   mat lam = exp(loglam);
// //   vec m_al = zeros(n_al), diag_al = ones(n_al)/100;
// //   Rprintf("set up initial values\n");
// //   
// //   
// //   // Start MCMC
// //   for(int s = 0; s < niter; s++) {
// //     
// //     // update alpha
// //     if( updateCOV ){
// //       if( (s+1 >= adaptInterval) && (s+1 - (adaptInterval * trunc((s+1) / adaptInterval)) == 0) ){
// //         dummyaccprob = accprob.col(0);
// //         rhat[0] = sum( dummyaccprob.rows(s+1-adaptInterval, s-1) ) / (adaptInterval-1);
// //         gamma1[0] = 1 / pow(adapIter[0], c1);
// //         gamma2[0] = c0 * gamma1[0];
// //         sigma2[0] = exp( log(sigma2[0]) + gamma2[0] * (rhat[0] - ropt) );
// //         
// //         COValpha = COValpha + gamma1[0] * ( cov( postAlpha.rows(s+1-adaptInterval, s-1) ) - COValpha );
// //         cholCOValpha = trans( chol( sigma2[0] * ( COValpha + 0.0000000001 * diagmat(ones(n_al)) ) ) );
// //         
// //         adapIter[0] = adapIter[0] + 1;
// //       }
// //     }
// //     
// //     newalpha = alpha + cholCOValpha * randn(n_al);
// //     
// //     newlogZ = cppf_logZ_eta(Xzoop, newalpha, eta);
// //     newlogZobl = cppf_logZobl(newlogZ, indobl);
// //     newlogZsur = cppf_logZsur(newlogZ, indsur);
// //     newloglam = cppf_loglam_cZ_b0t_psit(Xwhale, newlogZ, beta, psi, indtwhale);
// //     newlam = exp(newloglam);
// //     
// //     logprob = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) -
// //       cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //       cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) -
// //       cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
// //       cppf_loglik_obl(logYobl, newlogZobl, sig2_obl) -
// //       cppf_loglik_obl(logYobl, logZobl, sig2_obl) +
// //       cppf_loglik_sur(logYsur, newlogZsur, lam0, lam1, sig2_sur) -
// //       cppf_loglik_sur(logYsur, logZsur, lam0, lam1, sig2_sur) +
// //       MVN_logh(newalpha, zeros(n_al), diagmat(ones(n_al))/100) -
// //       MVN_logh(alpha, zeros(n_al), diagmat(ones(n_al))/100);
// //     
// //     if (log(randu()) < logprob) {
// //       
// //       alpha = newalpha;
// //       logZ = newlogZ;
// //       logZobl = newlogZobl;
// //       logZsur = newlogZsur;
// //       loglam = newloglam;
// //       lam = newlam;
// //       
// //       accprob(s, 0) = 1;
// //     }
// //     postAlpha.row(s) = trans(alpha);
// //     if(s == 0){ Rprintf("sampled alpha\n"); }
// //     
// //     
// //     
// //     // update eta via Elliptical slice sampling
// //     llprev = cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //       cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
// //       cppf_loglik_obl(logYobl, logZobl, sig2_obl) +
// //       cppf_loglik_sur(logYsur, logZsur, lam0, lam1, sig2_sur);
// //     prioreta = sqrt(kappa_eta) * trans(cholReta) * randn(nG);
// //     
// //     thetamin = 0;
// //     thetamax = 2 * M_PI;
// //     theta = thetamin + randu() * (thetamax - thetamin);
// //     thetamin = theta - 2 * M_PI;
// //     thetamax = theta;
// //     
// //     llthred = llprev + log(randu());
// //     accept = false;
// //     count = 0;
// //     
// //     while(accept == false){
// //       count = count + 1;
// //       
// //       neweta = eta * cos(theta) + prioreta * sin(theta);
// //       neweta = neweta - mean(neweta);
// //       
// //       newlogZ = cppf_logZ_eta(Xzoop, alpha, neweta);
// //       newlogZobl = cppf_logZobl(newlogZ, indobl);
// //       newlogZsur = cppf_logZsur(newlogZ, indsur);
// //       newloglam = cppf_loglam_cZ_b0t_psit(Xwhale, newlogZ, beta, psi, indtwhale);
// //       newlam = exp(newloglam);
// //       
// //       llnew = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //         cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) +
// //         cppf_loglik_obl(logYobl, newlogZobl, sig2_obl) +
// //         cppf_loglik_sur(logYsur, newlogZsur, lam0, lam1, sig2_sur);
// //       
// //       if(llnew > llthred){
// //         llprev = llnew;
// //         accept = true;
// //       } else {
// //         if(theta < 0){ thetamin = theta; } else { thetamax = theta; }
// //         theta = thetamin + randu() * (thetamax - thetamin);
// //         if ( (count) % 1000 == 0 ) {
// //           Rprintf("ESS for eta: %d iterations...\n", count);
// //         }
// //         if (count > 20000) stop("ESS is not converging...\n");
// //       }
// //     }
// //     eta = neweta;
// //     logZ = newlogZ;
// //     logZobl = newlogZobl;
// //     logZsur = newlogZsur;
// //     loglam = newloglam;
// //     lam = newlam;
// //     
// //     postEta.row(s) = trans(eta);
// //     if(s == 0){ Rprintf("sampled eta\n"); }
// //     
// //     
// //     // update kappa_eta
// //     kappa_eta = cppf_gibbs_kappa_eta(eta, invReta, 2, 2);
// //     postKappa(s,0) = kappa_eta;
// //     if(s == 0){ Rprintf("sampled kappa_eta\n"); }
// //     
// //     
// //     // update lam0
// //     lam0 = cppf_gibbs_lam0(logYsur, logZsur, lam1, sig2_sur, 100);
// //     postLam(s, 0) = lam0;
// //     if(s == 0){ Rprintf("sampled lam0\n"); }
// //     
// //     
// //     // update lam1
// //     lam1 = cppf_gibbs_lam1(logYsur, logZsur, lam0, sig2_sur, 100);
// //     postLam(s, 1) = lam1;
// //     if(s == 0){ Rprintf("sampled lam1\n"); }
// //     
// //     
// //     // update sig2_sur
// //     sig2_sur = cppf_gibbs_sig2_sur(logYsur, logZsur, lam0, lam1, 2, 2);
// //     postSig2(s, 0) = sig2_sur;
// //     if(s == 0){ Rprintf("sampled sig2_sur\n"); }
// //     
// //     
// //     // update beta
// //     if( updateCOV ){
// //       if( (s+1 >= adaptInterval) && (s+1 - (adaptInterval * trunc((s+1) / adaptInterval)) == 0) ){
// //         dummyaccprob = accprob.col(1);
// //         rhat[1] = sum( dummyaccprob.rows(s+1-adaptInterval, s-1) ) / (adaptInterval-1);
// //         gamma1[1] = 1 / pow(adapIter[1], c1);
// //         gamma2[1] = c0 * gamma1[1];
// //         sigma2[1] = exp( log(sigma2[1]) + gamma2[1] * (rhat[1] - ropt) );
// //         
// //         COVbeta = COVbeta + gamma1[1] * ( cov( postBeta.rows(s+1-adaptInterval, s-1)) - COVbeta );
// //         cholCOVbeta = trans( chol( sigma2[1] * ( COVbeta + 0.0000000001 * diagmat(ones(n_be)) ) ) );
// //         
// //         adapIter[1] = adapIter[1] + 1;
// //       }
// //     }
// //     
// //     newbeta = beta + cholCOVbeta * randn(n_be);
// //     
// //     newloglam = cppf_loglam_cZ_b0t_psit(Xwhale, logZ, newbeta, psi, indtwhale);
// //     newlam = exp(newloglam);
// //     
// //     logprob = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) -
// //       cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //       cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) -
// //       cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
// //       MVN_logh(newbeta, zeros(n_be), diagmat(ones(n_be))/100) -
// //       MVN_logh(beta, zeros(n_be), diagmat(ones(n_be))/100);
// //     
// //     if (log(randu()) < logprob) {
// //       
// //       beta = newbeta;
// //       loglam = newloglam;
// //       lam = newlam;
// //       
// //       accprob(s, 1) = 1;
// //     }
// //     
// //     postBeta.row(s) = trans(beta);
// //     if(s == 0){ Rprintf("sampled beta\n"); }
// //     
// //     
// //     
// //     // update psi via Elliptical slice sampling
// //     for(int t = 0; t < ntwhale; t ++){
// //       newpsi = psi;
// //       
// //       llprev = cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //         cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea);
// //       priorpsit = sqrt(kappa_psi) * trans(cholRpsi) * randn(nG);
// //       
// //       thetamin = 0;
// //       thetamax = 2 * M_PI;
// //       theta = thetamin + randu() * (thetamax - thetamin);
// //       thetamin = theta - 2 * M_PI;
// //       thetamax = theta;
// //       
// //       llthred = llprev + log(randu());
// //       accept = false;
// //       count = 0;
// //       
// //       while(accept == false){
// //         count = count + 1;
// //         
// //         newpsi.col(t) = psi.col(t) * cos(theta) + priorpsit * sin(theta);
// //         newpsi.col(t) = newpsi.col(t) - mean(newpsi.col(t));
// //         
// //         newloglam = cppf_loglam_cZ_b0t_psit(Xwhale, logZ, beta, newpsi, indtwhale);
// //         newlam = exp(newloglam);
// //         
// //         llnew = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //           cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea);
// //         
// //         if(llnew > llthred){
// //           llprev = llnew;
// //           accept = true;
// //         } else {
// //           if(theta < 0){ thetamin = theta; } else { thetamax = theta; }
// //           theta = thetamin + randu() * (thetamax - thetamin);
// //           if ( (count) % 1000 == 0 ) {
// //             Rprintf("ESS for psi: %d iterations...\n", count);
// //           }
// //           if (count > 20000) stop("ESS is not converging...\n");
// //         }
// //       }
// //       psi = newpsi;
// //       loglam = newloglam;
// //       lam = newlam;
// //       
// //       postPsit = as<mat>(postPsi[t]);
// //       postPsit.row(s) = trans(psi.col(t));
// //       postPsi[t] = postPsit;
// //     }
// //     if(s == 0){ Rprintf("sampled psi\n"); }
// //     
// //     
// //     // update kappa_psi
// //     kappa_psi = cppf_gibbs_kappa_psit(psi, invRpsi, 2, 2);
// //     postKappa(s,1) = kappa_psi;
// //     if(s == 0){ Rprintf("sampled kappa_psi\n"); }
// //     
// //     
// //     // update c
// //     c = cppf_gibbs_ct(Wpam, lam, dur_dist, p_pam, cellarea, shape_c, scale_c);
// //     postC.row(s) = trans(c);
// //     if(s == 0){ Rprintf("sampled c\n"); }
// //     
// //     
// //     if ( (s+1) % 10 == 0 ) {
// //       Rprintf("Generated %d samples...\n", s+1);
// //     }
// //   }
// //   
// //   
// //   return Rcpp::List::create(Rcpp::Named("postAlpha") = postAlpha,
// //                             Rcpp::Named("postEta") = postEta,
// //                             Rcpp::Named("postLam") = postLam,
// //                             Rcpp::Named("postKappa") = postKappa,
// //                             Rcpp::Named("postSig2") = postSig2,
// //                             Rcpp::Named("postBeta") = postBeta,
// //                             Rcpp::Named("postC") = postC,
// //                             Rcpp::Named("postPsi") = postPsi,
// //                             Rcpp::Named("Accprob") = accprob,
// //                             Rcpp::Named("sigma2") = sigma2,
// //                             Rcpp::Named("adapIter") = adapIter,
// //                             Rcpp::Named("COValpha") = COValpha,
// //                             Rcpp::Named("COVbeta") = COVbeta);
// //   
// // }
// // 
// // 
// // 
// // 
// // 
// // // ============================================================================-
// // // Fix kappa_eta ----
// // // ============================================================================-
// // 
// // 
// // // [[Rcpp::export]]
// // List cppf_fitJoint_cZ_eta_b0t_psi_ct_fixKeta(
// //     int niter, mat distmat,
// //     List Xzoop, vec logYobl, vec logYsur,
// //     List Xwhale, mat Sdist, List Wpam, vec dur_dist,
// //     vec f_dist_obs, List f_dist, List p_pam,
// //     mat indobl, mat indsur, vec indtwhale, mat inddist, vec indmdist, vec indmwhale,
// //     vec alpha, vec eta, double lam0, double lam1,
// //     double kappa_eta, double phi_eta, double sig2_obl, double sig2_sur,
// //     vec beta, vec psi, vec pii, vec c, double shape_c, double scale_c,
// //     double kappa_psi, double phi_psi,
// //     double cellarea,
// //     vec sigma2, mat COValpha, mat COVbeta,
// //     bool updateCOV, int adaptInterval, double adaptFactorExponent, vec adapIter){
// //   
// //   double negativeInf = -std::numeric_limits<float>::infinity();;
// //   
// //   indobl = indobl-1, indsur = indsur-1, indtwhale = indtwhale-1;
// //   inddist = inddist-1, indmdist = indmdist-1, indmwhale = indmwhale-1;
// //   
// //   int nG = distmat.n_rows, n_al = alpha.size(), n_be = beta.size(), ntzoop = eta.n_cols, n_sigma2 = sigma2.size();
// //   
// //   mat postAlpha = zeros(niter, n_al), postLam = zeros(niter, 2);
// //   mat postSig2 = zeros(niter, 1), postKappa = zeros(niter, 1);
// //   mat postBeta = zeros(niter, n_be), postEta = zeros(niter, nG);
// //   mat postPsi = zeros(niter, nG), postC = zeros(niter, c.size());
// //   Rprintf("set up matrices for posterior samples\n");
// //   
// //   mat accprob = zeros(niter, n_sigma2);
// //   vec rhat = zeros(n_sigma2), gamma1 = zeros(n_sigma2), gamma2 = zeros(n_sigma2), dummyaccprob;
// //   double c0 = 1, c1 = adaptFactorExponent, ropt = 0.234, logprob;
// //   mat cholCOValpha = trans( chol( sigma2[0] * ( COValpha + 0.0000000001 * diagmat(ones(n_al)) ) ) );
// //   mat cholCOVbeta = trans( chol( sigma2[1] * ( COVbeta + 0.0000000001 * diagmat(ones(n_be)) ) ) );
// //   Rprintf("set up adaptive proposals\n");
// //   
// //   vec newalpha, newbeta, newlogZobl, newlogZsur;
// //   mat newlogZ, newloglam, newlam;
// //   vec dummyvec;
// //   
// //   int count;
// //   vec xvals;
// //   double llprev, llthred, llnew, thetamin, thetamax, theta;
// //   vec prioreta, neweta, priorpsi, newpsi;
// //   bool accept;
// //   
// //   mat Reta = exp(-distmat / phi_eta), Rpsi = exp(-distmat / phi_psi);
// //   mat invReta = inv(Reta), invRpsi= inv(Rpsi), cholReta = chol(Reta), cholRpsi= chol(Rpsi);
// //   rowvec oneinvReta = trans(ones(nG)) * invReta, oneinvRpsi = trans(ones(nG)) * invRpsi;
// //   double oneinvRetaone = sum( oneinvReta * ones(nG) ), oneinvRpsione = sum( oneinvRpsi * ones(nG) );
// //   Rprintf("set up GP covariance matrices\n");
// //   
// //   mat logZ = cppf_logZ_eta(Xzoop, alpha, eta);
// //   vec logZobl = cppf_logZobl(logZ, indobl);
// //   vec logZsur = cppf_logZsur(logZ, indsur);
// //   mat loglam = cppf_loglam_cZ_b0t_psi(Xwhale, logZ, beta, psi, indtwhale);
// //   mat lam = exp(loglam);
// //   vec m_al = zeros(n_al), diag_al = ones(n_al)/100;
// //   Rprintf("set up initial values\n");
// //   
// //   
// //   // Start MCMC
// //   for(int s = 0; s < niter; s++) {
// //     
// //     // update alpha
// //     if( updateCOV ){
// //       if( (s+1 >= adaptInterval) && (s+1 - (adaptInterval * trunc((s+1) / adaptInterval)) == 0) ){
// //         dummyaccprob = accprob.col(0);
// //         rhat[0] = sum( dummyaccprob.rows(s+1-adaptInterval, s-1) ) / (adaptInterval-1);
// //         gamma1[0] = 1 / pow(adapIter[0], c1);
// //         gamma2[0] = c0 * gamma1[0];
// //         sigma2[0] = exp( log(sigma2[0]) + gamma2[0] * (rhat[0] - ropt) );
// //         
// //         COValpha = COValpha + gamma1[0] * ( cov( postAlpha.rows(s+1-adaptInterval, s-1) ) - COValpha );
// //         cholCOValpha = trans( chol( sigma2[0] * ( COValpha + 0.0000000001 * diagmat(ones(n_al)) ) ) );
// //         
// //         adapIter[0] = adapIter[0] + 1;
// //       }
// //     }
// //     
// //     newalpha = alpha + cholCOValpha * randn(n_al);
// //     
// //     newlogZ = cppf_logZ_eta(Xzoop, newalpha, eta);
// //     newlogZobl = cppf_logZobl(newlogZ, indobl);
// //     newlogZsur = cppf_logZsur(newlogZ, indsur);
// //     newloglam = cppf_loglam_cZ_b0t_psi(Xwhale, newlogZ, beta, psi, indtwhale);
// //     newlam = exp(newloglam);
// //     
// //     logprob = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) -
// //       cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //       cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) -
// //       cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
// //       cppf_loglik_obl(logYobl, newlogZobl, sig2_obl) -
// //       cppf_loglik_obl(logYobl, logZobl, sig2_obl) +
// //       cppf_loglik_sur(logYsur, newlogZsur, lam0, lam1, sig2_sur) -
// //       cppf_loglik_sur(logYsur, logZsur, lam0, lam1, sig2_sur) +
// //       MVN_logh(newalpha, zeros(n_al), diagmat(ones(n_al))/100) -
// //       MVN_logh(alpha, zeros(n_al), diagmat(ones(n_al))/100);
// //     
// //     if (log(randu()) < logprob) {
// //       
// //       alpha = newalpha;
// //       logZ = newlogZ;
// //       logZobl = newlogZobl;
// //       logZsur = newlogZsur;
// //       loglam = newloglam;
// //       lam = newlam;
// //       
// //       accprob(s, 0) = 1;
// //     }
// //     postAlpha.row(s) = trans(alpha);
// //     if(s == 0){ Rprintf("sampled alpha\n"); }
// //     
// //     
// //     
// //     // update eta via Elliptical slice sampling
// //     llprev = cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //       cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
// //       cppf_loglik_obl(logYobl, logZobl, sig2_obl) +
// //       cppf_loglik_sur(logYsur, logZsur, lam0, lam1, sig2_sur);
// //     prioreta = sqrt(kappa_eta) * trans(cholReta) * randn(nG);
// //     
// //     thetamin = 0;
// //     thetamax = 2 * M_PI;
// //     theta = thetamin + randu() * (thetamax - thetamin);
// //     thetamin = theta - 2 * M_PI;
// //     thetamax = theta;
// //     
// //     llthred = llprev + log(randu());
// //     accept = false;
// //     count = 0;
// //     
// //     while(accept == false){
// //       count = count + 1;
// //       
// //       neweta = eta * cos(theta) + prioreta * sin(theta);
// //       neweta = neweta - mean(neweta);
// //       
// //       newlogZ = cppf_logZ_eta(Xzoop, alpha, neweta);
// //       newlogZobl = cppf_logZobl(newlogZ, indobl);
// //       newlogZsur = cppf_logZsur(newlogZ, indsur);
// //       newloglam = cppf_loglam_cZ_b0t_psi(Xwhale, newlogZ, beta, psi, indtwhale);
// //       newlam = exp(newloglam);
// //       
// //       llnew = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //         cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) +
// //         cppf_loglik_obl(logYobl, newlogZobl, sig2_obl) +
// //         cppf_loglik_sur(logYsur, newlogZsur, lam0, lam1, sig2_sur);
// //       
// //       if(llnew > llthred){
// //         llprev = llnew;
// //         accept = true;
// //       } else {
// //         if(theta < 0){ thetamin = theta; } else { thetamax = theta; }
// //         theta = thetamin + randu() * (thetamax - thetamin);
// //         if ( (count) % 1000 == 0 ) {
// //           Rprintf("ESS for eta: %d iterations...\n", count);
// //         }
// //         if (count > 20000) stop("ESS is not converging...\n");
// //       }
// //     }
// //     eta = neweta;
// //     logZ = newlogZ;
// //     logZobl = newlogZobl;
// //     logZsur = newlogZsur;
// //     loglam = newloglam;
// //     lam = newlam;
// //     
// //     postEta.row(s) = trans(eta);
// //     if(s == 0){ Rprintf("sampled eta\n"); }
// //     
// //     
// //     // update lam0
// //     lam0 = cppf_gibbs_lam0(logYsur, logZsur, lam1, sig2_sur, 100);
// //     postLam(s, 0) = lam0;
// //     if(s == 0){ Rprintf("sampled lam0\n"); }
// //     
// //     
// //     // update lam1
// //     lam1 = cppf_gibbs_lam1(logYsur, logZsur, lam0, sig2_sur, 100);
// //     postLam(s, 1) = lam1;
// //     if(s == 0){ Rprintf("sampled lam1\n"); }
// //     
// //     
// //     // update sig2_sur
// //     sig2_sur = cppf_gibbs_sig2_sur(logYsur, logZsur, lam0, lam1, 2, 2);
// //     postSig2(s, 0) = sig2_sur;
// //     if(s == 0){ Rprintf("sampled sig2_sur\n"); }
// //     
// //     
// //     // update beta
// //     if( updateCOV ){
// //       if( (s+1 >= adaptInterval) && (s+1 - (adaptInterval * trunc((s+1) / adaptInterval)) == 0) ){
// //         dummyaccprob = accprob.col(1);
// //         rhat[1] = sum( dummyaccprob.rows(s+1-adaptInterval, s-1) ) / (adaptInterval-1);
// //         gamma1[1] = 1 / pow(adapIter[1], c1);
// //         gamma2[1] = c0 * gamma1[1];
// //         sigma2[1] = exp( log(sigma2[1]) + gamma2[1] * (rhat[1] - ropt) );
// //         
// //         COVbeta = COVbeta + gamma1[1] * ( cov( postBeta.rows(s+1-adaptInterval, s-1)) - COVbeta );
// //         cholCOVbeta = trans( chol( sigma2[1] * ( COVbeta + 0.0000000001 * diagmat(ones(n_be)) ) ) );
// //         
// //         adapIter[1] = adapIter[1] + 1;
// //       }
// //     }
// //     
// //     newbeta = beta + cholCOVbeta * randn(n_be);
// //     
// //     newloglam = cppf_loglam_cZ_b0t_psi(Xwhale, logZ, newbeta, psi, indtwhale);
// //     newlam = exp(newloglam);
// //     
// //     logprob = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) -
// //       cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //       cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) -
// //       cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
// //       MVN_logh(newbeta, zeros(n_be), diagmat(ones(n_be))/100) -
// //       MVN_logh(beta, zeros(n_be), diagmat(ones(n_be))/100);
// //     
// //     if (log(randu()) < logprob) {
// //       
// //       beta = newbeta;
// //       loglam = newloglam;
// //       lam = newlam;
// //       
// //       accprob(s, 1) = 1;
// //     }
// //     
// //     postBeta.row(s) = trans(beta);
// //     if(s == 0){ Rprintf("sampled beta\n"); }
// //     
// //     
// //     
// //     // update psi via elliptical slice sampling
// //     llprev = cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //       cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea);
// //     priorpsi = sqrt(kappa_psi) * trans(cholRpsi) * randn(nG);
// //     
// //     thetamin = 0;
// //     thetamax = 2 * M_PI;
// //     theta = thetamin + randu() * (thetamax - thetamin);
// //     thetamin = theta - 2 * M_PI;
// //     thetamax = theta;
// //     
// //     llthred = llprev + log(randu());
// //     accept = false;
// //     count = 0;
// //     
// //     while(accept == false){
// //       count = count + 1;
// //       
// //       newpsi = psi * cos(theta) + priorpsi * sin(theta);
// //       newpsi = newpsi - mean(newpsi);
// //       
// //       newloglam = cppf_loglam_cZ_b0t_psi(Xwhale, logZ, beta, newpsi, indtwhale);
// //       newlam = exp(newloglam);
// //       
// //       llnew = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //         cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea);
// //       
// //       if(llnew > llthred){
// //         llprev = llnew;
// //         accept = true;
// //       } else {
// //         if(theta < 0){ thetamin = theta; } else { thetamax = theta; }
// //         theta = thetamin + randu() * (thetamax - thetamin);
// //         if ( (count) % 1000 == 0 ) {
// //           Rprintf("ESS for psi: %d iterations...\n", count);
// //         }
// //         if (count > 20000) stop("ESS is not converging...\n");
// //       }
// //     }
// //     psi = newpsi;
// //     loglam = newloglam;
// //     lam = newlam;
// //     
// //     postPsi.row(s) = trans(psi);
// //     if(s == 0){ Rprintf("sampled psi\n"); }
// //     
// //     
// //     // update kappa_psi
// //     kappa_psi = cppf_gibbs_kappa_psi(psi, invRpsi, 2, 2);
// //     postKappa(s,0) = kappa_psi;
// //     if(s == 0){ Rprintf("sampled kappa_psi\n"); }
// //     
// //     
// //     // update c
// //     c = cppf_gibbs_ct(Wpam, lam, dur_dist, p_pam, cellarea, shape_c, scale_c);
// //     postC.row(s) = trans(c);
// //     if(s == 0){ Rprintf("sampled c\n"); }
// //     
// //     
// //     if ( (s+1) % 10 == 0 ) {
// //       Rprintf("Generated %d samples...\n", s+1);
// //     }
// //   }
// //   
// //   
// //   return Rcpp::List::create(Rcpp::Named("postAlpha") = postAlpha,
// //                             Rcpp::Named("postEta") = postEta,
// //                             Rcpp::Named("postLam") = postLam,
// //                             Rcpp::Named("postKappa") = postKappa,
// //                             Rcpp::Named("postSig2") = postSig2,
// //                             Rcpp::Named("postBeta") = postBeta,
// //                             Rcpp::Named("postC") = postC,
// //                             Rcpp::Named("postPsi") = postPsi,
// //                             Rcpp::Named("Accprob") = accprob,
// //                             Rcpp::Named("sigma2") = sigma2,
// //                             Rcpp::Named("adapIter") = adapIter,
// //                             Rcpp::Named("COValpha") = COValpha,
// //                             Rcpp::Named("COVbeta") = COVbeta);
// //   
// // }
// // 
// // 
// // 
// // // [[Rcpp::export]]
// // List cppf_fitJoint_cZ_eta_b0t_psit_ct_fixKeta(
// //     int niter, mat distmat,
// //     List Xzoop, vec logYobl, vec logYsur,
// //     List Xwhale, mat Sdist, List Wpam, vec dur_dist,
// //     vec f_dist_obs, List f_dist, List p_pam,
// //     mat indobl, mat indsur, vec indtwhale, mat inddist, vec indmdist, vec indmwhale,
// //     vec alpha, vec eta, double lam0, double lam1,
// //     double kappa_eta, double phi_eta, double sig2_obl, double sig2_sur,
// //     vec beta, mat psi, vec pii, vec c, double shape_c, double scale_c,
// //     double kappa_psi, double phi_psi,
// //     double cellarea,
// //     vec sigma2, mat COValpha, mat COVbeta,
// //     bool updateCOV, int adaptInterval, double adaptFactorExponent, vec adapIter){
// //   
// //   double negativeInf = -std::numeric_limits<float>::infinity();;
// //   
// //   indobl = indobl-1, indsur = indsur-1, indtwhale = indtwhale-1;
// //   inddist = inddist-1, indmdist = indmdist-1, indmwhale = indmwhale-1;
// //   
// //   int nG = distmat.n_rows, n_al = alpha.size(), n_be = beta.size(), ntzoop = eta.n_cols, ntwhale = psi.n_cols, n_sigma2 = sigma2.size();
// //   
// //   mat postAlpha = zeros(niter, n_al), postLam = zeros(niter, 2);
// //   mat postSig2 = zeros(niter, 1), postKappa = zeros(niter, 1);
// //   mat postBeta = zeros(niter, n_be), postEta = zeros(niter, nG);
// //   mat postPsit = zeros(niter, nG), postC = zeros(niter, c.size());
// //   List postPsi(ntwhale);
// //   for(int i = 0; i < ntwhale; i ++){ postPsi[i] = zeros(niter, nG); }
// //   Rprintf("set up matrices for posterior samples\n");
// //   
// //   mat accprob = zeros(niter, n_sigma2);
// //   vec rhat = zeros(n_sigma2), gamma1 = zeros(n_sigma2), gamma2 = zeros(n_sigma2), dummyaccprob;
// //   double c0 = 1, c1 = adaptFactorExponent, ropt = 0.234, logprob;
// //   mat cholCOValpha = trans( chol( sigma2[0] * ( COValpha + 0.0000000001 * diagmat(ones(n_al)) ) ) );
// //   mat cholCOVbeta = trans( chol( sigma2[1] * ( COVbeta + 0.0000000001 * diagmat(ones(n_be)) ) ) );
// //   Rprintf("set up adaptive proposals\n");
// //   
// //   vec newalpha, newbeta, newlogZobl, newlogZsur;
// //   mat newlogZ, newloglam, newlam;
// //   vec dummyvec;
// //   
// //   int count;
// //   vec xvals;
// //   double llprev, llthred, llnew, thetamin, thetamax, theta;
// //   vec prioreta, neweta, priorpsit;
// //   mat newpsi;
// //   bool accept;
// //   
// //   mat Reta = exp(-distmat / phi_eta), Rpsi = exp(-distmat / phi_psi);
// //   mat invReta = inv(Reta), invRpsi= inv(Rpsi), cholReta = chol(Reta), cholRpsi= chol(Rpsi);
// //   rowvec oneinvReta = trans(ones(nG)) * invReta, oneinvRpsi = trans(ones(nG)) * invRpsi;
// //   double oneinvRetaone = sum( oneinvReta * ones(nG) ), oneinvRpsione = sum( oneinvRpsi * ones(nG) );
// //   Rprintf("set up GP covariance matrices\n");
// //   
// //   mat logZ = cppf_logZ_eta(Xzoop, alpha, eta);
// //   vec logZobl = cppf_logZobl(logZ, indobl);
// //   vec logZsur = cppf_logZsur(logZ, indsur);
// //   mat loglam = cppf_loglam_cZ_b0t_psit(Xwhale, logZ, beta, psi, indtwhale);
// //   mat lam = exp(loglam);
// //   vec m_al = zeros(n_al), diag_al = ones(n_al)/100;
// //   Rprintf("set up initial values\n");
// //   
// //   
// //   // Start MCMC
// //   for(int s = 0; s < niter; s++) {
// //     
// //     // update alpha
// //     if( updateCOV ){
// //       if( (s+1 >= adaptInterval) && (s+1 - (adaptInterval * trunc((s+1) / adaptInterval)) == 0) ){
// //         dummyaccprob = accprob.col(0);
// //         rhat[0] = sum( dummyaccprob.rows(s+1-adaptInterval, s-1) ) / (adaptInterval-1);
// //         gamma1[0] = 1 / pow(adapIter[0], c1);
// //         gamma2[0] = c0 * gamma1[0];
// //         sigma2[0] = exp( log(sigma2[0]) + gamma2[0] * (rhat[0] - ropt) );
// //         
// //         COValpha = COValpha + gamma1[0] * ( cov( postAlpha.rows(s+1-adaptInterval, s-1) ) - COValpha );
// //         cholCOValpha = trans( chol( sigma2[0] * ( COValpha + 0.0000000001 * diagmat(ones(n_al)) ) ) );
// //         
// //         adapIter[0] = adapIter[0] + 1;
// //       }
// //     }
// //     
// //     newalpha = alpha + cholCOValpha * randn(n_al);
// //     
// //     newlogZ = cppf_logZ_eta(Xzoop, newalpha, eta);
// //     newlogZobl = cppf_logZobl(newlogZ, indobl);
// //     newlogZsur = cppf_logZsur(newlogZ, indsur);
// //     newloglam = cppf_loglam_cZ_b0t_psit(Xwhale, newlogZ, beta, psi, indtwhale);
// //     newlam = exp(newloglam);
// //     
// //     logprob = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) -
// //       cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //       cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) -
// //       cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
// //       cppf_loglik_obl(logYobl, newlogZobl, sig2_obl) -
// //       cppf_loglik_obl(logYobl, logZobl, sig2_obl) +
// //       cppf_loglik_sur(logYsur, newlogZsur, lam0, lam1, sig2_sur) -
// //       cppf_loglik_sur(logYsur, logZsur, lam0, lam1, sig2_sur) +
// //       MVN_logh(newalpha, zeros(n_al), diagmat(ones(n_al))/100) -
// //       MVN_logh(alpha, zeros(n_al), diagmat(ones(n_al))/100);
// //     
// //     if (log(randu()) < logprob) {
// //       
// //       alpha = newalpha;
// //       logZ = newlogZ;
// //       logZobl = newlogZobl;
// //       logZsur = newlogZsur;
// //       loglam = newloglam;
// //       lam = newlam;
// //       
// //       accprob(s, 0) = 1;
// //     }
// //     postAlpha.row(s) = trans(alpha);
// //     if(s == 0){ Rprintf("sampled alpha\n"); }
// //     
// //     
// //     
// //     // update eta via Elliptical slice sampling
// //     llprev = cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //       cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
// //       cppf_loglik_obl(logYobl, logZobl, sig2_obl) +
// //       cppf_loglik_sur(logYsur, logZsur, lam0, lam1, sig2_sur);
// //     prioreta = sqrt(kappa_eta) * trans(cholReta) * randn(nG);
// //     
// //     thetamin = 0;
// //     thetamax = 2 * M_PI;
// //     theta = thetamin + randu() * (thetamax - thetamin);
// //     thetamin = theta - 2 * M_PI;
// //     thetamax = theta;
// //     
// //     llthred = llprev + log(randu());
// //     accept = false;
// //     count = 0;
// //     
// //     while(accept == false){
// //       count = count + 1;
// //       
// //       neweta = eta * cos(theta) + prioreta * sin(theta);
// //       neweta = neweta - mean(neweta);
// //       
// //       newlogZ = cppf_logZ_eta(Xzoop, alpha, neweta);
// //       newlogZobl = cppf_logZobl(newlogZ, indobl);
// //       newlogZsur = cppf_logZsur(newlogZ, indsur);
// //       newloglam = cppf_loglam_cZ_b0t_psit(Xwhale, newlogZ, beta, psi, indtwhale);
// //       newlam = exp(newloglam);
// //       
// //       llnew = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //         cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) +
// //         cppf_loglik_obl(logYobl, newlogZobl, sig2_obl) +
// //         cppf_loglik_sur(logYsur, newlogZsur, lam0, lam1, sig2_sur);
// //       
// //       if(llnew > llthred){
// //         llprev = llnew;
// //         accept = true;
// //       } else {
// //         if(theta < 0){ thetamin = theta; } else { thetamax = theta; }
// //         theta = thetamin + randu() * (thetamax - thetamin);
// //         if ( (count) % 1000 == 0 ) {
// //           Rprintf("ESS for eta: %d iterations...\n", count);
// //         }
// //         if (count > 20000) stop("ESS is not converging...\n");
// //       }
// //     }
// //     eta = neweta;
// //     logZ = newlogZ;
// //     logZobl = newlogZobl;
// //     logZsur = newlogZsur;
// //     loglam = newloglam;
// //     lam = newlam;
// //     
// //     postEta.row(s) = trans(eta);
// //     if(s == 0){ Rprintf("sampled eta\n"); }
// //     
// //     
// //     // update lam0
// //     lam0 = cppf_gibbs_lam0(logYsur, logZsur, lam1, sig2_sur, 100);
// //     postLam(s, 0) = lam0;
// //     if(s == 0){ Rprintf("sampled lam0\n"); }
// //     
// //     
// //     // update lam1
// //     lam1 = cppf_gibbs_lam1(logYsur, logZsur, lam0, sig2_sur, 100);
// //     postLam(s, 1) = lam1;
// //     if(s == 0){ Rprintf("sampled lam1\n"); }
// //     
// //     
// //     // update sig2_sur
// //     sig2_sur = cppf_gibbs_sig2_sur(logYsur, logZsur, lam0, lam1, 2, 2);
// //     postSig2(s, 0) = sig2_sur;
// //     if(s == 0){ Rprintf("sampled sig2_sur\n"); }
// //     
// //     
// //     // update beta
// //     if( updateCOV ){
// //       if( (s+1 >= adaptInterval) && (s+1 - (adaptInterval * trunc((s+1) / adaptInterval)) == 0) ){
// //         dummyaccprob = accprob.col(1);
// //         rhat[1] = sum( dummyaccprob.rows(s+1-adaptInterval, s-1) ) / (adaptInterval-1);
// //         gamma1[1] = 1 / pow(adapIter[1], c1);
// //         gamma2[1] = c0 * gamma1[1];
// //         sigma2[1] = exp( log(sigma2[1]) + gamma2[1] * (rhat[1] - ropt) );
// //         
// //         COVbeta = COVbeta + gamma1[1] * ( cov( postBeta.rows(s+1-adaptInterval, s-1)) - COVbeta );
// //         cholCOVbeta = trans( chol( sigma2[1] * ( COVbeta + 0.0000000001 * diagmat(ones(n_be)) ) ) );
// //         
// //         adapIter[1] = adapIter[1] + 1;
// //       }
// //     }
// //     
// //     newbeta = beta + cholCOVbeta * randn(n_be);
// //     
// //     newloglam = cppf_loglam_cZ_b0t_psit(Xwhale, logZ, newbeta, psi, indtwhale);
// //     newlam = exp(newloglam);
// //     
// //     logprob = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) -
// //       cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //       cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) -
// //       cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
// //       MVN_logh(newbeta, zeros(n_be), diagmat(ones(n_be))/100) -
// //       MVN_logh(beta, zeros(n_be), diagmat(ones(n_be))/100);
// //     
// //     if (log(randu()) < logprob) {
// //       
// //       beta = newbeta;
// //       loglam = newloglam;
// //       lam = newlam;
// //       
// //       accprob(s, 1) = 1;
// //     }
// //     
// //     postBeta.row(s) = trans(beta);
// //     if(s == 0){ Rprintf("sampled beta\n"); }
// //     
// //     
// //     
// //     // update psi via Elliptical slice sampling
// //     for(int t = 0; t < ntwhale; t ++){
// //       newpsi = psi;
// //       
// //       llprev = cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //         cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea);
// //       priorpsit = sqrt(kappa_psi) * trans(cholRpsi) * randn(nG);
// //       
// //       thetamin = 0;
// //       thetamax = 2 * M_PI;
// //       theta = thetamin + randu() * (thetamax - thetamin);
// //       thetamin = theta - 2 * M_PI;
// //       thetamax = theta;
// //       
// //       llthred = llprev + log(randu());
// //       accept = false;
// //       count = 0;
// //       
// //       while(accept == false){
// //         count = count + 1;
// //         
// //         newpsi.col(t) = psi.col(t) * cos(theta) + priorpsit * sin(theta);
// //         newpsi.col(t) = newpsi.col(t) - mean(newpsi.col(t));
// //         
// //         newloglam = cppf_loglam_cZ_b0t_psit(Xwhale, logZ, beta, newpsi, indtwhale);
// //         newlam = exp(newloglam);
// //         
// //         llnew = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //           cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea);
// //         
// //         if(llnew > llthred){
// //           llprev = llnew;
// //           accept = true;
// //         } else {
// //           if(theta < 0){ thetamin = theta; } else { thetamax = theta; }
// //           theta = thetamin + randu() * (thetamax - thetamin);
// //           if ( (count) % 1000 == 0 ) {
// //             Rprintf("ESS for psi: %d iterations...\n", count);
// //           }
// //           if (count > 20000) stop("ESS is not converging...\n");
// //         }
// //       }
// //       psi = newpsi;
// //       loglam = newloglam;
// //       lam = newlam;
// //       
// //       postPsit = as<mat>(postPsi[t]);
// //       postPsit.row(s) = trans(psi.col(t));
// //       postPsi[t] = postPsit;
// //     }
// //     if(s == 0){ Rprintf("sampled psi\n"); }
// //     
// //     
// //     // update kappa_psi
// //     kappa_psi = cppf_gibbs_kappa_psit(psi, invRpsi, 2, 2);
// //     postKappa(s,0) = kappa_psi;
// //     if(s == 0){ Rprintf("sampled kappa_psi\n"); }
// //     
// //     
// //     // update c
// //     c = cppf_gibbs_ct(Wpam, lam, dur_dist, p_pam, cellarea, shape_c, scale_c);
// //     postC.row(s) = trans(c);
// //     if(s == 0){ Rprintf("sampled c\n"); }
// //     
// //     
// //     if ( (s+1) % 10 == 0 ) {
// //       Rprintf("Generated %d samples...\n", s+1);
// //     }
// //   }
// //   
// //   
// //   return Rcpp::List::create(Rcpp::Named("postAlpha") = postAlpha,
// //                             Rcpp::Named("postEta") = postEta,
// //                             Rcpp::Named("postLam") = postLam,
// //                             Rcpp::Named("postKappa") = postKappa,
// //                             Rcpp::Named("postSig2") = postSig2,
// //                             Rcpp::Named("postBeta") = postBeta,
// //                             Rcpp::Named("postC") = postC,
// //                             Rcpp::Named("postPsi") = postPsi,
// //                             Rcpp::Named("Accprob") = accprob,
// //                             Rcpp::Named("sigma2") = sigma2,
// //                             Rcpp::Named("adapIter") = adapIter,
// //                             Rcpp::Named("COValpha") = COValpha,
// //                             Rcpp::Named("COVbeta") = COVbeta);
// //   
// // }
// // 
// // 
// // 
// // // [[Rcpp::export]]
// // List cppf_fitJoint_cZ_a0t_etat_b0t_psi_ct_ja_fixKeta(
// //     int niter, mat distmat,
// //     List Xzoop, vec logYobl, vec logYsur,
// //     List Xwhale, mat Sdist, List Wpam, vec dur_dist,
// //     vec f_dist_obs, List f_dist, List p_pam,
// //     mat indobl, mat indsur, vec indtwhale, mat inddist, vec indmdist, vec indmwhale,
// //     double alpha0tilde, vec alpha, mat eta, double lam0, double lam1,
// //     double tau2, double kappa_eta, double phi_eta, double sig2_obl, double sig2_sur,
// //     vec beta, vec psi, vec pii, vec c, double shape_c, double scale_c,
// //     double kappa_psi, double phi_psi,
// //     double cellarea,
// //     vec sigma2, mat COValpha, mat COVbeta,
// //     bool updateCOV, int adaptInterval, double adaptFactorExponent, vec adapIter){
// //   
// //   double negativeInf = -std::numeric_limits<float>::infinity();;
// //   
// //   indobl = indobl-1, indsur = indsur-1, indtwhale = indtwhale-1;
// //   inddist = inddist-1, indmdist = indmdist-1, indmwhale = indmwhale-1;
// //   
// //   int nG = distmat.n_rows, n_al = alpha.size(), n_be = beta.size(), ntzoop = eta.n_cols, n_sigma2 = sigma2.size();
// //   
// //   mat postAlpha0tilde = zeros(niter, 1), postAlpha = zeros(niter, n_al), postLam = zeros(niter, 2);
// //   mat postSig2 = zeros(niter, 1), postKappa = zeros(niter, 1);
// //   mat postBeta = zeros(niter, n_be), postC = zeros(niter, c.size());
// //   mat postPsi(niter, nG), postEtat = zeros(niter, nG);
// //   List postEta(ntzoop);
// //   for(int i = 0; i < ntzoop; i ++){
// //     postEta[i] = zeros(niter, nG);
// //   }
// //   Rprintf("set up matrices for posterior samples\n");
// //   
// //   mat accprob = zeros(niter, n_sigma2);
// //   vec rhat = zeros(n_sigma2), gamma1 = zeros(n_sigma2), gamma2 = zeros(n_sigma2), dummyaccprob;
// //   double c0 = 1, c1 = adaptFactorExponent, ropt = 0.234, logprob;
// //   mat cholCOValpha = trans( chol( sigma2[0] * ( COValpha + 0.0000000001 * diagmat(ones(n_al)) ) ) );
// //   mat cholCOVbeta = trans( chol( sigma2[1] * ( COVbeta + 0.0000000001 * diagmat(ones(n_be)) ) ) );
// //   Rprintf("set up adaptive proposals\n");
// //   
// //   vec newalpha, newbeta, newlogZobl, newlogZsur;
// //   mat newlogZ, newloglam, newlam;
// //   vec dummyvec;
// //   
// //   int count;
// //   vec xvals;
// //   double llprev, llthred, llnew, thetamin, thetamax, theta;
// //   vec prioretat, priorpsi, newpsi;
// //   mat neweta;
// //   bool accept;
// //   
// //   mat Reta = exp(-distmat / phi_eta), Rpsi = exp(-distmat / phi_psi);
// //   mat invReta = inv(Reta), invRpsi= inv(Rpsi), cholReta = chol(Reta), cholRpsi= chol(Rpsi);
// //   rowvec oneinvReta = trans(ones(nG)) * invReta, oneinvRpsi = trans(ones(nG)) * invRpsi;
// //   double oneinvRetaone = sum( oneinvReta * ones(nG) ), oneinvRpsione = sum( oneinvRpsi * ones(nG) );
// //   Rprintf("set up GP covariance matrices\n");
// //   
// //   mat logZ = cppf_logZ_a0t_etat_ja(Xzoop, alpha, eta);
// //   vec logZobl = cppf_logZobl(logZ, indobl);
// //   vec logZsur = cppf_logZsur(logZ, indsur);
// //   mat loglam = cppf_loglam_cZ_b0t_psi(Xwhale, logZ, beta, psi, indtwhale);
// //   mat lam = exp(loglam);
// //   Rprintf("set up initial values\n");
// //   
// //   
// //   // Start MCMC
// //   for(int s = 0; s < niter; s++) {
// //     
// //     // update alpha
// //     if( updateCOV ){
// //       if( (s+1 >= adaptInterval) && (s+1 - (adaptInterval * trunc((s+1) / adaptInterval)) == 0) ){
// //         dummyaccprob = accprob.col(0);
// //         rhat[0] = sum( dummyaccprob.rows(s+1-adaptInterval, s-1) ) / (adaptInterval-1);
// //         gamma1[0] = 1 / pow(adapIter[0], c1);
// //         gamma2[0] = c0 * gamma1[0];
// //         sigma2[0] = exp( log(sigma2[0]) + gamma2[0] * (rhat[0] - ropt) );
// //         
// //         COValpha = COValpha + gamma1[0] * ( cov( postAlpha.rows(s+1-adaptInterval, s-1) ) - COValpha );
// //         cholCOValpha = trans( chol( sigma2[0] * ( COValpha + 0.0000000001 * diagmat(ones(n_al)) ) ) );
// //         
// //         adapIter[0] = adapIter[0] + 1;
// //       }
// //     }
// //     
// //     newalpha = alpha + cholCOValpha * randn(n_al);
// //     
// //     newlogZ = cppf_logZ_a0t_etat_ja(Xzoop, newalpha, eta);
// //     newlogZobl = cppf_logZobl(newlogZ, indobl);
// //     newlogZsur = cppf_logZsur(newlogZ, indsur);
// //     newloglam = cppf_loglam_cZ_b0t_psi(Xwhale, newlogZ, beta, psi, indtwhale);
// //     newlam = exp(newloglam);
// //     
// //     logprob = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) -
// //       cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //       cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) -
// //       cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
// //       cppf_loglik_obl(logYobl, newlogZobl, sig2_obl) -
// //       cppf_loglik_obl(logYobl, logZobl, sig2_obl) +
// //       cppf_loglik_sur(logYsur, newlogZsur, lam0, lam1, sig2_sur) -
// //       cppf_loglik_sur(logYsur, logZsur, lam0, lam1, sig2_sur) +
// //       MVN_logh(newalpha, zeros(n_al), diagmat(ones(n_al))/100) -
// //       MVN_logh(alpha, zeros(n_al), diagmat(ones(n_al))/100);
// //     
// //     if (log(randu()) < logprob) {
// //       
// //       alpha = newalpha;
// //       logZ = newlogZ;
// //       logZobl = newlogZobl;
// //       logZsur = newlogZsur;
// //       loglam = newloglam;
// //       lam = newlam;
// //       
// //       accprob(s, 0) = 1;
// //     }
// //     postAlpha.row(s) = trans(alpha);
// //     if(s == 0){ Rprintf("sampled alpha\n"); }
// //     
// //     
// //     // Gibbs update for alpha0tilde
// //     alpha0tilde = cppf_gibbs_alpha0tilde(alpha.rows(0, ntzoop-1), tau2);
// //     postAlpha0tilde[s] = alpha0tilde;
// //     if(s == 0){ Rprintf("sampled alpha0tilde\n"); }
// //     
// //     
// //     // update eta via Elliptical slice sampling
// //     for(int t = 0; t < ntzoop; t ++){
// //       neweta = eta;
// //       
// //       llprev = cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //         cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
// //         cppf_loglik_obl(logYobl, logZobl, sig2_obl) +
// //         cppf_loglik_sur(logYsur, logZsur, lam0, lam1, sig2_sur);
// //       prioretat = sqrt(kappa_eta) * trans(cholReta) * randn(nG);
// //       
// //       thetamin = 0;
// //       thetamax = 2 * M_PI;
// //       theta = thetamin + randu() * (thetamax - thetamin);
// //       thetamin = theta - 2 * M_PI;
// //       thetamax = theta;
// //       
// //       llthred = llprev + log(randu());
// //       accept = false;
// //       count = 0;
// //       
// //       while(accept == false){
// //         count = count + 1;
// //         
// //         neweta.col(t) = eta.col(t) * cos(theta) + prioretat * sin(theta);
// //         neweta.col(t) = neweta.col(t) - mean(neweta.col(t));
// //         
// //         newlogZ = cppf_logZ_a0t_etat_ja(Xzoop, alpha, neweta);
// //         newlogZobl = cppf_logZobl(newlogZ, indobl);
// //         newlogZsur = cppf_logZsur(newlogZ, indsur);
// //         newloglam = cppf_loglam_cZ_b0t_psi(Xwhale, newlogZ, beta, psi, indtwhale);
// //         newlam = exp(newloglam);
// //         
// //         llnew = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //           cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) +
// //           cppf_loglik_obl(logYobl, newlogZobl, sig2_obl) +
// //           cppf_loglik_sur(logYsur, newlogZsur, lam0, lam1, sig2_sur);
// //         
// //         if(llnew > llthred){
// //           llprev = llnew;
// //           accept = true;
// //         } else {
// //           if(theta < 0){ thetamin = theta; } else { thetamax = theta; }
// //           theta = thetamin + randu() * (thetamax - thetamin);
// //           if ( (count) % 1000 == 0 ) {
// //             Rprintf("ESS for eta: %d iterations...\n", count);
// //           }
// //           if (count > 20000) stop("ESS is not converging...\n");
// //         }
// //       }
// //       eta = neweta;
// //       logZ = newlogZ;
// //       logZobl = newlogZobl;
// //       logZsur = newlogZsur;
// //       loglam = newloglam;
// //       lam = newlam;
// //       
// //       postEtat = as<mat>(postEta[t]);
// //       postEtat.row(s) = trans(eta.col(t));
// //       postEta[t] = postEtat;
// //     }
// //     if(s == 0){ Rprintf("sampled eta\n"); }
// //     
// //     
// //     
// //     // update lam0
// //     lam0 = cppf_gibbs_lam0(logYsur, logZsur, lam1, sig2_sur, 100);
// //     postLam(s, 0) = lam0;
// //     if(s == 0){ Rprintf("sampled lam0\n"); }
// //     
// //     
// //     // update lam1
// //     lam1 = cppf_gibbs_lam1(logYsur, logZsur, lam0, sig2_sur, 100);
// //     postLam(s, 1) = lam1;
// //     if(s == 0){ Rprintf("sampled lam1\n"); }
// //     
// //     
// //     // update sig2_sur
// //     sig2_sur = cppf_gibbs_sig2_sur(logYsur, logZsur, lam0, lam1, 2, 2);
// //     postSig2(s, 0) = sig2_sur;
// //     if(s == 0){ Rprintf("sampled sig2_sur\n"); }
// //     
// //     
// //     // update beta
// //     if( updateCOV ){
// //       if( (s+1 >= adaptInterval) && (s+1 - (adaptInterval * trunc((s+1) / adaptInterval)) == 0) ){
// //         dummyaccprob = accprob.col(1);
// //         rhat[1] = sum( dummyaccprob.rows(s+1-adaptInterval, s-1) ) / (adaptInterval-1);
// //         gamma1[1] = 1 / pow(adapIter[1], c1);
// //         gamma2[1] = c0 * gamma1[1];
// //         sigma2[1] = exp( log(sigma2[1]) + gamma2[1] * (rhat[1] - ropt) );
// //         
// //         COVbeta = COVbeta + gamma1[1] * ( cov( postBeta.rows(s+1-adaptInterval, s-1)) - COVbeta );
// //         cholCOVbeta = trans( chol( sigma2[1] * ( COVbeta + 0.0000000001 * diagmat(ones(n_be)) ) ) );
// //         
// //         adapIter[1] = adapIter[1] + 1;
// //       }
// //     }
// //     
// //     newbeta = beta + cholCOVbeta * randn(n_be);
// //     
// //     newloglam = cppf_loglam_cZ_b0t_psi(Xwhale, logZ, newbeta, psi, indtwhale);
// //     newlam = exp(newloglam);
// //     
// //     logprob = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) -
// //       cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //       cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) -
// //       cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
// //       MVN_logh(newbeta, zeros(n_be), diagmat(ones(n_be))/100) -
// //       MVN_logh(beta, zeros(n_be), diagmat(ones(n_be))/100);
// //     
// //     if (log(randu()) < logprob) {
// //       
// //       beta = newbeta;
// //       loglam = newloglam;
// //       lam = newlam;
// //       
// //       accprob(s, 1) = 1;
// //     }
// //     
// //     postBeta.row(s) = trans(beta);
// //     if(s == 0){ Rprintf("sampled beta\n"); }
// //     
// //     
// //     // update psi via elliptical slice sampling
// //     llprev = cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //       cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea);
// //     priorpsi = sqrt(kappa_psi) * trans(cholRpsi) * randn(nG);
// //     
// //     thetamin = 0;
// //     thetamax = 2 * M_PI;
// //     theta = thetamin + randu() * (thetamax - thetamin);
// //     thetamin = theta - 2 * M_PI;
// //     thetamax = theta;
// //     
// //     llthred = llprev + log(randu());
// //     accept = false;
// //     count = 0;
// //     
// //     while(accept == false){
// //       count = count + 1;
// //       
// //       newpsi = psi * cos(theta) + priorpsi * sin(theta);
// //       newpsi = newpsi - mean(newpsi);
// //       
// //       newloglam = cppf_loglam_cZ_b0t_psi(Xwhale, logZ, beta, newpsi, indtwhale);
// //       newlam = exp(newloglam);
// //       
// //       llnew = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //         cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea);
// //       
// //       if(llnew > llthred){
// //         llprev = llnew;
// //         accept = true;
// //       } else {
// //         if(theta < 0){ thetamin = theta; } else { thetamax = theta; }
// //         theta = thetamin + randu() * (thetamax - thetamin);
// //         if ( (count) % 1000 == 0 ) {
// //           Rprintf("ESS for psi: %d iterations...\n", count);
// //         }
// //         if (count > 20000) stop("ESS is not converging...\n");
// //       }
// //     }
// //     psi = newpsi;
// //     loglam = newloglam;
// //     lam = newlam;
// //     
// //     postPsi.row(s) = trans(psi);
// //     if(s == 0){ Rprintf("sampled psi\n"); }
// //     
// //     
// //     // update kappa_psi
// //     kappa_psi = cppf_gibbs_kappa_psi(psi, invRpsi, 2, 2);
// //     postKappa(s,0) = kappa_psi;
// //     if(s == 0){ Rprintf("sampled kappa_psi\n"); }
// //     
// //     
// //     // update c
// //     c = cppf_gibbs_ct(Wpam, lam, dur_dist, p_pam, cellarea, shape_c, scale_c);
// //     postC.row(s) = trans(c);
// //     if(s == 0){ Rprintf("sampled c\n"); }
// //     
// //     
// //     if ( (s+1) % 10 == 0 ) {
// //       Rprintf("Generated %d samples...\n", s+1);
// //     }
// //   }
// //   
// //   
// //   return Rcpp::List::create(Rcpp::Named("postAlpha0tilde") = postAlpha0tilde,
// //                             Rcpp::Named("postAlpha") = postAlpha,
// //                             Rcpp::Named("postEta") = postEta,
// //                             Rcpp::Named("postLam") = postLam,
// //                             Rcpp::Named("postKappa") = postKappa,
// //                             Rcpp::Named("postSig2") = postSig2,
// //                             Rcpp::Named("postBeta") = postBeta,
// //                             Rcpp::Named("postC") = postC,
// //                             Rcpp::Named("postPsi") = postPsi,
// //                             Rcpp::Named("Accprob") = accprob,
// //                             Rcpp::Named("sigma2") = sigma2,
// //                             Rcpp::Named("adapIter") = adapIter,
// //                             Rcpp::Named("COValpha") = COValpha,
// //                             Rcpp::Named("COVbeta") = COVbeta);
// //   
// // }
// // 
// // 
// 
// // 
// // // [[Rcpp::export]]
// // List cppf_fitJoint_a0t_eta_b0t_psit_ja_fixC(
// //     int niter, mat distmat,
// //     List Xzoop, vec logYobl, vec logYsur,
// //     List Xwhale, mat Sdist, List Wpam, vec dur_dist,
// //     vec f_dist_obs, List f_dist, List p_pam,
// //     mat indobl, mat indsur, vec indtwhale, mat inddist, vec indmdist, vec indmwhale,
// //     double alpha0tilde, vec alpha, vec eta, double lam0, double lam1,
// //     double tau2, double kappa_eta, double phi_eta, double sig2_obl, double sig2_sur,
// //     vec beta, mat psi, vec pii, vec c,
// //     double kappa_psi, double phi_psi,
// //     double cellarea,
// //     vec sigma2, mat COValpha, mat COVbeta,
// //     bool updateCOV, int adaptInterval, double adaptFactorExponent, vec adapIter){
// //   
// //   double negativeInf = -std::numeric_limits<float>::infinity();;
// //   
// //   indobl = indobl-1, indsur = indsur-1, indtwhale = indtwhale-1;
// //   inddist = inddist-1, indmdist = indmdist-1, indmwhale = indmwhale-1;
// //   
// //   int nG = distmat.n_rows, n_al = alpha.size(), n_be = beta.size(), ntzoop = eta.n_cols, ntwhale = psi.n_cols, n_sigma2 = sigma2.size();
// //   
// //   mat postAlpha0tilde = zeros(niter, 1), postAlpha = zeros(niter, n_al), postLam = zeros(niter, 2);
// //   mat postSig2 = zeros(niter, 1), postKappa = zeros(niter, 2);
// //   mat postBeta = zeros(niter, n_be), postEta = zeros(niter, nG);
// //   mat postPsit = zeros(niter, nG);
// //   List postPsi(ntwhale);
// //   for(int i = 0; i < ntwhale; i ++){ postPsi[i] = zeros(niter, nG); }
// //   Rprintf("set up matrices for posterior samples\n");
// //   
// //   mat accprob = zeros(niter, n_sigma2);
// //   vec rhat = zeros(n_sigma2), gamma1 = zeros(n_sigma2), gamma2 = zeros(n_sigma2), dummyaccprob;
// //   double c0 = 1, c1 = adaptFactorExponent, ropt = 0.234, logprob;
// //   mat cholCOValpha = trans( chol( sigma2[0] * ( COValpha + 0.0000000001 * diagmat(ones(n_al)) ) ) );
// //   mat cholCOVbeta = trans( chol( sigma2[1] * ( COVbeta + 0.0000000001 * diagmat(ones(n_be)) ) ) );
// //   Rprintf("set up adaptive proposals\n");
// //   
// //   vec newalpha, newbeta, newlogZobl, newlogZsur;
// //   mat newlogZ, newloglam, newlam;
// //   vec dummyvec;
// //   
// //   int count;
// //   vec xvals;
// //   double llprev, llthred, llnew, thetamin, thetamax, theta;
// //   vec prioreta, neweta, priorpsit;
// //   mat newpsi;
// //   bool accept;
// //   
// //   mat Reta = exp(-distmat / phi_eta), Rpsi = exp(-distmat / phi_psi);
// //   mat invReta = inv(Reta), invRpsi= inv(Rpsi), cholReta = chol(Reta), cholRpsi= chol(Rpsi);
// //   rowvec oneinvReta = trans(ones(nG)) * invReta, oneinvRpsi = trans(ones(nG)) * invRpsi;
// //   double oneinvRetaone = sum( oneinvReta * ones(nG) ), oneinvRpsione = sum( oneinvRpsi * ones(nG) );
// //   Rprintf("set up GP covariance matrices\n");
// //   
// //   mat logZ = cppf_logZ_a0t_eta_ja(Xzoop, alpha, eta);
// //   vec logZobl = cppf_logZobl(logZ, indobl);
// //   vec logZsur = cppf_logZsur(logZ, indsur);
// //   mat loglam = cppf_loglam_b0t_psit(Xwhale, logZ, beta, psi, indtwhale);
// //   mat lam = exp(loglam);
// //   vec m_al = zeros(n_al), diag_al = ones(n_al)/100;
// //   Rprintf("set up initial values\n");
// //   
// //   
// //   // Start MCMC
// //   for(int s = 0; s < niter; s++) {
// //     
// //     // update alpha
// //     if( updateCOV ){
// //       if( (s+1 >= adaptInterval) && (s+1 - (adaptInterval * trunc((s+1) / adaptInterval)) == 0) ){
// //         dummyaccprob = accprob.col(0);
// //         rhat[0] = sum( dummyaccprob.rows(s+1-adaptInterval, s-1) ) / (adaptInterval-1);
// //         gamma1[0] = 1 / pow(adapIter[0], c1);
// //         gamma2[0] = c0 * gamma1[0];
// //         sigma2[0] = exp( log(sigma2[0]) + gamma2[0] * (rhat[0] - ropt) );
// //         
// //         COValpha = COValpha + gamma1[0] * ( cov( postAlpha.rows(s+1-adaptInterval, s-1) ) - COValpha );
// //         cholCOValpha = trans( chol( sigma2[0] * ( COValpha + 0.0000000001 * diagmat(ones(n_al)) ) ) );
// //         
// //         adapIter[0] = adapIter[0] + 1;
// //       }
// //     }
// //     
// //     newalpha = alpha + cholCOValpha * randn(n_al);
// //     
// //     newlogZ = cppf_logZ_a0t_eta_ja(Xzoop, newalpha, eta);
// //     newlogZobl = cppf_logZobl(newlogZ, indobl);
// //     newlogZsur = cppf_logZsur(newlogZ, indsur);
// //     newloglam = cppf_loglam_b0t_psit(Xwhale, newlogZ, beta, psi, indtwhale);
// //     newlam = exp(newloglam);
// //     
// //     for(int t = 0; t < ntzoop; t ++){
// //       m_al[t] = alpha0tilde;
// //       diag_al[t] = 1/tau2;
// //     }
// //     
// //     logprob = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) -
// //       cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //       cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) -
// //       cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
// //       cppf_loglik_obl(logYobl, newlogZobl, sig2_obl) -
// //       cppf_loglik_obl(logYobl, logZobl, sig2_obl) +
// //       cppf_loglik_sur(logYsur, newlogZsur, lam0, lam1, sig2_sur) -
// //       cppf_loglik_sur(logYsur, logZsur, lam0, lam1, sig2_sur) +
// //       MVN_logh(newalpha, m_al, diagmat(diag_al)) -
// //       MVN_logh(alpha, m_al, diagmat(diag_al));
// //     
// //     if (log(randu()) < logprob) {
// //       
// //       alpha = newalpha;
// //       logZ = newlogZ;
// //       logZobl = newlogZobl;
// //       logZsur = newlogZsur;
// //       loglam = newloglam;
// //       lam = newlam;
// //       
// //       accprob(s, 0) = 1;
// //     }
// //     postAlpha.row(s) = trans(alpha);
// //     if(s == 0){ Rprintf("sampled alpha\n"); }
// //     
// //     
// //     // Gibbs update for alpha0tilde
// //     alpha0tilde = cppf_gibbs_alpha0tilde(alpha.rows(0, ntzoop-1), tau2);
// //     postAlpha0tilde[s] = alpha0tilde;
// //     if(s == 0){ Rprintf("sampled alpha0tilde\n"); }
// //     
// //     
// //     
// //     // update eta via Elliptical slice sampling
// //     llprev = cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //       cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
// //       cppf_loglik_obl(logYobl, logZobl, sig2_obl) +
// //       cppf_loglik_sur(logYsur, logZsur, lam0, lam1, sig2_sur);
// //     prioreta = sqrt(kappa_eta) * trans(cholReta) * randn(nG);
// //     
// //     thetamin = 0;
// //     thetamax = 2 * M_PI;
// //     theta = thetamin + randu() * (thetamax - thetamin);
// //     thetamin = theta - 2 * M_PI;
// //     thetamax = theta;
// //     
// //     llthred = llprev + log(randu());
// //     accept = false;
// //     count = 0;
// //     
// //     while(accept == false){
// //       count = count + 1;
// //       
// //       neweta = eta * cos(theta) + prioreta * sin(theta);
// //       neweta = neweta - mean(neweta);
// //       
// //       newlogZ = cppf_logZ_a0t_eta_ja(Xzoop, alpha, neweta);
// //       newlogZobl = cppf_logZobl(newlogZ, indobl);
// //       newlogZsur = cppf_logZsur(newlogZ, indsur);
// //       newloglam = cppf_loglam_b0t_psit(Xwhale, newlogZ, beta, psi, indtwhale);
// //       newlam = exp(newloglam);
// //       
// //       llnew = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //         cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) +
// //         cppf_loglik_obl(logYobl, newlogZobl, sig2_obl) +
// //         cppf_loglik_sur(logYsur, newlogZsur, lam0, lam1, sig2_sur);
// //       
// //       if(llnew > llthred){
// //         llprev = llnew;
// //         accept = true;
// //       } else {
// //         if(theta < 0){ thetamin = theta; } else { thetamax = theta; }
// //         theta = thetamin + randu() * (thetamax - thetamin);
// //         if ( (count) % 1000 == 0 ) {
// //           Rprintf("ESS for eta: %d iterations...\n", count);
// //         }
// //         if (count > 20000) stop("ESS is not converging...\n");
// //       }
// //     }
// //     eta = neweta;
// //     logZ = newlogZ;
// //     logZobl = newlogZobl;
// //     logZsur = newlogZsur;
// //     loglam = newloglam;
// //     lam = newlam;
// //     
// //     postEta.row(s) = trans(eta);
// //     if(s == 0){ Rprintf("sampled eta\n"); }
// //     
// //     
// //     // update kappa_eta
// //     kappa_eta = cppf_gibbs_kappa_eta(eta, invReta, 2, 2);
// //     postKappa(s,0) = kappa_eta;
// //     if(s == 0){ Rprintf("sampled kappa_eta\n"); }
// //     
// //     
// //     // update lam0
// //     lam0 = cppf_gibbs_lam0(logYsur, logZsur, lam1, sig2_sur, 100);
// //     postLam(s, 0) = lam0;
// //     if(s == 0){ Rprintf("sampled lam0\n"); }
// //     
// //     
// //     // update lam1
// //     lam1 = cppf_gibbs_lam1(logYsur, logZsur, lam0, sig2_sur, 100);
// //     postLam(s, 1) = lam1;
// //     if(s == 0){ Rprintf("sampled lam1\n"); }
// //     
// //     
// //     // update sig2_sur
// //     sig2_sur = cppf_gibbs_sig2_sur(logYsur, logZsur, lam0, lam1, 2, 2);
// //     postSig2(s, 0) = sig2_sur;
// //     if(s == 0){ Rprintf("sampled sig2_sur\n"); }
// //     
// //     
// //     // update beta
// //     if( updateCOV ){
// //       if( (s+1 >= adaptInterval) && (s+1 - (adaptInterval * trunc((s+1) / adaptInterval)) == 0) ){
// //         dummyaccprob = accprob.col(1);
// //         rhat[1] = sum( dummyaccprob.rows(s+1-adaptInterval, s-1) ) / (adaptInterval-1);
// //         gamma1[1] = 1 / pow(adapIter[1], c1);
// //         gamma2[1] = c0 * gamma1[1];
// //         sigma2[1] = exp( log(sigma2[1]) + gamma2[1] * (rhat[1] - ropt) );
// //         
// //         COVbeta = COVbeta + gamma1[1] * ( cov( postBeta.rows(s+1-adaptInterval, s-1)) - COVbeta );
// //         cholCOVbeta = trans( chol( sigma2[1] * ( COVbeta + 0.0000000001 * diagmat(ones(n_be)) ) ) );
// //         
// //         adapIter[1] = adapIter[1] + 1;
// //       }
// //     }
// //     
// //     newbeta = beta + cholCOVbeta * randn(n_be);
// //     
// //     newloglam = cppf_loglam_b0t_psit(Xwhale, logZ, newbeta, psi, indtwhale);
// //     newlam = exp(newloglam);
// //     
// //     logprob = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) -
// //       cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //       cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) -
// //       cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
// //       MVN_logh(newbeta, zeros(n_be), diagmat(ones(n_be))/100) -
// //       MVN_logh(beta, zeros(n_be), diagmat(ones(n_be))/100);
// //     
// //     if (log(randu()) < logprob) {
// //       
// //       beta = newbeta;
// //       loglam = newloglam;
// //       lam = newlam;
// //       
// //       accprob(s, 1) = 1;
// //     }
// //     
// //     postBeta.row(s) = trans(beta);
// //     if(s == 0){ Rprintf("sampled beta\n"); }
// //     
// //     
// //     
// //     // update psi via Elliptical slice sampling
// //     for(int t = 0; t < ntwhale; t ++){
// //       newpsi = psi;
// //       
// //       llprev = cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //         cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea);
// //       priorpsit = sqrt(kappa_psi) * trans(cholRpsi) * randn(nG);
// //       
// //       thetamin = 0;
// //       thetamax = 2 * M_PI;
// //       theta = thetamin + randu() * (thetamax - thetamin);
// //       thetamin = theta - 2 * M_PI;
// //       thetamax = theta;
// //       
// //       llthred = llprev + log(randu());
// //       accept = false;
// //       count = 0;
// //       
// //       while(accept == false){
// //         count = count + 1;
// //         
// //         newpsi.col(t) = psi.col(t) * cos(theta) + priorpsit * sin(theta);
// //         newpsi.col(t) = newpsi.col(t) - mean(newpsi.col(t));
// //         
// //         newloglam = cppf_loglam_b0t_psit(Xwhale, logZ, beta, newpsi, indtwhale);
// //         newlam = exp(newloglam);
// //         
// //         llnew = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //           cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea);
// //         
// //         if(llnew > llthred){
// //           llprev = llnew;
// //           accept = true;
// //         } else {
// //           if(theta < 0){ thetamin = theta; } else { thetamax = theta; }
// //           theta = thetamin + randu() * (thetamax - thetamin);
// //           if ( (count) % 1000 == 0 ) {
// //             Rprintf("ESS for psi: %d iterations...\n", count);
// //           }
// //           if (count > 20000) stop("ESS is not converging...\n");
// //         }
// //       }
// //       psi = newpsi;
// //       loglam = newloglam;
// //       lam = newlam;
// //       
// //       postPsit = as<mat>(postPsi[t]);
// //       postPsit.row(s) = trans(psi.col(t));
// //       postPsi[t] = postPsit;
// //     }
// //     if(s == 0){ Rprintf("sampled psi\n"); }
// //     
// //     
// //     // update kappa_psi
// //     kappa_psi = cppf_gibbs_kappa_psit(psi, invRpsi, 2, 2);
// //     postKappa(s,1) = kappa_psi;
// //     if(s == 0){ Rprintf("sampled kappa_psi\n"); }
// //     
// //     
// //     if ( (s+1) % 10 == 0 ) {
// //       Rprintf("Generated %d samples...\n", s+1);
// //     }
// //   }
// //   
// //   
// //   return Rcpp::List::create(Rcpp::Named("postAlpha0tilde") = postAlpha0tilde,
// //                             Rcpp::Named("postAlpha") = postAlpha,
// //                             Rcpp::Named("postEta") = postEta,
// //                             Rcpp::Named("postLam") = postLam,
// //                             Rcpp::Named("postKappa") = postKappa,
// //                             Rcpp::Named("postSig2") = postSig2,
// //                             Rcpp::Named("postBeta") = postBeta,
// //                             Rcpp::Named("postPsi") = postPsi,
// //                             Rcpp::Named("Accprob") = accprob,
// //                             Rcpp::Named("sigma2") = sigma2,
// //                             Rcpp::Named("adapIter") = adapIter,
// //                             Rcpp::Named("COValpha") = COValpha,
// //                             Rcpp::Named("COVbeta") = COVbeta);
// //   
// // }
// // 
// // 
// // 
// // // sig2_obl is known
// // // [[Rcpp::export]]
// // List cppf_fitJoint_b0(
// //     int niter, mat distmat,
// //     List Xzoop, vec logYobl, vec logYsur,
// //     List Xwhale, mat Sdist, mat Wpam, vec dur_dist,
// //     vec f_dist_obs, List f_dist, List p_pam,
// //     mat indobl, mat indsur, vec indtwhale, mat inddist, vec indmdist, vec indmwhale,
// //     double alpha0tilde, vec alpha0, vec alpha, mat eta, double lam0, double lam1,
// //     double tau2, double kappa_eta, double phi_eta, double sig2_obl, double sig2_sur,
// //     vec beta, vec pii, vec c, double shape_c, double scale_c,
// //     double cellarea,
// //     vec sigma2, vec COValpha0, mat COValpha, mat COVbeta,
// //     bool updateCOV, int adaptInterval, double adaptFactorExponent, vec adapIter){
// //   
// //   double negativeInf = -std::numeric_limits<float>::infinity();;
// //   
// //   indobl = indobl-1, indsur = indsur-1, indtwhale = indtwhale-1;
// //   inddist = inddist-1, indmdist = indmdist-1, indmwhale = indmwhale-1;
// //   
// //   int nG = distmat.n_rows, n_al0 = alpha0.size(), n_al = alpha.size(), n_be = beta.size(), ntzoop = eta.n_cols, n_sigma2 = sigma2.size();
// //   
// //   mat postAlpha0tilde = zeros(niter, 1), postAlpha0 = zeros(niter, n_al0), postAlpha = zeros(niter, n_al), postLam = zeros(niter, 2);
// //   mat postSig2 = zeros(niter, 1), postKappa = zeros(niter, 1);
// //   mat postBeta = zeros(niter, n_be), postC = zeros(niter, c.size()), postEtat = zeros(niter, nG);
// //   List postEta(ntzoop);
// //   for(int i = 0; i < ntzoop; i ++){
// //     postEta[i] = zeros(niter, nG);
// //   }
// //   Rprintf("set up matrices for posterior samples\n");
// //   
// //   mat accprob = zeros(niter, n_sigma2);
// //   vec rhat = zeros(n_sigma2), gamma1 = zeros(n_sigma2), gamma2 = zeros(n_sigma2), dummyaccprob;
// //   double c0 = 1, c1 = adaptFactorExponent, ropt = 0.234, logprob;
// //   vec cholCOValpha0(n_al0);
// //   for(int i = 0; i < n_al0; i ++){ cholCOValpha0[i] = sqrt(sigma2[i] * COValpha0[i]); }
// //   mat cholCOValpha = trans( chol( sigma2[n_al0] * ( COValpha + 0.0000000001 * diagmat(ones(n_al)) ) ) );
// //   mat cholCOVbeta = trans( chol( sigma2[n_al0+1] * ( COVbeta + 0.0000000001 * diagmat(ones(n_be)) ) ) );
// //   Rprintf("set up adaptive proposals\n");
// //   
// //   double alpha0i, newalpha0i;
// //   vec newalpha0, newalpha, newbeta, newlogZobl, newlogZsur;
// //   mat newlogZ, newloglam, newlam;
// //   vec dummyvec;
// //   
// //   int count;
// //   vec xvals;
// //   double llprev, llthred, llnew, thetamin, thetamax, theta;
// //   vec prioretat;
// //   mat neweta;
// //   bool accept;
// //   
// //   mat Reta = exp(-distmat / phi_eta);
// //   mat invReta = inv(Reta), cholReta = chol(Reta);
// //   rowvec oneinvReta = trans(ones(nG)) * invReta;
// //   double oneinvRetaone = sum( oneinvReta * ones(nG) );
// //   Rprintf("set up GP covariance matrices\n");
// //   
// //   mat logZ = cppf_logZ_a0t_etat(Xzoop, alpha0, alpha, eta);
// //   vec logZobl = cppf_logZobl(logZ, indobl);
// //   vec logZsur = cppf_logZsur(logZ, indsur);
// //   mat loglam = cppf_loglam_b0(Xwhale, logZ, beta, indtwhale);
// //   mat lam = exp(loglam);
// //   Rprintf("set up initial values\n");
// //   
// //   
// //   // Start MCMC
// //   for(int s = 0; s < niter; s++) {
// //     
// //     // update alpha0
// //     for(int i = 0; i < n_al0; i ++){
// //       
// //       if( updateCOV ){
// //         if( (s+1 >= adaptInterval) && (s+1 - (adaptInterval * trunc((s+1) / adaptInterval)) == 0) ){
// //           dummyaccprob = accprob.col(i);
// //           rhat[i] = sum( dummyaccprob.rows(s+1-adaptInterval, s-1) ) / (adaptInterval-1);
// //           gamma1[i] = 1 / pow(adapIter[i], c1);
// //           gamma2[i] = c0 * gamma1[i];
// //           sigma2[i] = exp( log(sigma2[i]) + gamma2[i] * (rhat[i] - ropt) );
// //           
// //           
// //           dummyvec = postAlpha0.col(i);
// //           COValpha0[i] = COValpha0[i] + gamma1[i] * ( var( dummyvec.rows(s+1-adaptInterval, s-1) ) - COValpha0[i] );
// //           cholCOValpha0[i] = sqrt(  sigma2[i] * COValpha0[i] );
// //           
// //           adapIter[i] = adapIter[i] + 1;
// //         }
// //       }
// //       
// //       alpha0i = alpha0[i];
// //       newalpha0i = alpha0i + cholCOValpha0[i] * randn();
// //       newalpha0 = alpha0;
// //       newalpha0[i] = newalpha0i;
// //       
// //       newlogZ = cppf_logZ_a0t_etat(Xzoop, newalpha0, alpha, eta);
// //       newlogZobl = cppf_logZobl(newlogZ, indobl);
// //       newlogZsur = cppf_logZsur(newlogZ, indsur);
// //       newloglam = cppf_loglam_b0(Xwhale, newlogZ, beta, indtwhale);
// //       newlam = exp(newloglam);
// //       
// //       logprob = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) -
// //         cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //         cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) -
// //         cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
// //         cppf_loglik_obl(logYobl, newlogZobl, sig2_obl) -
// //         cppf_loglik_obl(logYobl, logZobl, sig2_obl) +
// //         cppf_loglik_sur(logYsur, newlogZsur, lam0, lam1, sig2_sur) -
// //         cppf_loglik_sur(logYsur, logZsur, lam0, lam1, sig2_sur) +
// //         normal_logh(newalpha0i, alpha0tilde, tau2) -
// //         normal_logh(alpha0i, alpha0tilde, tau2);
// //       
// //       if (log(randu()) < logprob) {
// //         
// //         alpha0 = newalpha0;
// //         logZ = newlogZ;
// //         logZobl = newlogZobl;
// //         logZsur = newlogZsur;
// //         loglam = newloglam;
// //         lam = newlam;
// //         
// //         accprob(s, i) = 1;
// //       }
// //     }
// //     postAlpha0.row(s) = trans(alpha0);
// //     if(s == 0){ Rprintf("sampled alpha0\n"); }
// //     
// //     
// //     // Gibbs update for alpha0tilde
// //     alpha0tilde = cppf_gibbs_alpha0tilde(alpha0, tau2);
// //     postAlpha0tilde[s] = alpha0tilde;
// //     if(s == 0){ Rprintf("sampled alpha0tilde\n"); }
// //     
// //     
// //     // update alpha
// //     if( updateCOV ){
// //       if( (s+1 >= adaptInterval) && (s+1 - (adaptInterval * trunc((s+1) / adaptInterval)) == 0) ){
// //         dummyaccprob = accprob.col(n_al0);
// //         rhat[n_al0] = sum( dummyaccprob.rows(s+1-adaptInterval, s-1) ) / (adaptInterval-1);
// //         gamma1[n_al0] = 1 / pow(adapIter[n_al0], c1);
// //         gamma2[n_al0] = c0 * gamma1[n_al0];
// //         sigma2[n_al0] = exp( log(sigma2[n_al0]) + gamma2[n_al0] * (rhat[n_al0] - ropt) );
// //         
// //         COValpha = COValpha + gamma1[n_al0] * ( cov( postAlpha.rows(s+1-adaptInterval, s-1) ) - COValpha );
// //         cholCOValpha = trans( chol( sigma2[n_al0] * ( COValpha + 0.0000000001 * diagmat(ones(n_al)) ) ) );
// //         
// //         adapIter[n_al0] = adapIter[n_al0] + 1;
// //       }
// //     }
// //     
// //     newalpha = alpha + cholCOValpha * randn(n_al);
// //     
// //     newlogZ = cppf_logZ_a0t_etat(Xzoop, alpha0, newalpha, eta);
// //     newlogZobl = cppf_logZobl(newlogZ, indobl);
// //     newlogZsur = cppf_logZsur(newlogZ, indsur);
// //     newloglam = cppf_loglam_b0(Xwhale, newlogZ, beta, indtwhale);
// //     newlam = exp(newloglam);
// //     
// //     logprob = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) -
// //       cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //       cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) -
// //       cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
// //       cppf_loglik_obl(logYobl, newlogZobl, sig2_obl) -
// //       cppf_loglik_obl(logYobl, logZobl, sig2_obl) +
// //       cppf_loglik_sur(logYsur, newlogZsur, lam0, lam1, sig2_sur) -
// //       cppf_loglik_sur(logYsur, logZsur, lam0, lam1, sig2_sur) +
// //       MVN_logh(newalpha, zeros(n_al), diagmat(ones(n_al))/100) -
// //       MVN_logh(alpha, zeros(n_al), diagmat(ones(n_al))/100);
// //     
// //     if (log(randu()) < logprob) {
// //       
// //       alpha = newalpha;
// //       logZ = newlogZ;
// //       logZobl = newlogZobl;
// //       logZsur = newlogZsur;
// //       loglam = newloglam;
// //       lam = newlam;
// //       
// //       accprob(s, n_al0) = 1;
// //     }
// //     
// //     postAlpha.row(s) = trans(alpha);
// //     if(s == 0){ Rprintf("sampled alpha\n"); }
// //     
// //     
// //     // update eta via Elliptical slice sampling
// //     for(int t = 0; t < ntzoop; t ++){
// //       neweta = eta;
// //       
// //       llprev = cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //         cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
// //         cppf_loglik_obl(logYobl, logZobl, sig2_obl) +
// //         cppf_loglik_sur(logYsur, logZsur, lam0, lam1, sig2_sur);
// //       prioretat = sqrt(kappa_eta) * trans(cholReta) * randn(nG);
// //       
// //       thetamin = 0;
// //       thetamax = 2 * M_PI;
// //       theta = thetamin + randu() * (thetamax - thetamin);
// //       thetamin = theta - 2 * M_PI;
// //       thetamax = theta;
// //       
// //       llthred = llprev + log(randu());
// //       accept = false;
// //       count = 0;
// //       
// //       while(accept == false){
// //         count = count + 1;
// //         
// //         neweta.col(t) = eta.col(t) * cos(theta) + prioretat * sin(theta);
// //         neweta.col(t) = neweta.col(t) - mean(neweta.col(t));
// //         
// //         newlogZ = cppf_logZ_a0t_etat(Xzoop, alpha0, alpha, neweta);
// //         newlogZobl = cppf_logZobl(newlogZ, indobl);
// //         newlogZsur = cppf_logZsur(newlogZ, indsur);
// //         newloglam = cppf_loglam_b0(Xwhale, newlogZ, beta, indtwhale);
// //         newlam = exp(newloglam);
// //         
// //         llnew = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //           cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) +
// //           cppf_loglik_obl(logYobl, newlogZobl, sig2_obl) +
// //           cppf_loglik_sur(logYsur, newlogZsur, lam0, lam1, sig2_sur);
// //         
// //         if(llnew > llthred){
// //           llprev = llnew;
// //           accept = true;
// //         } else {
// //           if(theta < 0){ thetamin = theta; } else { thetamax = theta; }
// //           theta = thetamin + randu() * (thetamax - thetamin);
// //           if ( (count) % 1000 == 0 ) {
// //             Rprintf("ESS for eta: %d iterations...\n", count);
// //           }
// //           if (count > 20000) stop("ESS is not converging...\n");
// //         }
// //       }
// //       eta = neweta;
// //       logZ = newlogZ;
// //       logZobl = newlogZobl;
// //       logZsur = newlogZsur;
// //       loglam = newloglam;
// //       lam = newlam;
// //       
// //       postEtat = as<mat>(postEta[t]);
// //       postEtat.row(s) = trans(eta.col(t));
// //       postEta[t] = postEtat;
// //     }
// //     if(s == 0){ Rprintf("sampled eta\n"); }
// //     
// //     
// //     // update kappa_eta
// //     kappa_eta = cppf_gibbs_kappa_eta(eta, invReta, 2, 2);
// //     postKappa(s,0) = kappa_eta;
// //     if(s == 0){ Rprintf("sampled kappa_eta\n"); }
// //     
// //     
// //     // update lam0
// //     lam0 = cppf_gibbs_lam0(logYsur, logZsur, lam1, sig2_sur, 100);
// //     postLam(s, 0) = lam0;
// //     if(s == 0){ Rprintf("sampled lam0\n"); }
// //     
// //     
// //     // update lam1
// //     lam1 = cppf_gibbs_lam1(logYsur, logZsur, lam0, sig2_sur, 100);
// //     postLam(s, 1) = lam1;
// //     if(s == 0){ Rprintf("sampled lam1\n"); }
// //     
// //     
// //     // update sig2_sur
// //     sig2_sur = cppf_gibbs_sig2_sur(logYsur, logZsur, lam0, lam1, 2, 2);
// //     postSig2(s, 0) = sig2_sur;
// //     if(s == 0){ Rprintf("sampled sig2_sur\n"); }
// //     
// //     
// //     
// //     // update beta
// //     if( updateCOV ){
// //       if( (s+1 >= adaptInterval) && (s+1 - (adaptInterval * trunc((s+1) / adaptInterval)) == 0) ){
// //         dummyaccprob = accprob.col(n_al0+1);
// //         rhat[n_al0+1] = sum( dummyaccprob.rows(s+1-adaptInterval, s-1) ) / (adaptInterval-1);
// //         gamma1[n_al0+1] = 1 / pow(adapIter[n_al0+1], c1);
// //         gamma2[n_al0+1] = c0 * gamma1[n_al0+1];
// //         sigma2[n_al0+1] = exp( log(sigma2[n_al0+1]) + gamma2[n_al0+1] * (rhat[n_al0+1] - ropt) );
// //         
// //         COVbeta = COVbeta + gamma1[n_al0+1] * ( cov( postBeta.rows(s+1-adaptInterval, s-1)) - COVbeta );
// //         cholCOVbeta = trans( chol( sigma2[n_al0+1] * ( COVbeta + 0.0000000001 * diagmat(ones(n_be)) ) ) );
// //         
// //         adapIter[n_al0+1] = adapIter[n_al0+1] + 1;
// //       }
// //     }
// //     
// //     newbeta = beta + cholCOVbeta * randn(n_be);
// //     
// //     newloglam = cppf_loglam_b0(Xwhale, logZ, newbeta, indtwhale);
// //     newlam = exp(newloglam);
// //     
// //     logprob = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) -
// //       cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //       cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) -
// //       cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
// //       MVN_logh(newbeta, zeros(n_be), diagmat(ones(n_be))/100) -
// //       MVN_logh(beta, zeros(n_be), diagmat(ones(n_be))/100);
// //     
// //     if (log(randu()) < logprob) {
// //       
// //       beta = newbeta;
// //       loglam = newloglam;
// //       lam = newlam;
// //       
// //       accprob(s, n_al0+1) = 1;
// //     }
// //     
// //     postBeta.row(s) = trans(beta);
// //     if(s == 0){ Rprintf("sampled beta\n"); }
// //     
// //     
// //     // update c
// //     c = cppf_gibbs_ct(Wpam, lam, dur_dist, p_pam, cellarea, shape_c, scale_c);
// //     postC.row(s) = trans(c);
// //     if(s == 0){ Rprintf("sampled c\n"); }
// //     
// //     
// //     if ( (s+1) % 10 == 0 ) {
// //       Rprintf("Generated %d samples...\n", s+1);
// //     }
// //   }
// //   
// //   
// //   return Rcpp::List::create(Rcpp::Named("postAlpha0") = postAlpha0,
// //                             Rcpp::Named("postAlpha0tilde") = postAlpha0tilde,
// //                             Rcpp::Named("postAlpha") = postAlpha,
// //                             Rcpp::Named("postEta") = postEta,
// //                             Rcpp::Named("postLam") = postLam,
// //                             Rcpp::Named("postKappa") = postKappa,
// //                             Rcpp::Named("postSig2") = postSig2,
// //                             Rcpp::Named("postBeta") = postBeta,
// //                             Rcpp::Named("postC") = postC,
// //                             Rcpp::Named("Accprob") = accprob,
// //                             Rcpp::Named("sigma2") = sigma2,
// //                             Rcpp::Named("adapIter") = adapIter,
// //                             Rcpp::Named("COValpha0") = COValpha0,
// //                             Rcpp::Named("COValpha") = COValpha,
// //                             Rcpp::Named("COVbeta") = COVbeta);
// //   
// // }
// // 
// // 
// // 
// // // sig2_obl is known
// // // annual GP for whale component
// // // [[Rcpp::export]]
// // List cppf_fitJoint_b0_gp(
// //     int niter, mat distmat,
// //     List Xzoop, vec logYobl, vec logYsur,
// //     List Xwhale, mat Sdist, mat Wpam, vec dur_dist,
// //     vec f_dist_obs, List f_dist, List p_pam,
// //     mat indobl, mat indsur, vec indtwhale, mat inddist, vec indmdist, vec indmwhale,
// //     double alpha0tilde, vec alpha0, vec alpha, mat eta, double lam0, double lam1,
// //     double tau2, double kappa_eta, double phi_eta, double sig2_obl, double sig2_sur,
// //     vec beta, vec psi, vec pii, vec c, double shape_c, double scale_c,
// //     double kappa_psi, double phi_psi,
// //     double cellarea,
// //     vec sigma2, vec COValpha0, mat COValpha, mat COVbeta,
// //     bool updateCOV, int adaptInterval, double adaptFactorExponent, vec adapIter){
// //   
// //   double negativeInf = -std::numeric_limits<float>::infinity();;
// //   
// //   indobl = indobl-1, indsur = indsur-1, indtwhale = indtwhale-1;
// //   inddist = inddist-1, indmdist = indmdist-1, indmwhale = indmwhale-1;
// //   
// //   int nG = distmat.n_rows, n_al0 = alpha0.size(), n_al = alpha.size(), n_be = beta.size(), ntzoop = eta.n_cols, n_sigma2 = sigma2.size(); // ntwhale = psi.n_cols
// //   
// //   mat postAlpha0tilde = zeros(niter, 1), postAlpha0 = zeros(niter, n_al0), postAlpha = zeros(niter, n_al), postLam = zeros(niter, 2);
// //   mat postSig2 = zeros(niter, 1), postKappa = zeros(niter, 2);
// //   mat postBeta = zeros(niter, n_be), postPsi(niter, nG), postC = zeros(niter, c.size()), postEtat = zeros(niter, nG);
// //   List postEta(ntzoop);
// //   for(int i = 0; i < ntzoop; i ++){
// //     postEta[i] = zeros(niter, nG);
// //   }
// //   Rprintf("set up matrices for posterior samples\n");
// //   
// //   mat accprob = zeros(niter, n_sigma2);
// //   vec rhat = zeros(n_sigma2), gamma1 = zeros(n_sigma2), gamma2 = zeros(n_sigma2), dummyaccprob;
// //   double c0 = 1, c1 = adaptFactorExponent, ropt = 0.234, logprob;
// //   vec cholCOValpha0(n_al0);
// //   for(int i = 0; i < n_al0; i ++){ cholCOValpha0[i] = sqrt(sigma2[i] * COValpha0[i]); }
// //   mat cholCOValpha = trans( chol( sigma2[n_al0] * ( COValpha + 0.0000000001 * diagmat(ones(n_al)) ) ) );
// //   mat cholCOVbeta = trans( chol( sigma2[n_al0+1] * ( COVbeta + 0.0000000001 * diagmat(ones(n_be)) ) ) );
// //   Rprintf("set up adaptive proposals\n");
// //   
// //   double alpha0i, newalpha0i;
// //   vec newalpha0, newalpha, newbeta, newlogZobl, newlogZsur;
// //   mat newlogZ, newloglam, newlam;
// //   vec dummyvec;
// //   
// //   int count;
// //   vec xvals;
// //   double llprev, llthred, llnew, thetamin, thetamax, theta;
// //   vec prioretat, priorpsi, newpsi;
// //   mat neweta;
// //   bool accept;
// //   
// //   mat Reta = exp(-distmat / phi_eta), Rpsi = exp(-distmat / phi_psi);
// //   mat invReta = inv(Reta), invRpsi= inv(Rpsi), cholReta = chol(Reta), cholRpsi= chol(Rpsi);
// //   rowvec oneinvReta = trans(ones(nG)) * invReta, oneinvRpsi = trans(ones(nG)) * invRpsi;
// //   double oneinvRetaone = sum( oneinvReta * ones(nG) ), oneinvRpsione = sum( oneinvRpsi * ones(nG) );
// //   Rprintf("set up GP covariance matrices\n");
// //   
// //   mat logZ = cppf_logZ_a0t_etat(Xzoop, alpha0, alpha, eta);
// //   vec logZobl = cppf_logZobl(logZ, indobl);
// //   vec logZsur = cppf_logZsur(logZ, indsur);
// //   mat loglam = cppf_loglam_b0_gp(Xwhale, logZ, beta, psi, indtwhale);
// //   mat lam = exp(loglam);
// //   Rprintf("set up initial values\n");
// //   
// //   
// //   // Start MCMC
// //   for(int s = 0; s < niter; s++) {
// //     
// //     // update alpha0
// //     for(int i = 0; i < n_al0; i ++){
// //       
// //       if( updateCOV ){
// //         if( (s+1 >= adaptInterval) && (s+1 - (adaptInterval * trunc((s+1) / adaptInterval)) == 0) ){
// //           dummyaccprob = accprob.col(i);
// //           rhat[i] = sum( dummyaccprob.rows(s+1-adaptInterval, s-1) ) / (adaptInterval-1);
// //           gamma1[i] = 1 / pow(adapIter[i], c1);
// //           gamma2[i] = c0 * gamma1[i];
// //           sigma2[i] = exp( log(sigma2[i]) + gamma2[i] * (rhat[i] - ropt) );
// //           
// //           
// //           dummyvec = postAlpha0.col(i);
// //           COValpha0[i] = COValpha0[i] + gamma1[i] * ( var( dummyvec.rows(s+1-adaptInterval, s-1) ) - COValpha0[i] );
// //           cholCOValpha0[i] = sqrt(  sigma2[i] * COValpha0[i] );
// //           
// //           adapIter[i] = adapIter[i] + 1;
// //         }
// //       }
// //       
// //       alpha0i = alpha0[i];
// //       newalpha0i = alpha0i + cholCOValpha0[i] * randn();
// //       newalpha0 = alpha0;
// //       newalpha0[i] = newalpha0i;
// //       
// //       newlogZ = cppf_logZ_a0t_etat(Xzoop, newalpha0, alpha, eta);
// //       newlogZobl = cppf_logZobl(newlogZ, indobl);
// //       newlogZsur = cppf_logZsur(newlogZ, indsur);
// //       newloglam = cppf_loglam_b0(Xwhale, newlogZ, beta, indtwhale);
// //       newlam = exp(newloglam);
// //       
// //       logprob = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) -
// //         cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //         cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) -
// //         cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
// //         cppf_loglik_obl(logYobl, newlogZobl, sig2_obl) -
// //         cppf_loglik_obl(logYobl, logZobl, sig2_obl) +
// //         cppf_loglik_sur(logYsur, newlogZsur, lam0, lam1, sig2_sur) -
// //         cppf_loglik_sur(logYsur, logZsur, lam0, lam1, sig2_sur) +
// //         normal_logh(newalpha0i, alpha0tilde, tau2) -
// //         normal_logh(alpha0i, alpha0tilde, tau2);
// //       
// //       if (log(randu()) < logprob) {
// //         
// //         alpha0 = newalpha0;
// //         logZ = newlogZ;
// //         logZobl = newlogZobl;
// //         logZsur = newlogZsur;
// //         loglam = newloglam;
// //         lam = newlam;
// //         
// //         accprob(s, i) = 1;
// //       }
// //     }
// //     postAlpha0.row(s) = trans(alpha0);
// //     if(s == 0){ Rprintf("sampled alpha0\n"); }
// //     
// //     
// //     // Gibbs update for alpha0tilde
// //     alpha0tilde = cppf_gibbs_alpha0tilde(alpha0, tau2);
// //     postAlpha0tilde[s] = alpha0tilde;
// //     if(s == 0){ Rprintf("sampled alpha0tilde\n"); }
// //     
// //     
// //     // update alpha
// //     if( updateCOV ){
// //       if( (s+1 >= adaptInterval) && (s+1 - (adaptInterval * trunc((s+1) / adaptInterval)) == 0) ){
// //         dummyaccprob = accprob.col(n_al0);
// //         rhat[n_al0] = sum( dummyaccprob.rows(s+1-adaptInterval, s-1) ) / (adaptInterval-1);
// //         gamma1[n_al0] = 1 / pow(adapIter[n_al0], c1);
// //         gamma2[n_al0] = c0 * gamma1[n_al0];
// //         sigma2[n_al0] = exp( log(sigma2[n_al0]) + gamma2[n_al0] * (rhat[n_al0] - ropt) );
// //         
// //         COValpha = COValpha + gamma1[n_al0] * ( cov( postAlpha.rows(s+1-adaptInterval, s-1) ) - COValpha );
// //         cholCOValpha = trans( chol( sigma2[n_al0] * ( COValpha + 0.0000000001 * diagmat(ones(n_al)) ) ) );
// //         
// //         adapIter[n_al0] = adapIter[n_al0] + 1;
// //       }
// //     }
// //     
// //     newalpha = alpha + cholCOValpha * randn(n_al);
// //     
// //     newlogZ = cppf_logZ_a0t_etat(Xzoop, alpha0, newalpha, eta);
// //     newlogZobl = cppf_logZobl(newlogZ, indobl);
// //     newlogZsur = cppf_logZsur(newlogZ, indsur);
// //     newloglam = cppf_loglam_b0_gp(Xwhale, newlogZ, beta, psi, indtwhale);
// //     newlam = exp(newloglam);
// //     
// //     logprob = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) -
// //       cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //       cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) -
// //       cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
// //       cppf_loglik_obl(logYobl, newlogZobl, sig2_obl) -
// //       cppf_loglik_obl(logYobl, logZobl, sig2_obl) +
// //       cppf_loglik_sur(logYsur, newlogZsur, lam0, lam1, sig2_sur) -
// //       cppf_loglik_sur(logYsur, logZsur, lam0, lam1, sig2_sur) +
// //       MVN_logh(newalpha, zeros(n_al), diagmat(ones(n_al))/100) -
// //       MVN_logh(alpha, zeros(n_al), diagmat(ones(n_al))/100);
// //     
// //     if (log(randu()) < logprob) {
// //       
// //       alpha = newalpha;
// //       logZ = newlogZ;
// //       logZobl = newlogZobl;
// //       logZsur = newlogZsur;
// //       loglam = newloglam;
// //       lam = newlam;
// //       
// //       accprob(s, n_al0) = 1;
// //     }
// //     
// //     postAlpha.row(s) = trans(alpha);
// //     if(s == 0){ Rprintf("sampled alpha\n"); }
// //     
// //     
// //     // update eta via Elliptical slice sampling
// //     for(int t = 0; t < ntzoop; t ++){
// //       neweta = eta;
// //       
// //       llprev = cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //         cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
// //         cppf_loglik_obl(logYobl, logZobl, sig2_obl) +
// //         cppf_loglik_sur(logYsur, logZsur, lam0, lam1, sig2_sur);
// //       prioretat = sqrt(kappa_eta) * trans(cholReta) * randn(nG);
// //       
// //       thetamin = 0;
// //       thetamax = 2 * M_PI;
// //       theta = thetamin + randu() * (thetamax - thetamin);
// //       thetamin = theta - 2 * M_PI;
// //       thetamax = theta;
// //       
// //       llthred = llprev + log(randu());
// //       accept = false;
// //       count = 0;
// //       
// //       while(accept == false){
// //         count = count + 1;
// //         
// //         neweta.col(t) = eta.col(t) * cos(theta) + prioretat * sin(theta);
// //         neweta.col(t) = neweta.col(t) - mean(neweta.col(t));
// //         
// //         newlogZ = cppf_logZ_a0t_etat(Xzoop, alpha0, alpha, neweta);
// //         newlogZobl = cppf_logZobl(newlogZ, indobl);
// //         newlogZsur = cppf_logZsur(newlogZ, indsur);
// //         newloglam = cppf_loglam_b0_gp(Xwhale, newlogZ, beta, psi, indtwhale);
// //         newlam = exp(newloglam);
// //         
// //         llnew = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //           cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) +
// //           cppf_loglik_obl(logYobl, newlogZobl, sig2_obl) +
// //           cppf_loglik_sur(logYsur, newlogZsur, lam0, lam1, sig2_sur);
// //         
// //         if(llnew > llthred){
// //           llprev = llnew;
// //           accept = true;
// //         } else {
// //           if(theta < 0){ thetamin = theta; } else { thetamax = theta; }
// //           theta = thetamin + randu() * (thetamax - thetamin);
// //           if ( (count) % 1000 == 0 ) {
// //             Rprintf("ESS for eta: %d iterations...\n", count);
// //           }
// //           if (count > 20000) stop("ESS is not converging...\n");
// //         }
// //       }
// //       eta = neweta;
// //       logZ = newlogZ;
// //       logZobl = newlogZobl;
// //       logZsur = newlogZsur;
// //       loglam = newloglam;
// //       lam = newlam;
// //       
// //       postEtat = as<mat>(postEta[t]);
// //       postEtat.row(s) = trans(eta.col(t));
// //       postEta[t] = postEtat;
// //     }
// //     if(s == 0){ Rprintf("sampled eta\n"); }
// //     
// //     
// //     // update kappa_eta
// //     kappa_eta = cppf_gibbs_kappa_eta(eta, invReta, 2, 2);
// //     postKappa(s,0) = kappa_eta;
// //     if(s == 0){ Rprintf("sampled kappa_eta\n"); }
// //     
// //     
// //     // update lam0
// //     lam0 = cppf_gibbs_lam0(logYsur, logZsur, lam1, sig2_sur, 100);
// //     postLam(s, 0) = lam0;
// //     if(s == 0){ Rprintf("sampled lam0\n"); }
// //     
// //     
// //     // update lam1
// //     lam1 = cppf_gibbs_lam1(logYsur, logZsur, lam0, sig2_sur, 100);
// //     postLam(s, 1) = lam1;
// //     if(s == 0){ Rprintf("sampled lam1\n"); }
// //     
// //     
// //     // update sig2_sur
// //     sig2_sur = cppf_gibbs_sig2_sur(logYsur, logZsur, lam0, lam1, 2, 2);
// //     postSig2(s, 0) = sig2_sur;
// //     if(s == 0){ Rprintf("sampled sig2_sur\n"); }
// //     
// //     
// //     // update beta
// //     if( updateCOV ){
// //       if( (s+1 >= adaptInterval) && (s+1 - (adaptInterval * trunc((s+1) / adaptInterval)) == 0) ){
// //         dummyaccprob = accprob.col(n_al0+1);
// //         rhat[n_al0+1] = sum( dummyaccprob.rows(s+1-adaptInterval, s-1) ) / (adaptInterval-1);
// //         gamma1[n_al0+1] = 1 / pow(adapIter[n_al0+1], c1);
// //         gamma2[n_al0+1] = c0 * gamma1[n_al0+1];
// //         sigma2[n_al0+1] = exp( log(sigma2[n_al0+1]) + gamma2[n_al0+1] * (rhat[n_al0+1] - ropt) );
// //         
// //         COVbeta = COVbeta + gamma1[n_al0+1] * ( cov( postBeta.rows(s+1-adaptInterval, s-1)) - COVbeta );
// //         cholCOVbeta = trans( chol( sigma2[n_al0+1] * ( COVbeta + 0.0000000001 * diagmat(ones(n_be)) ) ) );
// //         
// //         adapIter[n_al0+1] = adapIter[n_al0+1] + 1;
// //       }
// //     }
// //     
// //     newbeta = beta + cholCOVbeta * randn(n_be);
// //     
// //     newloglam = cppf_loglam_b0_gp(Xwhale, logZ, newbeta, psi, indtwhale);
// //     newlam = exp(newloglam);
// //     
// //     logprob = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) -
// //       cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //       cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) -
// //       cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
// //       MVN_logh(newbeta, zeros(n_be), diagmat(ones(n_be))/100) -
// //       MVN_logh(beta, zeros(n_be), diagmat(ones(n_be))/100);
// //     
// //     if (log(randu()) < logprob) {
// //       
// //       beta = newbeta;
// //       loglam = newloglam;
// //       lam = newlam;
// //       
// //       accprob(s, n_al0+1) = 1;
// //     }
// //     
// //     postBeta.row(s) = trans(beta);
// //     if(s == 0){ Rprintf("sampled beta\n"); }
// //     
// //     
// //     // update psi via elliptical slice sampling
// //     llprev = cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //       cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea);
// //     priorpsi = sqrt(kappa_psi) * trans(cholRpsi) * randn(nG);
// //     
// //     thetamin = 0;
// //     thetamax = 2 * M_PI;
// //     theta = thetamin + randu() * (thetamax - thetamin);
// //     thetamin = theta - 2 * M_PI;
// //     thetamax = theta;
// //     
// //     llthred = llprev + log(randu());
// //     accept = false;
// //     count = 0;
// //     
// //     while(accept == false){
// //       count = count + 1;
// //       
// //       newpsi = psi * cos(theta) + priorpsi * sin(theta);
// //       newpsi = newpsi - mean(newpsi);
// //       
// //       newloglam = cppf_loglam_b0_gp(Xwhale, logZ, beta, newpsi, indtwhale);
// //       newlam = exp(newloglam);
// //       
// //       llnew = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //         cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea);
// //       
// //       if(llnew > llthred){
// //         llprev = llnew;
// //         accept = true;
// //       } else {
// //         if(theta < 0){ thetamin = theta; } else { thetamax = theta; }
// //         theta = thetamin + randu() * (thetamax - thetamin);
// //         if ( (count) % 1000 == 0 ) {
// //           Rprintf("ESS for psi: %d iterations...\n", count);
// //         }
// //         if (count > 20000) stop("ESS is not converging...\n");
// //       }
// //     }
// //     psi = newpsi;
// //     loglam = newloglam;
// //     lam = newlam;
// //     
// //     postPsi.row(s) = trans(psi);
// //     if(s == 0){ Rprintf("sampled psi\n"); }
// //     
// //     
// //     // update kappa_psi
// //     kappa_psi = cppf_gibbs_kappa_psi(psi, invRpsi, 2, 2);
// //     postKappa(s,1) = kappa_psi;
// //     if(s == 0){ Rprintf("sampled kappa_psi\n"); }
// //     
// //     
// //     // update c
// //     c = cppf_gibbs_ct(Wpam, lam, dur_dist, p_pam, cellarea, shape_c, scale_c);
// //     postC.row(s) = trans(c);
// //     if(s == 0){ Rprintf("sampled c\n"); }
// //     
// //     
// //     if ( (s+1) % 10 == 0 ) {
// //       Rprintf("Generated %d samples...\n", s+1);
// //     }
// //   }
// //   
// //   
// //   return Rcpp::List::create(Rcpp::Named("postAlpha0") = postAlpha0,
// //                             Rcpp::Named("postAlpha0tilde") = postAlpha0tilde,
// //                             Rcpp::Named("postAlpha") = postAlpha,
// //                             Rcpp::Named("postEta") = postEta,
// //                             Rcpp::Named("postLam") = postLam,
// //                             Rcpp::Named("postKappa") = postKappa,
// //                             Rcpp::Named("postSig2") = postSig2,
// //                             Rcpp::Named("postBeta") = postBeta,
// //                             Rcpp::Named("postPsi") = postPsi,
// //                             Rcpp::Named("postC") = postC,
// //                             Rcpp::Named("Accprob") = accprob,
// //                             Rcpp::Named("sigma2") = sigma2,
// //                             Rcpp::Named("adapIter") = adapIter,
// //                             Rcpp::Named("COValpha0") = COValpha0,
// //                             Rcpp::Named("COValpha") = COValpha,
// //                             Rcpp::Named("COVbeta") = COVbeta);
// //   
// // }
// // 
// // 
// // 
// // // sig2_obl is known
// // // time-varying intercept for whale component
// // // [[Rcpp::export]]
// // List cppf_fitJoint_b0t(
// //     int niter, mat distmat,
// //     List Xzoop, vec logYobl, vec logYsur,
// //     List Xwhale, mat Sdist, mat Wpam, vec dur_dist,
// //     vec f_dist_obs, List f_dist, List p_pam,
// //     mat indobl, mat indsur, vec indtwhale, mat inddist, vec indmdist, vec indmwhale,
// //     double alpha0tilde, vec alpha0, vec alpha, mat eta, double lam0, double lam1,
// //     double tau2, double kappa_eta, double phi_eta, double sig2_obl, double sig2_sur, 
// //     vec beta, vec pii, vec c, double shape_c, double scale_c,
// //     double cellarea,
// //     vec sigma2, vec COValpha0, mat COValpha, mat COVbeta,
// //     bool updateCOV, int adaptInterval, double adaptFactorExponent, vec adapIter){
// //   
// //   double negativeInf = -std::numeric_limits<float>::infinity();;
// //   
// //   indobl = indobl-1, indsur = indsur-1, indtwhale = indtwhale-1;
// //   inddist = inddist-1, indmdist = indmdist-1, indmwhale = indmwhale-1;
// //   
// //   int nG = distmat.n_rows, n_al0 = alpha0.size(), n_al = alpha.size(), n_be = beta.size(), ntzoop = eta.n_cols, n_sigma2 = sigma2.size();
// //   vec dummyvec;
// //   
// //   mat postAlpha0tilde = zeros(niter, 1), postAlpha0 = zeros(niter, n_al0), postAlpha = zeros(niter, n_al), postLam = zeros(niter, 2);
// //   mat postSig2 = zeros(niter, 1), postKappa = zeros(niter, 1);
// //   mat postBeta = zeros(niter, n_be);
// //   mat postC = zeros(niter, c.size()), postEtat = zeros(niter, nG);
// //   List postEta(ntzoop);
// //   for(int i = 0; i < ntzoop; i ++){
// //     postEta[i] = zeros(niter, nG);
// //   }
// //   Rprintf("set up matrices for posterior samples\n");
// //   
// //   mat accprob = zeros(niter, n_sigma2);
// //   vec rhat = zeros(n_sigma2), gamma1 = zeros(n_sigma2), gamma2 = zeros(n_sigma2), dummyaccprob;
// //   double c0 = 1, c1 = adaptFactorExponent, ropt = 0.234, logprob;
// //   vec cholCOValpha0(n_al0);
// //   for(int i = 0; i < n_al0; i ++){ cholCOValpha0[i] = sqrt(sigma2[i] * COValpha0[i]); }
// //   mat cholCOValpha = trans( chol( sigma2[n_al0] * ( COValpha + 0.0000000001 * diagmat(ones(n_al)) ) ) );
// //   mat cholCOVbeta = trans( chol( sigma2[n_al0+1] * ( COVbeta + 0.0000000001 * diagmat(ones(n_be)) ) ) );
// //   Rprintf("set up adaptive proposals\n");
// //   
// //   double alpha0i, newalpha0i;
// //   vec newalpha0, newalpha, newbeta, newlogZobl, newlogZsur;
// //   mat newlogZ, newloglam, newlam;
// //   
// //   int count;
// //   vec xvals;
// //   double llprev, llthred, llnew, thetamin, thetamax, theta;
// //   vec prioretat;
// //   mat neweta;
// //   bool accept;
// //   
// //   mat Reta = exp(-distmat / phi_eta);
// //   mat invReta = inv(Reta),  cholReta = chol(Reta);
// //   rowvec oneinvReta = trans(ones(nG)) * invReta;
// //   double oneinvRetaone = sum( oneinvReta * ones(nG) );
// //   Rprintf("set up GP covariance matrices\n");
// //   
// //   mat logZ = cppf_logZ_a0t_etat(Xzoop, alpha0, alpha, eta);
// //   vec logZobl = cppf_logZobl(logZ, indobl);
// //   vec logZsur = cppf_logZsur(logZ, indsur);
// //   mat loglam = cppf_loglam_b0t(Xwhale, logZ, beta, indtwhale);
// //   mat lam = exp(loglam);
// //   Rprintf("set up initial values\n");
// //   
// //   
// //   // Start MCMC
// //   for(int s = 0; s < niter; s++) {
// //     
// //     // update alpha0
// //     for(int i = 0; i < n_al0; i ++){
// //       
// //       if( updateCOV ){
// //         if( (s+1 >= adaptInterval) && (s+1 - (adaptInterval * trunc((s+1) / adaptInterval)) == 0) ){
// //           dummyaccprob = accprob.col(i);
// //           rhat[i] = sum( dummyaccprob.rows(s+1-adaptInterval, s-1) ) / (adaptInterval-1);
// //           gamma1[i] = 1 / pow(adapIter[i], c1);
// //           gamma2[i] = c0 * gamma1[i];
// //           sigma2[i] = exp( log(sigma2[i]) + gamma2[i] * (rhat[i] - ropt) );
// //           
// //           
// //           dummyvec = postAlpha0.col(i);
// //           COValpha0[i] = COValpha0[i] + gamma1[i] * ( var( dummyvec.rows(s+1-adaptInterval, s-1) ) - COValpha0[i] );
// //           cholCOValpha0[i] = sqrt(  sigma2[i] * COValpha0[i] );
// //           
// //           adapIter[i] = adapIter[i] + 1;
// //         }
// //       }
// //       
// //       alpha0i = alpha0[i];
// //       newalpha0i = alpha0i + cholCOValpha0[i] * randn();
// //       newalpha0 = alpha0;
// //       newalpha0[i] = newalpha0i;
// //       
// //       newlogZ = cppf_logZ_a0t_etat(Xzoop, newalpha0, alpha, eta);
// //       newlogZobl = cppf_logZobl(newlogZ, indobl);
// //       newlogZsur = cppf_logZsur(newlogZ, indsur);
// //       newloglam = cppf_loglam_b0t(Xwhale, newlogZ, beta, indtwhale);
// //       newlam = exp(newloglam);
// //       
// //       logprob = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) -
// //         cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //         cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) -
// //         cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
// //         cppf_loglik_obl(logYobl, newlogZobl, sig2_obl) -
// //         cppf_loglik_obl(logYobl, logZobl, sig2_obl) +
// //         cppf_loglik_sur(logYsur, newlogZsur, lam0, lam1, sig2_sur) -
// //         cppf_loglik_sur(logYsur, logZsur, lam0, lam1, sig2_sur) +
// //         normal_logh(newalpha0i, alpha0tilde, tau2) -
// //         normal_logh(alpha0i, alpha0tilde, tau2);
// //       
// //       if (log(randu()) < logprob) {
// //         
// //         alpha0 = newalpha0;
// //         logZ = newlogZ;
// //         logZobl = newlogZobl;
// //         logZsur = newlogZsur;
// //         loglam = newloglam;
// //         lam = newlam;
// //         
// //         accprob(s, i) = 1;
// //       }
// //     }
// //     postAlpha0.row(s) = trans(alpha0);
// //     if(s == 0){ Rprintf("sampled alpha0\n"); }
// //     
// //     
// //     // Gibbs update for alpha0tilde
// //     alpha0tilde = cppf_gibbs_alpha0tilde(alpha0, tau2);
// //     postAlpha0tilde[s] = alpha0tilde;
// //     if(s == 0){ Rprintf("sampled alpha0tilde\n"); }
// //     
// //     
// //     // update alpha
// //     if( updateCOV ){
// //       if( (s+1 >= adaptInterval) && (s+1 - (adaptInterval * trunc((s+1) / adaptInterval)) == 0) ){
// //         dummyaccprob = accprob.col(n_al0);
// //         rhat[n_al0] = sum( dummyaccprob.rows(s+1-adaptInterval, s-1) ) / (adaptInterval-1);
// //         gamma1[n_al0] = 1 / pow(adapIter[n_al0], c1);
// //         gamma2[n_al0] = c0 * gamma1[n_al0];
// //         sigma2[n_al0] = exp( log(sigma2[n_al0]) + gamma2[n_al0] * (rhat[n_al0] - ropt) );
// //         
// //         COValpha = COValpha + gamma1[n_al0] * ( cov( postAlpha.rows(s+1-adaptInterval, s-1) ) - COValpha );
// //         cholCOValpha = trans( chol( sigma2[n_al0] * ( COValpha + 0.0000000001 * diagmat(ones(n_al)) ) ) );
// //         
// //         adapIter[n_al0] = adapIter[n_al0] + 1;
// //       }
// //     }
// //     
// //     newalpha = alpha + cholCOValpha * randn(n_al);
// //     
// //     newlogZ = cppf_logZ_a0t_etat(Xzoop, alpha0, newalpha, eta);
// //     newlogZobl = cppf_logZobl(newlogZ, indobl);
// //     newlogZsur = cppf_logZsur(newlogZ, indsur);
// //     newloglam = cppf_loglam_b0t(Xwhale, newlogZ, beta, indtwhale);
// //     newlam = exp(newloglam);
// //     
// //     logprob = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) -
// //       cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //       cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) -
// //       cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
// //       cppf_loglik_obl(logYobl, newlogZobl, sig2_obl) -
// //       cppf_loglik_obl(logYobl, logZobl, sig2_obl) +
// //       cppf_loglik_sur(logYsur, newlogZsur, lam0, lam1, sig2_sur) -
// //       cppf_loglik_sur(logYsur, logZsur, lam0, lam1, sig2_sur) +
// //       MVN_logh(newalpha, zeros(n_al), diagmat(ones(n_al))/100) -
// //       MVN_logh(alpha, zeros(n_al), diagmat(ones(n_al))/100);
// //     
// //     if (log(randu()) < logprob) {
// //       
// //       alpha = newalpha;
// //       logZ = newlogZ;
// //       logZobl = newlogZobl;
// //       logZsur = newlogZsur;
// //       loglam = newloglam;
// //       lam = newlam;
// //       
// //       accprob(s, n_al0) = 1;
// //     }
// //     
// //     postAlpha.row(s) = trans(alpha);
// //     if(s == 0){ Rprintf("sampled alpha\n"); }
// //     
// //     
// //     // update eta via Elliptical slice sampling
// //     for(int t = 0; t < ntzoop; t ++){
// //       neweta = eta;
// //       
// //       llprev = cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //         cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
// //         cppf_loglik_obl(logYobl, logZobl, sig2_obl) +
// //         cppf_loglik_sur(logYsur, logZsur, lam0, lam1, sig2_sur);
// //       prioretat = sqrt(kappa_eta) * trans(cholReta) * randn(nG);
// //       
// //       thetamin = 0;
// //       thetamax = 2 * M_PI;
// //       theta = thetamin + randu() * (thetamax - thetamin);
// //       thetamin = theta - 2 * M_PI;
// //       thetamax = theta;
// //       
// //       llthred = llprev + log(randu());
// //       accept = false;
// //       count = 0;
// //       
// //       while(accept == false){
// //         count = count + 1;
// //         
// //         neweta.col(t) = eta.col(t) * cos(theta) + prioretat * sin(theta);
// //         neweta.col(t) = neweta.col(t) - mean(neweta.col(t));
// //         
// //         newlogZ = cppf_logZ_a0t_etat(Xzoop, alpha0, alpha, neweta);
// //         newlogZobl = cppf_logZobl(newlogZ, indobl);
// //         newlogZsur = cppf_logZsur(newlogZ, indsur);
// //         newloglam = cppf_loglam_b0t(Xwhale, newlogZ, beta, indtwhale);
// //         newlam = exp(newloglam);
// //         
// //         llnew = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //           cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) +
// //           cppf_loglik_obl(logYobl, newlogZobl, sig2_obl) +
// //           cppf_loglik_sur(logYsur, newlogZsur, lam0, lam1, sig2_sur);
// //         
// //         if(llnew > llthred){
// //           llprev = llnew;
// //           accept = true;
// //         } else {
// //           if(theta < 0){ thetamin = theta; } else { thetamax = theta; }
// //           theta = thetamin + randu() * (thetamax - thetamin);
// //           if ( (count) % 1000 == 0 ) {
// //             Rprintf("ESS for eta: %d iterations...\n", count);
// //           }
// //           if (count > 20000) stop("ESS is not converging...\n");
// //         }
// //       }
// //       eta = neweta;
// //       logZ = newlogZ;
// //       logZobl = newlogZobl;
// //       logZsur = newlogZsur;
// //       loglam = newloglam;
// //       lam = newlam;
// //       
// //       postEtat = as<mat>(postEta[t]);
// //       postEtat.row(s) = trans(eta.col(t));
// //       postEta[t] = postEtat;
// //     }
// //     if(s == 0){ Rprintf("sampled eta\n"); }
// //     
// //     
// //     // update kappa_eta
// //     kappa_eta = cppf_gibbs_kappa_eta(eta, invReta, 2, 2);
// //     postKappa(s,0) = kappa_eta;
// //     if(s == 0){ Rprintf("sampled kappa_eta\n"); }
// //     
// //     
// //     // update lam0
// //     lam0 = cppf_gibbs_lam0(logYsur, logZsur, lam1, sig2_sur, 100);
// //     postLam(s, 0) = lam0;
// //     if(s == 0){ Rprintf("sampled lam0\n"); }
// //     
// //     
// //     // update lam1
// //     lam1 = cppf_gibbs_lam1(logYsur, logZsur, lam0, sig2_sur, 100);
// //     postLam(s, 1) = lam1;
// //     if(s == 0){ Rprintf("sampled lam1\n"); }
// //     
// //     
// //     // update sig2_sur
// //     sig2_sur = cppf_gibbs_sig2_sur(logYsur, logZsur, lam0, lam1, 2, 2);
// //     postSig2(s, 0) = sig2_sur;
// //     if(s == 0){ Rprintf("sampled sig2_sur\n"); }
// //     
// //     
// //     // update beta
// //     if( updateCOV ){
// //       if( (s+1 >= adaptInterval) && (s+1 - (adaptInterval * trunc((s+1) / adaptInterval)) == 0) ){
// //         dummyaccprob = accprob.col(n_al0+1);
// //         rhat[n_al0+1] = sum( dummyaccprob.rows(s+1-adaptInterval, s-1) ) / (adaptInterval-1);
// //         gamma1[n_al0+1] = 1 / pow(adapIter[n_al0+1], c1);
// //         gamma2[n_al0+1] = c0 * gamma1[n_al0+1];
// //         sigma2[n_al0+1] = exp( log(sigma2[n_al0+1]) + gamma2[n_al0+1] * (rhat[n_al0+1] - ropt) );
// //         
// //         COVbeta = COVbeta + gamma1[n_al0+1] * ( cov( postBeta.rows(s+1-adaptInterval, s-1)) - COVbeta );
// //         cholCOVbeta = trans( chol( sigma2[n_al0+1] * ( COVbeta + 0.0000000001 * diagmat(ones(n_be)) ) ) );
// //         
// //         adapIter[n_al0+1] = adapIter[n_al0+1] + 1;
// //       }
// //     }
// //     
// //     newbeta = beta + cholCOVbeta * randn(n_be);
// //     
// //     newloglam = cppf_loglam_b0t(Xwhale, logZ, newbeta, indtwhale);
// //     newlam = exp(newloglam);
// //     
// //     logprob = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) -
// //       cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //       cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) -
// //       cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
// //       MVN_logh(newbeta, zeros(n_be), diagmat(ones(n_be))/100) -
// //       MVN_logh(beta, zeros(n_be), diagmat(ones(n_be))/100);
// //     
// //     if (log(randu()) < logprob) {
// //       
// //       beta = newbeta;
// //       loglam = newloglam;
// //       lam = newlam;
// //       
// //       accprob(s, n_al0+1) = 1;
// //     }
// //     
// //     postBeta.row(s) = trans(beta);
// //     if(s == 0){ Rprintf("sampled beta\n"); }
// //     
// //     
// //     // update c
// //     c = cppf_gibbs_ct(Wpam, lam, dur_dist, p_pam, cellarea, shape_c, scale_c);
// //     postC.row(s) = trans(c);
// //     if(s == 0){ Rprintf("sampled c\n"); }
// //     
// //     
// //     if ( (s+1) % 10 == 0 ) {
// //       Rprintf("Generated %d samples...\n", s+1);
// //     }
// //   }
// //   
// //   
// //   return Rcpp::List::create(Rcpp::Named("postAlpha0") = postAlpha0,
// //                             Rcpp::Named("postAlpha0tilde") = postAlpha0tilde,
// //                             Rcpp::Named("postAlpha") = postAlpha,
// //                             Rcpp::Named("postEta") = postEta,
// //                             Rcpp::Named("postLam") = postLam,
// //                             Rcpp::Named("postKappa") = postKappa,
// //                             Rcpp::Named("postSig2") = postSig2,
// //                             Rcpp::Named("postBeta") = postBeta,
// //                             Rcpp::Named("postC") = postC,
// //                             Rcpp::Named("Accprob") = accprob,
// //                             Rcpp::Named("sigma2") = sigma2,
// //                             Rcpp::Named("adapIter") = adapIter,
// //                             Rcpp::Named("COValpha0") = COValpha0,
// //                             Rcpp::Named("COValpha") = COValpha,
// //                             Rcpp::Named("COVbeta") = COVbeta);
// //   
// // }
// // 
// // 
// // 
// // // [[Rcpp::export]]
// // List cppf_fitJoint_b0t_pref(
// //     int niter, mat distmat,
// //     List Xzoop, vec logYobl, vec logYsur, mat Spref,
// //     List Xwhale, mat Sdist, mat Wpam, vec dur_dist,
// //     vec f_dist_obs, List f_dist, List p_pam,
// //     mat indobl, mat indsur, vec indpref, vec indtwhale, mat inddist, vec indmdist, vec indmwhale,
// //     double alpha0tilde, vec alpha0, vec alpha, mat eta, 
// //     double mu0, vec delta, double lam0, double lam1,
// //     double tau2, double kappa_eta, double phi_eta, double kappa_delta, double phi_delta,
// //     double sig2_obl, double sig2_sur, 
// //     vec beta, vec pii, vec c, double shape_c, double scale_c,
// //     double cellarea,
// //     vec sigma2, vec COValpha0, mat COValpha, mat COVbeta, double COVmu0,
// //     bool updateCOV, int adaptInterval, double adaptFactorExponent, vec adapIter){
// //   
// //   double negativeInf = -std::numeric_limits<float>::infinity();;
// //   
// //   indobl = indobl-1, indsur = indsur-1, indpref = indpref-1, indtwhale = indtwhale-1;
// //   inddist = inddist-1, indmdist = indmdist-1, indmwhale = indmwhale-1;
// //   
// //   int nG = distmat.n_rows, n_al0 = alpha0.size(), n_al = alpha.size(), n_be = beta.size(), ntzoop = eta.n_cols, n_sigma2 = sigma2.size();
// //   vec dummyvec, postj;
// //   
// //   mat postAlpha0tilde = zeros(niter, 1), postAlpha0 = zeros(niter, n_al0), postAlpha = zeros(niter, n_al);
// //   mat postDelta = zeros(niter, delta.size()), postMu0(niter, 1), postLam = zeros(niter, 2);
// //   mat postSig2 = zeros(niter, 1), postKappa = zeros(niter, 2);
// //   mat postBeta = zeros(niter, n_be);
// //   mat postC = zeros(niter, c.size()), postEtat = zeros(niter, nG);
// //   List postEta(ntzoop);
// //   for(int i = 0; i < ntzoop; i ++){
// //     postEta[i] = zeros(niter, nG);
// //   }
// //   Rprintf("set up matrices for posterior samples\n");
// //   
// //   mat accprob = zeros(niter, n_sigma2);
// //   vec rhat = zeros(n_sigma2), gamma1 = zeros(n_sigma2), gamma2 = zeros(n_sigma2), dummyaccprob;
// //   double c0 = 1, c1 = adaptFactorExponent, ropt = 0.234, logprob;
// //   vec cholCOValpha0(n_al0);
// //   for(int i = 0; i < n_al0; i ++){ cholCOValpha0[i] = sqrt(sigma2[i] * COValpha0[i]); }
// //   mat cholCOValpha = trans( chol( sigma2[n_al0] * ( COValpha + 0.0000000001 * diagmat(ones(n_al)) ) ) );
// //   mat cholCOVbeta = trans( chol( sigma2[n_al0+1] * ( COVbeta + 0.0000000001 * diagmat(ones(n_be)) ) ) );
// //   double cholCOVmu0 = sqrt( sigma2[n_al0+2] * COVmu0 );
// //   Rprintf("set up adaptive proposals\n");
// //   
// //   double alpha0i, newalpha0i, newmu0;
// //   vec newalpha0, newalpha, newbeta, newlogZobl, newlogZsur, newdelta, newloglam_pref, newlam_pref;
// //   mat newlogZ, newloglam, newlam;
// //   
// //   int count;
// //   vec xvals;
// //   double llprev, llthred, llnew, thetamin, thetamax, theta;
// //   vec prioretat, priordelta;
// //   mat neweta;
// //   bool accept;
// //   
// //   mat Reta = exp(-distmat / phi_eta), Rdelta = exp(-distmat / phi_delta);
// //   mat invReta = inv(Reta), invRdelta= inv(Rdelta), cholReta = chol(Reta), cholRdelta = chol(Rdelta);
// //   rowvec oneinvReta = trans(ones(nG)) * invReta, oneinvRdelta = trans(ones(nG)) * invRdelta;
// //   double oneinvRetaone = sum( oneinvReta * ones(nG) ), oneinvRdeltaone = sum( oneinvRdelta * ones(nG) );
// //   Rprintf("set up GP covariance matrices\n");
// //   
// //   mat logZ = cppf_logZ_a0t_etat_pref(Xzoop, alpha0, alpha, eta, delta);
// //   vec logZobl = cppf_logZobl(logZ, indobl);
// //   vec logZsur = cppf_logZsur(logZ, indsur);
// //   mat loglam = cppf_loglam_b0t(Xwhale, logZ, beta, indtwhale);
// //   mat lam = exp(loglam);
// //   vec loglam_pref = cppf_loglam_pref(mu0, delta);
// //   vec lam_pref = exp(loglam_pref);
// //   Rprintf("set up initial values\n");
// //   
// //   
// //   // Start MCMC
// //   for(int s = 0; s < niter; s++) {
// //     
// //     // update alpha0
// //     for(int i = 0; i < n_al0; i ++){
// //       
// //       if( updateCOV ){
// //         if( (s+1 >= adaptInterval) && (s+1 - (adaptInterval * trunc((s+1) / adaptInterval)) == 0) ){
// //           dummyaccprob = accprob.col(i);
// //           rhat[i] = sum( dummyaccprob.rows(s+1-adaptInterval, s-1) ) / (adaptInterval-1);
// //           gamma1[i] = 1 / pow(adapIter[i], c1);
// //           gamma2[i] = c0 * gamma1[i];
// //           sigma2[i] = exp( log(sigma2[i]) + gamma2[i] * (rhat[i] - ropt) );
// //           
// //           
// //           dummyvec = postAlpha0.col(i);
// //           COValpha0[i] = COValpha0[i] + gamma1[i] * ( var( dummyvec.rows(s+1-adaptInterval, s-1) ) - COValpha0[i] );
// //           cholCOValpha0[i] = sqrt(  sigma2[i] * COValpha0[i] );
// //           
// //           adapIter[i] = adapIter[i] + 1;
// //         }
// //       }
// //       
// //       alpha0i = alpha0[i];
// //       newalpha0i = alpha0i + cholCOValpha0[i] * randn();
// //       newalpha0 = alpha0;
// //       newalpha0[i] = newalpha0i;
// //       
// //       newlogZ = cppf_logZ_a0t_etat_pref(Xzoop, newalpha0, alpha, eta, delta);
// //       newlogZobl = cppf_logZobl(newlogZ, indobl);
// //       newlogZsur = cppf_logZsur(newlogZ, indsur);
// //       newloglam = cppf_loglam_b0t(Xwhale, newlogZ, beta, indtwhale);
// //       newlam = exp(newloglam);
// //       
// //       logprob = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) -
// //         cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //         cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) -
// //         cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
// //         cppf_loglik_obl(logYobl, newlogZobl, sig2_obl) -
// //         cppf_loglik_obl(logYobl, logZobl, sig2_obl) +
// //         cppf_loglik_sur(logYsur, newlogZsur, lam0, lam1, sig2_sur) -
// //         cppf_loglik_sur(logYsur, logZsur, lam0, lam1, sig2_sur) +
// //         normal_logh(newalpha0i, alpha0tilde, tau2) -
// //         normal_logh(alpha0i, alpha0tilde, tau2);
// //       
// //       if (log(randu()) < logprob) {
// //         
// //         alpha0 = newalpha0;
// //         logZ = newlogZ;
// //         logZobl = newlogZobl;
// //         logZsur = newlogZsur;
// //         loglam = newloglam;
// //         lam = newlam;
// //         
// //         accprob(s, i) = 1;
// //       }
// //     }
// //     postAlpha0.row(s) = trans(alpha0);
// //     if(s == 0){ Rprintf("sampled alpha0\n"); }
// //     
// //     
// //     // Gibbs update for alpha0tilde
// //     alpha0tilde = cppf_gibbs_alpha0tilde(alpha0, tau2);
// //     postAlpha0tilde[s] = alpha0tilde;
// //     if(s == 0){ Rprintf("sampled alpha0tilde\n"); }
// //     
// //     
// //     // update alpha
// //     if( updateCOV ){
// //       if( (s+1 >= adaptInterval) && (s+1 - (adaptInterval * trunc((s+1) / adaptInterval)) == 0) ){
// //         dummyaccprob = accprob.col(n_al0);
// //         rhat[n_al0] = sum( dummyaccprob.rows(s+1-adaptInterval, s-1) ) / (adaptInterval-1);
// //         gamma1[n_al0] = 1 / pow(adapIter[n_al0], c1);
// //         gamma2[n_al0] = c0 * gamma1[n_al0];
// //         sigma2[n_al0] = exp( log(sigma2[n_al0]) + gamma2[n_al0] * (rhat[n_al0] - ropt) );
// //         
// //         COValpha = COValpha + gamma1[n_al0] * ( cov( postAlpha.rows(s+1-adaptInterval, s-1) ) - COValpha );
// //         cholCOValpha = trans( chol( sigma2[n_al0] * ( COValpha + 0.0000000001 * diagmat(ones(n_al)) ) ) );
// //         
// //         adapIter[n_al0] = adapIter[n_al0] + 1;
// //       }
// //     }
// //     
// //     newalpha = alpha + cholCOValpha * randn(n_al);
// //     
// //     newlogZ = cppf_logZ_a0t_etat_pref(Xzoop, alpha0, newalpha, eta, delta);
// //     newlogZobl = cppf_logZobl(newlogZ, indobl);
// //     newlogZsur = cppf_logZsur(newlogZ, indsur);
// //     newloglam = cppf_loglam_b0t(Xwhale, newlogZ, beta, indtwhale);
// //     newlam = exp(newloglam);
// //     
// //     logprob = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) -
// //       cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //       cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) -
// //       cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
// //       cppf_loglik_obl(logYobl, newlogZobl, sig2_obl) -
// //       cppf_loglik_obl(logYobl, logZobl, sig2_obl) +
// //       cppf_loglik_sur(logYsur, newlogZsur, lam0, lam1, sig2_sur) -
// //       cppf_loglik_sur(logYsur, logZsur, lam0, lam1, sig2_sur) +
// //       MVN_logh(newalpha, zeros(n_al), diagmat(ones(n_al))/100) -
// //       MVN_logh(alpha, zeros(n_al), diagmat(ones(n_al))/100);
// //     
// //     if (log(randu()) < logprob) {
// //       
// //       alpha = newalpha;
// //       logZ = newlogZ;
// //       logZobl = newlogZobl;
// //       logZsur = newlogZsur;
// //       loglam = newloglam;
// //       lam = newlam;
// //       
// //       accprob(s, n_al0) = 1;
// //     }
// //     
// //     postAlpha.row(s) = trans(alpha);
// //     if(s == 0){ Rprintf("sampled alpha\n"); }
// //     
// //     
// //     // update eta via Elliptical slice sampling
// //     for(int t = 0; t < ntzoop; t ++){
// //       neweta = eta;
// //       
// //       llprev = cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //         cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
// //         cppf_loglik_obl(logYobl, logZobl, sig2_obl) +
// //         cppf_loglik_sur(logYsur, logZsur, lam0, lam1, sig2_sur);
// //       prioretat = sqrt(kappa_eta) * trans(cholReta) * randn(nG);
// //       
// //       thetamin = 0;
// //       thetamax = 2 * M_PI;
// //       theta = thetamin + randu() * (thetamax - thetamin);
// //       thetamin = theta - 2 * M_PI;
// //       thetamax = theta;
// //       
// //       llthred = llprev + log(randu());
// //       accept = false;
// //       count = 0;
// //       
// //       while(accept == false){
// //         count = count + 1;
// //         
// //         neweta.col(t) = eta.col(t) * cos(theta) + prioretat * sin(theta);
// //         neweta.col(t) = neweta.col(t) - mean(neweta.col(t));
// //         
// //         newlogZ = cppf_logZ_a0t_etat_pref(Xzoop, alpha0, alpha, neweta, delta);
// //         newlogZobl = cppf_logZobl(newlogZ, indobl);
// //         newlogZsur = cppf_logZsur(newlogZ, indsur);
// //         newloglam = cppf_loglam_b0t(Xwhale, newlogZ, beta, indtwhale);
// //         newlam = exp(newloglam);
// //         
// //         llnew = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //           cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) +
// //           cppf_loglik_obl(logYobl, newlogZobl, sig2_obl) +
// //           cppf_loglik_sur(logYsur, newlogZsur, lam0, lam1, sig2_sur);
// //         
// //         if(llnew > llthred){
// //           llprev = llnew;
// //           accept = true;
// //         } else {
// //           if(theta < 0){ thetamin = theta; } else { thetamax = theta; }
// //           theta = thetamin + randu() * (thetamax - thetamin);
// //           if ( (count) % 1000 == 0 ) {
// //             Rprintf("ESS for eta: %d iterations...\n", count);
// //           }
// //           if (count > 20000) stop("ESS is not converging...\n");
// //         }
// //       }
// //       eta = neweta;
// //       logZ = newlogZ;
// //       logZobl = newlogZobl;
// //       logZsur = newlogZsur;
// //       loglam = newloglam;
// //       lam = newlam;
// //       
// //       postEtat = as<mat>(postEta[t]);
// //       postEtat.row(s) = trans(eta.col(t));
// //       postEta[t] = postEtat;
// //     }
// //     if(s == 0){ Rprintf("sampled eta\n"); }
// //     
// //     
// //     // update kappa_eta
// //     kappa_eta = cppf_gibbs_kappa_eta(eta, invReta, 2, 2);
// //     postKappa(s,0) = kappa_eta;
// //     if(s == 0){ Rprintf("sampled kappa_eta\n"); }
// //     
// //     
// //     // update lam0
// //     lam0 = cppf_gibbs_lam0(logYsur, logZsur, lam1, sig2_sur, 100);
// //     postLam(s, 0) = lam0;
// //     if(s == 0){ Rprintf("sampled lam0\n"); }
// //     
// //     
// //     // update lam1
// //     lam1 = cppf_gibbs_lam1(logYsur, logZsur, lam0, sig2_sur, 100);
// //     postLam(s, 1) = lam1;
// //     if(s == 0){ Rprintf("sampled lam1\n"); }
// //     
// //     
// //     // update sig2_sur
// //     sig2_sur = cppf_gibbs_sig2_sur(logYsur, logZsur, lam0, lam1, 2, 2);
// //     postSig2(s, 0) = sig2_sur;
// //     if(s == 0){ Rprintf("sampled sig2_sur\n"); }
// //     
// //     
// //     // update mu0
// //     if( updateCOV ){
// //       if( (s+1 >= adaptInterval) && (s+1 - (adaptInterval * trunc((s+1) / adaptInterval)) == 0) ){
// //         dummyaccprob = accprob.col(n_al0+2);
// //         rhat[n_al0+2] = sum( dummyaccprob.rows(s+1-adaptInterval, s-1) ) / (adaptInterval-1);
// //         gamma1[n_al0+2] = 1 / pow(adapIter[n_al0+2], c1);
// //         gamma2[n_al0+2] = c0 * gamma1[n_al0+2];
// //         sigma2[n_al0+2] = exp( log(sigma2[n_al0+2]) + gamma2[n_al0+2] * (rhat[n_al0+2] - ropt) );
// //         
// //         postj = postMu0.col(0);
// //         COVmu0 = COVmu0 + gamma1[n_al0+2] * ( var( postj.rows(s+1-adaptInterval, s-1)) - COVmu0 );
// //         cholCOVmu0 = sqrt( sigma2[n_al0+2] * COVmu0 );
// //         
// //         adapIter[n_al0+2] = adapIter[n_al0+2] + 1;
// //       }
// //     }
// //     
// //     newmu0 = mu0 + cholCOVmu0 * randn();
// //     
// //     newloglam_pref = cppf_loglam_pref(newmu0, delta);
// //     newlam_pref = exp(newloglam_pref);
// //     
// //     logprob = cppf_loglik_pref(newloglam_pref, newlam_pref, indpref, cellarea) - 
// //       cppf_loglik_pref(loglam_pref, lam_pref, indpref, cellarea) +
// //       normal_logh(newmu0, 0, 100) - normal_logh(mu0, 0, 100);
// //     
// //     if (log(randu()) < logprob) {
// //       
// //       mu0 = newmu0;
// //       loglam_pref = newloglam_pref;
// //       lam_pref = newlam_pref;
// //       
// //       accprob(s, n_al0+2) = 1;
// //     }
// //     
// //     postMu0(s,0) = mu0;
// //     if(s == 0){ Rprintf("sampled mu0\n"); }
// //     
// //     
// //     // update delta
// //     llprev = cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //       cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
// //       cppf_loglik_obl(logYobl, logZobl, sig2_obl) +
// //       cppf_loglik_sur(logYsur, logZsur, lam0, lam1, sig2_sur) +
// //       cppf_loglik_pref(loglam_pref, lam_pref, indpref, cellarea);
// //     priordelta = sqrt(kappa_delta) * trans(cholRdelta) * randn(nG);
// //     
// //     thetamin = 0;
// //     thetamax = 2 * M_PI;
// //     theta = thetamin + randu() * (thetamax - thetamin);
// //     thetamin = theta - 2 * M_PI;
// //     thetamax = theta;
// //     
// //     llthred = llprev + log(randu());
// //     accept = false;
// //     count = 0;
// //     
// //     while(accept == false){
// //       count = count + 1;
// //       
// //       newdelta = delta * cos(theta) + priordelta * sin(theta);
// //       newdelta = newdelta - mean(newdelta);
// //       
// //       newlogZ = cppf_logZ_a0t_etat_pref(Xzoop, alpha0, alpha, eta, newdelta);
// //       newlogZobl = cppf_logZobl(newlogZ, indobl);
// //       newlogZsur = cppf_logZsur(newlogZ, indsur);
// //       newloglam = cppf_loglam_b0t(Xwhale, newlogZ, beta, indtwhale);
// //       newlam = exp(newloglam);
// //       newloglam_pref = cppf_loglam_pref(mu0, newdelta);
// //       newlam_pref = exp(newloglam_pref);
// //       
// //       llnew = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //         cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) +
// //         cppf_loglik_obl(logYobl, newlogZobl, sig2_obl) +
// //         cppf_loglik_sur(logYsur, newlogZsur, lam0, lam1, sig2_sur) +
// //         cppf_loglik_pref(newloglam_pref, newlam_pref, indpref, cellarea);
// //       
// //       if(llnew > llthred){
// //         llprev = llnew;
// //         accept = true;
// //       } else {
// //         if(theta < 0){ thetamin = theta; } else { thetamax = theta; }
// //         theta = thetamin + randu() * (thetamax - thetamin);
// //         if ( (count) % 1000 == 0 ) {
// //           Rprintf("ESS for delta: %d iterations...\n", count);
// //         }
// //         if (count > 20000) stop("ESS is not converging...\n");
// //       }
// //     }
// //     delta = newdelta;
// //     logZ = newlogZ;
// //     logZobl = newlogZobl;
// //     logZsur = newlogZsur;
// //     loglam = newloglam;
// //     lam = newlam;
// //     loglam_pref = newloglam_pref;
// //     lam_pref = newlam_pref;
// //     
// //     postDelta.row(s) = trans(delta);
// //     if(s == 0){ Rprintf("sampled delta\n"); }
// //     
// //     
// //     // update kappa_delta
// //     kappa_delta = cppf_gibbs_kappa_delta(delta, invRdelta, 2, 2);
// //     postKappa(s, 1) = kappa_delta;
// //     if(s == 0){ Rprintf("sampled kappa_delta\n"); }
// //     
// //     
// //     // update beta
// //     if( updateCOV ){
// //       if( (s+1 >= adaptInterval) && (s+1 - (adaptInterval * trunc((s+1) / adaptInterval)) == 0) ){
// //         dummyaccprob = accprob.col(n_al0+1);
// //         rhat[n_al0+1] = sum( dummyaccprob.rows(s+1-adaptInterval, s-1) ) / (adaptInterval-1);
// //         gamma1[n_al0+1] = 1 / pow(adapIter[n_al0+1], c1);
// //         gamma2[n_al0+1] = c0 * gamma1[n_al0+1];
// //         sigma2[n_al0+1] = exp( log(sigma2[n_al0+1]) + gamma2[n_al0+1] * (rhat[n_al0+1] - ropt) );
// //         
// //         COVbeta = COVbeta + gamma1[n_al0+1] * ( cov( postBeta.rows(s+1-adaptInterval, s-1)) - COVbeta );
// //         cholCOVbeta = trans( chol( sigma2[n_al0+1] * ( COVbeta + 0.0000000001 * diagmat(ones(n_be)) ) ) );
// //         
// //         adapIter[n_al0+1] = adapIter[n_al0+1] + 1;
// //       }
// //     }
// //     
// //     newbeta = beta + cholCOVbeta * randn(n_be);
// //     
// //     newloglam = cppf_loglam_b0t(Xwhale, logZ, newbeta, indtwhale);
// //     newlam = exp(newloglam);
// //     
// //     logprob = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) -
// //       cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //       cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) -
// //       cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
// //       MVN_logh(newbeta, zeros(n_be), diagmat(ones(n_be))/100) -
// //       MVN_logh(beta, zeros(n_be), diagmat(ones(n_be))/100);
// //     
// //     if (log(randu()) < logprob) {
// //       
// //       beta = newbeta;
// //       loglam = newloglam;
// //       lam = newlam;
// //       
// //       accprob(s, n_al0+1) = 1;
// //     }
// //     
// //     postBeta.row(s) = trans(beta);
// //     if(s == 0){ Rprintf("sampled beta\n"); }
// //     
// //     
// //     // update c
// //     c = cppf_gibbs_ct(Wpam, lam, dur_dist, p_pam, cellarea, shape_c, scale_c);
// //     postC.row(s) = trans(c);
// //     if(s == 0){ Rprintf("sampled c\n"); }
// //     
// //     
// //     if ( (s+1) % 10 == 0 ) {
// //       Rprintf("Generated %d samples...\n", s+1);
// //     }
// //   }
// //   
// //   
// //   return Rcpp::List::create(Rcpp::Named("postAlpha0") = postAlpha0,
// //                             Rcpp::Named("postAlpha0tilde") = postAlpha0tilde,
// //                             Rcpp::Named("postAlpha") = postAlpha,
// //                             Rcpp::Named("postEta") = postEta,
// //                             Rcpp::Named("postMu0") = postMu0,
// //                             Rcpp::Named("postDelta") = postDelta,
// //                             Rcpp::Named("postLam") = postLam,
// //                             Rcpp::Named("postKappa") = postKappa,
// //                             Rcpp::Named("postSig2") = postSig2,
// //                             Rcpp::Named("postBeta") = postBeta,
// //                             Rcpp::Named("postC") = postC,
// //                             Rcpp::Named("Accprob") = accprob,
// //                             Rcpp::Named("sigma2") = sigma2,
// //                             Rcpp::Named("adapIter") = adapIter,
// //                             Rcpp::Named("COValpha0") = COValpha0,
// //                             Rcpp::Named("COValpha") = COValpha,
// //                             Rcpp::Named("COVbeta") = COVbeta,
// //                             Rcpp::Named("COVmu0") = COVmu0);
// //   
// // }
// // 
// // 
// // 
// // // sig2_obl is known
// // // time-varying intercept for whale component
// // // annual GP for whale component
// // // [[Rcpp::export]]
// // List cppf_fitJoint_b0t_psi(
// //     int niter, mat distmat,
// //     List Xzoop, vec logYobl, vec logYsur,
// //     List Xwhale, mat Sdist, mat Wpam, vec dur_dist,
// //     vec f_dist_obs, List f_dist, List p_pam,
// //     mat indobl, mat indsur, vec indtwhale, mat inddist, vec indmdist, vec indmwhale,
// //     double alpha0tilde, vec alpha0, vec alpha, mat eta, double lam0, double lam1,
// //     double tau2, double kappa_eta, double phi_eta, double sig2_obl, double sig2_sur,
// //     vec beta, vec psi, vec pii, vec c, double shape_c, double scale_c,
// //     double kappa_psi, double phi_psi,
// //     double cellarea,
// //     vec sigma2, vec COValpha0, mat COValpha, mat COVbeta,
// //     bool updateCOV, int adaptInterval, double adaptFactorExponent, vec adapIter){
// //   
// //   double negativeInf = -std::numeric_limits<float>::infinity();;
// //   
// //   indobl = indobl-1, indsur = indsur-1, indtwhale = indtwhale-1;
// //   inddist = inddist-1, indmdist = indmdist-1, indmwhale = indmwhale-1;
// //   
// //   int nG = distmat.n_rows, n_al0 = alpha0.size(), n_al = alpha.size(), n_be = beta.size(), ntzoop = eta.n_cols, n_sigma2 = sigma2.size();
// //   
// //   mat postAlpha0tilde = zeros(niter, 1), postAlpha0 = zeros(niter, n_al0), postAlpha = zeros(niter, n_al), postLam = zeros(niter, 2);
// //   mat postSig2 = zeros(niter, 1), postKappa = zeros(niter, 2);
// //   mat postBeta = zeros(niter, n_be);
// //   mat postPsi(niter, nG), postC = zeros(niter, c.size()), postEtat = zeros(niter, nG);
// //   List postEta(ntzoop);
// //   for(int i = 0; i < ntzoop; i ++){
// //     postEta[i] = zeros(niter, nG);
// //   }
// //   Rprintf("set up matrices for posterior samples\n");
// //   
// //   mat accprob = zeros(niter, n_sigma2);
// //   vec rhat = zeros(n_sigma2), gamma1 = zeros(n_sigma2), gamma2 = zeros(n_sigma2), dummyaccprob;
// //   double c0 = 1, c1 = adaptFactorExponent, ropt = 0.234, logprob;
// //   vec cholCOValpha0(n_al0);
// //   for(int i = 0; i < n_al0; i ++){ cholCOValpha0[i] = sqrt(sigma2[i] * COValpha0[i]); }
// //   mat cholCOValpha = trans( chol( sigma2[n_al0] * ( COValpha + 0.0000000001 * diagmat(ones(n_al)) ) ) );
// //   mat cholCOVbeta = trans( chol( sigma2[n_al0+1] * ( COVbeta + 0.0000000001 * diagmat(ones(n_be)) ) ) );
// //   Rprintf("set up adaptive proposals\n");
// //   
// //   double alpha0i, newalpha0i, beta0i, newbeta0i;
// //   vec newalpha0, newalpha, newbeta0, newbeta, newlogZobl, newlogZsur;
// //   mat newlogZ, newloglam, newlam;
// //   vec dummyvec;
// //   
// //   int count;
// //   vec xvals;
// //   double llprev, llthred, llnew, thetamin, thetamax, theta;
// //   vec prioretat, priorpsi, newpsi;
// //   mat neweta;
// //   bool accept;
// //   
// //   mat Reta = exp(-distmat / phi_eta), Rpsi = exp(-distmat / phi_psi);
// //   mat invReta = inv(Reta), invRpsi= inv(Rpsi), cholReta = chol(Reta), cholRpsi= chol(Rpsi);
// //   rowvec oneinvReta = trans(ones(nG)) * invReta, oneinvRpsi = trans(ones(nG)) * invRpsi;
// //   double oneinvRetaone = sum( oneinvReta * ones(nG) ), oneinvRpsione = sum( oneinvRpsi * ones(nG) );
// //   Rprintf("set up GP covariance matrices\n");
// //   
// //   mat logZ = cppf_logZ_a0t_etat(Xzoop, alpha0, alpha, eta);
// //   vec logZobl = cppf_logZobl(logZ, indobl);
// //   vec logZsur = cppf_logZsur(logZ, indsur);
// //   mat loglam = cppf_loglam_b0t_psi(Xwhale, logZ, beta, psi, indtwhale);
// //   mat lam = exp(loglam);
// //   Rprintf("set up initial values\n");
// //   
// //   
// //   // Start MCMC
// //   for(int s = 0; s < niter; s++) {
// //     
// //     // update alpha0
// //     for(int i = 0; i < n_al0; i ++){
// //       
// //       if( updateCOV ){
// //         if( (s+1 >= adaptInterval) && (s+1 - (adaptInterval * trunc((s+1) / adaptInterval)) == 0) ){
// //           dummyaccprob = accprob.col(i);
// //           rhat[i] = sum( dummyaccprob.rows(s+1-adaptInterval, s-1) ) / (adaptInterval-1);
// //           gamma1[i] = 1 / pow(adapIter[i], c1);
// //           gamma2[i] = c0 * gamma1[i];
// //           sigma2[i] = exp( log(sigma2[i]) + gamma2[i] * (rhat[i] - ropt) );
// //           
// //           
// //           dummyvec = postAlpha0.col(i);
// //           COValpha0[i] = COValpha0[i] + gamma1[i] * ( var( dummyvec.rows(s+1-adaptInterval, s-1) ) - COValpha0[i] );
// //           cholCOValpha0[i] = sqrt(  sigma2[i] * COValpha0[i] );
// //           
// //           adapIter[i] = adapIter[i] + 1;
// //         }
// //       }
// //       
// //       alpha0i = alpha0[i];
// //       newalpha0i = alpha0i + cholCOValpha0[i] * randn();
// //       newalpha0 = alpha0;
// //       newalpha0[i] = newalpha0i;
// //       
// //       newlogZ = cppf_logZ_a0t_etat(Xzoop, newalpha0, alpha, eta);
// //       newlogZobl = cppf_logZobl(newlogZ, indobl);
// //       newlogZsur = cppf_logZsur(newlogZ, indsur);
// //       newloglam = cppf_loglam_b0t(Xwhale, newlogZ, beta, indtwhale);
// //       newlam = exp(newloglam);
// //       
// //       logprob = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) -
// //         cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //         cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) -
// //         cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
// //         cppf_loglik_obl(logYobl, newlogZobl, sig2_obl) -
// //         cppf_loglik_obl(logYobl, logZobl, sig2_obl) +
// //         cppf_loglik_sur(logYsur, newlogZsur, lam0, lam1, sig2_sur) -
// //         cppf_loglik_sur(logYsur, logZsur, lam0, lam1, sig2_sur) +
// //         normal_logh(newalpha0i, alpha0tilde, tau2) -
// //         normal_logh(alpha0i, alpha0tilde, tau2);
// //       
// //       if (log(randu()) < logprob) {
// //         
// //         alpha0 = newalpha0;
// //         logZ = newlogZ;
// //         logZobl = newlogZobl;
// //         logZsur = newlogZsur;
// //         loglam = newloglam;
// //         lam = newlam;
// //         
// //         accprob(s, i) = 1;
// //       }
// //     }
// //     postAlpha0.row(s) = trans(alpha0);
// //     if(s == 0){ Rprintf("sampled alpha0\n"); }
// //     
// //     
// //     // Gibbs update for alpha0tilde
// //     alpha0tilde = cppf_gibbs_alpha0tilde(alpha0, tau2);
// //     postAlpha0tilde[s] = alpha0tilde;
// //     if(s == 0){ Rprintf("sampled alpha0tilde\n"); }
// //     
// //     
// //     // update alpha
// //     if( updateCOV ){
// //       if( (s+1 >= adaptInterval) && (s+1 - (adaptInterval * trunc((s+1) / adaptInterval)) == 0) ){
// //         dummyaccprob = accprob.col(n_al0);
// //         rhat[n_al0] = sum( dummyaccprob.rows(s+1-adaptInterval, s-1) ) / (adaptInterval-1);
// //         gamma1[n_al0] = 1 / pow(adapIter[n_al0], c1);
// //         gamma2[n_al0] = c0 * gamma1[n_al0];
// //         sigma2[n_al0] = exp( log(sigma2[n_al0]) + gamma2[n_al0] * (rhat[n_al0] - ropt) );
// //         
// //         COValpha = COValpha + gamma1[n_al0] * ( cov( postAlpha.rows(s+1-adaptInterval, s-1) ) - COValpha );
// //         cholCOValpha = trans( chol( sigma2[n_al0] * ( COValpha + 0.0000000001 * diagmat(ones(n_al)) ) ) );
// //         
// //         adapIter[n_al0] = adapIter[n_al0] + 1;
// //       }
// //     }
// //     
// //     newalpha = alpha + cholCOValpha * randn(n_al);
// //     
// //     newlogZ = cppf_logZ_a0t_etat(Xzoop, alpha0, newalpha, eta);
// //     newlogZobl = cppf_logZobl(newlogZ, indobl);
// //     newlogZsur = cppf_logZsur(newlogZ, indsur);
// //     newloglam = cppf_loglam_b0t_psi(Xwhale, newlogZ, beta, psi, indtwhale);
// //     newlam = exp(newloglam);
// //     
// //     logprob = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) -
// //       cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //       cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) -
// //       cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
// //       cppf_loglik_obl(logYobl, newlogZobl, sig2_obl) -
// //       cppf_loglik_obl(logYobl, logZobl, sig2_obl) +
// //       cppf_loglik_sur(logYsur, newlogZsur, lam0, lam1, sig2_sur) -
// //       cppf_loglik_sur(logYsur, logZsur, lam0, lam1, sig2_sur) +
// //       MVN_logh(newalpha, zeros(n_al), diagmat(ones(n_al))/100) -
// //       MVN_logh(alpha, zeros(n_al), diagmat(ones(n_al))/100);
// //     
// //     if (log(randu()) < logprob) {
// //       
// //       alpha = newalpha;
// //       logZ = newlogZ;
// //       logZobl = newlogZobl;
// //       logZsur = newlogZsur;
// //       loglam = newloglam;
// //       lam = newlam;
// //       
// //       accprob(s, n_al0) = 1;
// //     }
// //     
// //     postAlpha.row(s) = trans(alpha);
// //     if(s == 0){ Rprintf("sampled alpha\n"); }
// //     
// //     
// //     // update eta via Elliptical slice sampling
// //     for(int t = 0; t < ntzoop; t ++){
// //       neweta = eta;
// //       
// //       llprev = cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //         cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
// //         cppf_loglik_obl(logYobl, logZobl, sig2_obl) +
// //         cppf_loglik_sur(logYsur, logZsur, lam0, lam1, sig2_sur);
// //       prioretat = sqrt(kappa_eta) * trans(cholReta) * randn(nG);
// //       
// //       thetamin = 0;
// //       thetamax = 2 * M_PI;
// //       theta = thetamin + randu() * (thetamax - thetamin);
// //       thetamin = theta - 2 * M_PI;
// //       thetamax = theta;
// //       
// //       llthred = llprev + log(randu());
// //       accept = false;
// //       count = 0;
// //       
// //       while(accept == false){
// //         count = count + 1;
// //         
// //         neweta.col(t) = eta.col(t) * cos(theta) + prioretat * sin(theta);
// //         neweta.col(t) = neweta.col(t) - mean(neweta.col(t));
// //         
// //         newlogZ = cppf_logZ_a0t_etat(Xzoop, alpha0, alpha, neweta);
// //         newlogZobl = cppf_logZobl(newlogZ, indobl);
// //         newlogZsur = cppf_logZsur(newlogZ, indsur);
// //         newloglam = cppf_loglam_b0t_psi(Xwhale, newlogZ, beta, psi, indtwhale);
// //         newlam = exp(newloglam);
// //         
// //         llnew = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //           cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) +
// //           cppf_loglik_obl(logYobl, newlogZobl, sig2_obl) +
// //           cppf_loglik_sur(logYsur, newlogZsur, lam0, lam1, sig2_sur);
// //         
// //         if(llnew > llthred){
// //           llprev = llnew;
// //           accept = true;
// //         } else {
// //           if(theta < 0){ thetamin = theta; } else { thetamax = theta; }
// //           theta = thetamin + randu() * (thetamax - thetamin);
// //           if ( (count) % 1000 == 0 ) {
// //             Rprintf("ESS for eta: %d iterations...\n", count);
// //           }
// //           if (count > 20000) stop("ESS is not converging...\n");
// //         }
// //       }
// //       eta = neweta;
// //       logZ = newlogZ;
// //       logZobl = newlogZobl;
// //       logZsur = newlogZsur;
// //       loglam = newloglam;
// //       lam = newlam;
// //       
// //       postEtat = as<mat>(postEta[t]);
// //       postEtat.row(s) = trans(eta.col(t));
// //       postEta[t] = postEtat;
// //     }
// //     if(s == 0){ Rprintf("sampled eta\n"); }
// //     
// //     
// //     // update kappa_eta
// //     kappa_eta = cppf_gibbs_kappa_eta(eta, invReta, 2, 2);
// //     postKappa(s,0) = kappa_eta;
// //     if(s == 0){ Rprintf("sampled kappa_eta\n"); }
// //     
// //     
// //     // update lam0
// //     lam0 = cppf_gibbs_lam0(logYsur, logZsur, lam1, sig2_sur, 100);
// //     postLam(s, 0) = lam0;
// //     if(s == 0){ Rprintf("sampled lam0\n"); }
// //     
// //     
// //     // update lam1
// //     lam1 = cppf_gibbs_lam1(logYsur, logZsur, lam0, sig2_sur, 100);
// //     postLam(s, 1) = lam1;
// //     if(s == 0){ Rprintf("sampled lam1\n"); }
// //     
// //     
// //     // update sig2_sur
// //     sig2_sur = cppf_gibbs_sig2_sur(logYsur, logZsur, lam0, lam1, 2, 2);
// //     postSig2(s, 0) = sig2_sur;
// //     if(s == 0){ Rprintf("sampled sig2_sur\n"); }
// //     
// //     
// //     // update beta
// //     if( updateCOV ){
// //       if( (s+1 >= adaptInterval) && (s+1 - (adaptInterval * trunc((s+1) / adaptInterval)) == 0) ){
// //         dummyaccprob = accprob.col(n_al0+1);
// //         rhat[n_al0+1] = sum( dummyaccprob.rows(s+1-adaptInterval, s-1) ) / (adaptInterval-1);
// //         gamma1[n_al0+1] = 1 / pow(adapIter[n_al0+1], c1);
// //         gamma2[n_al0+1] = c0 * gamma1[n_al0+1];
// //         sigma2[n_al0+1] = exp( log(sigma2[n_al0+1]) + gamma2[n_al0+1] * (rhat[n_al0+1] - ropt) );
// //         
// //         COVbeta = COVbeta + gamma1[n_al0+1] * ( cov( postBeta.rows(s+1-adaptInterval, s-1)) - COVbeta );
// //         cholCOVbeta = trans( chol( sigma2[n_al0+1] * ( COVbeta + 0.0000000001 * diagmat(ones(n_be)) ) ) );
// //         
// //         adapIter[n_al0+1] = adapIter[n_al0+1] + 1;
// //       }
// //     }
// //     
// //     newbeta = beta + cholCOVbeta * randn(n_be);
// //     
// //     newloglam = cppf_loglam_b0t_psi(Xwhale, logZ, newbeta, psi, indtwhale);
// //     newlam = exp(newloglam);
// //     
// //     logprob = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) -
// //       cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //       cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) -
// //       cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
// //       MVN_logh(newbeta, zeros(n_be), diagmat(ones(n_be))/100) -
// //       MVN_logh(beta, zeros(n_be), diagmat(ones(n_be))/100);
// //     
// //     if (log(randu()) < logprob) {
// //       
// //       beta = newbeta;
// //       loglam = newloglam;
// //       lam = newlam;
// //       
// //       accprob(s, n_al0+1) = 1;
// //     }
// //     
// //     postBeta.row(s) = trans(beta);
// //     if(s == 0){ Rprintf("sampled beta\n"); }
// //     
// //     
// //     // update psi via elliptical slice sampling
// //     llprev = cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //       cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea);
// //     priorpsi = sqrt(kappa_psi) * trans(cholRpsi) * randn(nG);
// //     
// //     thetamin = 0;
// //     thetamax = 2 * M_PI;
// //     theta = thetamin + randu() * (thetamax - thetamin);
// //     thetamin = theta - 2 * M_PI;
// //     thetamax = theta;
// //     
// //     llthred = llprev + log(randu());
// //     accept = false;
// //     count = 0;
// //     
// //     while(accept == false){
// //       count = count + 1;
// //       
// //       newpsi = psi * cos(theta) + priorpsi * sin(theta);
// //       newpsi = newpsi - mean(newpsi);
// //       
// //       newloglam = cppf_loglam_b0t_psi(Xwhale, logZ, beta, newpsi, indtwhale);
// //       newlam = exp(newloglam);
// //       
// //       llnew = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //         cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea);
// //       
// //       if(llnew > llthred){
// //         llprev = llnew;
// //         accept = true;
// //       } else {
// //         if(theta < 0){ thetamin = theta; } else { thetamax = theta; }
// //         theta = thetamin + randu() * (thetamax - thetamin);
// //         if ( (count) % 1000 == 0 ) {
// //           Rprintf("ESS for psi: %d iterations...\n", count);
// //         }
// //         if (count > 20000) stop("ESS is not converging...\n");
// //       }
// //     }
// //     psi = newpsi;
// //     loglam = newloglam;
// //     lam = newlam;
// //     
// //     postPsi.row(s) = trans(psi);
// //     if(s == 0){ Rprintf("sampled psi\n"); }
// //     
// //     
// //     // update kappa_psi
// //     kappa_psi = cppf_gibbs_kappa_psi(psi, invRpsi, 2, 2);
// //     postKappa(s,1) = kappa_psi;
// //     if(s == 0){ Rprintf("sampled kappa_psi\n"); }
// //     
// //     
// //     // update c
// //     c = cppf_gibbs_ct(Wpam, lam, dur_dist, p_pam, cellarea, shape_c, scale_c);
// //     postC.row(s) = trans(c);
// //     if(s == 0){ Rprintf("sampled c\n"); }
// //     
// //     
// //     if ( (s+1) % 10 == 0 ) {
// //       Rprintf("Generated %d samples...\n", s+1);
// //     }
// //   }
// //   
// //   
// //   return Rcpp::List::create(Rcpp::Named("postAlpha0") = postAlpha0,
// //                             Rcpp::Named("postAlpha0tilde") = postAlpha0tilde,
// //                             Rcpp::Named("postAlpha") = postAlpha,
// //                             Rcpp::Named("postEta") = postEta,
// //                             Rcpp::Named("postLam") = postLam,
// //                             Rcpp::Named("postKappa") = postKappa,
// //                             Rcpp::Named("postSig2") = postSig2,
// //                             Rcpp::Named("postBeta") = postBeta,
// //                             Rcpp::Named("postPsi") = postPsi,
// //                             Rcpp::Named("postC") = postC,
// //                             Rcpp::Named("Accprob") = accprob,
// //                             Rcpp::Named("sigma2") = sigma2,
// //                             Rcpp::Named("adapIter") = adapIter,
// //                             Rcpp::Named("COValpha0") = COValpha0,
// //                             Rcpp::Named("COValpha") = COValpha,
// //                             Rcpp::Named("COVbeta") = COVbeta);
// //   
// // }
// // 
// // 
// // 
// // // sig2_obl is known
// // // time-varying intercept for whale component
// // // annual GP for whale component
// // // [[Rcpp::export]]
// // List cppf_fitJoint_b0t_psit(
// //     int niter, mat distmat,
// //     List Xzoop, vec logYobl, vec logYsur,
// //     List Xwhale, mat Sdist, mat Wpam, vec dur_dist,
// //     vec f_dist_obs, List f_dist, List p_pam,
// //     mat indobl, mat indsur, vec indtwhale, mat inddist, vec indmdist, vec indmwhale,
// //     double alpha0tilde, vec alpha0, vec alpha, mat eta, double lam0, double lam1,
// //     double tau2, double kappa_eta, double phi_eta, double sig2_obl, double sig2_sur,
// //     vec beta, mat psi, vec pii, vec c, double shape_c, double scale_c,
// //     double kappa_psi, double phi_psi,
// //     double cellarea,
// //     vec sigma2, vec COValpha0, mat COValpha, mat COVbeta,
// //     bool updateCOV, int adaptInterval, double adaptFactorExponent, vec adapIter){
// //   
// //   double negativeInf = -std::numeric_limits<float>::infinity();;
// //   
// //   indobl = indobl-1, indsur = indsur-1, indtwhale = indtwhale-1;
// //   inddist = inddist-1, indmdist = indmdist-1, indmwhale = indmwhale-1;
// //   
// //   int nG = distmat.n_rows, n_al0 = alpha0.size(), n_al = alpha.size(), n_be = beta.size(), ntzoop = eta.n_cols, ntwhale = psi.n_cols, n_sigma2 = sigma2.size();
// //   
// //   mat postAlpha0tilde = zeros(niter, 1), postAlpha0 = zeros(niter, n_al0), postAlpha = zeros(niter, n_al), postLam = zeros(niter, 2);
// //   mat postSig2 = zeros(niter, 1), postKappa = zeros(niter, 2);
// //   mat postBeta = zeros(niter, n_be);
// //   mat postC = zeros(niter, c.size()), postEtat = zeros(niter, nG), postPsit = zeros(niter, nG);
// //   List postEta(ntzoop), postPsi(ntwhale);
// //   for(int i = 0; i < ntzoop; i ++){ postEta[i] = zeros(niter, nG); }
// //   for(int i = 0; i < ntwhale; i ++){ postPsi[i] = zeros(niter, nG); }
// //   Rprintf("set up matrices for posterior samples\n");
// //   
// //   mat accprob = zeros(niter, n_sigma2);
// //   vec rhat = zeros(n_sigma2), gamma1 = zeros(n_sigma2), gamma2 = zeros(n_sigma2), dummyaccprob;
// //   double c0 = 1, c1 = adaptFactorExponent, ropt = 0.234, logprob;
// //   vec cholCOValpha0(n_al0);
// //   for(int i = 0; i < n_al0; i ++){ cholCOValpha0[i] = sqrt(sigma2[i] * COValpha0[i]); }
// //   mat cholCOValpha = trans( chol( sigma2[n_al0] * ( COValpha + 0.0000000001 * diagmat(ones(n_al)) ) ) );
// //   mat cholCOVbeta = trans( chol( sigma2[n_al0+1] * ( COVbeta + 0.0000000001 * diagmat(ones(n_be)) ) ) );
// //   Rprintf("set up adaptive proposals\n");
// //   
// //   double alpha0i, newalpha0i, beta0i, newbeta0i;
// //   vec newalpha0, newalpha, newbeta0, newbeta, newlogZobl, newlogZsur;
// //   mat newlogZ, newloglam, newlam;
// //   vec dummyvec;
// //   
// //   int count;
// //   vec xvals;
// //   double llprev, llthred, llnew, thetamin, thetamax, theta;
// //   vec prioretat, priorpsit;
// //   mat neweta, newpsi;
// //   bool accept;
// //   
// //   mat Reta = exp(-distmat / phi_eta), Rpsi = exp(-distmat / phi_psi);
// //   mat invReta = inv(Reta), invRpsi= inv(Rpsi), cholReta = chol(Reta), cholRpsi= chol(Rpsi);
// //   rowvec oneinvReta = trans(ones(nG)) * invReta, oneinvRpsi = trans(ones(nG)) * invRpsi;
// //   double oneinvRetaone = sum( oneinvReta * ones(nG) ), oneinvRpsione = sum( oneinvRpsi * ones(nG) );
// //   Rprintf("set up GP covariance matrices\n");
// //   
// //   mat logZ = cppf_logZ_a0t_etat(Xzoop, alpha0, alpha, eta);
// //   vec logZobl = cppf_logZobl(logZ, indobl);
// //   vec logZsur = cppf_logZsur(logZ, indsur);
// //   mat loglam = cppf_loglam_b0t_psit(Xwhale, logZ, beta, psi, indtwhale);
// //   mat lam = exp(loglam);
// //   Rprintf("set up initial values\n");
// //   
// //   
// //   // Start MCMC
// //   for(int s = 0; s < niter; s++) {
// //     
// //     // update alpha0
// //     for(int i = 0; i < n_al0; i ++){
// //       
// //       if( updateCOV ){
// //         if( (s+1 >= adaptInterval) && (s+1 - (adaptInterval * trunc((s+1) / adaptInterval)) == 0) ){
// //           dummyaccprob = accprob.col(i);
// //           rhat[i] = sum( dummyaccprob.rows(s+1-adaptInterval, s-1) ) / (adaptInterval-1);
// //           gamma1[i] = 1 / pow(adapIter[i], c1);
// //           gamma2[i] = c0 * gamma1[i];
// //           sigma2[i] = exp( log(sigma2[i]) + gamma2[i] * (rhat[i] - ropt) );
// //           
// //           
// //           dummyvec = postAlpha0.col(i);
// //           COValpha0[i] = COValpha0[i] + gamma1[i] * ( var( dummyvec.rows(s+1-adaptInterval, s-1) ) - COValpha0[i] );
// //           cholCOValpha0[i] = sqrt(  sigma2[i] * COValpha0[i] );
// //           
// //           adapIter[i] = adapIter[i] + 1;
// //         }
// //       }
// //       
// //       alpha0i = alpha0[i];
// //       newalpha0i = alpha0i + cholCOValpha0[i] * randn();
// //       newalpha0 = alpha0;
// //       newalpha0[i] = newalpha0i;
// //       
// //       newlogZ = cppf_logZ_a0t_etat(Xzoop, newalpha0, alpha, eta);
// //       newlogZobl = cppf_logZobl(newlogZ, indobl);
// //       newlogZsur = cppf_logZsur(newlogZ, indsur);
// //       newloglam = cppf_loglam_b0t(Xwhale, newlogZ, beta, indtwhale);
// //       newlam = exp(newloglam);
// //       
// //       logprob = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) -
// //         cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //         cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) -
// //         cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
// //         cppf_loglik_obl(logYobl, newlogZobl, sig2_obl) -
// //         cppf_loglik_obl(logYobl, logZobl, sig2_obl) +
// //         cppf_loglik_sur(logYsur, newlogZsur, lam0, lam1, sig2_sur) -
// //         cppf_loglik_sur(logYsur, logZsur, lam0, lam1, sig2_sur) +
// //         normal_logh(newalpha0i, alpha0tilde, tau2) -
// //         normal_logh(alpha0i, alpha0tilde, tau2);
// //       
// //       if (log(randu()) < logprob) {
// //         
// //         alpha0 = newalpha0;
// //         logZ = newlogZ;
// //         logZobl = newlogZobl;
// //         logZsur = newlogZsur;
// //         loglam = newloglam;
// //         lam = newlam;
// //         
// //         accprob(s, i) = 1;
// //       }
// //     }
// //     postAlpha0.row(s) = trans(alpha0);
// //     if(s == 0){ Rprintf("sampled alpha0\n"); }
// //     
// //     
// //     // Gibbs update for alpha0tilde
// //     alpha0tilde = cppf_gibbs_alpha0tilde(alpha0, tau2);
// //     postAlpha0tilde[s] = alpha0tilde;
// //     if(s == 0){ Rprintf("sampled alpha0tilde\n"); }
// //     
// //     
// //     // update alpha
// //     if( updateCOV ){
// //       if( (s+1 >= adaptInterval) && (s+1 - (adaptInterval * trunc((s+1) / adaptInterval)) == 0) ){
// //         dummyaccprob = accprob.col(n_al0);
// //         rhat[n_al0] = sum( dummyaccprob.rows(s+1-adaptInterval, s-1) ) / (adaptInterval-1);
// //         gamma1[n_al0] = 1 / pow(adapIter[n_al0], c1);
// //         gamma2[n_al0] = c0 * gamma1[n_al0];
// //         sigma2[n_al0] = exp( log(sigma2[n_al0]) + gamma2[n_al0] * (rhat[n_al0] - ropt) );
// //         
// //         COValpha = COValpha + gamma1[n_al0] * ( cov( postAlpha.rows(s+1-adaptInterval, s-1) ) - COValpha );
// //         cholCOValpha = trans( chol( sigma2[n_al0] * ( COValpha + 0.0000000001 * diagmat(ones(n_al)) ) ) );
// //         
// //         adapIter[n_al0] = adapIter[n_al0] + 1;
// //       }
// //     }
// //     
// //     newalpha = alpha + cholCOValpha * randn(n_al);
// //     
// //     newlogZ = cppf_logZ_a0t_etat(Xzoop, alpha0, newalpha, eta);
// //     newlogZobl = cppf_logZobl(newlogZ, indobl);
// //     newlogZsur = cppf_logZsur(newlogZ, indsur);
// //     newloglam = cppf_loglam_b0t_psit(Xwhale, newlogZ, beta, psi, indtwhale);
// //     newlam = exp(newloglam);
// //     
// //     logprob = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) -
// //       cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //       cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) -
// //       cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
// //       cppf_loglik_obl(logYobl, newlogZobl, sig2_obl) -
// //       cppf_loglik_obl(logYobl, logZobl, sig2_obl) +
// //       cppf_loglik_sur(logYsur, newlogZsur, lam0, lam1, sig2_sur) -
// //       cppf_loglik_sur(logYsur, logZsur, lam0, lam1, sig2_sur) +
// //       MVN_logh(newalpha, zeros(n_al), diagmat(ones(n_al))/100) -
// //       MVN_logh(alpha, zeros(n_al), diagmat(ones(n_al))/100);
// //     
// //     if (log(randu()) < logprob) {
// //       
// //       alpha = newalpha;
// //       logZ = newlogZ;
// //       logZobl = newlogZobl;
// //       logZsur = newlogZsur;
// //       loglam = newloglam;
// //       lam = newlam;
// //       
// //       accprob(s, n_al0) = 1;
// //     }
// //     
// //     postAlpha.row(s) = trans(alpha);
// //     if(s == 0){ Rprintf("sampled alpha\n"); }
// //     
// //     
// //     // update eta via Elliptical slice sampling
// //     for(int t = 0; t < ntzoop; t ++){
// //       neweta = eta;
// //       
// //       llprev = cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //         cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
// //         cppf_loglik_obl(logYobl, logZobl, sig2_obl) +
// //         cppf_loglik_sur(logYsur, logZsur, lam0, lam1, sig2_sur);
// //       prioretat = sqrt(kappa_eta) * trans(cholReta) * randn(nG);
// //       
// //       thetamin = 0;
// //       thetamax = 2 * M_PI;
// //       theta = thetamin + randu() * (thetamax - thetamin);
// //       thetamin = theta - 2 * M_PI;
// //       thetamax = theta;
// //       
// //       llthred = llprev + log(randu());
// //       accept = false;
// //       count = 0;
// //       
// //       while(accept == false){
// //         count = count + 1;
// //         
// //         neweta.col(t) = eta.col(t) * cos(theta) + prioretat * sin(theta);
// //         neweta.col(t) = neweta.col(t) - mean(neweta.col(t));
// //         
// //         newlogZ = cppf_logZ_a0t_etat(Xzoop, alpha0, alpha, neweta);
// //         newlogZobl = cppf_logZobl(newlogZ, indobl);
// //         newlogZsur = cppf_logZsur(newlogZ, indsur);
// //         newloglam = cppf_loglam_b0t_psit(Xwhale, newlogZ, beta, psi, indtwhale);
// //         newlam = exp(newloglam);
// //         
// //         llnew = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //           cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) +
// //           cppf_loglik_obl(logYobl, newlogZobl, sig2_obl) +
// //           cppf_loglik_sur(logYsur, newlogZsur, lam0, lam1, sig2_sur);
// //         
// //         if(llnew > llthred){
// //           llprev = llnew;
// //           accept = true;
// //         } else {
// //           if(theta < 0){ thetamin = theta; } else { thetamax = theta; }
// //           theta = thetamin + randu() * (thetamax - thetamin);
// //           if ( (count) % 1000 == 0 ) {
// //             Rprintf("ESS for eta: %d iterations...\n", count);
// //           }
// //           if (count > 20000) stop("ESS is not converging...\n");
// //         }
// //       }
// //       eta = neweta;
// //       logZ = newlogZ;
// //       logZobl = newlogZobl;
// //       logZsur = newlogZsur;
// //       loglam = newloglam;
// //       lam = newlam;
// //       
// //       postEtat = as<mat>(postEta[t]);
// //       postEtat.row(s) = trans(eta.col(t));
// //       postEta[t] = postEtat;
// //     }
// //     if(s == 0){ Rprintf("sampled eta\n"); }
// //     
// //     
// //     // update kappa_eta
// //     kappa_eta = cppf_gibbs_kappa_eta(eta, invReta, 2, 2);
// //     postKappa(s,0) = kappa_eta;
// //     if(s == 0){ Rprintf("sampled kappa_eta\n"); }
// //     
// //     
// //     // update lam0
// //     lam0 = cppf_gibbs_lam0(logYsur, logZsur, lam1, sig2_sur, 100);
// //     postLam(s, 0) = lam0;
// //     if(s == 0){ Rprintf("sampled lam0\n"); }
// //     
// //     
// //     // update lam1
// //     lam1 = cppf_gibbs_lam1(logYsur, logZsur, lam0, sig2_sur, 100);
// //     postLam(s, 1) = lam1;
// //     if(s == 0){ Rprintf("sampled lam1\n"); }
// //     
// //     
// //     // update sig2_sur
// //     sig2_sur = cppf_gibbs_sig2_sur(logYsur, logZsur, lam0, lam1, 2, 2);
// //     postSig2(s, 0) = sig2_sur;
// //     if(s == 0){ Rprintf("sampled sig2_sur\n"); }
// //     
// //     
// //     // update beta
// //     if( updateCOV ){
// //       if( (s+1 >= adaptInterval) && (s+1 - (adaptInterval * trunc((s+1) / adaptInterval)) == 0) ){
// //         dummyaccprob = accprob.col(n_al0+1);
// //         rhat[n_al0+1] = sum( dummyaccprob.rows(s+1-adaptInterval, s-1) ) / (adaptInterval-1);
// //         gamma1[n_al0+1] = 1 / pow(adapIter[n_al0+1], c1);
// //         gamma2[n_al0+1] = c0 * gamma1[n_al0+1];
// //         sigma2[n_al0+1] = exp( log(sigma2[n_al0+1]) + gamma2[n_al0+1] * (rhat[n_al0+1] - ropt) );
// //         
// //         COVbeta = COVbeta + gamma1[n_al0+1] * ( cov( postBeta.rows(s+1-adaptInterval, s-1)) - COVbeta );
// //         cholCOVbeta = trans( chol( sigma2[n_al0+1] * ( COVbeta + 0.0000000001 * diagmat(ones(n_be)) ) ) );
// //         
// //         adapIter[n_al0+1] = adapIter[n_al0+1] + 1;
// //       }
// //     }
// //     
// //     newbeta = beta + cholCOVbeta * randn(n_be);
// //     
// //     newloglam = cppf_loglam_b0t_psit(Xwhale, logZ, newbeta, psi, indtwhale);
// //     newlam = exp(newloglam);
// //     
// //     logprob = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) -
// //       cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //       cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) -
// //       cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
// //       MVN_logh(newbeta, zeros(n_be), diagmat(ones(n_be))/100) -
// //       MVN_logh(beta, zeros(n_be), diagmat(ones(n_be))/100);
// //     
// //     if (log(randu()) < logprob) {
// //       
// //       beta = newbeta;
// //       loglam = newloglam;
// //       lam = newlam;
// //       
// //       accprob(s, n_al0+1) = 1;
// //     }
// //     
// //     postBeta.row(s) = trans(beta);
// //     if(s == 0){ Rprintf("sampled beta\n"); }
// //     
// //     
// //     
// //     // update psi via Elliptical slice sampling
// //     for(int t = 0; t < ntwhale; t ++){
// //       newpsi = psi;
// //       
// //       llprev = cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //         cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea);
// //       priorpsit = sqrt(kappa_psi) * trans(cholRpsi) * randn(nG);
// //       
// //       thetamin = 0;
// //       thetamax = 2 * M_PI;
// //       theta = thetamin + randu() * (thetamax - thetamin);
// //       thetamin = theta - 2 * M_PI;
// //       thetamax = theta;
// //       
// //       llthred = llprev + log(randu());
// //       accept = false;
// //       count = 0;
// //       
// //       while(accept == false){
// //         count = count + 1;
// //         
// //         newpsi.col(t) = psi.col(t) * cos(theta) + priorpsit * sin(theta);
// //         newpsi.col(t) = newpsi.col(t) - mean(newpsi.col(t));
// //         
// //         newloglam = cppf_loglam_b0t_psit(Xwhale, logZ, beta, newpsi, indtwhale);
// //         newlam = exp(newloglam);
// //         
// //         llnew = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //           cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea);
// //         
// //         if(llnew > llthred){
// //           llprev = llnew;
// //           accept = true;
// //         } else {
// //           if(theta < 0){ thetamin = theta; } else { thetamax = theta; }
// //           theta = thetamin + randu() * (thetamax - thetamin);
// //           if ( (count) % 1000 == 0 ) {
// //             Rprintf("ESS for psi: %d iterations...\n", count);
// //           }
// //           if (count > 20000) stop("ESS is not converging...\n");
// //         }
// //       }
// //       psi = newpsi;
// //       loglam = newloglam;
// //       lam = newlam;
// //       
// //       postPsit = as<mat>(postPsi[t]);
// //       postPsit.row(s) = trans(psi.col(t));
// //       postPsi[t] = postPsit;
// //     }
// //     if(s == 0){ Rprintf("sampled psi\n"); }
// //     
// //     
// //     // update kappa_psi
// //     kappa_psi = cppf_gibbs_kappa_psi(psi, invRpsi, 2, 2);
// //     postKappa(s,1) = kappa_psi;
// //     if(s == 0){ Rprintf("sampled kappa_psi\n"); }
// //     
// //     
// //     // update c
// //     c = cppf_gibbs_ct(Wpam, lam, dur_dist, p_pam, cellarea, shape_c, scale_c);
// //     postC.row(s) = trans(c);
// //     if(s == 0){ Rprintf("sampled c\n"); }
// //     
// //     
// //     if ( (s+1) % 10 == 0 ) {
// //       Rprintf("Generated %d samples...\n", s+1);
// //     }
// //   }
// //   
// //   
// //   return Rcpp::List::create(Rcpp::Named("postAlpha0") = postAlpha0,
// //                             Rcpp::Named("postAlpha0tilde") = postAlpha0tilde,
// //                             Rcpp::Named("postAlpha") = postAlpha,
// //                             Rcpp::Named("postEta") = postEta,
// //                             Rcpp::Named("postLam") = postLam,
// //                             Rcpp::Named("postKappa") = postKappa,
// //                             Rcpp::Named("postSig2") = postSig2,
// //                             Rcpp::Named("postBeta") = postBeta,
// //                             Rcpp::Named("postPsi") = postPsi,
// //                             Rcpp::Named("postC") = postC,
// //                             Rcpp::Named("Accprob") = accprob,
// //                             Rcpp::Named("sigma2") = sigma2,
// //                             Rcpp::Named("adapIter") = adapIter,
// //                             Rcpp::Named("COValpha0") = COValpha0,
// //                             Rcpp::Named("COValpha") = COValpha,
// //                             Rcpp::Named("COVbeta") = COVbeta);
// //   
// // }
// // 
// // 
// // 
// // // sig2_obl is known
// // // time-varying intercept for whale component
// // // annual GP for whale component
// // // [[Rcpp::export]]
// // List cppf_fitJoint_b0t_psit_1day(
// //     int niter, mat distmat,
// //     List Xzoop, vec logYobl, vec logYsur,
// //     List Xwhale, mat Sdist, mat Wpam, vec dur_dist,
// //     vec f_dist_obs, List f_dist, List p_pam,
// //     mat indobl, mat indsur, vec indtwhale, mat inddist, vec indmdist, vec indmwhale,
// //     vec alpha, mat eta, double lam0, double lam1,
// //     double tau2, double kappa_eta, double phi_eta, double sig2_obl, double sig2_sur,
// //     vec beta, mat psi, vec pii, vec c, double shape_c, double scale_c,
// //     double kappa_psi, double phi_psi,
// //     double cellarea,
// //     vec sigma2, mat COValpha, mat COVbeta,
// //     bool updateCOV, int adaptInterval, double adaptFactorExponent, vec adapIter){
// //   
// //   double negativeInf = -std::numeric_limits<float>::infinity();;
// //   
// //   indobl = indobl-1, indsur = indsur-1, indtwhale = indtwhale-1;
// //   inddist = inddist-1, indmdist = indmdist-1, indmwhale = indmwhale-1;
// //   
// //   int nG = distmat.n_rows, n_al = alpha.size(), n_be = beta.size(), ntzoop = eta.n_cols, ntwhale = psi.n_cols, n_sigma2 = sigma2.size();
// //   
// //   mat postAlpha = zeros(niter, n_al), postLam = zeros(niter, 2);
// //   mat postSig2 = zeros(niter, 1), postKappa = zeros(niter, 2);
// //   mat postBeta = zeros(niter, n_be);
// //   mat postC = zeros(niter, c.size()), postEtat = zeros(niter, nG), postPsit = zeros(niter, nG);
// //   List postEta(ntzoop), postPsi(ntwhale);
// //   for(int i = 0; i < ntzoop; i ++){ postEta[i] = zeros(niter, nG); }
// //   for(int i = 0; i < ntwhale; i ++){ postPsi[i] = zeros(niter, nG); }
// //   Rprintf("set up matrices for posterior samples\n");
// //   
// //   mat accprob = zeros(niter, n_sigma2);
// //   vec rhat = zeros(n_sigma2), gamma1 = zeros(n_sigma2), gamma2 = zeros(n_sigma2), dummyaccprob;
// //   double c0 = 1, c1 = adaptFactorExponent, ropt = 0.234, logprob;
// //   mat cholCOValpha = trans( chol( sigma2[0] * ( COValpha + 0.0000000001 * diagmat(ones(n_al)) ) ) );
// //   mat cholCOVbeta = trans( chol( sigma2[1] * ( COVbeta + 0.0000000001 * diagmat(ones(n_be)) ) ) );
// //   Rprintf("set up adaptive proposals\n");
// //   
// //   double beta0i, newbeta0i;
// //   vec newalpha, newbeta0, newbeta, newlogZobl, newlogZsur;
// //   mat newlogZ, newloglam, newlam;
// //   vec dummyvec;
// //   
// //   int count;
// //   vec xvals;
// //   double llprev, llthred, llnew, thetamin, thetamax, theta;
// //   vec prioretat, priorpsit;
// //   mat neweta, newpsi;
// //   bool accept;
// //   
// //   mat Reta = exp(-distmat / phi_eta), Rpsi = exp(-distmat / phi_psi);
// //   mat invReta = inv(Reta), invRpsi= inv(Rpsi), cholReta = chol(Reta), cholRpsi= chol(Rpsi);
// //   rowvec oneinvReta = trans(ones(nG)) * invReta, oneinvRpsi = trans(ones(nG)) * invRpsi;
// //   double oneinvRetaone = sum( oneinvReta * ones(nG) ), oneinvRpsione = sum( oneinvRpsi * ones(nG) );
// //   Rprintf("set up GP covariance matrices\n");
// //   
// //   mat logZ = cppf_logZ_a0t_etat_ja(Xzoop, alpha, eta);
// //   vec logZobl = cppf_logZobl(logZ, indobl);
// //   vec logZsur = cppf_logZsur(logZ, indsur);
// //   mat loglam = cppf_loglam_b0t_psit(Xwhale, logZ, beta, psi, indtwhale);
// //   mat lam = exp(loglam);
// //   Rprintf("set up initial values\n");
// //   
// //   
// //   // Start MCMC
// //   for(int s = 0; s < niter; s++) {
// //     
// //     // update alpha
// //     if( updateCOV ){
// //       if( (s+1 >= adaptInterval) && (s+1 - (adaptInterval * trunc((s+1) / adaptInterval)) == 0) ){
// //         dummyaccprob = accprob.col(0);
// //         rhat[0] = sum( dummyaccprob.rows(s+1-adaptInterval, s-1) ) / (adaptInterval-1);
// //         gamma1[0] = 1 / pow(adapIter[0], c1);
// //         gamma2[0] = c0 * gamma1[0];
// //         sigma2[0] = exp( log(sigma2[0]) + gamma2[0] * (rhat[0] - ropt) );
// //         
// //         COValpha = COValpha + gamma1[0] * ( cov( postAlpha.rows(s+1-adaptInterval, s-1) ) - COValpha );
// //         cholCOValpha = trans( chol( sigma2[0] * ( COValpha + 0.0000000001 * diagmat(ones(n_al)) ) ) );
// //         
// //         adapIter[0] = adapIter[0] + 1;
// //       }
// //     }
// //     
// //     newalpha = alpha + cholCOValpha * randn(n_al);
// //     
// //     newlogZ = cppf_logZ_a0t_etat_ja(Xzoop, newalpha, eta);
// //     newlogZobl = cppf_logZobl(newlogZ, indobl);
// //     newlogZsur = cppf_logZsur(newlogZ, indsur);
// //     newloglam = cppf_loglam_b0t_psit(Xwhale, newlogZ, beta, psi, indtwhale);
// //     newlam = exp(newloglam);
// //     
// //     logprob = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) -
// //       cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //       cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) -
// //       cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
// //       cppf_loglik_obl(logYobl, newlogZobl, sig2_obl) -
// //       cppf_loglik_obl(logYobl, logZobl, sig2_obl) +
// //       cppf_loglik_sur(logYsur, newlogZsur, lam0, lam1, sig2_sur) -
// //       cppf_loglik_sur(logYsur, logZsur, lam0, lam1, sig2_sur) +
// //       MVN_logh(newalpha, zeros(n_al), diagmat(ones(n_al))/100) -
// //       MVN_logh(alpha, zeros(n_al), diagmat(ones(n_al))/100);
// //     
// //     if (log(randu()) < logprob) {
// //       
// //       alpha = newalpha;
// //       logZ = newlogZ;
// //       logZobl = newlogZobl;
// //       logZsur = newlogZsur;
// //       loglam = newloglam;
// //       lam = newlam;
// //       
// //       accprob(s, 0) = 1;
// //     }
// //     
// //     postAlpha.row(s) = trans(alpha);
// //     if(s == 0){ Rprintf("sampled alpha\n"); }
// //     
// //     
// //     // update eta via Elliptical slice sampling
// //     for(int t = 0; t < ntzoop; t ++){
// //       neweta = eta;
// //       
// //       llprev = cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //         cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
// //         cppf_loglik_obl(logYobl, logZobl, sig2_obl) +
// //         cppf_loglik_sur(logYsur, logZsur, lam0, lam1, sig2_sur);
// //       prioretat = sqrt(kappa_eta) * trans(cholReta) * randn(nG);
// //       
// //       thetamin = 0;
// //       thetamax = 2 * M_PI;
// //       theta = thetamin + randu() * (thetamax - thetamin);
// //       thetamin = theta - 2 * M_PI;
// //       thetamax = theta;
// //       
// //       llthred = llprev + log(randu());
// //       accept = false;
// //       count = 0;
// //       
// //       while(accept == false){
// //         count = count + 1;
// //         
// //         neweta.col(t) = eta.col(t) * cos(theta) + prioretat * sin(theta);
// //         neweta.col(t) = neweta.col(t) - mean(neweta.col(t));
// //         
// //         newlogZ = cppf_logZ_a0t_etat_ja(Xzoop, alpha, neweta);
// //         newlogZobl = cppf_logZobl(newlogZ, indobl);
// //         newlogZsur = cppf_logZsur(newlogZ, indsur);
// //         newloglam = cppf_loglam_b0t_psit(Xwhale, newlogZ, beta, psi, indtwhale);
// //         newlam = exp(newloglam);
// //         
// //         llnew = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //           cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) +
// //           cppf_loglik_obl(logYobl, newlogZobl, sig2_obl) +
// //           cppf_loglik_sur(logYsur, newlogZsur, lam0, lam1, sig2_sur);
// //         
// //         if(llnew > llthred){
// //           llprev = llnew;
// //           accept = true;
// //         } else {
// //           if(theta < 0){ thetamin = theta; } else { thetamax = theta; }
// //           theta = thetamin + randu() * (thetamax - thetamin);
// //           if ( (count) % 1000 == 0 ) {
// //             Rprintf("ESS for eta: %d iterations...\n", count);
// //           }
// //           if (count > 20000) stop("ESS is not converging...\n");
// //         }
// //       }
// //       eta = neweta;
// //       logZ = newlogZ;
// //       logZobl = newlogZobl;
// //       logZsur = newlogZsur;
// //       loglam = newloglam;
// //       lam = newlam;
// //       
// //       postEtat = as<mat>(postEta[t]);
// //       postEtat.row(s) = trans(eta.col(t));
// //       postEta[t] = postEtat;
// //     }
// //     if(s == 0){ Rprintf("sampled eta\n"); }
// //     
// //     
// //     // update kappa_eta
// //     kappa_eta = cppf_gibbs_kappa_eta(eta, invReta, 2, 2);
// //     postKappa(s,0) = kappa_eta;
// //     if(s == 0){ Rprintf("sampled kappa_eta\n"); }
// //     
// //     
// //     // update lam0
// //     lam0 = cppf_gibbs_lam0(logYsur, logZsur, lam1, sig2_sur, 100);
// //     postLam(s, 0) = lam0;
// //     if(s == 0){ Rprintf("sampled lam0\n"); }
// //     
// //     
// //     // update lam1
// //     lam1 = cppf_gibbs_lam1(logYsur, logZsur, lam0, sig2_sur, 100);
// //     postLam(s, 1) = lam1;
// //     if(s == 0){ Rprintf("sampled lam1\n"); }
// //     
// //     
// //     // update sig2_sur
// //     sig2_sur = cppf_gibbs_sig2_sur(logYsur, logZsur, lam0, lam1, 2, 2);
// //     postSig2(s, 0) = sig2_sur;
// //     if(s == 0){ Rprintf("sampled sig2_sur\n"); }
// //     
// //     
// //     // update beta
// //     if( updateCOV ){
// //       if( (s+1 >= adaptInterval) && (s+1 - (adaptInterval * trunc((s+1) / adaptInterval)) == 0) ){
// //         dummyaccprob = accprob.col(1);
// //         rhat[1] = sum( dummyaccprob.rows(s+1-adaptInterval, s-1) ) / (adaptInterval-1);
// //         gamma1[1] = 1 / pow(adapIter[1], c1);
// //         gamma2[1] = c0 * gamma1[1];
// //         sigma2[1] = exp( log(sigma2[1]) + gamma2[1] * (rhat[1] - ropt) );
// //         
// //         COVbeta = COVbeta + gamma1[1] * ( cov( postBeta.rows(s+1-adaptInterval, s-1)) - COVbeta );
// //         cholCOVbeta = trans( chol( sigma2[1] * ( COVbeta + 0.0000000001 * diagmat(ones(n_be)) ) ) );
// //         
// //         adapIter[1] = adapIter[1] + 1;
// //       }
// //     }
// //     
// //     newbeta = beta + cholCOVbeta * randn(n_be);
// //     
// //     newloglam = cppf_loglam_b0t_psit(Xwhale, logZ, newbeta, psi, indtwhale);
// //     newlam = exp(newloglam);
// //     
// //     logprob = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) -
// //       cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //       cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) -
// //       cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
// //       MVN_logh(newbeta, zeros(n_be), diagmat(ones(n_be))/100) -
// //       MVN_logh(beta, zeros(n_be), diagmat(ones(n_be))/100);
// //     
// //     if (log(randu()) < logprob) {
// //       
// //       beta = newbeta;
// //       loglam = newloglam;
// //       lam = newlam;
// //       
// //       accprob(s, 1) = 1;
// //     }
// //     
// //     postBeta.row(s) = trans(beta);
// //     if(s == 0){ Rprintf("sampled beta\n"); }
// //     
// //     
// //     
// //     // update psi via Elliptical slice sampling
// //     for(int t = 0; t < ntwhale; t ++){
// //       newpsi = psi;
// //       
// //       llprev = cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //         cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea);
// //       priorpsit = sqrt(kappa_psi) * trans(cholRpsi) * randn(nG);
// //       
// //       thetamin = 0;
// //       thetamax = 2 * M_PI;
// //       theta = thetamin + randu() * (thetamax - thetamin);
// //       thetamin = theta - 2 * M_PI;
// //       thetamax = theta;
// //       
// //       llthred = llprev + log(randu());
// //       accept = false;
// //       count = 0;
// //       
// //       while(accept == false){
// //         count = count + 1;
// //         
// //         newpsi.col(t) = psi.col(t) * cos(theta) + priorpsit * sin(theta);
// //         newpsi.col(t) = newpsi.col(t) - mean(newpsi.col(t));
// //         
// //         newloglam = cppf_loglam_b0t_psit(Xwhale, logZ, beta, newpsi, indtwhale);
// //         newlam = exp(newloglam);
// //         
// //         llnew = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //           cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea);
// //         
// //         if(llnew > llthred){
// //           llprev = llnew;
// //           accept = true;
// //         } else {
// //           if(theta < 0){ thetamin = theta; } else { thetamax = theta; }
// //           theta = thetamin + randu() * (thetamax - thetamin);
// //           if ( (count) % 1000 == 0 ) {
// //             Rprintf("ESS for psi: %d iterations...\n", count);
// //           }
// //           if (count > 20000) stop("ESS is not converging...\n");
// //         }
// //       }
// //       psi = newpsi;
// //       loglam = newloglam;
// //       lam = newlam;
// //       
// //       postPsit = as<mat>(postPsi[t]);
// //       postPsit.row(s) = trans(psi.col(t));
// //       postPsi[t] = postPsit;
// //     }
// //     if(s == 0){ Rprintf("sampled psi\n"); }
// //     
// //     
// //     // update kappa_psi
// //     kappa_psi = cppf_gibbs_kappa_psi(psi, invRpsi, 2, 2);
// //     postKappa(s,1) = kappa_psi;
// //     if(s == 0){ Rprintf("sampled kappa_psi\n"); }
// //     
// //     
// //     // update c
// //     c = cppf_gibbs_ct(Wpam, lam, dur_dist, p_pam, cellarea, shape_c, scale_c);
// //     postC.row(s) = trans(c);
// //     if(s == 0){ Rprintf("sampled c\n"); }
// //     
// //     
// //     if ( (s+1) % 10 == 0 ) {
// //       Rprintf("Generated %d samples...\n", s+1);
// //     }
// //   }
// //   
// //   
// //   return Rcpp::List::create(Rcpp::Named("postAlpha") = postAlpha,
// //                             Rcpp::Named("postEta") = postEta,
// //                             Rcpp::Named("postLam") = postLam,
// //                             Rcpp::Named("postKappa") = postKappa,
// //                             Rcpp::Named("postSig2") = postSig2,
// //                             Rcpp::Named("postBeta") = postBeta,
// //                             Rcpp::Named("postPsi") = postPsi,
// //                             Rcpp::Named("postC") = postC,
// //                             Rcpp::Named("Accprob") = accprob,
// //                             Rcpp::Named("sigma2") = sigma2,
// //                             Rcpp::Named("adapIter") = adapIter,
// //                             Rcpp::Named("COValpha") = COValpha,
// //                             Rcpp::Named("COVbeta") = COVbeta);
// //   
// // }
// // 
// // 
// // 
// // 
// // 
// // 
// // // sig2_obl is known
// // // GP eta
// // // [[Rcpp::export]]
// // List cppf_fitZoop(
// //     int niter, mat distmat,
// //     List Xzoop, vec logYobl, vec logYsur,
// //     mat indobl, mat indsur,
// //     double alpha0tilde, vec alpha0, vec alpha, mat eta, double lam0, double lam1,
// //     double tau2, double kappa_eta, double phi_eta, double sig2_obl, double sig2_sur,
// //     vec sigma2, vec COValpha0, mat COValpha,
// //     bool updateCOV, int adaptInterval, double adaptFactorExponent, vec adapIter){
// //   
// //   double negativeInf = -std::numeric_limits<float>::infinity();;
// //   
// //   indobl = indobl-1, indsur = indsur-1;
// //   
// //   int nG = distmat.n_rows, n_al0 = alpha0.size(), n_al = alpha.size(), ntzoop = eta.n_cols, n_sigma2 = sigma2.size();
// //   
// //   mat postAlpha0tilde = zeros(niter, 1), postAlpha0 = zeros(niter, n_al0), postAlpha = zeros(niter, n_al), postLam = zeros(niter, 2);
// //   mat postSig2 = zeros(niter, 1), postKappa = zeros(niter, 1);
// //   mat postEtat = zeros(niter, nG);
// //   List postEta(ntzoop);
// //   for(int i = 0; i < ntzoop; i ++){
// //     postEta[i] = zeros(niter, nG);
// //   }
// //   Rprintf("set up matrices for posterior samples\n");
// //   
// //   mat accprob = zeros(niter, n_sigma2);
// //   vec rhat = zeros(n_sigma2), gamma1 = zeros(n_sigma2), gamma2 = zeros(n_sigma2), dummyaccprob;
// //   double c0 = 1, c1 = adaptFactorExponent, ropt = 0.234, logprob;
// //   vec cholCOValpha0(n_al0);
// //   for(int i = 0; i < n_al0; i ++){ cholCOValpha0[i] = sqrt(sigma2[i] * COValpha0[i]); }
// //   mat cholCOValpha = trans( chol( sigma2[n_al0] * ( COValpha + 0.0000000001 * diagmat(ones(n_al)) ) ) );
// //   Rprintf("set up adaptive proposals\n");
// //   
// //   double alpha0i, newalpha0i;
// //   vec newalpha0, newalpha, newlogZobl, newlogZsur;
// //   mat newlogZ;
// //   vec dummyvec;
// //   
// //   int count;
// //   vec xvals;
// //   double llprev, llthred, llnew, thetamin, thetamax, theta;
// //   vec prioretat;
// //   mat neweta;
// //   bool accept;
// //   
// //   mat Reta = exp(-distmat / phi_eta);
// //   mat invReta = inv(Reta), cholReta = chol(Reta);
// //   rowvec oneinvReta = trans(ones(nG)) * invReta;
// //   double oneinvRetaone = sum( oneinvReta * ones(nG) );
// //   Rprintf("set up GP covariance matrices\n");
// //   
// //   mat logZ = cppf_logZ_a0t_etat(Xzoop, alpha0, alpha, eta);
// //   vec logZobl = cppf_logZobl(logZ, indobl);
// //   vec logZsur = cppf_logZsur(logZ, indsur);
// //   Rprintf("set up initial values\n");
// //   
// //   
// //   // Start MCMC
// //   for(int s = 0; s < niter; s++) {
// //     
// //     // update alpha0
// //     for(int i = 0; i < n_al0; i ++){
// //       
// //       if( updateCOV ){
// //         if( (s+1 >= adaptInterval) && (s+1 - (adaptInterval * trunc((s+1) / adaptInterval)) == 0) ){
// //           dummyaccprob = accprob.col(i);
// //           rhat[i] = sum( dummyaccprob.rows(s+1-adaptInterval, s-1) ) / (adaptInterval-1);
// //           gamma1[i] = 1 / pow(adapIter[i], c1);
// //           gamma2[i] = c0 * gamma1[i];
// //           sigma2[i] = exp( log(sigma2[i]) + gamma2[i] * (rhat[i] - ropt) );
// //           
// //           
// //           dummyvec = postAlpha0.col(i);
// //           COValpha0[i] = COValpha0[i] + gamma1[i] * ( var( dummyvec.rows(s+1-adaptInterval, s-1) ) - COValpha0[i] );
// //           cholCOValpha0[i] = sqrt(  sigma2[i] * COValpha0[i] );
// //           
// //           adapIter[i] = adapIter[i] + 1;
// //         }
// //       }
// //       
// //       alpha0i = alpha0[i];
// //       newalpha0i = alpha0i + cholCOValpha0[i] * randn();
// //       newalpha0 = alpha0;
// //       newalpha0[i] = newalpha0i;
// //       
// //       newlogZ = cppf_logZ_a0t_etat(Xzoop, newalpha0, alpha, eta);
// //       newlogZobl = cppf_logZobl(newlogZ, indobl);
// //       newlogZsur = cppf_logZsur(newlogZ, indsur);
// //       
// //       logprob = cppf_loglik_obl(logYobl, newlogZobl, sig2_obl) -
// //         cppf_loglik_obl(logYobl, logZobl, sig2_obl) +
// //         cppf_loglik_sur(logYsur, newlogZsur, lam0, lam1, sig2_sur) -
// //         cppf_loglik_sur(logYsur, logZsur, lam0, lam1, sig2_sur) +
// //         normal_logh(newalpha0i, alpha0tilde, tau2) -
// //         normal_logh(alpha0i, alpha0tilde, tau2);
// //       
// //       if (log(randu()) < logprob) {
// //         
// //         alpha0 = newalpha0;
// //         logZ = newlogZ;
// //         logZobl = newlogZobl;
// //         logZsur = newlogZsur;
// //         
// //         accprob(s, i) = 1;
// //       }
// //     }
// //     postAlpha0.row(s) = trans(alpha0);
// //     if(s == 0){ Rprintf("sampled alpha0\n"); }
// //     
// //     
// //     // Gibbs update for alpha0tilde
// //     alpha0tilde = cppf_gibbs_alpha0tilde(alpha0, tau2);
// //     postAlpha0tilde[s] = alpha0tilde;
// //     if(s == 0){ Rprintf("sampled alpha0tilde\n"); }
// //     
// //     
// //     // update alpha
// //     if( updateCOV ){
// //       if( (s+1 >= adaptInterval) && (s+1 - (adaptInterval * trunc((s+1) / adaptInterval)) == 0) ){
// //         dummyaccprob = accprob.col(n_al0);
// //         rhat[n_al0] = sum( dummyaccprob.rows(s+1-adaptInterval, s-1) ) / (adaptInterval-1);
// //         gamma1[n_al0] = 1 / pow(adapIter[n_al0], c1);
// //         gamma2[n_al0] = c0 * gamma1[n_al0];
// //         sigma2[n_al0] = exp( log(sigma2[n_al0]) + gamma2[n_al0] * (rhat[n_al0] - ropt) );
// //         
// //         COValpha = COValpha + gamma1[n_al0] * ( cov( postAlpha.rows(s+1-adaptInterval, s-1) ) - COValpha );
// //         cholCOValpha = trans( chol( sigma2[n_al0] * ( COValpha + 0.0000000001 * diagmat(ones(n_al)) ) ) );
// //         
// //         adapIter[n_al0] = adapIter[n_al0] + 1;
// //       }
// //     }
// //     
// //     newalpha = alpha + cholCOValpha * randn(n_al);
// //     
// //     newlogZ = cppf_logZ_a0t_etat(Xzoop, alpha0, newalpha, eta);
// //     newlogZobl = cppf_logZobl(newlogZ, indobl);
// //     newlogZsur = cppf_logZsur(newlogZ, indsur);
// //     
// //     logprob = cppf_loglik_obl(logYobl, newlogZobl, sig2_obl) -
// //       cppf_loglik_obl(logYobl, logZobl, sig2_obl) +
// //       cppf_loglik_sur(logYsur, newlogZsur, lam0, lam1, sig2_sur) -
// //       cppf_loglik_sur(logYsur, logZsur, lam0, lam1, sig2_sur) +
// //       MVN_logh(newalpha, zeros(n_al), diagmat(ones(n_al))/100) -
// //       MVN_logh(alpha, zeros(n_al), diagmat(ones(n_al))/100);
// //     
// //     if (log(randu()) < logprob) {
// //       
// //       alpha = newalpha;
// //       logZ = newlogZ;
// //       logZobl = newlogZobl;
// //       logZsur = newlogZsur;
// //       
// //       accprob(s, n_al0) = 1;
// //     }
// //     
// //     postAlpha.row(s) = trans(alpha);
// //     if(s == 0){ Rprintf("sampled alpha\n"); }
// //     
// //     
// //     // update eta via Elliptical slice sampling
// //     for(int t = 0; t < ntzoop; t ++){
// //       neweta = eta;
// //       
// //       llprev = cppf_loglik_obl(logYobl, logZobl, sig2_obl) +
// //         cppf_loglik_sur(logYsur, logZsur, lam0, lam1, sig2_sur);
// //       prioretat = sqrt(kappa_eta) * trans(cholReta) * randn(nG);
// //       
// //       thetamin = 0;
// //       thetamax = 2 * M_PI;
// //       theta = thetamin + randu() * (thetamax - thetamin);
// //       thetamin = theta - 2 * M_PI;
// //       thetamax = theta;
// //       
// //       llthred = llprev + log(randu());
// //       accept = false;
// //       count = 0;
// //       
// //       while(accept == false){
// //         count = count + 1;
// //         
// //         neweta.col(t) = eta.col(t) * cos(theta) + prioretat * sin(theta);
// //         neweta.col(t) = neweta.col(t) - mean(neweta.col(t));
// //         
// //         newlogZ = cppf_logZ_a0t_etat(Xzoop, alpha0, alpha, neweta);
// //         newlogZobl = cppf_logZobl(newlogZ, indobl);
// //         newlogZsur = cppf_logZsur(newlogZ, indsur);
// //         
// //         llnew = cppf_loglik_obl(logYobl, newlogZobl, sig2_obl) +
// //           cppf_loglik_sur(logYsur, newlogZsur, lam0, lam1, sig2_sur);
// //         
// //         if(llnew > llthred){
// //           llprev = llnew;
// //           accept = true;
// //         } else {
// //           if(theta < 0){ thetamin = theta; } else { thetamax = theta; }
// //           theta = thetamin + randu() * (thetamax - thetamin);
// //           if ( (count) % 1000 == 0 ) {
// //             Rprintf("ESS for eta: %d iterations...\n", count);
// //           }
// //           if (count > 20000) stop("ESS is not converging...\n");
// //         }
// //       }
// //       eta = neweta;
// //       logZ = newlogZ;
// //       logZobl = newlogZobl;
// //       logZsur = newlogZsur;
// //       
// //       postEtat = as<mat>(postEta[t]);
// //       postEtat.row(s) = trans(eta.col(t));
// //       postEta[t] = postEtat;
// //     }
// //     if(s == 0){ Rprintf("sampled eta\n"); }
// //     
// //     
// //     // update kappa_eta
// //     kappa_eta = cppf_gibbs_kappa_eta(eta, invReta, 2, 2);
// //     postKappa(s,0) = kappa_eta;
// //     if(s == 0){ Rprintf("sampled kappa_eta\n"); }
// //     
// //     
// //     // update lam0
// //     lam0 = cppf_gibbs_lam0(logYsur, logZsur, lam1, sig2_sur, 100);
// //     postLam(s, 0) = lam0;
// //     if(s == 0){ Rprintf("sampled lam0\n"); }
// //     
// //     
// //     // update lam1
// //     lam1 = cppf_gibbs_lam1(logYsur, logZsur, lam0, sig2_sur, 100);
// //     postLam(s, 1) = lam1;
// //     if(s == 0){ Rprintf("sampled lam1\n"); }
// //     
// //     
// //     // update sig2_sur
// //     sig2_sur = cppf_gibbs_sig2_sur(logYsur, logZsur, lam0, lam1, 2, 2);
// //     postSig2(s, 0) = sig2_sur;
// //     if(s == 0){ Rprintf("sampled sig2_sur\n"); }
// //     
// //     
// //     if ( (s+1) % 10 == 0 ) {
// //       Rprintf("Generated %d samples...\n", s+1);
// //     }
// //   }
// //   
// //   
// //   return Rcpp::List::create(Rcpp::Named("postAlpha0") = postAlpha0,
// //                             Rcpp::Named("postAlpha0tilde") = postAlpha0tilde,
// //                             Rcpp::Named("postAlpha") = postAlpha,
// //                             Rcpp::Named("postEta") = postEta,
// //                             Rcpp::Named("postLam") = postLam,
// //                             Rcpp::Named("postSig2") = postSig2,
// //                             Rcpp::Named("postKappa") = postKappa,
// //                             Rcpp::Named("Accprob") = accprob,
// //                             Rcpp::Named("sigma2") = sigma2,
// //                             Rcpp::Named("adapIter") = adapIter,
// //                             Rcpp::Named("COValpha0") = COValpha0,
// //                             Rcpp::Named("COValpha") = COValpha);
// //   
// // }
// // 
// // 
// // 
// // // [[Rcpp::export]]
// // List cppf_fitWhale_b0t_psi(
// //     int niter, mat distmat,
// //     List Xwhale, mat Sdist, mat Wpam, vec dur_dist,
// //     vec f_dist_obs, List f_dist, List p_pam,
// //     mat inddist, vec indmdist, vec indmwhale,
// //     vec beta, mat psi, vec pii, vec c, double shape_c, double scale_c,
// //     double kappa_psi, double phi_psi,
// //     double cellarea,
// //     vec sigma2, mat COVbeta,
// //     bool updateCOV, int adaptInterval, double adaptFactorExponent, vec adapIter){
// //   
// //   double negativeInf = -std::numeric_limits<float>::infinity();;
// //   
// //   inddist = inddist-1, indmdist = indmdist-1, indmwhale = indmwhale-1;
// //   
// //   int nG = distmat.n_rows, n_be = beta.size(), ntwhale = psi.n_cols, n_sigma2 = sigma2.size();
// //   
// //   mat postKappa = zeros(niter, 1), postBeta = zeros(niter, n_be);
// //   mat postC = zeros(niter, c.size()), postPsit = zeros(niter, nG);
// //   List postPsi(ntwhale);
// //   for(int i = 0; i < ntwhale; i ++){
// //     postPsi[i] = zeros(niter, nG);
// //   }
// //   Rprintf("set up matrices for posterior samples\n");
// //   
// //   mat accprob = zeros(niter, n_sigma2);
// //   vec rhat = zeros(n_sigma2), gamma1 = zeros(n_sigma2), gamma2 = zeros(n_sigma2), dummyaccprob;
// //   double c0 = 1, c1 = adaptFactorExponent, ropt = 0.234, logprob;
// //   mat cholCOVbeta = trans( chol( sigma2[0] * ( COVbeta + 0.0000000001 * diagmat(ones(n_be)) ) ) );
// //   Rprintf("set up adaptive proposals\n");
// //   
// //   vec newbeta;
// //   mat newloglam, newlam;
// //   vec dummyvec;
// //   
// //   int count;
// //   vec xvals;
// //   double llprev, llthred, llnew, thetamin, thetamax, theta;
// //   vec prioreta, priorpsit, neweta;
// //   mat newpsi;
// //   bool accept;
// //   
// //   mat Rpsi = exp(-distmat / phi_psi);
// //   mat invRpsi= inv(Rpsi), cholRpsi= chol(Rpsi);
// //   rowvec oneinvRpsi = trans(ones(nG)) * invRpsi;
// //   double oneinvRpsione = sum( oneinvRpsi * ones(nG) );
// //   Rprintf("set up GP covariance matrices\n");
// //   
// //   mat loglam = cppf_loglam_b0t_woZ(Xwhale, beta, psi);
// //   mat lam = exp(loglam);
// //   Rprintf("set up initial values\n");
// //   
// //   
// //   // Start MCMC
// //   for(int s = 0; s < niter; s++) {
// //     
// //     // update beta
// //     if( updateCOV ){
// //       if( (s+1 >= adaptInterval) && (s+1 - (adaptInterval * trunc((s+1) / adaptInterval)) == 0) ){
// //         dummyaccprob = accprob.col(0);
// //         rhat[0] = sum( dummyaccprob.rows(s+1-adaptInterval, s-1) ) / (adaptInterval-1);
// //         gamma1[0] = 1 / pow(adapIter[0], c1);
// //         gamma2[0] = c0 * gamma1[0];
// //         sigma2[0] = exp( log(sigma2[0]) + gamma2[0] * (rhat[0] - ropt) );
// //         
// //         COVbeta = COVbeta + gamma1[0] * ( cov( postBeta.rows(s+1-adaptInterval, s-1)) - COVbeta );
// //         cholCOVbeta = trans( chol( sigma2[0] * ( COVbeta + 0.0000000001 * diagmat(ones(n_be)) ) ) );
// //         
// //         adapIter[0] = adapIter[0] + 1;
// //       }
// //     }
// //     
// //     newbeta = beta + cholCOVbeta * randn(n_be);
// //     
// //     newloglam = cppf_loglam_b0t_woZ(Xwhale, newbeta, psi);
// //     newlam = exp(newloglam);
// //     
// //     logprob = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) -
// //       cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //       cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) -
// //       cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
// //       MVN_logh(newbeta, zeros(n_be), diagmat(ones(n_be))/100) -
// //       MVN_logh(beta, zeros(n_be), diagmat(ones(n_be))/100);
// //     
// //     if (log(randu()) < logprob) {
// //       
// //       beta = newbeta;
// //       loglam = newloglam;
// //       lam = newlam;
// //       
// //       accprob(s, 0) = 1;
// //     }
// //     
// //     postBeta.row(s) = trans(beta);
// //     if(s == 0){ Rprintf("sampled beta\n"); }
// //     
// //     
// //     // update psi via elliptical slice sampling
// //     for(int t = 0; t < ntwhale; t ++){
// //       newpsi = psi;
// //       
// //       llprev = cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //         cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea);
// //       priorpsit = sqrt(kappa_psi) * trans(cholRpsi) * randn(nG);
// //       
// //       thetamin = 0;
// //       thetamax = 2 * M_PI;
// //       theta = thetamin + randu() * (thetamax - thetamin);
// //       thetamin = theta - 2 * M_PI;
// //       thetamax = theta;
// //       
// //       llthred = llprev + log(randu());
// //       accept = false;
// //       count = 0;
// //       
// //       while(accept == false){
// //         count = count + 1;
// //         
// //         newpsi.col(t) = psi.col(t) * cos(theta) + priorpsit * sin(theta);
// //         newpsi.col(t) = newpsi.col(t) - mean(newpsi.col(t));
// //         
// //         newloglam = cppf_loglam_b0t_woZ(Xwhale, beta, newpsi);
// //         newlam = exp(newloglam);
// //         
// //         llnew = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //           cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea);
// //         
// //         if(llnew > llthred){
// //           llprev = llnew;
// //           accept = true;
// //         } else {
// //           if(theta < 0){ thetamin = theta; } else { thetamax = theta; }
// //           theta = thetamin + randu() * (thetamax - thetamin);
// //           if ( (count) % 1000 == 0 ) {
// //             Rprintf("ESS for psi: %d iterations...\n", count);
// //           }
// //           if (count > 20000) stop("ESS is not converging...\n");
// //         }
// //       }
// //       psi = newpsi;
// //       loglam = newloglam;
// //       lam = newlam;
// //       
// //       postPsit = as<mat>(postPsi[t]);
// //       postPsit.row(s) = trans(psi.col(t));
// //       postPsi[t] = postPsit;
// //     }
// //     if(s == 0){ Rprintf("sampled psi\n"); }
// //     
// //     
// //     // update kappa_psi
// //     kappa_psi = cppf_gibbs_kappa_psis(psi, invRpsi, 2, 2);
// //     postKappa(s,0) = kappa_psi;
// //     if(s == 0){ Rprintf("sampled kappa_psi\n"); }
// //     
// //     
// //     // update c
// //     c = cppf_gibbs_ct(Wpam, lam, dur_dist, p_pam, cellarea, shape_c, scale_c);
// //     postC.row(s) = trans(c);
// //     if(s == 0){ Rprintf("sampled c\n"); }
// //     
// //     
// //     if ( (s+1) % 10 == 0 ) {
// //       Rprintf("Generated %d samples...\n", s+1);
// //     }
// //   }
// //   
// //   
// //   return Rcpp::List::create(Rcpp::Named("postKappa") = postKappa,
// //                             Rcpp::Named("postBeta") = postBeta,
// //                             Rcpp::Named("postPsi") = postPsi,
// //                             Rcpp::Named("postC") = postC,
// //                             Rcpp::Named("Accprob") = accprob,
// //                             Rcpp::Named("sigma2") = sigma2,
// //                             Rcpp::Named("adapIter") = adapIter,
// //                             Rcpp::Named("COVbeta") = COVbeta);
// //   
// // }
// // 
// // 
// // 
// // // =============================================================================
// // // fit the joint model after some number of iterations ----
// // // =============================================================================
// // 
// // // sig2_obl is known
// // // time-varying intercept for whale component
// // // [[Rcpp::export]]
// // List cppf_fitSeparJoint_b0t(
// //     int niter, mat distmat,
// //     List Xzoop, vec logYobl, vec logYsur,
// //     List Xwhale, mat Sdist, mat Wpam, vec dur_dist,
// //     vec f_dist_obs, List f_dist, List p_pam,
// //     mat indobl, mat indsur, vec indtwhale, mat inddist, vec indmdist, vec indmwhale,
// //     double alpha0tilde, vec alpha0, vec alpha, mat eta, double lam0, double lam1,
// //     double tau2, double kappa_eta, double phi_eta, double sig2_obl, double sig2_sur, 
// //     vec beta, vec pii, vec c, double shape_c, double scale_c,
// //     double cellarea,
// //     vec sigma2, vec COValpha0, mat COValpha, mat COVbeta,
// //     bool updateCOV, int adaptInterval, double adaptFactorExponent, vec adapIter,
// //     int curniter, int joiniter){
// //   
// //   double negativeInf = -std::numeric_limits<float>::infinity();;
// //   
// //   indobl = indobl-1, indsur = indsur-1, indtwhale = indtwhale-1;
// //   inddist = inddist-1, indmdist = indmdist-1, indmwhale = indmwhale-1;
// //   
// //   int nG = distmat.n_rows, n_al0 = alpha0.size(), n_al = alpha.size(), n_be = beta.size(), ntzoop = eta.n_cols, n_sigma2 = sigma2.size();
// //   vec dummyvec;
// //   
// //   mat postAlpha0tilde = zeros(niter, 1), postAlpha0 = zeros(niter, n_al0), postAlpha = zeros(niter, n_al), postLam = zeros(niter, 2);
// //   mat postSig2 = zeros(niter, 1), postKappa = zeros(niter, 1);
// //   mat postBeta = zeros(niter, n_be);
// //   mat postC = zeros(niter, c.size()), postEtat = zeros(niter, nG);
// //   List postEta(ntzoop);
// //   for(int i = 0; i < ntzoop; i ++){
// //     postEta[i] = zeros(niter, nG);
// //   }
// //   Rprintf("set up matrices for posterior samples\n");
// //   
// //   mat accprob = zeros(niter, n_sigma2);
// //   vec rhat = zeros(n_sigma2), gamma1 = zeros(n_sigma2), gamma2 = zeros(n_sigma2), dummyaccprob;
// //   double c0 = 1, c1 = adaptFactorExponent, ropt = 0.234, logprob;
// //   vec cholCOValpha0(n_al0);
// //   for(int i = 0; i < n_al0; i ++){ cholCOValpha0[i] = sqrt(sigma2[i] * COValpha0[i]); }
// //   mat cholCOValpha = trans( chol( sigma2[n_al0] * ( COValpha + 0.0000000001 * diagmat(ones(n_al)) ) ) );
// //   mat cholCOVbeta = trans( chol( sigma2[n_al0+1] * ( COVbeta + 0.0000000001 * diagmat(ones(n_be)) ) ) );
// //   Rprintf("set up adaptive proposals\n");
// //   
// //   double alpha0i, newalpha0i;
// //   vec newalpha0, newalpha, newbeta, newlogZobl, newlogZsur;
// //   mat newlogZ, newloglam, newlam;
// //   
// //   int count;
// //   vec xvals;
// //   double llprev, llthred, llnew, thetamin, thetamax, theta;
// //   vec prioretat;
// //   mat neweta;
// //   bool accept;
// //   
// //   mat Reta = exp(-distmat / phi_eta);
// //   mat invReta = inv(Reta),  cholReta = chol(Reta);
// //   rowvec oneinvReta = trans(ones(nG)) * invReta;
// //   double oneinvRetaone = sum( oneinvReta * ones(nG) );
// //   Rprintf("set up GP covariance matrices\n");
// //   
// //   mat logZ = cppf_logZ_a0t_etat(Xzoop, alpha0, alpha, eta);
// //   vec logZobl = cppf_logZobl(logZ, indobl);
// //   vec logZsur = cppf_logZsur(logZ, indsur);
// //   mat loglam = cppf_loglam_b0t(Xwhale, logZ, beta, indtwhale);
// //   mat lam = exp(loglam);
// //   Rprintf("set up initial values\n");
// //   
// //   
// //   // Start MCMC
// //   for(int s = 0; s < niter; s++) {
// //     
// //     // update alpha0
// //     for(int i = 0; i < n_al0; i ++){
// //       
// //       if( updateCOV ){
// //         if( (s+1 >= adaptInterval) && (s+1 - (adaptInterval * trunc((s+1) / adaptInterval)) == 0) ){
// //           dummyaccprob = accprob.col(i);
// //           rhat[i] = sum( dummyaccprob.rows(s+1-adaptInterval, s-1) ) / (adaptInterval-1);
// //           gamma1[i] = 1 / pow(adapIter[i], c1);
// //           gamma2[i] = c0 * gamma1[i];
// //           sigma2[i] = exp( log(sigma2[i]) + gamma2[i] * (rhat[i] - ropt) );
// //           
// //           
// //           dummyvec = postAlpha0.col(i);
// //           COValpha0[i] = COValpha0[i] + gamma1[i] * ( var( dummyvec.rows(s+1-adaptInterval, s-1) ) - COValpha0[i] );
// //           cholCOValpha0[i] = sqrt(  sigma2[i] * COValpha0[i] );
// //           
// //           adapIter[i] = adapIter[i] + 1;
// //         }
// //       }
// //       
// //       alpha0i = alpha0[i];
// //       newalpha0i = alpha0i + cholCOValpha0[i] * randn();
// //       newalpha0 = alpha0;
// //       newalpha0[i] = newalpha0i;
// //       
// //       newlogZ = cppf_logZ_a0t_etat(Xzoop, newalpha0, alpha, eta);
// //       newlogZobl = cppf_logZobl(newlogZ, indobl);
// //       newlogZsur = cppf_logZsur(newlogZ, indsur);
// //       newloglam = cppf_loglam_b0t(Xwhale, newlogZ, beta, indtwhale);
// //       newlam = exp(newloglam);
// //       
// //       if(curniter > joiniter){
// //         logprob = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) -
// //           cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //           cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) -
// //           cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
// //           cppf_loglik_obl(logYobl, newlogZobl, sig2_obl) -
// //           cppf_loglik_obl(logYobl, logZobl, sig2_obl) +
// //           cppf_loglik_sur(logYsur, newlogZsur, lam0, lam1, sig2_sur) -
// //           cppf_loglik_sur(logYsur, logZsur, lam0, lam1, sig2_sur) +
// //           normal_logh(newalpha0i, alpha0tilde, tau2) -
// //           normal_logh(alpha0i, alpha0tilde, tau2);  
// //         
// //       } else {
// //         logprob = cppf_loglik_obl(logYobl, newlogZobl, sig2_obl) -
// //           cppf_loglik_obl(logYobl, logZobl, sig2_obl) +
// //           cppf_loglik_sur(logYsur, newlogZsur, lam0, lam1, sig2_sur) -
// //           cppf_loglik_sur(logYsur, logZsur, lam0, lam1, sig2_sur) +
// //           normal_logh(newalpha0i, alpha0tilde, tau2) -
// //           normal_logh(alpha0i, alpha0tilde, tau2);
// //       }
// //       
// //       if (log(randu()) < logprob) {
// //         
// //         alpha0 = newalpha0;
// //         logZ = newlogZ;
// //         logZobl = newlogZobl;
// //         logZsur = newlogZsur;
// //         loglam = newloglam;
// //         lam = newlam;
// //         
// //         accprob(s, i) = 1;
// //       }
// //     }
// //     postAlpha0.row(s) = trans(alpha0);
// //     if(s == 0){ Rprintf("sampled alpha0\n"); }
// //     
// //     
// //     // Gibbs update for alpha0tilde
// //     alpha0tilde = cppf_gibbs_alpha0tilde(alpha0, tau2);
// //     postAlpha0tilde[s] = alpha0tilde;
// //     if(s == 0){ Rprintf("sampled alpha0tilde\n"); }
// //     
// //     
// //     // update alpha
// //     if( updateCOV ){
// //       if( (s+1 >= adaptInterval) && (s+1 - (adaptInterval * trunc((s+1) / adaptInterval)) == 0) ){
// //         dummyaccprob = accprob.col(n_al0);
// //         rhat[n_al0] = sum( dummyaccprob.rows(s+1-adaptInterval, s-1) ) / (adaptInterval-1);
// //         gamma1[n_al0] = 1 / pow(adapIter[n_al0], c1);
// //         gamma2[n_al0] = c0 * gamma1[n_al0];
// //         sigma2[n_al0] = exp( log(sigma2[n_al0]) + gamma2[n_al0] * (rhat[n_al0] - ropt) );
// //         
// //         COValpha = COValpha + gamma1[n_al0] * ( cov( postAlpha.rows(s+1-adaptInterval, s-1) ) - COValpha );
// //         cholCOValpha = trans( chol( sigma2[n_al0] * ( COValpha + 0.0000000001 * diagmat(ones(n_al)) ) ) );
// //         
// //         adapIter[n_al0] = adapIter[n_al0] + 1;
// //       }
// //     }
// //     
// //     newalpha = alpha + cholCOValpha * randn(n_al);
// //     
// //     newlogZ = cppf_logZ_a0t_etat(Xzoop, alpha0, newalpha, eta);
// //     newlogZobl = cppf_logZobl(newlogZ, indobl);
// //     newlogZsur = cppf_logZsur(newlogZ, indsur);
// //     newloglam = cppf_loglam_b0t(Xwhale, newlogZ, beta, indtwhale);
// //     newlam = exp(newloglam);
// //     
// //     if(curniter > joiniter){
// //       logprob = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) -
// //         cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //         cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) -
// //         cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
// //         cppf_loglik_obl(logYobl, newlogZobl, sig2_obl) -
// //         cppf_loglik_obl(logYobl, logZobl, sig2_obl) +
// //         cppf_loglik_sur(logYsur, newlogZsur, lam0, lam1, sig2_sur) -
// //         cppf_loglik_sur(logYsur, logZsur, lam0, lam1, sig2_sur) +
// //         MVN_logh(newalpha, zeros(n_al), diagmat(ones(n_al))/100) -
// //         MVN_logh(alpha, zeros(n_al), diagmat(ones(n_al))/100);
// //       
// //     } else {
// //       logprob = cppf_loglik_obl(logYobl, newlogZobl, sig2_obl) -
// //         cppf_loglik_obl(logYobl, logZobl, sig2_obl) +
// //         cppf_loglik_sur(logYsur, newlogZsur, lam0, lam1, sig2_sur) -
// //         cppf_loglik_sur(logYsur, logZsur, lam0, lam1, sig2_sur) +
// //         MVN_logh(newalpha, zeros(n_al), diagmat(ones(n_al))/100) -
// //         MVN_logh(alpha, zeros(n_al), diagmat(ones(n_al))/100);
// //     }
// //     
// //     if (log(randu()) < logprob) {
// //       
// //       alpha = newalpha;
// //       logZ = newlogZ;
// //       logZobl = newlogZobl;
// //       logZsur = newlogZsur;
// //       loglam = newloglam;
// //       lam = newlam;
// //       
// //       accprob(s, n_al0) = 1;
// //     }
// //     
// //     postAlpha.row(s) = trans(alpha);
// //     if(s == 0){ Rprintf("sampled alpha\n"); }
// //     
// //     
// //     // update eta via Elliptical slice sampling
// //     for(int t = 0; t < ntzoop; t ++){
// //       neweta = eta;
// //       
// //       if(curniter > joiniter){
// //         llprev = cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //           cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
// //           cppf_loglik_obl(logYobl, logZobl, sig2_obl) +
// //           cppf_loglik_sur(logYsur, logZsur, lam0, lam1, sig2_sur);
// //         prioretat = sqrt(kappa_eta) * trans(cholReta) * randn(nG);
// //         
// //       } else {
// //         llprev = cppf_loglik_obl(logYobl, logZobl, sig2_obl) +
// //           cppf_loglik_sur(logYsur, logZsur, lam0, lam1, sig2_sur);
// //         prioretat = sqrt(kappa_eta) * trans(cholReta) * randn(nG);
// //       }
// //       
// //       thetamin = 0;
// //       thetamax = 2 * M_PI;
// //       theta = thetamin + randu() * (thetamax - thetamin);
// //       thetamin = theta - 2 * M_PI;
// //       thetamax = theta;
// //       
// //       llthred = llprev + log(randu());
// //       accept = false;
// //       count = 0;
// //       
// //       while(accept == false){
// //         count = count + 1;
// //         
// //         neweta.col(t) = eta.col(t) * cos(theta) + prioretat * sin(theta);
// //         neweta.col(t) = neweta.col(t) - mean(neweta.col(t));
// //         
// //         newlogZ = cppf_logZ_a0t_etat(Xzoop, alpha0, alpha, neweta);
// //         newlogZobl = cppf_logZobl(newlogZ, indobl);
// //         newlogZsur = cppf_logZsur(newlogZ, indsur);
// //         newloglam = cppf_loglam_b0t(Xwhale, newlogZ, beta, indtwhale);
// //         newlam = exp(newloglam);
// //         
// //         if(curniter > joiniter){
// //           llnew = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //             cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) +
// //             cppf_loglik_obl(logYobl, newlogZobl, sig2_obl) +
// //             cppf_loglik_sur(logYsur, newlogZsur, lam0, lam1, sig2_sur);
// //           
// //         } else {
// //           llnew = cppf_loglik_obl(logYobl, newlogZobl, sig2_obl) +
// //             cppf_loglik_sur(logYsur, newlogZsur, lam0, lam1, sig2_sur);
// //         }
// //         
// //         if(llnew > llthred){
// //           llprev = llnew;
// //           accept = true;
// //         } else {
// //           if(theta < 0){ thetamin = theta; } else { thetamax = theta; }
// //           theta = thetamin + randu() * (thetamax - thetamin);
// //           if ( (count) % 1000 == 0 ) {
// //             Rprintf("ESS for eta: %d iterations...\n", count);
// //           }
// //           if (count > 20000) stop("ESS is not converging...\n");
// //         }
// //       }
// //       eta = neweta;
// //       logZ = newlogZ;
// //       logZobl = newlogZobl;
// //       logZsur = newlogZsur;
// //       loglam = newloglam;
// //       lam = newlam;
// //       
// //       postEtat = as<mat>(postEta[t]);
// //       postEtat.row(s) = trans(eta.col(t));
// //       postEta[t] = postEtat;
// //     }
// //     if(s == 0){ Rprintf("sampled eta\n"); }
// //     
// //     
// //     // update kappa_eta
// //     kappa_eta = cppf_gibbs_kappa_eta(eta, invReta, 2, 2);
// //     postKappa(s,0) = kappa_eta;
// //     if(s == 0){ Rprintf("sampled kappa_eta\n"); }
// //     
// //     
// //     // update lam0
// //     lam0 = cppf_gibbs_lam0(logYsur, logZsur, lam1, sig2_sur, 100);
// //     postLam(s, 0) = lam0;
// //     if(s == 0){ Rprintf("sampled lam0\n"); }
// //     
// //     
// //     // update lam1
// //     lam1 = cppf_gibbs_lam1(logYsur, logZsur, lam0, sig2_sur, 100);
// //     postLam(s, 1) = lam1;
// //     if(s == 0){ Rprintf("sampled lam1\n"); }
// //     
// //     
// //     // update sig2_sur
// //     sig2_sur = cppf_gibbs_sig2_sur(logYsur, logZsur, lam0, lam1, 2, 2);
// //     postSig2(s, 0) = sig2_sur;
// //     if(s == 0){ Rprintf("sampled sig2_sur\n"); }
// //     
// //     
// //     // update beta
// //     if( updateCOV ){
// //       if( (s+1 >= adaptInterval) && (s+1 - (adaptInterval * trunc((s+1) / adaptInterval)) == 0) ){
// //         dummyaccprob = accprob.col(n_al0+1);
// //         rhat[n_al0+1] = sum( dummyaccprob.rows(s+1-adaptInterval, s-1) ) / (adaptInterval-1);
// //         gamma1[n_al0+1] = 1 / pow(adapIter[n_al0+1], c1);
// //         gamma2[n_al0+1] = c0 * gamma1[n_al0+1];
// //         sigma2[n_al0+1] = exp( log(sigma2[n_al0+1]) + gamma2[n_al0+1] * (rhat[n_al0+1] - ropt) );
// //         
// //         COVbeta = COVbeta + gamma1[n_al0+1] * ( cov( postBeta.rows(s+1-adaptInterval, s-1)) - COVbeta );
// //         cholCOVbeta = trans( chol( sigma2[n_al0+1] * ( COVbeta + 0.0000000001 * diagmat(ones(n_be)) ) ) );
// //         
// //         adapIter[n_al0+1] = adapIter[n_al0+1] + 1;
// //       }
// //     }
// //     
// //     newbeta = beta + cholCOVbeta * randn(n_be);
// //     
// //     newloglam = cppf_loglam_b0t(Xwhale, logZ, newbeta, indtwhale);
// //     newlam = exp(newloglam);
// //     
// //     logprob = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) -
// //       cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //       cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) -
// //       cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
// //       MVN_logh(newbeta, zeros(n_be), diagmat(ones(n_be))/100) -
// //       MVN_logh(beta, zeros(n_be), diagmat(ones(n_be))/100);
// //     
// //     if (log(randu()) < logprob) {
// //       
// //       beta = newbeta;
// //       loglam = newloglam;
// //       lam = newlam;
// //       
// //       accprob(s, n_al0+1) = 1;
// //     }
// //     
// //     postBeta.row(s) = trans(beta);
// //     if(s == 0){ Rprintf("sampled beta\n"); }
// //     
// //     
// //     // update c
// //     c = cppf_gibbs_ct(Wpam, lam, dur_dist, p_pam, cellarea, shape_c, scale_c);
// //     postC.row(s) = trans(c);
// //     if(s == 0){ Rprintf("sampled c\n"); }
// //     
// //     
// //     if ( (s+1) % 10 == 0 ) {
// //       Rprintf("Generated %d samples...\n", s+1);
// //     }
// //   }
// //   
// //   return Rcpp::List::create(Rcpp::Named("postAlpha0") = postAlpha0,
// //                             Rcpp::Named("postAlpha0tilde") = postAlpha0tilde,
// //                             Rcpp::Named("postAlpha") = postAlpha,
// //                             Rcpp::Named("postEta") = postEta,
// //                             Rcpp::Named("postLam") = postLam,
// //                             Rcpp::Named("postKappa") = postKappa,
// //                             Rcpp::Named("postSig2") = postSig2,
// //                             Rcpp::Named("postBeta") = postBeta,
// //                             Rcpp::Named("postC") = postC,
// //                             Rcpp::Named("Accprob") = accprob,
// //                             Rcpp::Named("sigma2") = sigma2,
// //                             Rcpp::Named("adapIter") = adapIter,
// //                             Rcpp::Named("COValpha0") = COValpha0,
// //                             Rcpp::Named("COValpha") = COValpha,
// //                             Rcpp::Named("COVbeta") = COVbeta);
// //   
// // }
// // 
// // 
// // 
// // 
// // 
// // 
// // // =============================================================================
// // // fix call rates at truth ----
// // // =============================================================================
// // 
// // 
// // // sig2_obl is known
// // // time-varying intercept for whale component
// // // [[Rcpp::export]]
// // List cppf_fitJoint_b0t_fixC(
// //     int niter, mat distmat,
// //     List Xzoop, vec logYobl, vec logYsur,
// //     List Xwhale, mat Sdist, mat Wpam, vec dur_dist,
// //     vec f_dist_obs, List f_dist, List p_pam,
// //     mat indobl, mat indsur, vec indtwhale, mat inddist, vec indmdist, vec indmwhale,
// //     double alpha0tilde, vec alpha0, vec alpha, mat eta, double lam0, double lam1,
// //     double tau2, double kappa_eta, double phi_eta, double sig2_obl, double sig2_sur, 
// //     vec beta, vec pii, vec c,
// //     double cellarea,
// //     vec sigma2, vec COValpha0, mat COValpha, mat COVbeta,
// //     bool updateCOV, int adaptInterval, double adaptFactorExponent, vec adapIter){
// //   
// //   double negativeInf = -std::numeric_limits<float>::infinity();;
// //   
// //   indobl = indobl-1, indsur = indsur-1, indtwhale = indtwhale-1;
// //   inddist = inddist-1, indmdist = indmdist-1, indmwhale = indmwhale-1;
// //   
// //   int nG = distmat.n_rows, n_al0 = alpha0.size(), n_al = alpha.size(), n_be = beta.size(), ntzoop = eta.n_cols, n_sigma2 = sigma2.size();
// //   vec dummyvec;
// //   
// //   mat postAlpha0tilde = zeros(niter, 1), postAlpha0 = zeros(niter, n_al0), postAlpha = zeros(niter, n_al), postLam = zeros(niter, 2);
// //   mat postSig2 = zeros(niter, 1), postKappa = zeros(niter, 1);
// //   mat postBeta = zeros(niter, n_be);
// //   mat postEtat = zeros(niter, nG);
// //   List postEta(ntzoop);
// //   for(int i = 0; i < ntzoop; i ++){
// //     postEta[i] = zeros(niter, nG);
// //   }
// //   Rprintf("set up matrices for posterior samples\n");
// //   
// //   mat accprob = zeros(niter, n_sigma2);
// //   vec rhat = zeros(n_sigma2), gamma1 = zeros(n_sigma2), gamma2 = zeros(n_sigma2), dummyaccprob;
// //   double c0 = 1, c1 = adaptFactorExponent, ropt = 0.234, logprob;
// //   vec cholCOValpha0(n_al0);
// //   for(int i = 0; i < n_al0; i ++){ cholCOValpha0[i] = sqrt(sigma2[i] * COValpha0[i]); }
// //   mat cholCOValpha = trans( chol( sigma2[n_al0] * ( COValpha + 0.0000000001 * diagmat(ones(n_al)) ) ) );
// //   mat cholCOVbeta = trans( chol( sigma2[n_al0+1] * ( COVbeta + 0.0000000001 * diagmat(ones(n_be)) ) ) );
// //   Rprintf("set up adaptive proposals\n");
// //   
// //   double alpha0i, newalpha0i;
// //   vec newalpha0, newalpha, newbeta, newlogZobl, newlogZsur;
// //   mat newlogZ, newloglam, newlam;
// //   
// //   int count;
// //   vec xvals;
// //   double llprev, llthred, llnew, thetamin, thetamax, theta;
// //   vec prioretat;
// //   mat neweta;
// //   bool accept;
// //   
// //   mat Reta = exp(-distmat / phi_eta);
// //   mat invReta = inv(Reta),  cholReta = chol(Reta);
// //   rowvec oneinvReta = trans(ones(nG)) * invReta;
// //   double oneinvRetaone = sum( oneinvReta * ones(nG) );
// //   Rprintf("set up GP covariance matrices\n");
// //   
// //   mat logZ = cppf_logZ_a0t_etat(Xzoop, alpha0, alpha, eta);
// //   vec logZobl = cppf_logZobl(logZ, indobl);
// //   vec logZsur = cppf_logZsur(logZ, indsur);
// //   Rprintf("set up initial values for zoop\n");
// //   mat loglam = cppf_loglam_b0t(Xwhale, logZ, beta, indtwhale);
// //   mat lam = exp(loglam);
// //   Rprintf("set up initial values\n");
// //   
// //   
// //   // Start MCMC
// //   for(int s = 0; s < niter; s++) {
// //     
// //     // update alpha0
// //     for(int i = 0; i < n_al0; i ++){
// //       
// //       if( updateCOV ){
// //         if( (s+1 >= adaptInterval) && (s+1 - (adaptInterval * trunc((s+1) / adaptInterval)) == 0) ){
// //           dummyaccprob = accprob.col(i);
// //           rhat[i] = sum( dummyaccprob.rows(s+1-adaptInterval, s-1) ) / (adaptInterval-1);
// //           gamma1[i] = 1 / pow(adapIter[i], c1);
// //           gamma2[i] = c0 * gamma1[i];
// //           sigma2[i] = exp( log(sigma2[i]) + gamma2[i] * (rhat[i] - ropt) );
// //           
// //           
// //           dummyvec = postAlpha0.col(i);
// //           COValpha0[i] = COValpha0[i] + gamma1[i] * ( var( dummyvec.rows(s+1-adaptInterval, s-1) ) - COValpha0[i] );
// //           cholCOValpha0[i] = sqrt(  sigma2[i] * COValpha0[i] );
// //           
// //           adapIter[i] = adapIter[i] + 1;
// //         }
// //       }
// //       
// //       alpha0i = alpha0[i];
// //       newalpha0i = alpha0i + cholCOValpha0[i] * randn();
// //       newalpha0 = alpha0;
// //       newalpha0[i] = newalpha0i;
// //       
// //       newlogZ = cppf_logZ_a0t_etat(Xzoop, newalpha0, alpha, eta);
// //       newlogZobl = cppf_logZobl(newlogZ, indobl);
// //       newlogZsur = cppf_logZsur(newlogZ, indsur);
// //       newloglam = cppf_loglam_b0t(Xwhale, newlogZ, beta, indtwhale);
// //       newlam = exp(newloglam);
// //       
// //       logprob = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) -
// //         cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //         cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) -
// //         cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
// //         cppf_loglik_obl(logYobl, newlogZobl, sig2_obl) -
// //         cppf_loglik_obl(logYobl, logZobl, sig2_obl) +
// //         cppf_loglik_sur(logYsur, newlogZsur, lam0, lam1, sig2_sur) -
// //         cppf_loglik_sur(logYsur, logZsur, lam0, lam1, sig2_sur) +
// //         normal_logh(newalpha0i, alpha0tilde, tau2) -
// //         normal_logh(alpha0i, alpha0tilde, tau2);
// //       
// //       if (log(randu()) < logprob) {
// //         
// //         alpha0 = newalpha0;
// //         logZ = newlogZ;
// //         logZobl = newlogZobl;
// //         logZsur = newlogZsur;
// //         loglam = newloglam;
// //         lam = newlam;
// //         
// //         accprob(s, i) = 1;
// //       }
// //     }
// //     postAlpha0.row(s) = trans(alpha0);
// //     if(s == 0){ Rprintf("sampled alpha0\n"); }
// //     
// //     
// //     // Gibbs update for alpha0tilde
// //     alpha0tilde = cppf_gibbs_alpha0tilde(alpha0, tau2);
// //     postAlpha0tilde[s] = alpha0tilde;
// //     if(s == 0){ Rprintf("sampled alpha0tilde\n"); }
// //     
// //     
// //     // update alpha
// //     if( updateCOV ){
// //       if( (s+1 >= adaptInterval) && (s+1 - (adaptInterval * trunc((s+1) / adaptInterval)) == 0) ){
// //         dummyaccprob = accprob.col(n_al0);
// //         rhat[n_al0] = sum( dummyaccprob.rows(s+1-adaptInterval, s-1) ) / (adaptInterval-1);
// //         gamma1[n_al0] = 1 / pow(adapIter[n_al0], c1);
// //         gamma2[n_al0] = c0 * gamma1[n_al0];
// //         sigma2[n_al0] = exp( log(sigma2[n_al0]) + gamma2[n_al0] * (rhat[n_al0] - ropt) );
// //         
// //         COValpha = COValpha + gamma1[n_al0] * ( cov( postAlpha.rows(s+1-adaptInterval, s-1) ) - COValpha );
// //         cholCOValpha = trans( chol( sigma2[n_al0] * ( COValpha + 0.0000000001 * diagmat(ones(n_al)) ) ) );
// //         
// //         adapIter[n_al0] = adapIter[n_al0] + 1;
// //       }
// //     }
// //     
// //     newalpha = alpha + cholCOValpha * randn(n_al);
// //     
// //     newlogZ = cppf_logZ_a0t_etat(Xzoop, alpha0, newalpha, eta);
// //     newlogZobl = cppf_logZobl(newlogZ, indobl);
// //     newlogZsur = cppf_logZsur(newlogZ, indsur);
// //     newloglam = cppf_loglam_b0t(Xwhale, newlogZ, beta, indtwhale);
// //     newlam = exp(newloglam);
// //     
// //     logprob = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) -
// //       cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //       cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) -
// //       cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
// //       cppf_loglik_obl(logYobl, newlogZobl, sig2_obl) -
// //       cppf_loglik_obl(logYobl, logZobl, sig2_obl) +
// //       cppf_loglik_sur(logYsur, newlogZsur, lam0, lam1, sig2_sur) -
// //       cppf_loglik_sur(logYsur, logZsur, lam0, lam1, sig2_sur) +
// //       MVN_logh(newalpha, zeros(n_al), diagmat(ones(n_al))/100) -
// //       MVN_logh(alpha, zeros(n_al), diagmat(ones(n_al))/100);
// //     
// //     if (log(randu()) < logprob) {
// //       
// //       alpha = newalpha;
// //       logZ = newlogZ;
// //       logZobl = newlogZobl;
// //       logZsur = newlogZsur;
// //       loglam = newloglam;
// //       lam = newlam;
// //       
// //       accprob(s, n_al0) = 1;
// //     }
// //     
// //     postAlpha.row(s) = trans(alpha);
// //     if(s == 0){ Rprintf("sampled alpha\n"); }
// //     
// //     
// //     // update eta via Elliptical slice sampling
// //     for(int t = 0; t < ntzoop; t ++){
// //       neweta = eta;
// //       
// //       llprev = cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //         cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
// //         cppf_loglik_obl(logYobl, logZobl, sig2_obl) +
// //         cppf_loglik_sur(logYsur, logZsur, lam0, lam1, sig2_sur);
// //       prioretat = sqrt(kappa_eta) * trans(cholReta) * randn(nG);
// //       
// //       thetamin = 0;
// //       thetamax = 2 * M_PI;
// //       theta = thetamin + randu() * (thetamax - thetamin);
// //       thetamin = theta - 2 * M_PI;
// //       thetamax = theta;
// //       
// //       llthred = llprev + log(randu());
// //       accept = false;
// //       count = 0;
// //       
// //       while(accept == false){
// //         count = count + 1;
// //         
// //         neweta.col(t) = eta.col(t) * cos(theta) + prioretat * sin(theta);
// //         neweta.col(t) = neweta.col(t) - mean(neweta.col(t));
// //         
// //         newlogZ = cppf_logZ_a0t_etat(Xzoop, alpha0, alpha, neweta);
// //         newlogZobl = cppf_logZobl(newlogZ, indobl);
// //         newlogZsur = cppf_logZsur(newlogZ, indsur);
// //         newloglam = cppf_loglam_b0t(Xwhale, newlogZ, beta, indtwhale);
// //         newlam = exp(newloglam);
// //         
// //         llnew = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //           cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) +
// //           cppf_loglik_obl(logYobl, newlogZobl, sig2_obl) +
// //           cppf_loglik_sur(logYsur, newlogZsur, lam0, lam1, sig2_sur);
// //         
// //         if(llnew > llthred){
// //           llprev = llnew;
// //           accept = true;
// //         } else {
// //           if(theta < 0){ thetamin = theta; } else { thetamax = theta; }
// //           theta = thetamin + randu() * (thetamax - thetamin);
// //           if ( (count) % 1000 == 0 ) {
// //             Rprintf("ESS for eta: %d iterations...\n", count);
// //           }
// //           if (count > 20000) stop("ESS is not converging...\n");
// //         }
// //       }
// //       eta = neweta;
// //       logZ = newlogZ;
// //       logZobl = newlogZobl;
// //       logZsur = newlogZsur;
// //       loglam = newloglam;
// //       lam = newlam;
// //       
// //       postEtat = as<mat>(postEta[t]);
// //       postEtat.row(s) = trans(eta.col(t));
// //       postEta[t] = postEtat;
// //     }
// //     if(s == 0){ Rprintf("sampled eta\n"); }
// //     
// //     
// //     // update kappa_eta
// //     kappa_eta = cppf_gibbs_kappa_eta(eta, invReta, 2, 2);
// //     postKappa(s,0) = kappa_eta;
// //     if(s == 0){ Rprintf("sampled kappa_eta\n"); }
// //     
// //     
// //     // update lam0
// //     lam0 = cppf_gibbs_lam0(logYsur, logZsur, lam1, sig2_sur, 100);
// //     postLam(s, 0) = lam0;
// //     if(s == 0){ Rprintf("sampled lam0\n"); }
// //     
// //     
// //     // update lam1
// //     lam1 = cppf_gibbs_lam1(logYsur, logZsur, lam0, sig2_sur, 100);
// //     postLam(s, 1) = lam1;
// //     if(s == 0){ Rprintf("sampled lam1\n"); }
// //     
// //     
// //     // update sig2_sur
// //     sig2_sur = cppf_gibbs_sig2_sur(logYsur, logZsur, lam0, lam1, 2, 2);
// //     postSig2(s, 0) = sig2_sur;
// //     if(s == 0){ Rprintf("sampled sig2_sur\n"); }
// //     
// //     
// //     // update beta
// //     if( updateCOV ){
// //       if( (s+1 >= adaptInterval) && (s+1 - (adaptInterval * trunc((s+1) / adaptInterval)) == 0) ){
// //         dummyaccprob = accprob.col(n_al0+1);
// //         rhat[n_al0+1] = sum( dummyaccprob.rows(s+1-adaptInterval, s-1) ) / (adaptInterval-1);
// //         gamma1[n_al0+1] = 1 / pow(adapIter[n_al0+1], c1);
// //         gamma2[n_al0+1] = c0 * gamma1[n_al0+1];
// //         sigma2[n_al0+1] = exp( log(sigma2[n_al0+1]) + gamma2[n_al0+1] * (rhat[n_al0+1] - ropt) );
// //         
// //         COVbeta = COVbeta + gamma1[n_al0+1] * ( cov( postBeta.rows(s+1-adaptInterval, s-1)) - COVbeta );
// //         cholCOVbeta = trans( chol( sigma2[n_al0+1] * ( COVbeta + 0.0000000001 * diagmat(ones(n_be)) ) ) );
// //         
// //         adapIter[n_al0+1] = adapIter[n_al0+1] + 1;
// //       }
// //     }
// //     
// //     newbeta = beta + cholCOVbeta * randn(n_be);
// //     
// //     newloglam = cppf_loglam_b0t(Xwhale, logZ, newbeta, indtwhale);
// //     newlam = exp(newloglam);
// //     
// //     logprob = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) -
// //       cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //       cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) -
// //       cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
// //       MVN_logh(newbeta, zeros(n_be), diagmat(ones(n_be))/100) -
// //       MVN_logh(beta, zeros(n_be), diagmat(ones(n_be))/100);
// //     
// //     if (log(randu()) < logprob) {
// //       
// //       beta = newbeta;
// //       loglam = newloglam;
// //       lam = newlam;
// //       
// //       accprob(s, n_al0+1) = 1;
// //     }
// //     
// //     postBeta.row(s) = trans(beta);
// //     if(s == 0){ Rprintf("sampled beta\n"); }
// //     
// //     
// //     if ( (s+1) % 10 == 0 ) {
// //       Rprintf("Generated %d samples...\n", s+1);
// //     }
// //   }
// //   
// //   
// //   return Rcpp::List::create(Rcpp::Named("postAlpha0") = postAlpha0,
// //                             Rcpp::Named("postAlpha0tilde") = postAlpha0tilde,
// //                             Rcpp::Named("postAlpha") = postAlpha,
// //                             Rcpp::Named("postEta") = postEta,
// //                             Rcpp::Named("postLam") = postLam,
// //                             Rcpp::Named("postKappa") = postKappa,
// //                             Rcpp::Named("postSig2") = postSig2,
// //                             Rcpp::Named("postBeta") = postBeta,
// //                             Rcpp::Named("Accprob") = accprob,
// //                             Rcpp::Named("sigma2") = sigma2,
// //                             Rcpp::Named("adapIter") = adapIter,
// //                             Rcpp::Named("COValpha0") = COValpha0,
// //                             Rcpp::Named("COValpha") = COValpha,
// //                             Rcpp::Named("COVbeta") = COVbeta);
// //   
// // }
// // 
// // 
// // // [[Rcpp::export]]
// // List cppf_fitJoint_b0t_psi_fixC(
// //     int niter, mat distmat,
// //     List Xzoop, vec logYobl, vec logYsur,
// //     List Xwhale, mat Sdist, mat Wpam, vec dur_dist,
// //     vec f_dist_obs, List f_dist, List p_pam,
// //     mat indobl, mat indsur, vec indtwhale, mat inddist, vec indmdist, vec indmwhale,
// //     double alpha0tilde, vec alpha0, vec alpha, mat eta, double lam0, double lam1,
// //     double tau2, double kappa_eta, double phi_eta, double sig2_obl, double sig2_sur,
// //     vec beta, vec psi, vec pii, vec c, 
// //     double kappa_psi, double phi_psi,
// //     double cellarea,
// //     vec sigma2, vec COValpha0, mat COValpha, mat COVbeta,
// //     bool updateCOV, int adaptInterval, double adaptFactorExponent, vec adapIter){
// //   
// //   double negativeInf = -std::numeric_limits<float>::infinity();;
// //   
// //   indobl = indobl-1, indsur = indsur-1, indtwhale = indtwhale-1;
// //   inddist = inddist-1, indmdist = indmdist-1, indmwhale = indmwhale-1;
// //   
// //   int nG = distmat.n_rows, n_al0 = alpha0.size(), n_al = alpha.size(), n_be = beta.size(), ntzoop = eta.n_cols, n_sigma2 = sigma2.size();
// //   
// //   mat postAlpha0tilde = zeros(niter, 1), postAlpha0 = zeros(niter, n_al0), postAlpha = zeros(niter, n_al), postLam = zeros(niter, 2);
// //   mat postSig2 = zeros(niter, 1), postKappa = zeros(niter, 2);
// //   mat postBeta = zeros(niter, n_be);
// //   mat postPsi(niter, nG), postEtat = zeros(niter, nG);
// //   List postEta(ntzoop);
// //   for(int i = 0; i < ntzoop; i ++){
// //     postEta[i] = zeros(niter, nG);
// //   }
// //   Rprintf("set up matrices for posterior samples\n");
// //   
// //   mat accprob = zeros(niter, n_sigma2);
// //   vec rhat = zeros(n_sigma2), gamma1 = zeros(n_sigma2), gamma2 = zeros(n_sigma2), dummyaccprob;
// //   double c0 = 1, c1 = adaptFactorExponent, ropt = 0.234, logprob;
// //   vec cholCOValpha0(n_al0);
// //   for(int i = 0; i < n_al0; i ++){ cholCOValpha0[i] = sqrt(sigma2[i] * COValpha0[i]); }
// //   mat cholCOValpha = trans( chol( sigma2[n_al0] * ( COValpha + 0.0000000001 * diagmat(ones(n_al)) ) ) );
// //   mat cholCOVbeta = trans( chol( sigma2[n_al0+1] * ( COVbeta + 0.0000000001 * diagmat(ones(n_be)) ) ) );
// //   Rprintf("set up adaptive proposals\n");
// //   
// //   double alpha0i, newalpha0i, beta0i, newbeta0i;
// //   vec newalpha0, newalpha, newbeta0, newbeta, newlogZobl, newlogZsur;
// //   mat newlogZ, newloglam, newlam;
// //   vec dummyvec;
// //   
// //   int count;
// //   vec xvals;
// //   double llprev, llthred, llnew, thetamin, thetamax, theta;
// //   vec prioretat, priorpsi, newpsi;
// //   mat neweta;
// //   bool accept;
// //   
// //   mat Reta = exp(-distmat / phi_eta), Rpsi = exp(-distmat / phi_psi);
// //   mat invReta = inv(Reta), invRpsi= inv(Rpsi), cholReta = chol(Reta), cholRpsi= chol(Rpsi);
// //   rowvec oneinvReta = trans(ones(nG)) * invReta, oneinvRpsi = trans(ones(nG)) * invRpsi;
// //   double oneinvRetaone = sum( oneinvReta * ones(nG) ), oneinvRpsione = sum( oneinvRpsi * ones(nG) );
// //   Rprintf("set up GP covariance matrices\n");
// //   
// //   mat logZ = cppf_logZ_a0t_etat(Xzoop, alpha0, alpha, eta);
// //   vec logZobl = cppf_logZobl(logZ, indobl);
// //   vec logZsur = cppf_logZsur(logZ, indsur);
// //   mat loglam = cppf_loglam_b0t_psi(Xwhale, logZ, beta, psi, indtwhale);
// //   mat lam = exp(loglam);
// //   Rprintf("set up initial values\n");
// //   
// //   
// //   // Start MCMC
// //   for(int s = 0; s < niter; s++) {
// //     
// //     // update alpha0
// //     for(int i = 0; i < n_al0; i ++){
// //       
// //       if( updateCOV ){
// //         if( (s+1 >= adaptInterval) && (s+1 - (adaptInterval * trunc((s+1) / adaptInterval)) == 0) ){
// //           dummyaccprob = accprob.col(i);
// //           rhat[i] = sum( dummyaccprob.rows(s+1-adaptInterval, s-1) ) / (adaptInterval-1);
// //           gamma1[i] = 1 / pow(adapIter[i], c1);
// //           gamma2[i] = c0 * gamma1[i];
// //           sigma2[i] = exp( log(sigma2[i]) + gamma2[i] * (rhat[i] - ropt) );
// //           
// //           
// //           dummyvec = postAlpha0.col(i);
// //           COValpha0[i] = COValpha0[i] + gamma1[i] * ( var( dummyvec.rows(s+1-adaptInterval, s-1) ) - COValpha0[i] );
// //           cholCOValpha0[i] = sqrt(  sigma2[i] * COValpha0[i] );
// //           
// //           adapIter[i] = adapIter[i] + 1;
// //         }
// //       }
// //       
// //       alpha0i = alpha0[i];
// //       newalpha0i = alpha0i + cholCOValpha0[i] * randn();
// //       newalpha0 = alpha0;
// //       newalpha0[i] = newalpha0i;
// //       
// //       newlogZ = cppf_logZ_a0t_etat(Xzoop, newalpha0, alpha, eta);
// //       newlogZobl = cppf_logZobl(newlogZ, indobl);
// //       newlogZsur = cppf_logZsur(newlogZ, indsur);
// //       newloglam = cppf_loglam_b0t(Xwhale, newlogZ, beta, indtwhale);
// //       newlam = exp(newloglam);
// //       
// //       logprob = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) -
// //         cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //         cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) -
// //         cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
// //         cppf_loglik_obl(logYobl, newlogZobl, sig2_obl) -
// //         cppf_loglik_obl(logYobl, logZobl, sig2_obl) +
// //         cppf_loglik_sur(logYsur, newlogZsur, lam0, lam1, sig2_sur) -
// //         cppf_loglik_sur(logYsur, logZsur, lam0, lam1, sig2_sur) +
// //         normal_logh(newalpha0i, alpha0tilde, tau2) -
// //         normal_logh(alpha0i, alpha0tilde, tau2);
// //       
// //       if (log(randu()) < logprob) {
// //         
// //         alpha0 = newalpha0;
// //         logZ = newlogZ;
// //         logZobl = newlogZobl;
// //         logZsur = newlogZsur;
// //         loglam = newloglam;
// //         lam = newlam;
// //         
// //         accprob(s, i) = 1;
// //       }
// //     }
// //     postAlpha0.row(s) = trans(alpha0);
// //     if(s == 0){ Rprintf("sampled alpha0\n"); }
// //     
// //     
// //     // Gibbs update for alpha0tilde
// //     alpha0tilde = cppf_gibbs_alpha0tilde(alpha0, tau2);
// //     postAlpha0tilde[s] = alpha0tilde;
// //     if(s == 0){ Rprintf("sampled alpha0tilde\n"); }
// //     
// //     
// //     // update alpha
// //     if( updateCOV ){
// //       if( (s+1 >= adaptInterval) && (s+1 - (adaptInterval * trunc((s+1) / adaptInterval)) == 0) ){
// //         dummyaccprob = accprob.col(n_al0);
// //         rhat[n_al0] = sum( dummyaccprob.rows(s+1-adaptInterval, s-1) ) / (adaptInterval-1);
// //         gamma1[n_al0] = 1 / pow(adapIter[n_al0], c1);
// //         gamma2[n_al0] = c0 * gamma1[n_al0];
// //         sigma2[n_al0] = exp( log(sigma2[n_al0]) + gamma2[n_al0] * (rhat[n_al0] - ropt) );
// //         
// //         COValpha = COValpha + gamma1[n_al0] * ( cov( postAlpha.rows(s+1-adaptInterval, s-1) ) - COValpha );
// //         cholCOValpha = trans( chol( sigma2[n_al0] * ( COValpha + 0.0000000001 * diagmat(ones(n_al)) ) ) );
// //         
// //         adapIter[n_al0] = adapIter[n_al0] + 1;
// //       }
// //     }
// //     
// //     newalpha = alpha + cholCOValpha * randn(n_al);
// //     
// //     newlogZ = cppf_logZ_a0t_etat(Xzoop, alpha0, newalpha, eta);
// //     newlogZobl = cppf_logZobl(newlogZ, indobl);
// //     newlogZsur = cppf_logZsur(newlogZ, indsur);
// //     newloglam = cppf_loglam_b0t_psi(Xwhale, newlogZ, beta, psi, indtwhale);
// //     newlam = exp(newloglam);
// //     
// //     logprob = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) -
// //       cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //       cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) -
// //       cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
// //       cppf_loglik_obl(logYobl, newlogZobl, sig2_obl) -
// //       cppf_loglik_obl(logYobl, logZobl, sig2_obl) +
// //       cppf_loglik_sur(logYsur, newlogZsur, lam0, lam1, sig2_sur) -
// //       cppf_loglik_sur(logYsur, logZsur, lam0, lam1, sig2_sur) +
// //       MVN_logh(newalpha, zeros(n_al), diagmat(ones(n_al))/100) -
// //       MVN_logh(alpha, zeros(n_al), diagmat(ones(n_al))/100);
// //     
// //     if (log(randu()) < logprob) {
// //       
// //       alpha = newalpha;
// //       logZ = newlogZ;
// //       logZobl = newlogZobl;
// //       logZsur = newlogZsur;
// //       loglam = newloglam;
// //       lam = newlam;
// //       
// //       accprob(s, n_al0) = 1;
// //     }
// //     
// //     postAlpha.row(s) = trans(alpha);
// //     if(s == 0){ Rprintf("sampled alpha\n"); }
// //     
// //     
// //     // update eta via Elliptical slice sampling
// //     for(int t = 0; t < ntzoop; t ++){
// //       neweta = eta;
// //       
// //       llprev = cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //         cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
// //         cppf_loglik_obl(logYobl, logZobl, sig2_obl) +
// //         cppf_loglik_sur(logYsur, logZsur, lam0, lam1, sig2_sur);
// //       prioretat = sqrt(kappa_eta) * trans(cholReta) * randn(nG);
// //       
// //       thetamin = 0;
// //       thetamax = 2 * M_PI;
// //       theta = thetamin + randu() * (thetamax - thetamin);
// //       thetamin = theta - 2 * M_PI;
// //       thetamax = theta;
// //       
// //       llthred = llprev + log(randu());
// //       accept = false;
// //       count = 0;
// //       
// //       while(accept == false){
// //         count = count + 1;
// //         
// //         neweta.col(t) = eta.col(t) * cos(theta) + prioretat * sin(theta);
// //         neweta.col(t) = neweta.col(t) - mean(neweta.col(t));
// //         
// //         newlogZ = cppf_logZ_a0t_etat(Xzoop, alpha0, alpha, neweta);
// //         newlogZobl = cppf_logZobl(newlogZ, indobl);
// //         newlogZsur = cppf_logZsur(newlogZ, indsur);
// //         newloglam = cppf_loglam_b0t_psi(Xwhale, newlogZ, beta, psi, indtwhale);
// //         newlam = exp(newloglam);
// //         
// //         llnew = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //           cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) +
// //           cppf_loglik_obl(logYobl, newlogZobl, sig2_obl) +
// //           cppf_loglik_sur(logYsur, newlogZsur, lam0, lam1, sig2_sur);
// //         
// //         if(llnew > llthred){
// //           llprev = llnew;
// //           accept = true;
// //         } else {
// //           if(theta < 0){ thetamin = theta; } else { thetamax = theta; }
// //           theta = thetamin + randu() * (thetamax - thetamin);
// //           if ( (count) % 1000 == 0 ) {
// //             Rprintf("ESS for eta: %d iterations...\n", count);
// //           }
// //           if (count > 20000) stop("ESS is not converging...\n");
// //         }
// //       }
// //       eta = neweta;
// //       logZ = newlogZ;
// //       logZobl = newlogZobl;
// //       logZsur = newlogZsur;
// //       loglam = newloglam;
// //       lam = newlam;
// //       
// //       postEtat = as<mat>(postEta[t]);
// //       postEtat.row(s) = trans(eta.col(t));
// //       postEta[t] = postEtat;
// //     }
// //     if(s == 0){ Rprintf("sampled eta\n"); }
// //     
// //     
// //     // update kappa_eta
// //     kappa_eta = cppf_gibbs_kappa_eta(eta, invReta, 2, 2);
// //     postKappa(s,0) = kappa_eta;
// //     if(s == 0){ Rprintf("sampled kappa_eta\n"); }
// //     
// //     
// //     // update lam0
// //     lam0 = cppf_gibbs_lam0(logYsur, logZsur, lam1, sig2_sur, 100);
// //     postLam(s, 0) = lam0;
// //     if(s == 0){ Rprintf("sampled lam0\n"); }
// //     
// //     
// //     // update lam1
// //     lam1 = cppf_gibbs_lam1(logYsur, logZsur, lam0, sig2_sur, 100);
// //     postLam(s, 1) = lam1;
// //     if(s == 0){ Rprintf("sampled lam1\n"); }
// //     
// //     
// //     // update sig2_sur
// //     sig2_sur = cppf_gibbs_sig2_sur(logYsur, logZsur, lam0, lam1, 2, 2);
// //     postSig2(s, 0) = sig2_sur;
// //     if(s == 0){ Rprintf("sampled sig2_sur\n"); }
// //     
// //     
// //     // update beta
// //     if( updateCOV ){
// //       if( (s+1 >= adaptInterval) && (s+1 - (adaptInterval * trunc((s+1) / adaptInterval)) == 0) ){
// //         dummyaccprob = accprob.col(n_al0+1);
// //         rhat[n_al0+1] = sum( dummyaccprob.rows(s+1-adaptInterval, s-1) ) / (adaptInterval-1);
// //         gamma1[n_al0+1] = 1 / pow(adapIter[n_al0+1], c1);
// //         gamma2[n_al0+1] = c0 * gamma1[n_al0+1];
// //         sigma2[n_al0+1] = exp( log(sigma2[n_al0+1]) + gamma2[n_al0+1] * (rhat[n_al0+1] - ropt) );
// //         
// //         COVbeta = COVbeta + gamma1[n_al0+1] * ( cov( postBeta.rows(s+1-adaptInterval, s-1)) - COVbeta );
// //         cholCOVbeta = trans( chol( sigma2[n_al0+1] * ( COVbeta + 0.0000000001 * diagmat(ones(n_be)) ) ) );
// //         
// //         adapIter[n_al0+1] = adapIter[n_al0+1] + 1;
// //       }
// //     }
// //     
// //     newbeta = beta + cholCOVbeta * randn(n_be);
// //     
// //     newloglam = cppf_loglam_b0t_psi(Xwhale, logZ, newbeta, psi, indtwhale);
// //     newlam = exp(newloglam);
// //     
// //     logprob = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) -
// //       cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //       cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) -
// //       cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
// //       MVN_logh(newbeta, zeros(n_be), diagmat(ones(n_be))/100) -
// //       MVN_logh(beta, zeros(n_be), diagmat(ones(n_be))/100);
// //     
// //     if (log(randu()) < logprob) {
// //       
// //       beta = newbeta;
// //       loglam = newloglam;
// //       lam = newlam;
// //       
// //       accprob(s, n_al0+1) = 1;
// //     }
// //     
// //     postBeta.row(s) = trans(beta);
// //     if(s == 0){ Rprintf("sampled beta\n"); }
// //     
// //     
// //     // update psi via elliptical slice sampling
// //     llprev = cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //       cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea);
// //     priorpsi = sqrt(kappa_psi) * trans(cholRpsi) * randn(nG);
// //     
// //     thetamin = 0;
// //     thetamax = 2 * M_PI;
// //     theta = thetamin + randu() * (thetamax - thetamin);
// //     thetamin = theta - 2 * M_PI;
// //     thetamax = theta;
// //     
// //     llthred = llprev + log(randu());
// //     accept = false;
// //     count = 0;
// //     
// //     while(accept == false){
// //       count = count + 1;
// //       
// //       newpsi = psi * cos(theta) + priorpsi * sin(theta);
// //       newpsi = newpsi - mean(newpsi);
// //       
// //       newloglam = cppf_loglam_b0t_psi(Xwhale, logZ, beta, newpsi, indtwhale);
// //       newlam = exp(newloglam);
// //       
// //       llnew = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //         cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea);
// //       
// //       if(llnew > llthred){
// //         llprev = llnew;
// //         accept = true;
// //       } else {
// //         if(theta < 0){ thetamin = theta; } else { thetamax = theta; }
// //         theta = thetamin + randu() * (thetamax - thetamin);
// //         if ( (count) % 1000 == 0 ) {
// //           Rprintf("ESS for psi: %d iterations...\n", count);
// //         }
// //         if (count > 20000) stop("ESS is not converging...\n");
// //       }
// //     }
// //     psi = newpsi;
// //     loglam = newloglam;
// //     lam = newlam;
// //     
// //     postPsi.row(s) = trans(psi);
// //     if(s == 0){ Rprintf("sampled psi\n"); }
// //     
// //     
// //     // update kappa_psi
// //     kappa_psi = cppf_gibbs_kappa_psi(psi, invRpsi, 2, 2);
// //     postKappa(s,1) = kappa_psi;
// //     if(s == 0){ Rprintf("sampled kappa_psi\n"); }
// //     
// //     
// //     if ( (s+1) % 10 == 0 ) {
// //       Rprintf("Generated %d samples...\n", s+1);
// //     }
// //   }
// //   
// //   
// //   return Rcpp::List::create(Rcpp::Named("postAlpha0") = postAlpha0,
// //                             Rcpp::Named("postAlpha0tilde") = postAlpha0tilde,
// //                             Rcpp::Named("postAlpha") = postAlpha,
// //                             Rcpp::Named("postEta") = postEta,
// //                             Rcpp::Named("postLam") = postLam,
// //                             Rcpp::Named("postKappa") = postKappa,
// //                             Rcpp::Named("postSig2") = postSig2,
// //                             Rcpp::Named("postBeta") = postBeta,
// //                             Rcpp::Named("postPsi") = postPsi,
// //                             Rcpp::Named("Accprob") = accprob,
// //                             Rcpp::Named("sigma2") = sigma2,
// //                             Rcpp::Named("adapIter") = adapIter,
// //                             Rcpp::Named("COValpha0") = COValpha0,
// //                             Rcpp::Named("COValpha") = COValpha,
// //                             Rcpp::Named("COVbeta") = COVbeta);
// //   
// // }
// // 
// // 
// // 
// // // [[Rcpp::export]]
// // List cppf_fitJoint_b0t_psit_fixC(
// //     int niter, mat distmat,
// //     List Xzoop, vec logYobl, vec logYsur,
// //     List Xwhale, mat Sdist, mat Wpam, vec dur_dist,
// //     vec f_dist_obs, List f_dist, List p_pam,
// //     mat indobl, mat indsur, vec indtwhale, mat inddist, vec indmdist, vec indmwhale,
// //     double alpha0tilde, vec alpha0, vec alpha, mat eta, double lam0, double lam1,
// //     double tau2, double kappa_eta, double phi_eta, double sig2_obl, double sig2_sur,
// //     vec beta, mat psi, vec pii, vec c,
// //     double kappa_psi, double phi_psi,
// //     double cellarea,
// //     vec sigma2, vec COValpha0, mat COValpha, mat COVbeta,
// //     bool updateCOV, int adaptInterval, double adaptFactorExponent, vec adapIter){
// //   
// //   double negativeInf = -std::numeric_limits<float>::infinity();;
// //   
// //   indobl = indobl-1, indsur = indsur-1, indtwhale = indtwhale-1;
// //   inddist = inddist-1, indmdist = indmdist-1, indmwhale = indmwhale-1;
// //   
// //   int nG = distmat.n_rows, n_al0 = alpha0.size(), n_al = alpha.size(), n_be = beta.size(), ntzoop = eta.n_cols, ntwhale = psi.n_cols, n_sigma2 = sigma2.size();
// //   
// //   mat postAlpha0tilde = zeros(niter, 1), postAlpha0 = zeros(niter, n_al0), postAlpha = zeros(niter, n_al), postLam = zeros(niter, 2);
// //   mat postSig2 = zeros(niter, 1), postKappa = zeros(niter, 2);
// //   mat postBeta = zeros(niter, n_be);
// //   mat postEtat = zeros(niter, nG), postPsit = zeros(niter, nG);
// //   List postEta(ntzoop), postPsi(ntwhale);
// //   for(int i = 0; i < ntzoop; i ++){ postEta[i] = zeros(niter, nG); }
// //   for(int i = 0; i < ntwhale; i ++){ postPsi[i] = zeros(niter, nG); }
// //   Rprintf("set up matrices for posterior samples\n");
// //   
// //   mat accprob = zeros(niter, n_sigma2);
// //   vec rhat = zeros(n_sigma2), gamma1 = zeros(n_sigma2), gamma2 = zeros(n_sigma2), dummyaccprob;
// //   double c0 = 1, c1 = adaptFactorExponent, ropt = 0.234, logprob;
// //   vec cholCOValpha0(n_al0);
// //   for(int i = 0; i < n_al0; i ++){ cholCOValpha0[i] = sqrt(sigma2[i] * COValpha0[i]); }
// //   mat cholCOValpha = trans( chol( sigma2[n_al0] * ( COValpha + 0.0000000001 * diagmat(ones(n_al)) ) ) );
// //   mat cholCOVbeta = trans( chol( sigma2[n_al0+1] * ( COVbeta + 0.0000000001 * diagmat(ones(n_be)) ) ) );
// //   Rprintf("set up adaptive proposals\n");
// //   
// //   double alpha0i, newalpha0i, beta0i, newbeta0i;
// //   vec newalpha0, newalpha, newbeta0, newbeta, newlogZobl, newlogZsur;
// //   mat newlogZ, newloglam, newlam;
// //   vec dummyvec;
// //   
// //   int count;
// //   vec xvals;
// //   double llprev, llthred, llnew, thetamin, thetamax, theta;
// //   vec prioretat, priorpsit;
// //   mat neweta, newpsi;
// //   bool accept;
// //   
// //   mat Reta = exp(-distmat / phi_eta), Rpsi = exp(-distmat / phi_psi);
// //   mat invReta = inv(Reta), invRpsi= inv(Rpsi), cholReta = chol(Reta), cholRpsi= chol(Rpsi);
// //   rowvec oneinvReta = trans(ones(nG)) * invReta, oneinvRpsi = trans(ones(nG)) * invRpsi;
// //   double oneinvRetaone = sum( oneinvReta * ones(nG) ), oneinvRpsione = sum( oneinvRpsi * ones(nG) );
// //   Rprintf("set up GP covariance matrices\n");
// //   
// //   mat logZ = cppf_logZ_a0t_etat(Xzoop, alpha0, alpha, eta);
// //   vec logZobl = cppf_logZobl(logZ, indobl);
// //   vec logZsur = cppf_logZsur(logZ, indsur);
// //   mat loglam = cppf_loglam_b0t_psit(Xwhale, logZ, beta, psi, indtwhale);
// //   mat lam = exp(loglam);
// //   Rprintf("set up initial values\n");
// //   
// //   
// //   // Start MCMC
// //   for(int s = 0; s < niter; s++) {
// //     
// //     // update alpha0
// //     for(int i = 0; i < n_al0; i ++){
// //       
// //       if( updateCOV ){
// //         if( (s+1 >= adaptInterval) && (s+1 - (adaptInterval * trunc((s+1) / adaptInterval)) == 0) ){
// //           dummyaccprob = accprob.col(i);
// //           rhat[i] = sum( dummyaccprob.rows(s+1-adaptInterval, s-1) ) / (adaptInterval-1);
// //           gamma1[i] = 1 / pow(adapIter[i], c1);
// //           gamma2[i] = c0 * gamma1[i];
// //           sigma2[i] = exp( log(sigma2[i]) + gamma2[i] * (rhat[i] - ropt) );
// //           
// //           
// //           dummyvec = postAlpha0.col(i);
// //           COValpha0[i] = COValpha0[i] + gamma1[i] * ( var( dummyvec.rows(s+1-adaptInterval, s-1) ) - COValpha0[i] );
// //           cholCOValpha0[i] = sqrt(  sigma2[i] * COValpha0[i] );
// //           
// //           adapIter[i] = adapIter[i] + 1;
// //         }
// //       }
// //       
// //       alpha0i = alpha0[i];
// //       newalpha0i = alpha0i + cholCOValpha0[i] * randn();
// //       newalpha0 = alpha0;
// //       newalpha0[i] = newalpha0i;
// //       
// //       newlogZ = cppf_logZ_a0t_etat(Xzoop, newalpha0, alpha, eta);
// //       newlogZobl = cppf_logZobl(newlogZ, indobl);
// //       newlogZsur = cppf_logZsur(newlogZ, indsur);
// //       newloglam = cppf_loglam_b0t(Xwhale, newlogZ, beta, indtwhale);
// //       newlam = exp(newloglam);
// //       
// //       logprob = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) -
// //         cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //         cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) -
// //         cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
// //         cppf_loglik_obl(logYobl, newlogZobl, sig2_obl) -
// //         cppf_loglik_obl(logYobl, logZobl, sig2_obl) +
// //         cppf_loglik_sur(logYsur, newlogZsur, lam0, lam1, sig2_sur) -
// //         cppf_loglik_sur(logYsur, logZsur, lam0, lam1, sig2_sur) +
// //         normal_logh(newalpha0i, alpha0tilde, tau2) -
// //         normal_logh(alpha0i, alpha0tilde, tau2);
// //       
// //       if (log(randu()) < logprob) {
// //         
// //         alpha0 = newalpha0;
// //         logZ = newlogZ;
// //         logZobl = newlogZobl;
// //         logZsur = newlogZsur;
// //         loglam = newloglam;
// //         lam = newlam;
// //         
// //         accprob(s, i) = 1;
// //       }
// //     }
// //     postAlpha0.row(s) = trans(alpha0);
// //     if(s == 0){ Rprintf("sampled alpha0\n"); }
// //     
// //     
// //     // Gibbs update for alpha0tilde
// //     alpha0tilde = cppf_gibbs_alpha0tilde(alpha0, tau2);
// //     postAlpha0tilde[s] = alpha0tilde;
// //     if(s == 0){ Rprintf("sampled alpha0tilde\n"); }
// //     
// //     
// //     // update alpha
// //     if( updateCOV ){
// //       if( (s+1 >= adaptInterval) && (s+1 - (adaptInterval * trunc((s+1) / adaptInterval)) == 0) ){
// //         dummyaccprob = accprob.col(n_al0);
// //         rhat[n_al0] = sum( dummyaccprob.rows(s+1-adaptInterval, s-1) ) / (adaptInterval-1);
// //         gamma1[n_al0] = 1 / pow(adapIter[n_al0], c1);
// //         gamma2[n_al0] = c0 * gamma1[n_al0];
// //         sigma2[n_al0] = exp( log(sigma2[n_al0]) + gamma2[n_al0] * (rhat[n_al0] - ropt) );
// //         
// //         COValpha = COValpha + gamma1[n_al0] * ( cov( postAlpha.rows(s+1-adaptInterval, s-1) ) - COValpha );
// //         cholCOValpha = trans( chol( sigma2[n_al0] * ( COValpha + 0.0000000001 * diagmat(ones(n_al)) ) ) );
// //         
// //         adapIter[n_al0] = adapIter[n_al0] + 1;
// //       }
// //     }
// //     
// //     newalpha = alpha + cholCOValpha * randn(n_al);
// //     
// //     newlogZ = cppf_logZ_a0t_etat(Xzoop, alpha0, newalpha, eta);
// //     newlogZobl = cppf_logZobl(newlogZ, indobl);
// //     newlogZsur = cppf_logZsur(newlogZ, indsur);
// //     newloglam = cppf_loglam_b0t_psit(Xwhale, newlogZ, beta, psi, indtwhale);
// //     newlam = exp(newloglam);
// //     
// //     logprob = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) -
// //       cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //       cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) -
// //       cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
// //       cppf_loglik_obl(logYobl, newlogZobl, sig2_obl) -
// //       cppf_loglik_obl(logYobl, logZobl, sig2_obl) +
// //       cppf_loglik_sur(logYsur, newlogZsur, lam0, lam1, sig2_sur) -
// //       cppf_loglik_sur(logYsur, logZsur, lam0, lam1, sig2_sur) +
// //       MVN_logh(newalpha, zeros(n_al), diagmat(ones(n_al))/100) -
// //       MVN_logh(alpha, zeros(n_al), diagmat(ones(n_al))/100);
// //     
// //     if (log(randu()) < logprob) {
// //       
// //       alpha = newalpha;
// //       logZ = newlogZ;
// //       logZobl = newlogZobl;
// //       logZsur = newlogZsur;
// //       loglam = newloglam;
// //       lam = newlam;
// //       
// //       accprob(s, n_al0) = 1;
// //     }
// //     
// //     postAlpha.row(s) = trans(alpha);
// //     if(s == 0){ Rprintf("sampled alpha\n"); }
// //     
// //     
// //     // update eta via Elliptical slice sampling
// //     for(int t = 0; t < ntzoop; t ++){
// //       neweta = eta;
// //       
// //       llprev = cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //         cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
// //         cppf_loglik_obl(logYobl, logZobl, sig2_obl) +
// //         cppf_loglik_sur(logYsur, logZsur, lam0, lam1, sig2_sur);
// //       prioretat = sqrt(kappa_eta) * trans(cholReta) * randn(nG);
// //       
// //       thetamin = 0;
// //       thetamax = 2 * M_PI;
// //       theta = thetamin + randu() * (thetamax - thetamin);
// //       thetamin = theta - 2 * M_PI;
// //       thetamax = theta;
// //       
// //       llthred = llprev + log(randu());
// //       accept = false;
// //       count = 0;
// //       
// //       while(accept == false){
// //         count = count + 1;
// //         
// //         neweta.col(t) = eta.col(t) * cos(theta) + prioretat * sin(theta);
// //         neweta.col(t) = neweta.col(t) - mean(neweta.col(t));
// //         
// //         newlogZ = cppf_logZ_a0t_etat(Xzoop, alpha0, alpha, neweta);
// //         newlogZobl = cppf_logZobl(newlogZ, indobl);
// //         newlogZsur = cppf_logZsur(newlogZ, indsur);
// //         newloglam = cppf_loglam_b0t_psit(Xwhale, newlogZ, beta, psi, indtwhale);
// //         newlam = exp(newloglam);
// //         
// //         llnew = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //           cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) +
// //           cppf_loglik_obl(logYobl, newlogZobl, sig2_obl) +
// //           cppf_loglik_sur(logYsur, newlogZsur, lam0, lam1, sig2_sur);
// //         
// //         if(llnew > llthred){
// //           llprev = llnew;
// //           accept = true;
// //         } else {
// //           if(theta < 0){ thetamin = theta; } else { thetamax = theta; }
// //           theta = thetamin + randu() * (thetamax - thetamin);
// //           if ( (count) % 1000 == 0 ) {
// //             Rprintf("ESS for eta: %d iterations...\n", count);
// //           }
// //           if (count > 20000) stop("ESS is not converging...\n");
// //         }
// //       }
// //       eta = neweta;
// //       logZ = newlogZ;
// //       logZobl = newlogZobl;
// //       logZsur = newlogZsur;
// //       loglam = newloglam;
// //       lam = newlam;
// //       
// //       postEtat = as<mat>(postEta[t]);
// //       postEtat.row(s) = trans(eta.col(t));
// //       postEta[t] = postEtat;
// //     }
// //     if(s == 0){ Rprintf("sampled eta\n"); }
// //     
// //     
// //     // update kappa_eta
// //     kappa_eta = cppf_gibbs_kappa_eta(eta, invReta, 2, 2);
// //     postKappa(s,0) = kappa_eta;
// //     if(s == 0){ Rprintf("sampled kappa_eta\n"); }
// //     
// //     
// //     // update lam0
// //     lam0 = cppf_gibbs_lam0(logYsur, logZsur, lam1, sig2_sur, 100);
// //     postLam(s, 0) = lam0;
// //     if(s == 0){ Rprintf("sampled lam0\n"); }
// //     
// //     
// //     // update lam1
// //     lam1 = cppf_gibbs_lam1(logYsur, logZsur, lam0, sig2_sur, 100);
// //     postLam(s, 1) = lam1;
// //     if(s == 0){ Rprintf("sampled lam1\n"); }
// //     
// //     
// //     // update sig2_sur
// //     sig2_sur = cppf_gibbs_sig2_sur(logYsur, logZsur, lam0, lam1, 2, 2);
// //     postSig2(s, 0) = sig2_sur;
// //     if(s == 0){ Rprintf("sampled sig2_sur\n"); }
// //     
// //     
// //     // update beta
// //     if( updateCOV ){
// //       if( (s+1 >= adaptInterval) && (s+1 - (adaptInterval * trunc((s+1) / adaptInterval)) == 0) ){
// //         dummyaccprob = accprob.col(n_al0+1);
// //         rhat[n_al0+1] = sum( dummyaccprob.rows(s+1-adaptInterval, s-1) ) / (adaptInterval-1);
// //         gamma1[n_al0+1] = 1 / pow(adapIter[n_al0+1], c1);
// //         gamma2[n_al0+1] = c0 * gamma1[n_al0+1];
// //         sigma2[n_al0+1] = exp( log(sigma2[n_al0+1]) + gamma2[n_al0+1] * (rhat[n_al0+1] - ropt) );
// //         
// //         COVbeta = COVbeta + gamma1[n_al0+1] * ( cov( postBeta.rows(s+1-adaptInterval, s-1)) - COVbeta );
// //         cholCOVbeta = trans( chol( sigma2[n_al0+1] * ( COVbeta + 0.0000000001 * diagmat(ones(n_be)) ) ) );
// //         
// //         adapIter[n_al0+1] = adapIter[n_al0+1] + 1;
// //       }
// //     }
// //     
// //     newbeta = beta + cholCOVbeta * randn(n_be);
// //     
// //     newloglam = cppf_loglam_b0t_psit(Xwhale, logZ, newbeta, psi, indtwhale);
// //     newlam = exp(newloglam);
// //     
// //     logprob = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) -
// //       cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //       cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) -
// //       cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
// //       MVN_logh(newbeta, zeros(n_be), diagmat(ones(n_be))/100) -
// //       MVN_logh(beta, zeros(n_be), diagmat(ones(n_be))/100);
// //     
// //     if (log(randu()) < logprob) {
// //       
// //       beta = newbeta;
// //       loglam = newloglam;
// //       lam = newlam;
// //       
// //       accprob(s, n_al0+1) = 1;
// //     }
// //     
// //     postBeta.row(s) = trans(beta);
// //     if(s == 0){ Rprintf("sampled beta\n"); }
// //     
// //     
// //     
// //     // update psi via Elliptical slice sampling
// //     for(int t = 0; t < ntwhale; t ++){
// //       newpsi = psi;
// //       
// //       llprev = cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //         cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea);
// //       priorpsit = sqrt(kappa_psi) * trans(cholRpsi) * randn(nG);
// //       
// //       thetamin = 0;
// //       thetamax = 2 * M_PI;
// //       theta = thetamin + randu() * (thetamax - thetamin);
// //       thetamin = theta - 2 * M_PI;
// //       thetamax = theta;
// //       
// //       llthred = llprev + log(randu());
// //       accept = false;
// //       count = 0;
// //       
// //       while(accept == false){
// //         count = count + 1;
// //         
// //         newpsi.col(t) = psi.col(t) * cos(theta) + priorpsit * sin(theta);
// //         newpsi.col(t) = newpsi.col(t) - mean(newpsi.col(t));
// //         
// //         newloglam = cppf_loglam_b0t_psit(Xwhale, logZ, beta, newpsi, indtwhale);
// //         newlam = exp(newloglam);
// //         
// //         llnew = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //           cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea);
// //         
// //         if(llnew > llthred){
// //           llprev = llnew;
// //           accept = true;
// //         } else {
// //           if(theta < 0){ thetamin = theta; } else { thetamax = theta; }
// //           theta = thetamin + randu() * (thetamax - thetamin);
// //           if ( (count) % 1000 == 0 ) {
// //             Rprintf("ESS for psi: %d iterations...\n", count);
// //           }
// //           if (count > 20000) stop("ESS is not converging...\n");
// //         }
// //       }
// //       psi = newpsi;
// //       loglam = newloglam;
// //       lam = newlam;
// //       
// //       postPsit = as<mat>(postPsi[t]);
// //       postPsit.row(s) = trans(psi.col(t));
// //       postPsi[t] = postPsit;
// //     }
// //     if(s == 0){ Rprintf("sampled psi\n"); }
// //     
// //     
// //     // update kappa_psi
// //     kappa_psi = cppf_gibbs_kappa_psit(psi, invRpsi, 2, 2);
// //     postKappa(s,1) = kappa_psi;
// //     if(s == 0){ Rprintf("sampled kappa_psi\n"); }
// //     
// //     
// //     if ( (s+1) % 10 == 0 ) {
// //       Rprintf("Generated %d samples...\n", s+1);
// //     }
// //   }
// //   
// //   
// //   return Rcpp::List::create(Rcpp::Named("postAlpha0") = postAlpha0,
// //                             Rcpp::Named("postAlpha0tilde") = postAlpha0tilde,
// //                             Rcpp::Named("postAlpha") = postAlpha,
// //                             Rcpp::Named("postEta") = postEta,
// //                             Rcpp::Named("postLam") = postLam,
// //                             Rcpp::Named("postKappa") = postKappa,
// //                             Rcpp::Named("postSig2") = postSig2,
// //                             Rcpp::Named("postBeta") = postBeta,
// //                             Rcpp::Named("postPsi") = postPsi,
// //                             Rcpp::Named("Accprob") = accprob,
// //                             Rcpp::Named("sigma2") = sigma2,
// //                             Rcpp::Named("adapIter") = adapIter,
// //                             Rcpp::Named("COValpha0") = COValpha0,
// //                             Rcpp::Named("COValpha") = COValpha,
// //                             Rcpp::Named("COVbeta") = COVbeta);
// //   
// // }
// // 
// // 
// // // [[Rcpp::export]]
// // List cppf_fitSeparJoint_b0t_fixC(
// //     int niter, mat distmat,
// //     List Xzoop, vec logYobl, vec logYsur,
// //     List Xwhale, mat Sdist, mat Wpam, vec dur_dist,
// //     vec f_dist_obs, List f_dist, List p_pam,
// //     mat indobl, mat indsur, vec indtwhale, mat inddist, vec indmdist, vec indmwhale,
// //     double alpha0tilde, vec alpha0, vec alpha, mat eta, double lam0, double lam1,
// //     double tau2, double kappa_eta, double phi_eta, double sig2_obl, double sig2_sur, 
// //     vec beta, vec pii, vec c,
// //     double cellarea,
// //     vec sigma2, vec COValpha0, mat COValpha, mat COVbeta,
// //     bool updateCOV, int adaptInterval, double adaptFactorExponent, vec adapIter,
// //     int curniter, int joiniter){
// //   
// //   double negativeInf = -std::numeric_limits<float>::infinity();;
// //   
// //   indobl = indobl-1, indsur = indsur-1, indtwhale = indtwhale-1;
// //   inddist = inddist-1, indmdist = indmdist-1, indmwhale = indmwhale-1;
// //   
// //   int nG = distmat.n_rows, n_al0 = alpha0.size(), n_al = alpha.size(), n_be = beta.size(), ntzoop = eta.n_cols, n_sigma2 = sigma2.size();
// //   vec dummyvec;
// //   
// //   mat postAlpha0tilde = zeros(niter, 1), postAlpha0 = zeros(niter, n_al0), postAlpha = zeros(niter, n_al), postLam = zeros(niter, 2);
// //   mat postSig2 = zeros(niter, 1), postKappa = zeros(niter, 1);
// //   mat postBeta = zeros(niter, n_be);
// //   mat postEtat = zeros(niter, nG);
// //   List postEta(ntzoop);
// //   for(int i = 0; i < ntzoop; i ++){
// //     postEta[i] = zeros(niter, nG);
// //   }
// //   Rprintf("set up matrices for posterior samples\n");
// //   
// //   mat accprob = zeros(niter, n_sigma2);
// //   vec rhat = zeros(n_sigma2), gamma1 = zeros(n_sigma2), gamma2 = zeros(n_sigma2), dummyaccprob;
// //   double c0 = 1, c1 = adaptFactorExponent, ropt = 0.234, logprob;
// //   vec cholCOValpha0(n_al0);
// //   for(int i = 0; i < n_al0; i ++){ cholCOValpha0[i] = sqrt(sigma2[i] * COValpha0[i]); }
// //   mat cholCOValpha = trans( chol( sigma2[n_al0] * ( COValpha + 0.0000000001 * diagmat(ones(n_al)) ) ) );
// //   mat cholCOVbeta = trans( chol( sigma2[n_al0+1] * ( COVbeta + 0.0000000001 * diagmat(ones(n_be)) ) ) );
// //   Rprintf("set up adaptive proposals\n");
// //   
// //   double alpha0i, newalpha0i;
// //   vec newalpha0, newalpha, newbeta, newlogZobl, newlogZsur;
// //   mat newlogZ, newloglam, newlam;
// //   
// //   int count;
// //   vec xvals;
// //   double llprev, llthred, llnew, thetamin, thetamax, theta;
// //   vec prioretat;
// //   mat neweta;
// //   bool accept;
// //   
// //   mat Reta = exp(-distmat / phi_eta);
// //   mat invReta = inv(Reta),  cholReta = chol(Reta);
// //   rowvec oneinvReta = trans(ones(nG)) * invReta;
// //   double oneinvRetaone = sum( oneinvReta * ones(nG) );
// //   Rprintf("set up GP covariance matrices\n");
// //   
// //   mat logZ = cppf_logZ_a0t_etat(Xzoop, alpha0, alpha, eta);
// //   vec logZobl = cppf_logZobl(logZ, indobl);
// //   vec logZsur = cppf_logZsur(logZ, indsur);
// //   mat loglam = cppf_loglam_b0t(Xwhale, logZ, beta, indtwhale);
// //   mat lam = exp(loglam);
// //   Rprintf("set up initial values\n");
// //   
// //   
// //   // Start MCMC
// //   for(int s = 0; s < niter; s++) {
// //     
// //     // update alpha0
// //     for(int i = 0; i < n_al0; i ++){
// //       
// //       if( updateCOV ){
// //         if( (s+1 >= adaptInterval) && (s+1 - (adaptInterval * trunc((s+1) / adaptInterval)) == 0) ){
// //           dummyaccprob = accprob.col(i);
// //           rhat[i] = sum( dummyaccprob.rows(s+1-adaptInterval, s-1) ) / (adaptInterval-1);
// //           gamma1[i] = 1 / pow(adapIter[i], c1);
// //           gamma2[i] = c0 * gamma1[i];
// //           sigma2[i] = exp( log(sigma2[i]) + gamma2[i] * (rhat[i] - ropt) );
// //           
// //           
// //           dummyvec = postAlpha0.col(i);
// //           COValpha0[i] = COValpha0[i] + gamma1[i] * ( var( dummyvec.rows(s+1-adaptInterval, s-1) ) - COValpha0[i] );
// //           cholCOValpha0[i] = sqrt(  sigma2[i] * COValpha0[i] );
// //           
// //           adapIter[i] = adapIter[i] + 1;
// //         }
// //       }
// //       
// //       alpha0i = alpha0[i];
// //       newalpha0i = alpha0i + cholCOValpha0[i] * randn();
// //       newalpha0 = alpha0;
// //       newalpha0[i] = newalpha0i;
// //       
// //       newlogZ = cppf_logZ_a0t_etat(Xzoop, newalpha0, alpha, eta);
// //       newlogZobl = cppf_logZobl(newlogZ, indobl);
// //       newlogZsur = cppf_logZsur(newlogZ, indsur);
// //       newloglam = cppf_loglam_b0t(Xwhale, newlogZ, beta, indtwhale);
// //       newlam = exp(newloglam);
// //       
// //       if(curniter > joiniter){
// //         logprob = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) -
// //           cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //           cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) -
// //           cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
// //           cppf_loglik_obl(logYobl, newlogZobl, sig2_obl) -
// //           cppf_loglik_obl(logYobl, logZobl, sig2_obl) +
// //           cppf_loglik_sur(logYsur, newlogZsur, lam0, lam1, sig2_sur) -
// //           cppf_loglik_sur(logYsur, logZsur, lam0, lam1, sig2_sur) +
// //           normal_logh(newalpha0i, alpha0tilde, tau2) -
// //           normal_logh(alpha0i, alpha0tilde, tau2);  
// //         
// //       } else {
// //         logprob = cppf_loglik_obl(logYobl, newlogZobl, sig2_obl) -
// //           cppf_loglik_obl(logYobl, logZobl, sig2_obl) +
// //           cppf_loglik_sur(logYsur, newlogZsur, lam0, lam1, sig2_sur) -
// //           cppf_loglik_sur(logYsur, logZsur, lam0, lam1, sig2_sur) +
// //           normal_logh(newalpha0i, alpha0tilde, tau2) -
// //           normal_logh(alpha0i, alpha0tilde, tau2);
// //       }
// //       
// //       if (log(randu()) < logprob) {
// //         
// //         alpha0 = newalpha0;
// //         logZ = newlogZ;
// //         logZobl = newlogZobl;
// //         logZsur = newlogZsur;
// //         loglam = newloglam;
// //         lam = newlam;
// //         
// //         accprob(s, i) = 1;
// //       }
// //     }
// //     postAlpha0.row(s) = trans(alpha0);
// //     if(s == 0){ Rprintf("sampled alpha0\n"); }
// //     
// //     
// //     // Gibbs update for alpha0tilde
// //     alpha0tilde = cppf_gibbs_alpha0tilde(alpha0, tau2);
// //     postAlpha0tilde[s] = alpha0tilde;
// //     if(s == 0){ Rprintf("sampled alpha0tilde\n"); }
// //     
// //     
// //     // update alpha
// //     if( updateCOV ){
// //       if( (s+1 >= adaptInterval) && (s+1 - (adaptInterval * trunc((s+1) / adaptInterval)) == 0) ){
// //         dummyaccprob = accprob.col(n_al0);
// //         rhat[n_al0] = sum( dummyaccprob.rows(s+1-adaptInterval, s-1) ) / (adaptInterval-1);
// //         gamma1[n_al0] = 1 / pow(adapIter[n_al0], c1);
// //         gamma2[n_al0] = c0 * gamma1[n_al0];
// //         sigma2[n_al0] = exp( log(sigma2[n_al0]) + gamma2[n_al0] * (rhat[n_al0] - ropt) );
// //         
// //         COValpha = COValpha + gamma1[n_al0] * ( cov( postAlpha.rows(s+1-adaptInterval, s-1) ) - COValpha );
// //         cholCOValpha = trans( chol( sigma2[n_al0] * ( COValpha + 0.0000000001 * diagmat(ones(n_al)) ) ) );
// //         
// //         adapIter[n_al0] = adapIter[n_al0] + 1;
// //       }
// //     }
// //     
// //     newalpha = alpha + cholCOValpha * randn(n_al);
// //     
// //     newlogZ = cppf_logZ_a0t_etat(Xzoop, alpha0, newalpha, eta);
// //     newlogZobl = cppf_logZobl(newlogZ, indobl);
// //     newlogZsur = cppf_logZsur(newlogZ, indsur);
// //     newloglam = cppf_loglam_b0t(Xwhale, newlogZ, beta, indtwhale);
// //     newlam = exp(newloglam);
// //     
// //     if(curniter > joiniter){
// //       logprob = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) -
// //         cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //         cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) -
// //         cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
// //         cppf_loglik_obl(logYobl, newlogZobl, sig2_obl) -
// //         cppf_loglik_obl(logYobl, logZobl, sig2_obl) +
// //         cppf_loglik_sur(logYsur, newlogZsur, lam0, lam1, sig2_sur) -
// //         cppf_loglik_sur(logYsur, logZsur, lam0, lam1, sig2_sur) +
// //         MVN_logh(newalpha, zeros(n_al), diagmat(ones(n_al))/100) -
// //         MVN_logh(alpha, zeros(n_al), diagmat(ones(n_al))/100);
// //       
// //     } else {
// //       logprob = cppf_loglik_obl(logYobl, newlogZobl, sig2_obl) -
// //         cppf_loglik_obl(logYobl, logZobl, sig2_obl) +
// //         cppf_loglik_sur(logYsur, newlogZsur, lam0, lam1, sig2_sur) -
// //         cppf_loglik_sur(logYsur, logZsur, lam0, lam1, sig2_sur) +
// //         MVN_logh(newalpha, zeros(n_al), diagmat(ones(n_al))/100) -
// //         MVN_logh(alpha, zeros(n_al), diagmat(ones(n_al))/100);
// //     }
// //     
// //     if (log(randu()) < logprob) {
// //       
// //       alpha = newalpha;
// //       logZ = newlogZ;
// //       logZobl = newlogZobl;
// //       logZsur = newlogZsur;
// //       loglam = newloglam;
// //       lam = newlam;
// //       
// //       accprob(s, n_al0) = 1;
// //     }
// //     
// //     postAlpha.row(s) = trans(alpha);
// //     if(s == 0){ Rprintf("sampled alpha\n"); }
// //     
// //     
// //     // update eta via Elliptical slice sampling
// //     for(int t = 0; t < ntzoop; t ++){
// //       neweta = eta;
// //       
// //       if(curniter > joiniter){
// //         llprev = cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //           cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
// //           cppf_loglik_obl(logYobl, logZobl, sig2_obl) +
// //           cppf_loglik_sur(logYsur, logZsur, lam0, lam1, sig2_sur);
// //         prioretat = sqrt(kappa_eta) * trans(cholReta) * randn(nG);
// //         
// //       } else {
// //         llprev = cppf_loglik_obl(logYobl, logZobl, sig2_obl) +
// //           cppf_loglik_sur(logYsur, logZsur, lam0, lam1, sig2_sur);
// //         prioretat = sqrt(kappa_eta) * trans(cholReta) * randn(nG);
// //       }
// //       
// //       thetamin = 0;
// //       thetamax = 2 * M_PI;
// //       theta = thetamin + randu() * (thetamax - thetamin);
// //       thetamin = theta - 2 * M_PI;
// //       thetamax = theta;
// //       
// //       llthred = llprev + log(randu());
// //       accept = false;
// //       count = 0;
// //       
// //       while(accept == false){
// //         count = count + 1;
// //         
// //         neweta.col(t) = eta.col(t) * cos(theta) + prioretat * sin(theta);
// //         neweta.col(t) = neweta.col(t) - mean(neweta.col(t));
// //         
// //         newlogZ = cppf_logZ_a0t_etat(Xzoop, alpha0, alpha, neweta);
// //         newlogZobl = cppf_logZobl(newlogZ, indobl);
// //         newlogZsur = cppf_logZsur(newlogZ, indsur);
// //         newloglam = cppf_loglam_b0t(Xwhale, newlogZ, beta, indtwhale);
// //         newlam = exp(newloglam);
// //         
// //         if(curniter > joiniter){
// //           llnew = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //             cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) +
// //             cppf_loglik_obl(logYobl, newlogZobl, sig2_obl) +
// //             cppf_loglik_sur(logYsur, newlogZsur, lam0, lam1, sig2_sur);
// //           
// //         } else {
// //           llnew = cppf_loglik_obl(logYobl, newlogZobl, sig2_obl) +
// //             cppf_loglik_sur(logYsur, newlogZsur, lam0, lam1, sig2_sur);
// //         }
// //         
// //         if(llnew > llthred){
// //           llprev = llnew;
// //           accept = true;
// //         } else {
// //           if(theta < 0){ thetamin = theta; } else { thetamax = theta; }
// //           theta = thetamin + randu() * (thetamax - thetamin);
// //           if ( (count) % 1000 == 0 ) {
// //             Rprintf("ESS for eta: %d iterations...\n", count);
// //           }
// //           if (count > 20000) stop("ESS is not converging...\n");
// //         }
// //       }
// //       eta = neweta;
// //       logZ = newlogZ;
// //       logZobl = newlogZobl;
// //       logZsur = newlogZsur;
// //       loglam = newloglam;
// //       lam = newlam;
// //       
// //       postEtat = as<mat>(postEta[t]);
// //       postEtat.row(s) = trans(eta.col(t));
// //       postEta[t] = postEtat;
// //     }
// //     if(s == 0){ Rprintf("sampled eta\n"); }
// //     
// //     
// //     // update kappa_eta
// //     kappa_eta = cppf_gibbs_kappa_eta(eta, invReta, 2, 2);
// //     postKappa(s,0) = kappa_eta;
// //     if(s == 0){ Rprintf("sampled kappa_eta\n"); }
// //     
// //     
// //     // update lam0
// //     lam0 = cppf_gibbs_lam0(logYsur, logZsur, lam1, sig2_sur, 100);
// //     postLam(s, 0) = lam0;
// //     if(s == 0){ Rprintf("sampled lam0\n"); }
// //     
// //     
// //     // update lam1
// //     lam1 = cppf_gibbs_lam1(logYsur, logZsur, lam0, sig2_sur, 100);
// //     postLam(s, 1) = lam1;
// //     if(s == 0){ Rprintf("sampled lam1\n"); }
// //     
// //     
// //     // update sig2_sur
// //     sig2_sur = cppf_gibbs_sig2_sur(logYsur, logZsur, lam0, lam1, 2, 2);
// //     postSig2(s, 0) = sig2_sur;
// //     if(s == 0){ Rprintf("sampled sig2_sur\n"); }
// //     
// //     
// //     // update beta
// //     if( updateCOV ){
// //       if( (s+1 >= adaptInterval) && (s+1 - (adaptInterval * trunc((s+1) / adaptInterval)) == 0) ){
// //         dummyaccprob = accprob.col(n_al0+1);
// //         rhat[n_al0+1] = sum( dummyaccprob.rows(s+1-adaptInterval, s-1) ) / (adaptInterval-1);
// //         gamma1[n_al0+1] = 1 / pow(adapIter[n_al0+1], c1);
// //         gamma2[n_al0+1] = c0 * gamma1[n_al0+1];
// //         sigma2[n_al0+1] = exp( log(sigma2[n_al0+1]) + gamma2[n_al0+1] * (rhat[n_al0+1] - ropt) );
// //         
// //         COVbeta = COVbeta + gamma1[n_al0+1] * ( cov( postBeta.rows(s+1-adaptInterval, s-1)) - COVbeta );
// //         cholCOVbeta = trans( chol( sigma2[n_al0+1] * ( COVbeta + 0.0000000001 * diagmat(ones(n_be)) ) ) );
// //         
// //         adapIter[n_al0+1] = adapIter[n_al0+1] + 1;
// //       }
// //     }
// //     
// //     newbeta = beta + cholCOVbeta * randn(n_be);
// //     
// //     newloglam = cppf_loglam_b0t(Xwhale, logZ, newbeta, indtwhale);
// //     newlam = exp(newloglam);
// //     
// //     logprob = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) -
// //       cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //       cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) -
// //       cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
// //       MVN_logh(newbeta, zeros(n_be), diagmat(ones(n_be))/100) -
// //       MVN_logh(beta, zeros(n_be), diagmat(ones(n_be))/100);
// //     
// //     if (log(randu()) < logprob) {
// //       
// //       beta = newbeta;
// //       loglam = newloglam;
// //       lam = newlam;
// //       
// //       accprob(s, n_al0+1) = 1;
// //     }
// //     
// //     postBeta.row(s) = trans(beta);
// //     if(s == 0){ Rprintf("sampled beta\n"); }
// //     
// //     
// //     
// //     if ( (s+1) % 10 == 0 ) {
// //       Rprintf("Generated %d samples...\n", s+1);
// //     }
// //   }
// //   
// //   return Rcpp::List::create(Rcpp::Named("postAlpha0") = postAlpha0,
// //                             Rcpp::Named("postAlpha0tilde") = postAlpha0tilde,
// //                             Rcpp::Named("postAlpha") = postAlpha,
// //                             Rcpp::Named("postEta") = postEta,
// //                             Rcpp::Named("postLam") = postLam,
// //                             Rcpp::Named("postKappa") = postKappa,
// //                             Rcpp::Named("postSig2") = postSig2,
// //                             Rcpp::Named("postBeta") = postBeta,
// //                             Rcpp::Named("Accprob") = accprob,
// //                             Rcpp::Named("sigma2") = sigma2,
// //                             Rcpp::Named("adapIter") = adapIter,
// //                             Rcpp::Named("COValpha0") = COValpha0,
// //                             Rcpp::Named("COValpha") = COValpha,
// //                             Rcpp::Named("COVbeta") = COVbeta);
// //   
// // }
// // 
// // 
// // 
// // 
// // 
// // // ============================================================================-
// // // Fix daily call rates at truth
// // // Preferential sampling
// // // ============================================================================-
// // 
// // // pref: preferential sampling for zoop
// // // b0t: daily intercepts for whale
// // // psit: daily GPs for whale
// // // ja: joint updating for alpha_0 and alpha
// // // fixC: fix daily call rates at truth
// // 
// // // [[Rcpp::export]]
// // List cppf_fitJoint_pref_b0t_psit_ja_fixC(
// //     int niter, mat distmat,
// //     List Xzoop, vec logYobl, vec logYsur, mat Spref,
// //     List Xwhale, mat Sdist, mat Wpam, vec dur_dist,
// //     vec f_dist_obs, List f_dist, List p_pam,
// //     mat indobl, mat indsur, vec indpref, vec indtwhale, mat inddist, vec indmdist, vec indmwhale,
// //     double alpha0tilde, vec alpha, mat eta, 
// //     double mu0, vec delta, double lam0, double lam1,
// //     double tau2, double kappa_eta, double phi_eta, double kappa_delta, double phi_delta,
// //     double sig2_obl, double sig2_sur,
// //     vec beta, mat psi, vec pii, vec c,
// //     double kappa_psi, double phi_psi,
// //     double cellarea,
// //     vec sigma2, mat COValpha, mat COVbeta, double COVmu0,
// //     bool updateCOV, int adaptInterval, double adaptFactorExponent, vec adapIter){
// //   
// //   double negativeInf = -std::numeric_limits<float>::infinity();;
// //   
// //   indobl = indobl-1, indsur = indsur-1, indpref = indpref-1, indtwhale = indtwhale-1;
// //   inddist = inddist-1, indmdist = indmdist-1, indmwhale = indmwhale-1;
// //   
// //   int nG = distmat.n_rows, n_al = alpha.size(), n_be = beta.size(), ntzoop = eta.n_cols, ntwhale = psi.n_cols, n_sigma2 = sigma2.size();
// //   
// //   mat postAlpha0tilde = zeros(niter, 1), postAlpha = zeros(niter, n_al), postLam = zeros(niter, 2);
// //   mat postDelta = zeros(niter, delta.size()), postMu0(niter, 1);
// //   mat postSig2 = zeros(niter, 1), postKappa = zeros(niter, 3);
// //   mat postBeta = zeros(niter, n_be);
// //   mat postEtat = zeros(niter, nG), postPsit = zeros(niter, nG);
// //   List postEta(ntzoop), postPsi(ntwhale);
// //   for(int i = 0; i < ntzoop; i ++){ postEta[i] = zeros(niter, nG); }
// //   for(int i = 0; i < ntwhale; i ++){ postPsi[i] = zeros(niter, nG); }
// //   Rprintf("set up matrices for posterior samples\n");
// //   
// //   mat accprob = zeros(niter, n_sigma2);
// //   vec rhat = zeros(n_sigma2), gamma1 = zeros(n_sigma2), gamma2 = zeros(n_sigma2), dummyaccprob, postj;
// //   double c0 = 1, c1 = adaptFactorExponent, ropt = 0.234, logprob;
// //   mat cholCOValpha = trans( chol( sigma2[0] * ( COValpha + 0.0000000001 * diagmat(ones(n_al)) ) ) );
// //   mat cholCOVbeta = trans( chol( sigma2[1] * ( COVbeta + 0.0000000001 * diagmat(ones(n_be)) ) ) );
// //   double cholCOVmu0 = sqrt( sigma2[2] * COVmu0 );
// //   Rprintf("set up adaptive proposals\n");
// //   
// //   double newmu0;
// //   vec newalpha, newbeta, newlogZobl, newlogZsur, newloglam_pref, newlam_pref;
// //   mat newlogZ, newloglam, newlam;
// //   vec dummyvec;
// //   
// //   int count;
// //   vec xvals;
// //   double llprev, llthred, llnew, thetamin, thetamax, theta;
// //   vec prioretat, priorpsit, priordelta, newdelta;
// //   mat neweta, newpsi;
// //   bool accept;
// //   
// //   mat Reta = exp(-distmat / phi_eta), Rpsi = exp(-distmat / phi_psi), Rdelta = exp(-distmat / phi_delta);
// //   mat invReta = inv(Reta), invRpsi= inv(Rpsi), invRdelta = inv(Rdelta);
// //   mat cholReta = chol(Reta), cholRpsi = chol(Rpsi), cholRdelta = chol(Rdelta);
// //   rowvec oneinvReta = trans(ones(nG)) * invReta, oneinvRpsi = trans(ones(nG)) * invRpsi, oneinvRdelta = trans(ones(nG)) * invRdelta;
// //   double oneinvRetaone = sum( oneinvReta * ones(nG) ), oneinvRpsione = sum( oneinvRpsi * ones(nG) ), oneinvRdeltaone = sum( oneinvRdelta * ones(nG) );
// //   Rprintf("set up GP covariance matrices\n");
// //   
// //   mat logZ = cppf_logZ_a0t_etat_ja_pref(Xzoop, alpha, eta, delta);
// //   vec logZobl = cppf_logZobl(logZ, indobl);
// //   vec logZsur = cppf_logZsur(logZ, indsur);
// //   mat loglam = cppf_loglam_b0t_psit(Xwhale, logZ, beta, psi, indtwhale);
// //   mat lam = exp(loglam);
// //   vec loglam_pref = cppf_loglam_pref(mu0, delta);
// //   vec lam_pref = exp(loglam_pref);
// //   Rprintf("set up initial values\n");
// //   
// //   
// //   // Start MCMC
// //   for(int s = 0; s < niter; s++) {
// //     
// //     // update alpha
// //     if( updateCOV ){
// //       if( (s+1 >= adaptInterval) && (s+1 - (adaptInterval * trunc((s+1) / adaptInterval)) == 0) ){
// //         dummyaccprob = accprob.col(0);
// //         rhat[0] = sum( dummyaccprob.rows(s+1-adaptInterval, s-1) ) / (adaptInterval-1);
// //         gamma1[0] = 1 / pow(adapIter[0], c1);
// //         gamma2[0] = c0 * gamma1[0];
// //         sigma2[0] = exp( log(sigma2[0]) + gamma2[0] * (rhat[0] - ropt) );
// //         
// //         COValpha = COValpha + gamma1[0] * ( cov( postAlpha.rows(s+1-adaptInterval, s-1) ) - COValpha );
// //         cholCOValpha = trans( chol( sigma2[0] * ( COValpha + 0.0000000001 * diagmat(ones(n_al)) ) ) );
// //         
// //         adapIter[0] = adapIter[0] + 1;
// //       }
// //     }
// //     
// //     newalpha = alpha + cholCOValpha * randn(n_al);
// //     
// //     newlogZ = cppf_logZ_a0t_etat_ja_pref(Xzoop, newalpha, eta, delta);
// //     newlogZobl = cppf_logZobl(newlogZ, indobl);
// //     newlogZsur = cppf_logZsur(newlogZ, indsur);
// //     newloglam = cppf_loglam_b0t_psit(Xwhale, newlogZ, beta, psi, indtwhale);
// //     newlam = exp(newloglam);
// //     
// //     logprob = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) -
// //       cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //       cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) -
// //       cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
// //       cppf_loglik_obl(logYobl, newlogZobl, sig2_obl) -
// //       cppf_loglik_obl(logYobl, logZobl, sig2_obl) +
// //       cppf_loglik_sur(logYsur, newlogZsur, lam0, lam1, sig2_sur) -
// //       cppf_loglik_sur(logYsur, logZsur, lam0, lam1, sig2_sur) +
// //       MVN_logh(newalpha, zeros(n_al), diagmat(ones(n_al))/100) -
// //       MVN_logh(alpha, zeros(n_al), diagmat(ones(n_al))/100);
// //     
// //     if (log(randu()) < logprob) {
// //       
// //       alpha = newalpha;
// //       logZ = newlogZ;
// //       logZobl = newlogZobl;
// //       logZsur = newlogZsur;
// //       loglam = newloglam;
// //       lam = newlam;
// //       
// //       accprob(s, 0) = 1;
// //     }
// //     postAlpha.row(s) = trans(alpha);
// //     if(s == 0){ Rprintf("sampled alpha\n"); }
// //     
// //     
// //     // Gibbs update for alpha0tilde
// //     alpha0tilde = cppf_gibbs_alpha0tilde(alpha.rows(0, ntzoop-1), tau2);
// //     postAlpha0tilde[s] = alpha0tilde;
// //     if(s == 0){ Rprintf("sampled alpha0tilde\n"); }
// //     
// //     
// //     
// //     // update eta via Elliptical slice sampling
// //     for(int t = 0; t < ntzoop; t ++){
// //       neweta = eta;
// //       
// //       llprev = cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //         cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
// //         cppf_loglik_obl(logYobl, logZobl, sig2_obl) +
// //         cppf_loglik_sur(logYsur, logZsur, lam0, lam1, sig2_sur);
// //       prioretat = sqrt(kappa_eta) * trans(cholReta) * randn(nG);
// //       
// //       thetamin = 0;
// //       thetamax = 2 * M_PI;
// //       theta = thetamin + randu() * (thetamax - thetamin);
// //       thetamin = theta - 2 * M_PI;
// //       thetamax = theta;
// //       
// //       llthred = llprev + log(randu());
// //       accept = false;
// //       count = 0;
// //       
// //       while(accept == false){
// //         count = count + 1;
// //         
// //         neweta.col(t) = eta.col(t) * cos(theta) + prioretat * sin(theta);
// //         neweta.col(t) = neweta.col(t) - mean(neweta.col(t));
// //         
// //         newlogZ = cppf_logZ_a0t_etat_ja_pref(Xzoop, alpha, neweta, delta);
// //         newlogZobl = cppf_logZobl(newlogZ, indobl);
// //         newlogZsur = cppf_logZsur(newlogZ, indsur);
// //         newloglam = cppf_loglam_b0t_psit(Xwhale, newlogZ, beta, psi, indtwhale);
// //         newlam = exp(newloglam);
// //         
// //         llnew = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //           cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) +
// //           cppf_loglik_obl(logYobl, newlogZobl, sig2_obl) +
// //           cppf_loglik_sur(logYsur, newlogZsur, lam0, lam1, sig2_sur);
// //         
// //         if(llnew > llthred){
// //           llprev = llnew;
// //           accept = true;
// //         } else {
// //           if(theta < 0){ thetamin = theta; } else { thetamax = theta; }
// //           theta = thetamin + randu() * (thetamax - thetamin);
// //           if ( (count) % 1000 == 0 ) {
// //             Rprintf("ESS for eta: %d iterations...\n", count);
// //           }
// //           if (count > 20000) stop("ESS is not converging...\n");
// //         }
// //       }
// //       eta = neweta;
// //       logZ = newlogZ;
// //       logZobl = newlogZobl;
// //       logZsur = newlogZsur;
// //       loglam = newloglam;
// //       lam = newlam;
// //       
// //       postEtat = as<mat>(postEta[t]);
// //       postEtat.row(s) = trans(eta.col(t));
// //       postEta[t] = postEtat;
// //     }
// //     if(s == 0){ Rprintf("sampled eta\n"); }
// //     
// //     
// //     // update kappa_eta
// //     kappa_eta = cppf_gibbs_kappa_eta(eta, invReta, 2, 2);
// //     postKappa(s,0) = kappa_eta;
// //     if(s == 0){ Rprintf("sampled kappa_eta\n"); }
// //     
// //     
// //     // update lam0
// //     lam0 = cppf_gibbs_lam0(logYsur, logZsur, lam1, sig2_sur, 100);
// //     postLam(s, 0) = lam0;
// //     if(s == 0){ Rprintf("sampled lam0\n"); }
// //     
// //     
// //     // update lam1
// //     lam1 = cppf_gibbs_lam1(logYsur, logZsur, lam0, sig2_sur, 100);
// //     postLam(s, 1) = lam1;
// //     if(s == 0){ Rprintf("sampled lam1\n"); }
// //     
// //     
// //     // update sig2_sur
// //     sig2_sur = cppf_gibbs_sig2_sur(logYsur, logZsur, lam0, lam1, 2, 2);
// //     postSig2(s, 0) = sig2_sur;
// //     if(s == 0){ Rprintf("sampled sig2_sur\n"); }
// //     
// //     
// //     // update mu0
// //     if( updateCOV ){
// //       if( (s+1 >= adaptInterval) && (s+1 - (adaptInterval * trunc((s+1) / adaptInterval)) == 0) ){
// //         dummyaccprob = accprob.col(2);
// //         rhat[2] = sum( dummyaccprob.rows(s+1-adaptInterval, s-1) ) / (adaptInterval-1);
// //         gamma1[2] = 1 / pow(adapIter[2], c1);
// //         gamma2[2] = c0 * gamma1[2];
// //         sigma2[2] = exp( log(sigma2[2]) + gamma2[2] * (rhat[2] - ropt) );
// //         
// //         postj = postMu0.col(0);
// //         COVmu0 = COVmu0 + gamma1[2] * ( var( postj.rows(s+1-adaptInterval, s-1)) - COVmu0 );
// //         cholCOVmu0 = sqrt( sigma2[2] * COVmu0 );
// //         
// //         adapIter[2] = adapIter[2] + 1;
// //       }
// //     }
// //     
// //     newmu0 = mu0 + cholCOVmu0 * randn();
// //     
// //     newloglam_pref = cppf_loglam_pref(newmu0, delta);
// //     newlam_pref = exp(newloglam_pref);
// //     
// //     logprob = cppf_loglik_pref(newloglam_pref, newlam_pref, indpref, cellarea) - 
// //       cppf_loglik_pref(loglam_pref, lam_pref, indpref, cellarea) +
// //       normal_logh(newmu0, 0, 100) - normal_logh(mu0, 0, 100);
// //     
// //     if (log(randu()) < logprob) {
// //       
// //       mu0 = newmu0;
// //       loglam_pref = newloglam_pref;
// //       lam_pref = newlam_pref;
// //       
// //       accprob(s, 2) = 1;
// //     }
// //     postMu0(s,0) = mu0;
// //     if(s == 0){ Rprintf("sampled mu0\n"); }
// //     
// //     
// //     // update delta
// //     llprev = cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //       cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
// //       cppf_loglik_obl(logYobl, logZobl, sig2_obl) +
// //       cppf_loglik_sur(logYsur, logZsur, lam0, lam1, sig2_sur) +
// //       cppf_loglik_pref(loglam_pref, lam_pref, indpref, cellarea);
// //     priordelta = sqrt(kappa_delta) * trans(cholRdelta) * randn(nG);
// //     
// //     thetamin = 0;
// //     thetamax = 2 * M_PI;
// //     theta = thetamin + randu() * (thetamax - thetamin);
// //     thetamin = theta - 2 * M_PI;
// //     thetamax = theta;
// //     
// //     llthred = llprev + log(randu());
// //     accept = false;
// //     count = 0;
// //     
// //     while(accept == false){
// //       count = count + 1;
// //       
// //       newdelta = delta * cos(theta) + priordelta * sin(theta);
// //       newdelta = newdelta - mean(newdelta);
// //       
// //       newlogZ = cppf_logZ_a0t_etat_ja_pref(Xzoop, alpha, eta, newdelta);
// //       newlogZobl = cppf_logZobl(newlogZ, indobl);
// //       newlogZsur = cppf_logZsur(newlogZ, indsur);
// //       newloglam = cppf_loglam_b0t_psit(Xwhale, newlogZ, beta, psi, indtwhale);
// //       newlam = exp(newloglam);
// //       newloglam_pref = cppf_loglam_pref(mu0, newdelta);
// //       newlam_pref = exp(newloglam_pref);
// //       
// //       llnew = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //         cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) +
// //         cppf_loglik_obl(logYobl, newlogZobl, sig2_obl) +
// //         cppf_loglik_sur(logYsur, newlogZsur, lam0, lam1, sig2_sur) +
// //         cppf_loglik_pref(newloglam_pref, newlam_pref, indpref, cellarea);
// //       
// //       if(llnew > llthred){
// //         llprev = llnew;
// //         accept = true;
// //       } else {
// //         if(theta < 0){ thetamin = theta; } else { thetamax = theta; }
// //         theta = thetamin + randu() * (thetamax - thetamin);
// //         if ( (count) % 1000 == 0 ) {
// //           Rprintf("ESS for delta: %d iterations...\n", count);
// //         }
// //         if (count > 20000) stop("ESS is not converging...\n");
// //       }
// //     }
// //     delta = newdelta;
// //     logZ = newlogZ;
// //     logZobl = newlogZobl;
// //     logZsur = newlogZsur;
// //     loglam = newloglam;
// //     lam = newlam;
// //     loglam_pref = newloglam_pref;
// //     lam_pref = newlam_pref;
// //     
// //     postDelta.row(s) = trans(delta);
// //     if(s == 0){ Rprintf("sampled delta\n"); }
// //     
// //     
// //     // update kappa_delta
// //     kappa_delta = cppf_gibbs_kappa_delta(delta, invRdelta, 2, 2);
// //     postKappa(s, 2) = kappa_delta;
// //     if(s == 0){ Rprintf("sampled kappa_delta\n"); }
// //     
// //     
// //     // update beta
// //     if( updateCOV ){
// //       if( (s+1 >= adaptInterval) && (s+1 - (adaptInterval * trunc((s+1) / adaptInterval)) == 0) ){
// //         dummyaccprob = accprob.col(1);
// //         rhat[1] = sum( dummyaccprob.rows(s+1-adaptInterval, s-1) ) / (adaptInterval-1);
// //         gamma1[1] = 1 / pow(adapIter[1], c1);
// //         gamma2[1] = c0 * gamma1[1];
// //         sigma2[1] = exp( log(sigma2[1]) + gamma2[1] * (rhat[1] - ropt) );
// //         
// //         COVbeta = COVbeta + gamma1[1] * ( cov( postBeta.rows(s+1-adaptInterval, s-1)) - COVbeta );
// //         cholCOVbeta = trans( chol( sigma2[1] * ( COVbeta + 0.0000000001 * diagmat(ones(n_be)) ) ) );
// //         
// //         adapIter[1] = adapIter[1] + 1;
// //       }
// //     }
// //     
// //     newbeta = beta + cholCOVbeta * randn(n_be);
// //     
// //     newloglam = cppf_loglam_b0t_psit(Xwhale, logZ, newbeta, psi, indtwhale);
// //     newlam = exp(newloglam);
// //     
// //     logprob = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) -
// //       cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //       cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) -
// //       cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
// //       MVN_logh(newbeta, zeros(n_be), diagmat(ones(n_be))/100) -
// //       MVN_logh(beta, zeros(n_be), diagmat(ones(n_be))/100);
// //     
// //     if (log(randu()) < logprob) {
// //       
// //       beta = newbeta;
// //       loglam = newloglam;
// //       lam = newlam;
// //       
// //       accprob(s, 1) = 1;
// //     }
// //     
// //     postBeta.row(s) = trans(beta);
// //     if(s == 0){ Rprintf("sampled beta\n"); }
// //     
// //     
// //     
// //     // update psi via Elliptical slice sampling
// //     for(int t = 0; t < ntwhale; t ++){
// //       newpsi = psi;
// //       
// //       llprev = cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //         cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea);
// //       priorpsit = sqrt(kappa_psi) * trans(cholRpsi) * randn(nG);
// //       
// //       thetamin = 0;
// //       thetamax = 2 * M_PI;
// //       theta = thetamin + randu() * (thetamax - thetamin);
// //       thetamin = theta - 2 * M_PI;
// //       thetamax = theta;
// //       
// //       llthred = llprev + log(randu());
// //       accept = false;
// //       count = 0;
// //       
// //       while(accept == false){
// //         count = count + 1;
// //         
// //         newpsi.col(t) = psi.col(t) * cos(theta) + priorpsit * sin(theta);
// //         newpsi.col(t) = newpsi.col(t) - mean(newpsi.col(t));
// //         
// //         newloglam = cppf_loglam_b0t_psit(Xwhale, logZ, beta, newpsi, indtwhale);
// //         newlam = exp(newloglam);
// //         
// //         llnew = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //           cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea);
// //         
// //         if(llnew > llthred){
// //           llprev = llnew;
// //           accept = true;
// //         } else {
// //           if(theta < 0){ thetamin = theta; } else { thetamax = theta; }
// //           theta = thetamin + randu() * (thetamax - thetamin);
// //           if ( (count) % 1000 == 0 ) {
// //             Rprintf("ESS for psi: %d iterations...\n", count);
// //           }
// //           if (count > 20000) stop("ESS is not converging...\n");
// //         }
// //       }
// //       psi = newpsi;
// //       loglam = newloglam;
// //       lam = newlam;
// //       
// //       postPsit = as<mat>(postPsi[t]);
// //       postPsit.row(s) = trans(psi.col(t));
// //       postPsi[t] = postPsit;
// //     }
// //     if(s == 0){ Rprintf("sampled psi\n"); }
// //     
// //     
// //     // update kappa_psi
// //     kappa_psi = cppf_gibbs_kappa_psit(psi, invRpsi, 2, 2);
// //     postKappa(s,1) = kappa_psi;
// //     if(s == 0){ Rprintf("sampled kappa_psi\n"); }
// //     
// //     
// //     if ( (s+1) % 10 == 0 ) {
// //       Rprintf("Generated %d samples...\n", s+1);
// //     }
// //   }
// //   
// //   
// //   return Rcpp::List::create(Rcpp::Named("postAlpha0tilde") = postAlpha0tilde,
// //                             Rcpp::Named("postAlpha") = postAlpha,
// //                             Rcpp::Named("postEta") = postEta,
// //                             Rcpp::Named("postLam") = postLam,
// //                             Rcpp::Named("postKappa") = postKappa,
// //                             Rcpp::Named("postSig2") = postSig2,
// //                             Rcpp::Named("postMu0") = postMu0,
// //                             Rcpp::Named("postDelta") = postDelta,
// //                             Rcpp::Named("postBeta") = postBeta,
// //                             Rcpp::Named("postPsi") = postPsi,
// //                             Rcpp::Named("Accprob") = accprob,
// //                             Rcpp::Named("sigma2") = sigma2,
// //                             Rcpp::Named("adapIter") = adapIter,
// //                             Rcpp::Named("COValpha") = COValpha,
// //                             Rcpp::Named("COVmu0") = COVmu0,
// //                             Rcpp::Named("COVbeta") = COVbeta);
// //   
// // }
// // 
// // 
// // 
// // // [[Rcpp::export]]
// // List cppf_fitJoint_pref_b0t_psi_ja_fixC(
// //     int niter, mat distmat,
// //     List Xzoop, vec logYobl, vec logYsur, mat Spref,
// //     List Xwhale, mat Sdist, mat Wpam, vec dur_dist,
// //     vec f_dist_obs, List f_dist, List p_pam,
// //     mat indobl, mat indsur, vec indpref, vec indtwhale, mat inddist, vec indmdist, vec indmwhale,
// //     double alpha0tilde, vec alpha, mat eta, 
// //     double mu0, vec delta, double lam0, double lam1,
// //     double tau2, double kappa_eta, double phi_eta, double kappa_delta, double phi_delta,
// //     double sig2_obl, double sig2_sur,
// //     vec beta, vec psi, vec pii, vec c,
// //     double kappa_psi, double phi_psi,
// //     double cellarea,
// //     vec sigma2, mat COValpha, mat COVbeta, double COVmu0,
// //     bool updateCOV, int adaptInterval, double adaptFactorExponent, vec adapIter){
// //   
// //   double negativeInf = -std::numeric_limits<float>::infinity();;
// //   
// //   indobl = indobl-1, indsur = indsur-1, indpref = indpref-1, indtwhale = indtwhale-1;
// //   inddist = inddist-1, indmdist = indmdist-1, indmwhale = indmwhale-1;
// //   
// //   int nG = distmat.n_rows, n_al = alpha.size(), n_be = beta.size(), ntzoop = Xzoop.size(), ntwhale = Xwhale.size(), n_sigma2 = sigma2.size();
// //   
// //   mat postAlpha0tilde = zeros(niter, 1), postAlpha = zeros(niter, n_al), postLam = zeros(niter, 2);
// //   mat postDelta = zeros(niter, delta.size()), postMu0(niter, 1);
// //   mat postSig2 = zeros(niter, 1), postKappa = zeros(niter, 3);
// //   mat postBeta = zeros(niter, n_be);
// //   mat postEtat = zeros(niter, nG), postPsi = zeros(niter, nG);
// //   List postEta(ntzoop);
// //   for(int i = 0; i < ntzoop; i ++){ postEta[i] = zeros(niter, nG); }
// //   Rprintf("set up matrices for posterior samples\n");
// //   
// //   mat accprob = zeros(niter, n_sigma2);
// //   vec rhat = zeros(n_sigma2), gamma1 = zeros(n_sigma2), gamma2 = zeros(n_sigma2), dummyaccprob, postj;
// //   double c0 = 1, c1 = adaptFactorExponent, ropt = 0.234, logprob;
// //   mat cholCOValpha = trans( chol( sigma2[0] * ( COValpha + 0.0000000001 * diagmat(ones(n_al)) ) ) );
// //   mat cholCOVbeta = trans( chol( sigma2[1] * ( COVbeta + 0.0000000001 * diagmat(ones(n_be)) ) ) );
// //   double cholCOVmu0 = sqrt( sigma2[2] * COVmu0 );
// //   Rprintf("set up adaptive proposals\n");
// //   
// //   double newmu0;
// //   vec newalpha, newbeta, newlogZobl, newlogZsur, newloglam_pref, newlam_pref;
// //   mat newlogZ, newloglam, newlam;
// //   vec dummyvec;
// //   
// //   int count;
// //   vec xvals;
// //   double llprev, llthred, llnew, thetamin, thetamax, theta;
// //   vec prioretat, priorpsi, newpsi, priordelta, newdelta;
// //   mat neweta;
// //   bool accept;
// //   
// //   mat Reta = exp(-distmat / phi_eta), Rpsi = exp(-distmat / phi_psi), Rdelta = exp(-distmat / phi_delta);
// //   mat invReta = inv(Reta), invRpsi= inv(Rpsi), invRdelta = inv(Rdelta);
// //   mat cholReta = chol(Reta), cholRpsi = chol(Rpsi), cholRdelta = chol(Rdelta);
// //   rowvec oneinvReta = trans(ones(nG)) * invReta, oneinvRpsi = trans(ones(nG)) * invRpsi, oneinvRdelta = trans(ones(nG)) * invRdelta;
// //   double oneinvRetaone = sum( oneinvReta * ones(nG) ), oneinvRpsione = sum( oneinvRpsi * ones(nG) ), oneinvRdeltaone = sum( oneinvRdelta * ones(nG) );
// //   Rprintf("set up GP covariance matrices\n");
// //   
// //   mat logZ = cppf_logZ_a0t_etat_ja_pref(Xzoop, alpha, eta, delta);
// //   vec logZobl = cppf_logZobl(logZ, indobl);
// //   vec logZsur = cppf_logZsur(logZ, indsur);
// //   mat loglam = cppf_loglam_b0t_psi(Xwhale, logZ, beta, psi, indtwhale);
// //   mat lam = exp(loglam);
// //   vec loglam_pref = cppf_loglam_pref(mu0, delta);
// //   vec lam_pref = exp(loglam_pref);
// //   Rprintf("set up initial values\n");
// //   
// //   
// //   // Start MCMC
// //   for(int s = 0; s < niter; s++) {
// //     
// //     // update alpha
// //     if( updateCOV ){
// //       if( (s+1 >= adaptInterval) && (s+1 - (adaptInterval * trunc((s+1) / adaptInterval)) == 0) ){
// //         dummyaccprob = accprob.col(0);
// //         rhat[0] = sum( dummyaccprob.rows(s+1-adaptInterval, s-1) ) / (adaptInterval-1);
// //         gamma1[0] = 1 / pow(adapIter[0], c1);
// //         gamma2[0] = c0 * gamma1[0];
// //         sigma2[0] = exp( log(sigma2[0]) + gamma2[0] * (rhat[0] - ropt) );
// //         
// //         COValpha = COValpha + gamma1[0] * ( cov( postAlpha.rows(s+1-adaptInterval, s-1) ) - COValpha );
// //         cholCOValpha = trans( chol( sigma2[0] * ( COValpha + 0.0000000001 * diagmat(ones(n_al)) ) ) );
// //         
// //         adapIter[0] = adapIter[0] + 1;
// //       }
// //     }
// //     
// //     newalpha = alpha + cholCOValpha * randn(n_al);
// //     
// //     newlogZ = cppf_logZ_a0t_etat_ja_pref(Xzoop, newalpha, eta, delta);
// //     newlogZobl = cppf_logZobl(newlogZ, indobl);
// //     newlogZsur = cppf_logZsur(newlogZ, indsur);
// //     newloglam = cppf_loglam_b0t_psi(Xwhale, newlogZ, beta, psi, indtwhale);
// //     newlam = exp(newloglam);
// //     
// //     logprob = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) -
// //       cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //       cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) -
// //       cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
// //       cppf_loglik_obl(logYobl, newlogZobl, sig2_obl) -
// //       cppf_loglik_obl(logYobl, logZobl, sig2_obl) +
// //       cppf_loglik_sur(logYsur, newlogZsur, lam0, lam1, sig2_sur) -
// //       cppf_loglik_sur(logYsur, logZsur, lam0, lam1, sig2_sur) +
// //       MVN_logh(newalpha, zeros(n_al), diagmat(ones(n_al))/100) -
// //       MVN_logh(alpha, zeros(n_al), diagmat(ones(n_al))/100);
// //     
// //     if (log(randu()) < logprob) {
// //       
// //       alpha = newalpha;
// //       logZ = newlogZ;
// //       logZobl = newlogZobl;
// //       logZsur = newlogZsur;
// //       loglam = newloglam;
// //       lam = newlam;
// //       
// //       accprob(s, 0) = 1;
// //     }
// //     postAlpha.row(s) = trans(alpha);
// //     if(s == 0){ Rprintf("sampled alpha\n"); }
// //     
// //     
// //     // Gibbs update for alpha0tilde
// //     alpha0tilde = cppf_gibbs_alpha0tilde(alpha.rows(0, ntzoop-1), tau2);
// //     postAlpha0tilde[s] = alpha0tilde;
// //     if(s == 0){ Rprintf("sampled alpha0tilde\n"); }
// //     
// //     
// //     
// //     // update eta via Elliptical slice sampling
// //     for(int t = 0; t < ntzoop; t ++){
// //       neweta = eta;
// //       
// //       llprev = cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //         cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
// //         cppf_loglik_obl(logYobl, logZobl, sig2_obl) +
// //         cppf_loglik_sur(logYsur, logZsur, lam0, lam1, sig2_sur);
// //       prioretat = sqrt(kappa_eta) * trans(cholReta) * randn(nG);
// //       
// //       thetamin = 0;
// //       thetamax = 2 * M_PI;
// //       theta = thetamin + randu() * (thetamax - thetamin);
// //       thetamin = theta - 2 * M_PI;
// //       thetamax = theta;
// //       
// //       llthred = llprev + log(randu());
// //       accept = false;
// //       count = 0;
// //       
// //       while(accept == false){
// //         count = count + 1;
// //         
// //         neweta.col(t) = eta.col(t) * cos(theta) + prioretat * sin(theta);
// //         neweta.col(t) = neweta.col(t) - mean(neweta.col(t));
// //         
// //         newlogZ = cppf_logZ_a0t_etat_ja_pref(Xzoop, alpha, neweta, delta);
// //         newlogZobl = cppf_logZobl(newlogZ, indobl);
// //         newlogZsur = cppf_logZsur(newlogZ, indsur);
// //         newloglam = cppf_loglam_b0t_psi(Xwhale, newlogZ, beta, psi, indtwhale);
// //         newlam = exp(newloglam);
// //         
// //         llnew = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //           cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) +
// //           cppf_loglik_obl(logYobl, newlogZobl, sig2_obl) +
// //           cppf_loglik_sur(logYsur, newlogZsur, lam0, lam1, sig2_sur);
// //         
// //         if(llnew > llthred){
// //           llprev = llnew;
// //           accept = true;
// //         } else {
// //           if(theta < 0){ thetamin = theta; } else { thetamax = theta; }
// //           theta = thetamin + randu() * (thetamax - thetamin);
// //           if ( (count) % 1000 == 0 ) {
// //             Rprintf("ESS for eta: %d iterations...\n", count);
// //           }
// //           if (count > 20000) stop("ESS is not converging...\n");
// //         }
// //       }
// //       eta = neweta;
// //       logZ = newlogZ;
// //       logZobl = newlogZobl;
// //       logZsur = newlogZsur;
// //       loglam = newloglam;
// //       lam = newlam;
// //       
// //       postEtat = as<mat>(postEta[t]);
// //       postEtat.row(s) = trans(eta.col(t));
// //       postEta[t] = postEtat;
// //     }
// //     if(s == 0){ Rprintf("sampled eta\n"); }
// //     
// //     
// //     // update kappa_eta
// //     kappa_eta = cppf_gibbs_kappa_eta(eta, invReta, 2, 2);
// //     postKappa(s,0) = kappa_eta;
// //     if(s == 0){ Rprintf("sampled kappa_eta\n"); }
// //     
// //     
// //     // update lam0
// //     lam0 = cppf_gibbs_lam0(logYsur, logZsur, lam1, sig2_sur, 100);
// //     postLam(s, 0) = lam0;
// //     if(s == 0){ Rprintf("sampled lam0\n"); }
// //     
// //     
// //     // update lam1
// //     lam1 = cppf_gibbs_lam1(logYsur, logZsur, lam0, sig2_sur, 100);
// //     postLam(s, 1) = lam1;
// //     if(s == 0){ Rprintf("sampled lam1\n"); }
// //     
// //     
// //     // update sig2_sur
// //     sig2_sur = cppf_gibbs_sig2_sur(logYsur, logZsur, lam0, lam1, 2, 2);
// //     postSig2(s, 0) = sig2_sur;
// //     if(s == 0){ Rprintf("sampled sig2_sur\n"); }
// //     
// //     
// //     // update mu0
// //     if( updateCOV ){
// //       if( (s+1 >= adaptInterval) && (s+1 - (adaptInterval * trunc((s+1) / adaptInterval)) == 0) ){
// //         dummyaccprob = accprob.col(2);
// //         rhat[2] = sum( dummyaccprob.rows(s+1-adaptInterval, s-1) ) / (adaptInterval-1);
// //         gamma1[2] = 1 / pow(adapIter[2], c1);
// //         gamma2[2] = c0 * gamma1[2];
// //         sigma2[2] = exp( log(sigma2[2]) + gamma2[2] * (rhat[2] - ropt) );
// //         
// //         postj = postMu0.col(0);
// //         COVmu0 = COVmu0 + gamma1[2] * ( var( postj.rows(s+1-adaptInterval, s-1)) - COVmu0 );
// //         cholCOVmu0 = sqrt( sigma2[2] * COVmu0 );
// //         
// //         adapIter[2] = adapIter[2] + 1;
// //       }
// //     }
// //     
// //     newmu0 = mu0 + cholCOVmu0 * randn();
// //     
// //     newloglam_pref = cppf_loglam_pref(newmu0, delta);
// //     newlam_pref = exp(newloglam_pref);
// //     
// //     logprob = cppf_loglik_pref(newloglam_pref, newlam_pref, indpref, cellarea) - 
// //       cppf_loglik_pref(loglam_pref, lam_pref, indpref, cellarea) +
// //       normal_logh(newmu0, 0, 100) - normal_logh(mu0, 0, 100);
// //     
// //     if (log(randu()) < logprob) {
// //       
// //       mu0 = newmu0;
// //       loglam_pref = newloglam_pref;
// //       lam_pref = newlam_pref;
// //       
// //       accprob(s, 2) = 1;
// //     }
// //     postMu0(s,0) = mu0;
// //     if(s == 0){ Rprintf("sampled mu0\n"); }
// //     
// //     
// //     // update delta
// //     llprev = cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //       cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
// //       cppf_loglik_obl(logYobl, logZobl, sig2_obl) +
// //       cppf_loglik_sur(logYsur, logZsur, lam0, lam1, sig2_sur) +
// //       cppf_loglik_pref(loglam_pref, lam_pref, indpref, cellarea);
// //     priordelta = sqrt(kappa_delta) * trans(cholRdelta) * randn(nG);
// //     
// //     thetamin = 0;
// //     thetamax = 2 * M_PI;
// //     theta = thetamin + randu() * (thetamax - thetamin);
// //     thetamin = theta - 2 * M_PI;
// //     thetamax = theta;
// //     
// //     llthred = llprev + log(randu());
// //     accept = false;
// //     count = 0;
// //     
// //     while(accept == false){
// //       count = count + 1;
// //       
// //       newdelta = delta * cos(theta) + priordelta * sin(theta);
// //       newdelta = newdelta - mean(newdelta);
// //       
// //       newlogZ = cppf_logZ_a0t_etat_ja_pref(Xzoop, alpha, eta, newdelta);
// //       newlogZobl = cppf_logZobl(newlogZ, indobl);
// //       newlogZsur = cppf_logZsur(newlogZ, indsur);
// //       newloglam = cppf_loglam_b0t_psi(Xwhale, newlogZ, beta, psi, indtwhale);
// //       newlam = exp(newloglam);
// //       newloglam_pref = cppf_loglam_pref(mu0, newdelta);
// //       newlam_pref = exp(newloglam_pref);
// //       
// //       llnew = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //         cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) +
// //         cppf_loglik_obl(logYobl, newlogZobl, sig2_obl) +
// //         cppf_loglik_sur(logYsur, newlogZsur, lam0, lam1, sig2_sur) +
// //         cppf_loglik_pref(newloglam_pref, newlam_pref, indpref, cellarea);
// //       
// //       if(llnew > llthred){
// //         llprev = llnew;
// //         accept = true;
// //       } else {
// //         if(theta < 0){ thetamin = theta; } else { thetamax = theta; }
// //         theta = thetamin + randu() * (thetamax - thetamin);
// //         if ( (count) % 1000 == 0 ) {
// //           Rprintf("ESS for delta: %d iterations...\n", count);
// //         }
// //         if (count > 20000) stop("ESS is not converging...\n");
// //       }
// //     }
// //     delta = newdelta;
// //     logZ = newlogZ;
// //     logZobl = newlogZobl;
// //     logZsur = newlogZsur;
// //     loglam = newloglam;
// //     lam = newlam;
// //     loglam_pref = newloglam_pref;
// //     lam_pref = newlam_pref;
// //     
// //     postDelta.row(s) = trans(delta);
// //     if(s == 0){ Rprintf("sampled delta\n"); }
// //     
// //     
// //     // update kappa_delta
// //     kappa_delta = cppf_gibbs_kappa_delta(delta, invRdelta, 2, 2);
// //     postKappa(s, 2) = kappa_delta;
// //     if(s == 0){ Rprintf("sampled kappa_delta\n"); }
// //     
// //     
// //     // update beta
// //     if( updateCOV ){
// //       if( (s+1 >= adaptInterval) && (s+1 - (adaptInterval * trunc((s+1) / adaptInterval)) == 0) ){
// //         dummyaccprob = accprob.col(1);
// //         rhat[1] = sum( dummyaccprob.rows(s+1-adaptInterval, s-1) ) / (adaptInterval-1);
// //         gamma1[1] = 1 / pow(adapIter[1], c1);
// //         gamma2[1] = c0 * gamma1[1];
// //         sigma2[1] = exp( log(sigma2[1]) + gamma2[1] * (rhat[1] - ropt) );
// //         
// //         COVbeta = COVbeta + gamma1[1] * ( cov( postBeta.rows(s+1-adaptInterval, s-1)) - COVbeta );
// //         cholCOVbeta = trans( chol( sigma2[1] * ( COVbeta + 0.0000000001 * diagmat(ones(n_be)) ) ) );
// //         
// //         adapIter[1] = adapIter[1] + 1;
// //       }
// //     }
// //     
// //     newbeta = beta + cholCOVbeta * randn(n_be);
// //     
// //     newloglam = cppf_loglam_b0t_psi(Xwhale, logZ, newbeta, psi, indtwhale);
// //     newlam = exp(newloglam);
// //     
// //     logprob = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) -
// //       cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //       cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) -
// //       cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
// //       MVN_logh(newbeta, zeros(n_be), diagmat(ones(n_be))/100) -
// //       MVN_logh(beta, zeros(n_be), diagmat(ones(n_be))/100);
// //     
// //     if (log(randu()) < logprob) {
// //       
// //       beta = newbeta;
// //       loglam = newloglam;
// //       lam = newlam;
// //       
// //       accprob(s, 1) = 1;
// //     }
// //     
// //     postBeta.row(s) = trans(beta);
// //     if(s == 0){ Rprintf("sampled beta\n"); }
// //     
// //     
// //     
// //     // update psi via elliptical slice sampling
// //     llprev = cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //       cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea);
// //     priorpsi = sqrt(kappa_psi) * trans(cholRpsi) * randn(nG);
// //     
// //     thetamin = 0;
// //     thetamax = 2 * M_PI;
// //     theta = thetamin + randu() * (thetamax - thetamin);
// //     thetamin = theta - 2 * M_PI;
// //     thetamax = theta;
// //     
// //     llthred = llprev + log(randu());
// //     accept = false;
// //     count = 0;
// //     
// //     while(accept == false){
// //       count = count + 1;
// //       
// //       newpsi = psi * cos(theta) + priorpsi * sin(theta);
// //       newpsi = newpsi - mean(newpsi);
// //       
// //       newloglam = cppf_loglam_b0t_psi(Xwhale, logZ, beta, newpsi, indtwhale);
// //       newlam = exp(newloglam);
// //       
// //       llnew = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //         cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea);
// //       
// //       if(llnew > llthred){
// //         llprev = llnew;
// //         accept = true;
// //       } else {
// //         if(theta < 0){ thetamin = theta; } else { thetamax = theta; }
// //         theta = thetamin + randu() * (thetamax - thetamin);
// //         if ( (count) % 1000 == 0 ) {
// //           Rprintf("ESS for psi: %d iterations...\n", count);
// //         }
// //         if (count > 20000) stop("ESS is not converging...\n");
// //       }
// //     }
// //     psi = newpsi;
// //     loglam = newloglam;
// //     lam = newlam;
// //     
// //     postPsi.row(s) = trans(psi);
// //     if(s == 0){ Rprintf("sampled psi\n"); }
// //     
// //     
// //     // update kappa_psi
// //     kappa_psi = cppf_gibbs_kappa_psi(psi, invRpsi, 2, 2);
// //     postKappa(s,1) = kappa_psi;
// //     if(s == 0){ Rprintf("sampled kappa_psi\n"); }
// //     
// //     
// //     if ( (s+1) % 10 == 0 ) {
// //       Rprintf("Generated %d samples...\n", s+1);
// //     }
// //   }
// //   
// //   
// //   return Rcpp::List::create(Rcpp::Named("postAlpha0tilde") = postAlpha0tilde,
// //                             Rcpp::Named("postAlpha") = postAlpha,
// //                             Rcpp::Named("postEta") = postEta,
// //                             Rcpp::Named("postLam") = postLam,
// //                             Rcpp::Named("postKappa") = postKappa,
// //                             Rcpp::Named("postSig2") = postSig2,
// //                             Rcpp::Named("postMu0") = postMu0,
// //                             Rcpp::Named("postDelta") = postDelta,
// //                             Rcpp::Named("postBeta") = postBeta,
// //                             Rcpp::Named("postPsi") = postPsi,
// //                             Rcpp::Named("Accprob") = accprob,
// //                             Rcpp::Named("sigma2") = sigma2,
// //                             Rcpp::Named("adapIter") = adapIter,
// //                             Rcpp::Named("COValpha") = COValpha,
// //                             Rcpp::Named("COVmu0") = COVmu0,
// //                             Rcpp::Named("COVbeta") = COVbeta);
// //   
// // }
// // 
// // 
// // 
// // // [[Rcpp::export]]
// // List cppf_fitJoint_pref_b0t_ja_fixC(
// //     int niter, mat distmat,
// //     List Xzoop, vec logYobl, vec logYsur, mat Spref,
// //     List Xwhale, mat Sdist, mat Wpam, vec dur_dist,
// //     vec f_dist_obs, List f_dist, List p_pam,
// //     mat indobl, mat indsur, vec indpref, vec indtwhale, mat inddist, vec indmdist, vec indmwhale,
// //     double alpha0tilde, vec alpha, mat eta, 
// //     double mu0, vec delta, double lam0, double lam1,
// //     double tau2, double kappa_eta, double phi_eta, double kappa_delta, double phi_delta,
// //     double sig2_obl, double sig2_sur,
// //     vec beta, vec pii, vec c,
// //     double cellarea,
// //     vec sigma2, mat COValpha, mat COVbeta, double COVmu0,
// //     bool updateCOV, int adaptInterval, double adaptFactorExponent, vec adapIter){
// //   
// //   double negativeInf = -std::numeric_limits<float>::infinity();;
// //   
// //   indobl = indobl-1, indsur = indsur-1, indpref = indpref-1, indtwhale = indtwhale-1;
// //   inddist = inddist-1, indmdist = indmdist-1, indmwhale = indmwhale-1;
// //   
// //   int nG = distmat.n_rows, n_al = alpha.size(), n_be = beta.size(), ntzoop = Xzoop.size(), ntwhale = Xwhale.size(), n_sigma2 = sigma2.size();
// //   
// //   mat postAlpha0tilde = zeros(niter, 1), postAlpha = zeros(niter, n_al), postLam = zeros(niter, 2);
// //   mat postDelta = zeros(niter, delta.size()), postMu0(niter, 1);
// //   mat postSig2 = zeros(niter, 1), postKappa = zeros(niter, 2);
// //   mat postBeta = zeros(niter, n_be);
// //   mat postEtat = zeros(niter, nG);
// //   List postEta(ntzoop);
// //   for(int i = 0; i < ntzoop; i ++){ postEta[i] = zeros(niter, nG); }
// //   Rprintf("set up matrices for posterior samples\n");
// //   
// //   mat accprob = zeros(niter, n_sigma2);
// //   vec rhat = zeros(n_sigma2), gamma1 = zeros(n_sigma2), gamma2 = zeros(n_sigma2), dummyaccprob, postj;
// //   double c0 = 1, c1 = adaptFactorExponent, ropt = 0.234, logprob;
// //   mat cholCOValpha = trans( chol( sigma2[0] * ( COValpha + 0.0000000001 * diagmat(ones(n_al)) ) ) );
// //   mat cholCOVbeta = trans( chol( sigma2[1] * ( COVbeta + 0.0000000001 * diagmat(ones(n_be)) ) ) );
// //   double cholCOVmu0 = sqrt( sigma2[2] * COVmu0 );
// //   Rprintf("set up adaptive proposals\n");
// //   
// //   double newmu0;
// //   vec newalpha, newbeta, newlogZobl, newlogZsur, newloglam_pref, newlam_pref;
// //   mat newlogZ, newloglam, newlam;
// //   vec dummyvec;
// //   
// //   int count;
// //   vec xvals;
// //   double llprev, llthred, llnew, thetamin, thetamax, theta;
// //   vec prioretat, priordelta, newdelta;
// //   mat neweta;
// //   bool accept;
// //   
// //   mat Reta = exp(-distmat / phi_eta), Rdelta = exp(-distmat / phi_delta);
// //   mat invReta = inv(Reta), invRdelta = inv(Rdelta);
// //   mat cholReta = chol(Reta), cholRdelta = chol(Rdelta);
// //   rowvec oneinvReta = trans(ones(nG)) * invReta, oneinvRdelta = trans(ones(nG)) * invRdelta;
// //   double oneinvRetaone = sum( oneinvReta * ones(nG) ), oneinvRdeltaone = sum( oneinvRdelta * ones(nG) );
// //   Rprintf("set up GP covariance matrices\n");
// //   
// //   mat logZ = cppf_logZ_a0t_etat_ja_pref(Xzoop, alpha, eta, delta);
// //   vec logZobl = cppf_logZobl(logZ, indobl);
// //   vec logZsur = cppf_logZsur(logZ, indsur);
// //   mat loglam = cppf_loglam_b0t(Xwhale, logZ, beta, indtwhale);
// //   mat lam = exp(loglam);
// //   vec loglam_pref = cppf_loglam_pref(mu0, delta);
// //   vec lam_pref = exp(loglam_pref);
// //   Rprintf("set up initial values\n");
// //   
// //   
// //   // Start MCMC
// //   for(int s = 0; s < niter; s++) {
// //     
// //     // update alpha
// //     if( updateCOV ){
// //       if( (s+1 >= adaptInterval) && (s+1 - (adaptInterval * trunc((s+1) / adaptInterval)) == 0) ){
// //         dummyaccprob = accprob.col(0);
// //         rhat[0] = sum( dummyaccprob.rows(s+1-adaptInterval, s-1) ) / (adaptInterval-1);
// //         gamma1[0] = 1 / pow(adapIter[0], c1);
// //         gamma2[0] = c0 * gamma1[0];
// //         sigma2[0] = exp( log(sigma2[0]) + gamma2[0] * (rhat[0] - ropt) );
// //         
// //         COValpha = COValpha + gamma1[0] * ( cov( postAlpha.rows(s+1-adaptInterval, s-1) ) - COValpha );
// //         cholCOValpha = trans( chol( sigma2[0] * ( COValpha + 0.0000000001 * diagmat(ones(n_al)) ) ) );
// //         
// //         adapIter[0] = adapIter[0] + 1;
// //       }
// //     }
// //     
// //     newalpha = alpha + cholCOValpha * randn(n_al);
// //     
// //     newlogZ = cppf_logZ_a0t_etat_ja_pref(Xzoop, newalpha, eta, delta);
// //     newlogZobl = cppf_logZobl(newlogZ, indobl);
// //     newlogZsur = cppf_logZsur(newlogZ, indsur);
// //     newloglam = cppf_loglam_b0t(Xwhale, newlogZ, beta, indtwhale);
// //     newlam = exp(newloglam);
// //     
// //     logprob = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) -
// //       cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //       cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) -
// //       cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
// //       cppf_loglik_obl(logYobl, newlogZobl, sig2_obl) -
// //       cppf_loglik_obl(logYobl, logZobl, sig2_obl) +
// //       cppf_loglik_sur(logYsur, newlogZsur, lam0, lam1, sig2_sur) -
// //       cppf_loglik_sur(logYsur, logZsur, lam0, lam1, sig2_sur) +
// //       MVN_logh(newalpha, zeros(n_al), diagmat(ones(n_al))/100) -
// //       MVN_logh(alpha, zeros(n_al), diagmat(ones(n_al))/100);
// //     
// //     if (log(randu()) < logprob) {
// //       
// //       alpha = newalpha;
// //       logZ = newlogZ;
// //       logZobl = newlogZobl;
// //       logZsur = newlogZsur;
// //       loglam = newloglam;
// //       lam = newlam;
// //       
// //       accprob(s, 0) = 1;
// //     }
// //     postAlpha.row(s) = trans(alpha);
// //     if(s == 0){ Rprintf("sampled alpha\n"); }
// //     
// //     
// //     // Gibbs update for alpha0tilde
// //     alpha0tilde = cppf_gibbs_alpha0tilde(alpha.rows(0, ntzoop-1), tau2);
// //     postAlpha0tilde[s] = alpha0tilde;
// //     if(s == 0){ Rprintf("sampled alpha0tilde\n"); }
// //     
// //     
// //     
// //     // update eta via Elliptical slice sampling
// //     for(int t = 0; t < ntzoop; t ++){
// //       neweta = eta;
// //       
// //       llprev = cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //         cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
// //         cppf_loglik_obl(logYobl, logZobl, sig2_obl) +
// //         cppf_loglik_sur(logYsur, logZsur, lam0, lam1, sig2_sur);
// //       prioretat = sqrt(kappa_eta) * trans(cholReta) * randn(nG);
// //       
// //       thetamin = 0;
// //       thetamax = 2 * M_PI;
// //       theta = thetamin + randu() * (thetamax - thetamin);
// //       thetamin = theta - 2 * M_PI;
// //       thetamax = theta;
// //       
// //       llthred = llprev + log(randu());
// //       accept = false;
// //       count = 0;
// //       
// //       while(accept == false){
// //         count = count + 1;
// //         
// //         neweta.col(t) = eta.col(t) * cos(theta) + prioretat * sin(theta);
// //         neweta.col(t) = neweta.col(t) - mean(neweta.col(t));
// //         
// //         newlogZ = cppf_logZ_a0t_etat_ja_pref(Xzoop, alpha, neweta, delta);
// //         newlogZobl = cppf_logZobl(newlogZ, indobl);
// //         newlogZsur = cppf_logZsur(newlogZ, indsur);
// //         newloglam = cppf_loglam_b0t(Xwhale, newlogZ, beta, indtwhale);
// //         newlam = exp(newloglam);
// //         
// //         llnew = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //           cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) +
// //           cppf_loglik_obl(logYobl, newlogZobl, sig2_obl) +
// //           cppf_loglik_sur(logYsur, newlogZsur, lam0, lam1, sig2_sur);
// //         
// //         if(llnew > llthred){
// //           llprev = llnew;
// //           accept = true;
// //         } else {
// //           if(theta < 0){ thetamin = theta; } else { thetamax = theta; }
// //           theta = thetamin + randu() * (thetamax - thetamin);
// //           if ( (count) % 1000 == 0 ) {
// //             Rprintf("ESS for eta: %d iterations...\n", count);
// //           }
// //           if (count > 20000) stop("ESS is not converging...\n");
// //         }
// //       }
// //       eta = neweta;
// //       logZ = newlogZ;
// //       logZobl = newlogZobl;
// //       logZsur = newlogZsur;
// //       loglam = newloglam;
// //       lam = newlam;
// //       
// //       postEtat = as<mat>(postEta[t]);
// //       postEtat.row(s) = trans(eta.col(t));
// //       postEta[t] = postEtat;
// //     }
// //     if(s == 0){ Rprintf("sampled eta\n"); }
// //     
// //     
// //     // update kappa_eta
// //     kappa_eta = cppf_gibbs_kappa_eta(eta, invReta, 2, 2);
// //     postKappa(s,0) = kappa_eta;
// //     if(s == 0){ Rprintf("sampled kappa_eta\n"); }
// //     
// //     
// //     // update lam0
// //     lam0 = cppf_gibbs_lam0(logYsur, logZsur, lam1, sig2_sur, 100);
// //     postLam(s, 0) = lam0;
// //     if(s == 0){ Rprintf("sampled lam0\n"); }
// //     
// //     
// //     // update lam1
// //     lam1 = cppf_gibbs_lam1(logYsur, logZsur, lam0, sig2_sur, 100);
// //     postLam(s, 1) = lam1;
// //     if(s == 0){ Rprintf("sampled lam1\n"); }
// //     
// //     
// //     // update sig2_sur
// //     sig2_sur = cppf_gibbs_sig2_sur(logYsur, logZsur, lam0, lam1, 2, 2);
// //     postSig2(s, 0) = sig2_sur;
// //     if(s == 0){ Rprintf("sampled sig2_sur\n"); }
// //     
// //     
// //     // update mu0
// //     if( updateCOV ){
// //       if( (s+1 >= adaptInterval) && (s+1 - (adaptInterval * trunc((s+1) / adaptInterval)) == 0) ){
// //         dummyaccprob = accprob.col(2);
// //         rhat[2] = sum( dummyaccprob.rows(s+1-adaptInterval, s-1) ) / (adaptInterval-1);
// //         gamma1[2] = 1 / pow(adapIter[2], c1);
// //         gamma2[2] = c0 * gamma1[2];
// //         sigma2[2] = exp( log(sigma2[2]) + gamma2[2] * (rhat[2] - ropt) );
// //         
// //         postj = postMu0.col(0);
// //         COVmu0 = COVmu0 + gamma1[2] * ( var( postj.rows(s+1-adaptInterval, s-1)) - COVmu0 );
// //         cholCOVmu0 = sqrt( sigma2[2] * COVmu0 );
// //         
// //         adapIter[2] = adapIter[2] + 1;
// //       }
// //     }
// //     
// //     newmu0 = mu0 + cholCOVmu0 * randn();
// //     
// //     newloglam_pref = cppf_loglam_pref(newmu0, delta);
// //     newlam_pref = exp(newloglam_pref);
// //     
// //     logprob = cppf_loglik_pref(newloglam_pref, newlam_pref, indpref, cellarea) - 
// //       cppf_loglik_pref(loglam_pref, lam_pref, indpref, cellarea) +
// //       normal_logh(newmu0, 0, 100) - normal_logh(mu0, 0, 100);
// //     
// //     if (log(randu()) < logprob) {
// //       
// //       mu0 = newmu0;
// //       loglam_pref = newloglam_pref;
// //       lam_pref = newlam_pref;
// //       
// //       accprob(s, 2) = 1;
// //     }
// //     postMu0(s,0) = mu0;
// //     if(s == 0){ Rprintf("sampled mu0\n"); }
// //     
// //     
// //     // update delta
// //     llprev = cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //       cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
// //       cppf_loglik_obl(logYobl, logZobl, sig2_obl) +
// //       cppf_loglik_sur(logYsur, logZsur, lam0, lam1, sig2_sur) +
// //       cppf_loglik_pref(loglam_pref, lam_pref, indpref, cellarea);
// //     priordelta = sqrt(kappa_delta) * trans(cholRdelta) * randn(nG);
// //     
// //     thetamin = 0;
// //     thetamax = 2 * M_PI;
// //     theta = thetamin + randu() * (thetamax - thetamin);
// //     thetamin = theta - 2 * M_PI;
// //     thetamax = theta;
// //     
// //     llthred = llprev + log(randu());
// //     accept = false;
// //     count = 0;
// //     
// //     while(accept == false){
// //       count = count + 1;
// //       
// //       newdelta = delta * cos(theta) + priordelta * sin(theta);
// //       newdelta = newdelta - mean(newdelta);
// //       
// //       newlogZ = cppf_logZ_a0t_etat_ja_pref(Xzoop, alpha, eta, newdelta);
// //       newlogZobl = cppf_logZobl(newlogZ, indobl);
// //       newlogZsur = cppf_logZsur(newlogZ, indsur);
// //       newloglam = cppf_loglam_b0t(Xwhale, newlogZ, beta, indtwhale);
// //       newlam = exp(newloglam);
// //       newloglam_pref = cppf_loglam_pref(mu0, newdelta);
// //       newlam_pref = exp(newloglam_pref);
// //       
// //       llnew = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //         cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) +
// //         cppf_loglik_obl(logYobl, newlogZobl, sig2_obl) +
// //         cppf_loglik_sur(logYsur, newlogZsur, lam0, lam1, sig2_sur) +
// //         cppf_loglik_pref(newloglam_pref, newlam_pref, indpref, cellarea);
// //       
// //       if(llnew > llthred){
// //         llprev = llnew;
// //         accept = true;
// //       } else {
// //         if(theta < 0){ thetamin = theta; } else { thetamax = theta; }
// //         theta = thetamin + randu() * (thetamax - thetamin);
// //         if ( (count) % 1000 == 0 ) {
// //           Rprintf("ESS for delta: %d iterations...\n", count);
// //         }
// //         if (count > 20000) stop("ESS is not converging...\n");
// //       }
// //     }
// //     delta = newdelta;
// //     logZ = newlogZ;
// //     logZobl = newlogZobl;
// //     logZsur = newlogZsur;
// //     loglam = newloglam;
// //     lam = newlam;
// //     loglam_pref = newloglam_pref;
// //     lam_pref = newlam_pref;
// //     
// //     postDelta.row(s) = trans(delta);
// //     if(s == 0){ Rprintf("sampled delta\n"); }
// //     
// //     
// //     // update kappa_delta
// //     kappa_delta = cppf_gibbs_kappa_delta(delta, invRdelta, 2, 2);
// //     postKappa(s, 1) = kappa_delta;
// //     if(s == 0){ Rprintf("sampled kappa_delta\n"); }
// //     
// //     
// //     // update beta
// //     if( updateCOV ){
// //       if( (s+1 >= adaptInterval) && (s+1 - (adaptInterval * trunc((s+1) / adaptInterval)) == 0) ){
// //         dummyaccprob = accprob.col(1);
// //         rhat[1] = sum( dummyaccprob.rows(s+1-adaptInterval, s-1) ) / (adaptInterval-1);
// //         gamma1[1] = 1 / pow(adapIter[1], c1);
// //         gamma2[1] = c0 * gamma1[1];
// //         sigma2[1] = exp( log(sigma2[1]) + gamma2[1] * (rhat[1] - ropt) );
// //         
// //         COVbeta = COVbeta + gamma1[1] * ( cov( postBeta.rows(s+1-adaptInterval, s-1)) - COVbeta );
// //         cholCOVbeta = trans( chol( sigma2[1] * ( COVbeta + 0.0000000001 * diagmat(ones(n_be)) ) ) );
// //         
// //         adapIter[1] = adapIter[1] + 1;
// //       }
// //     }
// //     
// //     newbeta = beta + cholCOVbeta * randn(n_be);
// //     
// //     newloglam = cppf_loglam_b0t(Xwhale, logZ, newbeta, indtwhale);
// //     newlam = exp(newloglam);
// //     
// //     logprob = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) -
// //       cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //       cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) -
// //       cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
// //       MVN_logh(newbeta, zeros(n_be), diagmat(ones(n_be))/100) -
// //       MVN_logh(beta, zeros(n_be), diagmat(ones(n_be))/100);
// //     
// //     if (log(randu()) < logprob) {
// //       
// //       beta = newbeta;
// //       loglam = newloglam;
// //       lam = newlam;
// //       
// //       accprob(s, 1) = 1;
// //     }
// //     
// //     postBeta.row(s) = trans(beta);
// //     if(s == 0){ Rprintf("sampled beta\n"); }
// //     
// //     
// //     if ( (s+1) % 10 == 0 ) {
// //       Rprintf("Generated %d samples...\n", s+1);
// //     }
// //   }
// //   
// //   
// //   return Rcpp::List::create(Rcpp::Named("postAlpha0tilde") = postAlpha0tilde,
// //                             Rcpp::Named("postAlpha") = postAlpha,
// //                             Rcpp::Named("postEta") = postEta,
// //                             Rcpp::Named("postLam") = postLam,
// //                             Rcpp::Named("postKappa") = postKappa,
// //                             Rcpp::Named("postSig2") = postSig2,
// //                             Rcpp::Named("postMu0") = postMu0,
// //                             Rcpp::Named("postDelta") = postDelta,
// //                             Rcpp::Named("postBeta") = postBeta,
// //                             Rcpp::Named("Accprob") = accprob,
// //                             Rcpp::Named("sigma2") = sigma2,
// //                             Rcpp::Named("adapIter") = adapIter,
// //                             Rcpp::Named("COValpha") = COValpha,
// //                             Rcpp::Named("COVmu0") = COVmu0,
// //                             Rcpp::Named("COVbeta") = COVbeta);
// //   
// // }
// // 
// // 
// // 
// // // ============================================================================-
// // // Modularized model
// // // Fix daily call rates at truth
// // // Preferential sampling
// // // ============================================================================-
// // 
// // // pref: preferential sampling for zoop
// // // b0t: daily intercepts for whale
// // // psit: daily GPs for whale
// // // ja: joint updating for alpha_0 and alpha
// // // fixC: fix daily call rates at truth
// // 
// // // [[Rcpp::export]]
// // List cppf_fitMod_pref_b0t_psit_ja_fixC(
// //     int niter, mat distmat,
// //     List Xzoop, vec logYobl, vec logYsur, mat Spref,
// //     List Xwhale, mat Sdist, mat Wpam, vec dur_dist,
// //     vec f_dist_obs, List f_dist, List p_pam,
// //     mat indobl, mat indsur, vec indpref, vec indtwhale, mat inddist, vec indmdist, vec indmwhale,
// //     double alpha0tilde, vec alpha, mat eta, 
// //     double mu0, vec delta, double lam0, double lam1,
// //     double tau2, double kappa_eta, double phi_eta, double kappa_delta, double phi_delta,
// //     double sig2_obl, double sig2_sur,
// //     vec beta, mat psi, vec pii, vec c,
// //     double kappa_psi, double phi_psi,
// //     double cellarea,
// //     vec sigma2, mat COValpha, mat COVbeta, double COVmu0,
// //     bool updateCOV, int adaptInterval, double adaptFactorExponent, vec adapIter){
// //   
// //   double negativeInf = -std::numeric_limits<float>::infinity();;
// //   
// //   indobl = indobl-1, indsur = indsur-1, indpref = indpref-1, indtwhale = indtwhale-1;
// //   inddist = inddist-1, indmdist = indmdist-1, indmwhale = indmwhale-1;
// //   
// //   int nG = distmat.n_rows, n_al = alpha.size(), n_be = beta.size(), ntzoop = eta.n_cols, ntwhale = psi.n_cols, n_sigma2 = sigma2.size();
// //   
// //   mat postAlpha0tilde = zeros(niter, 1), postAlpha = zeros(niter, n_al), postLam = zeros(niter, 2);
// //   mat postDelta = zeros(niter, delta.size()), postMu0(niter, 1);
// //   mat postSig2 = zeros(niter, 1), postKappa = zeros(niter, 3);
// //   mat postBeta = zeros(niter, n_be);
// //   mat postEtat = zeros(niter, nG), postPsit = zeros(niter, nG);
// //   List postEta(ntzoop), postPsi(ntwhale);
// //   for(int i = 0; i < ntzoop; i ++){ postEta[i] = zeros(niter, nG); }
// //   for(int i = 0; i < ntwhale; i ++){ postPsi[i] = zeros(niter, nG); }
// //   Rprintf("set up matrices for posterior samples\n");
// //   
// //   mat accprob = zeros(niter, n_sigma2);
// //   vec rhat = zeros(n_sigma2), gamma1 = zeros(n_sigma2), gamma2 = zeros(n_sigma2), dummyaccprob, postj;
// //   double c0 = 1, c1 = adaptFactorExponent, ropt = 0.234, logprob;
// //   mat cholCOValpha = trans( chol( sigma2[0] * ( COValpha + 0.0000000001 * diagmat(ones(n_al)) ) ) );
// //   mat cholCOVbeta = trans( chol( sigma2[1] * ( COVbeta + 0.0000000001 * diagmat(ones(n_be)) ) ) );
// //   double cholCOVmu0 = sqrt( sigma2[2] * COVmu0 );
// //   Rprintf("set up adaptive proposals\n");
// //   
// //   double newmu0;
// //   vec newalpha, newbeta, newlogZobl, newlogZsur, newloglam_pref, newlam_pref;
// //   mat newlogZ, newloglam, newlam;
// //   vec dummyvec;
// //   
// //   int count;
// //   vec xvals;
// //   double llprev, llthred, llnew, thetamin, thetamax, theta;
// //   vec prioretat, priorpsit, priordelta, newdelta;
// //   mat neweta, newpsi;
// //   bool accept;
// //   
// //   mat Reta = exp(-distmat / phi_eta), Rpsi = exp(-distmat / phi_psi), Rdelta = exp(-distmat / phi_delta);
// //   mat invReta = inv(Reta), invRpsi= inv(Rpsi), invRdelta = inv(Rdelta);
// //   mat cholReta = chol(Reta), cholRpsi = chol(Rpsi), cholRdelta = chol(Rdelta);
// //   rowvec oneinvReta = trans(ones(nG)) * invReta, oneinvRpsi = trans(ones(nG)) * invRpsi, oneinvRdelta = trans(ones(nG)) * invRdelta;
// //   double oneinvRetaone = sum( oneinvReta * ones(nG) ), oneinvRpsione = sum( oneinvRpsi * ones(nG) ), oneinvRdeltaone = sum( oneinvRdelta * ones(nG) );
// //   Rprintf("set up GP covariance matrices\n");
// //   
// //   mat logZ = cppf_logZ_a0t_etat_ja_pref(Xzoop, alpha, eta, delta);
// //   vec logZobl = cppf_logZobl(logZ, indobl);
// //   vec logZsur = cppf_logZsur(logZ, indsur);
// //   mat loglam = cppf_loglam_b0t_psit(Xwhale, logZ, beta, psi, indtwhale);
// //   mat lam = exp(loglam);
// //   vec loglam_pref = cppf_loglam_pref(mu0, delta);
// //   vec lam_pref = exp(loglam_pref);
// //   Rprintf("set up initial values\n");
// //   
// //   
// //   // Start MCMC
// //   for(int s = 0; s < niter; s++) {
// //     
// //     // update alpha
// //     if( updateCOV ){
// //       if( (s+1 >= adaptInterval) && (s+1 - (adaptInterval * trunc((s+1) / adaptInterval)) == 0) ){
// //         dummyaccprob = accprob.col(0);
// //         rhat[0] = sum( dummyaccprob.rows(s+1-adaptInterval, s-1) ) / (adaptInterval-1);
// //         gamma1[0] = 1 / pow(adapIter[0], c1);
// //         gamma2[0] = c0 * gamma1[0];
// //         sigma2[0] = exp( log(sigma2[0]) + gamma2[0] * (rhat[0] - ropt) );
// //         
// //         COValpha = COValpha + gamma1[0] * ( cov( postAlpha.rows(s+1-adaptInterval, s-1) ) - COValpha );
// //         cholCOValpha = trans( chol( sigma2[0] * ( COValpha + 0.0000000001 * diagmat(ones(n_al)) ) ) );
// //         
// //         adapIter[0] = adapIter[0] + 1;
// //       }
// //     }
// //     
// //     newalpha = alpha + cholCOValpha * randn(n_al);
// //     
// //     newlogZ = cppf_logZ_a0t_etat_ja_pref(Xzoop, newalpha, eta, delta);
// //     newlogZobl = cppf_logZobl(newlogZ, indobl);
// //     newlogZsur = cppf_logZsur(newlogZ, indsur);
// //     newloglam = cppf_loglam_b0t_psit(Xwhale, newlogZ, beta, psi, indtwhale);
// //     newlam = exp(newloglam);
// //     
// //     logprob = cppf_loglik_obl(logYobl, newlogZobl, sig2_obl) -
// //       cppf_loglik_obl(logYobl, logZobl, sig2_obl) +
// //       cppf_loglik_sur(logYsur, newlogZsur, lam0, lam1, sig2_sur) -
// //       cppf_loglik_sur(logYsur, logZsur, lam0, lam1, sig2_sur) +
// //       MVN_logh(newalpha, zeros(n_al), diagmat(ones(n_al))/100) -
// //       MVN_logh(alpha, zeros(n_al), diagmat(ones(n_al))/100);
// //     
// //     if (log(randu()) < logprob) {
// //       
// //       alpha = newalpha;
// //       logZ = newlogZ;
// //       logZobl = newlogZobl;
// //       logZsur = newlogZsur;
// //       loglam = newloglam;
// //       lam = newlam;
// //       
// //       accprob(s, 0) = 1;
// //     }
// //     postAlpha.row(s) = trans(alpha);
// //     if(s == 0){ Rprintf("sampled alpha\n"); }
// //     
// //     
// //     // Gibbs update for alpha0tilde
// //     alpha0tilde = cppf_gibbs_alpha0tilde(alpha.rows(0, ntzoop-1), tau2);
// //     postAlpha0tilde[s] = alpha0tilde;
// //     if(s == 0){ Rprintf("sampled alpha0tilde\n"); }
// //     
// //     
// //     
// //     // update eta via Elliptical slice sampling
// //     for(int t = 0; t < ntzoop; t ++){
// //       neweta = eta;
// //       
// //       llprev = cppf_loglik_obl(logYobl, logZobl, sig2_obl) +
// //         cppf_loglik_sur(logYsur, logZsur, lam0, lam1, sig2_sur);
// //       prioretat = sqrt(kappa_eta) * trans(cholReta) * randn(nG);
// //       
// //       thetamin = 0;
// //       thetamax = 2 * M_PI;
// //       theta = thetamin + randu() * (thetamax - thetamin);
// //       thetamin = theta - 2 * M_PI;
// //       thetamax = theta;
// //       
// //       llthred = llprev + log(randu());
// //       accept = false;
// //       count = 0;
// //       
// //       while(accept == false){
// //         count = count + 1;
// //         
// //         neweta.col(t) = eta.col(t) * cos(theta) + prioretat * sin(theta);
// //         neweta.col(t) = neweta.col(t) - mean(neweta.col(t));
// //         
// //         newlogZ = cppf_logZ_a0t_etat_ja_pref(Xzoop, alpha, neweta, delta);
// //         newlogZobl = cppf_logZobl(newlogZ, indobl);
// //         newlogZsur = cppf_logZsur(newlogZ, indsur);
// //         newloglam = cppf_loglam_b0t_psit(Xwhale, newlogZ, beta, psi, indtwhale);
// //         newlam = exp(newloglam);
// //         
// //         llnew = cppf_loglik_obl(logYobl, newlogZobl, sig2_obl) +
// //           cppf_loglik_sur(logYsur, newlogZsur, lam0, lam1, sig2_sur);
// //         
// //         if(llnew > llthred){
// //           llprev = llnew;
// //           accept = true;
// //         } else {
// //           if(theta < 0){ thetamin = theta; } else { thetamax = theta; }
// //           theta = thetamin + randu() * (thetamax - thetamin);
// //           if ( (count) % 1000 == 0 ) {
// //             Rprintf("ESS for eta: %d iterations...\n", count);
// //           }
// //           if (count > 20000) stop("ESS is not converging...\n");
// //         }
// //       }
// //       eta = neweta;
// //       logZ = newlogZ;
// //       logZobl = newlogZobl;
// //       logZsur = newlogZsur;
// //       loglam = newloglam;
// //       lam = newlam;
// //       
// //       postEtat = as<mat>(postEta[t]);
// //       postEtat.row(s) = trans(eta.col(t));
// //       postEta[t] = postEtat;
// //     }
// //     if(s == 0){ Rprintf("sampled eta\n"); }
// //     
// //     
// //     // update kappa_eta
// //     kappa_eta = cppf_gibbs_kappa_eta(eta, invReta, 2, 2);
// //     postKappa(s,0) = kappa_eta;
// //     if(s == 0){ Rprintf("sampled kappa_eta\n"); }
// //     
// //     
// //     // update lam0
// //     lam0 = cppf_gibbs_lam0(logYsur, logZsur, lam1, sig2_sur, 100);
// //     postLam(s, 0) = lam0;
// //     if(s == 0){ Rprintf("sampled lam0\n"); }
// //     
// //     
// //     // update lam1
// //     lam1 = cppf_gibbs_lam1(logYsur, logZsur, lam0, sig2_sur, 100);
// //     postLam(s, 1) = lam1;
// //     if(s == 0){ Rprintf("sampled lam1\n"); }
// //     
// //     
// //     // update sig2_sur
// //     sig2_sur = cppf_gibbs_sig2_sur(logYsur, logZsur, lam0, lam1, 2, 2);
// //     postSig2(s, 0) = sig2_sur;
// //     if(s == 0){ Rprintf("sampled sig2_sur\n"); }
// //     
// //     
// //     // update mu0
// //     if( updateCOV ){
// //       if( (s+1 >= adaptInterval) && (s+1 - (adaptInterval * trunc((s+1) / adaptInterval)) == 0) ){
// //         dummyaccprob = accprob.col(2);
// //         rhat[2] = sum( dummyaccprob.rows(s+1-adaptInterval, s-1) ) / (adaptInterval-1);
// //         gamma1[2] = 1 / pow(adapIter[2], c1);
// //         gamma2[2] = c0 * gamma1[2];
// //         sigma2[2] = exp( log(sigma2[2]) + gamma2[2] * (rhat[2] - ropt) );
// //         
// //         postj = postMu0.col(0);
// //         COVmu0 = COVmu0 + gamma1[2] * ( var( postj.rows(s+1-adaptInterval, s-1)) - COVmu0 );
// //         cholCOVmu0 = sqrt( sigma2[2] * COVmu0 );
// //         
// //         adapIter[2] = adapIter[2] + 1;
// //       }
// //     }
// //     
// //     newmu0 = mu0 + cholCOVmu0 * randn();
// //     
// //     newloglam_pref = cppf_loglam_pref(newmu0, delta);
// //     newlam_pref = exp(newloglam_pref);
// //     
// //     logprob = cppf_loglik_pref(newloglam_pref, newlam_pref, indpref, cellarea) - 
// //       cppf_loglik_pref(loglam_pref, lam_pref, indpref, cellarea) +
// //       normal_logh(newmu0, 0, 100) - normal_logh(mu0, 0, 100);
// //     
// //     if (log(randu()) < logprob) {
// //       
// //       mu0 = newmu0;
// //       loglam_pref = newloglam_pref;
// //       lam_pref = newlam_pref;
// //       
// //       accprob(s, 2) = 1;
// //     }
// //     postMu0(s,0) = mu0;
// //     if(s == 0){ Rprintf("sampled mu0\n"); }
// //     
// //     
// //     // update delta
// //     llprev = cppf_loglik_obl(logYobl, logZobl, sig2_obl) +
// //       cppf_loglik_sur(logYsur, logZsur, lam0, lam1, sig2_sur) +
// //       cppf_loglik_pref(loglam_pref, lam_pref, indpref, cellarea);
// //     priordelta = sqrt(kappa_delta) * trans(cholRdelta) * randn(nG);
// //     
// //     thetamin = 0;
// //     thetamax = 2 * M_PI;
// //     theta = thetamin + randu() * (thetamax - thetamin);
// //     thetamin = theta - 2 * M_PI;
// //     thetamax = theta;
// //     
// //     llthred = llprev + log(randu());
// //     accept = false;
// //     count = 0;
// //     
// //     while(accept == false){
// //       count = count + 1;
// //       
// //       newdelta = delta * cos(theta) + priordelta * sin(theta);
// //       newdelta = newdelta - mean(newdelta);
// //       
// //       newlogZ = cppf_logZ_a0t_etat_ja_pref(Xzoop, alpha, eta, newdelta);
// //       newlogZobl = cppf_logZobl(newlogZ, indobl);
// //       newlogZsur = cppf_logZsur(newlogZ, indsur);
// //       newloglam = cppf_loglam_b0t_psit(Xwhale, newlogZ, beta, psi, indtwhale);
// //       newlam = exp(newloglam);
// //       newloglam_pref = cppf_loglam_pref(mu0, newdelta);
// //       newlam_pref = exp(newloglam_pref);
// //       
// //       llnew = cppf_loglik_obl(logYobl, newlogZobl, sig2_obl) +
// //         cppf_loglik_sur(logYsur, newlogZsur, lam0, lam1, sig2_sur) +
// //         cppf_loglik_pref(newloglam_pref, newlam_pref, indpref, cellarea);
// //       
// //       if(llnew > llthred){
// //         llprev = llnew;
// //         accept = true;
// //       } else {
// //         if(theta < 0){ thetamin = theta; } else { thetamax = theta; }
// //         theta = thetamin + randu() * (thetamax - thetamin);
// //         if ( (count) % 1000 == 0 ) {
// //           Rprintf("ESS for delta: %d iterations...\n", count);
// //         }
// //         if (count > 20000) stop("ESS is not converging...\n");
// //       }
// //     }
// //     delta = newdelta;
// //     logZ = newlogZ;
// //     logZobl = newlogZobl;
// //     logZsur = newlogZsur;
// //     loglam = newloglam;
// //     lam = newlam;
// //     loglam_pref = newloglam_pref;
// //     lam_pref = newlam_pref;
// //     
// //     postDelta.row(s) = trans(delta);
// //     if(s == 0){ Rprintf("sampled delta\n"); }
// //     
// //     
// //     // update kappa_delta
// //     kappa_delta = cppf_gibbs_kappa_delta(delta, invRdelta, 2, 2);
// //     postKappa(s, 2) = kappa_delta;
// //     if(s == 0){ Rprintf("sampled kappa_delta\n"); }
// //     
// //     
// //     // update beta
// //     if( updateCOV ){
// //       if( (s+1 >= adaptInterval) && (s+1 - (adaptInterval * trunc((s+1) / adaptInterval)) == 0) ){
// //         dummyaccprob = accprob.col(1);
// //         rhat[1] = sum( dummyaccprob.rows(s+1-adaptInterval, s-1) ) / (adaptInterval-1);
// //         gamma1[1] = 1 / pow(adapIter[1], c1);
// //         gamma2[1] = c0 * gamma1[1];
// //         sigma2[1] = exp( log(sigma2[1]) + gamma2[1] * (rhat[1] - ropt) );
// //         
// //         COVbeta = COVbeta + gamma1[1] * ( cov( postBeta.rows(s+1-adaptInterval, s-1)) - COVbeta );
// //         cholCOVbeta = trans( chol( sigma2[1] * ( COVbeta + 0.0000000001 * diagmat(ones(n_be)) ) ) );
// //         
// //         adapIter[1] = adapIter[1] + 1;
// //       }
// //     }
// //     
// //     newbeta = beta + cholCOVbeta * randn(n_be);
// //     
// //     newloglam = cppf_loglam_b0t_psit(Xwhale, logZ, newbeta, psi, indtwhale);
// //     newlam = exp(newloglam);
// //     
// //     logprob = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) -
// //       cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //       cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) -
// //       cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
// //       MVN_logh(newbeta, zeros(n_be), diagmat(ones(n_be))/100) -
// //       MVN_logh(beta, zeros(n_be), diagmat(ones(n_be))/100);
// //     
// //     if (log(randu()) < logprob) {
// //       
// //       beta = newbeta;
// //       loglam = newloglam;
// //       lam = newlam;
// //       
// //       accprob(s, 1) = 1;
// //     }
// //     
// //     postBeta.row(s) = trans(beta);
// //     if(s == 0){ Rprintf("sampled beta\n"); }
// //     
// //     
// //     
// //     // update psi via Elliptical slice sampling
// //     for(int t = 0; t < ntwhale; t ++){
// //       newpsi = psi;
// //       
// //       llprev = cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //         cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea);
// //       priorpsit = sqrt(kappa_psi) * trans(cholRpsi) * randn(nG);
// //       
// //       thetamin = 0;
// //       thetamax = 2 * M_PI;
// //       theta = thetamin + randu() * (thetamax - thetamin);
// //       thetamin = theta - 2 * M_PI;
// //       thetamax = theta;
// //       
// //       llthred = llprev + log(randu());
// //       accept = false;
// //       count = 0;
// //       
// //       while(accept == false){
// //         count = count + 1;
// //         
// //         newpsi.col(t) = psi.col(t) * cos(theta) + priorpsit * sin(theta);
// //         newpsi.col(t) = newpsi.col(t) - mean(newpsi.col(t));
// //         
// //         newloglam = cppf_loglam_b0t_psit(Xwhale, logZ, beta, newpsi, indtwhale);
// //         newlam = exp(newloglam);
// //         
// //         llnew = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //           cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea);
// //         
// //         if(llnew > llthred){
// //           llprev = llnew;
// //           accept = true;
// //         } else {
// //           if(theta < 0){ thetamin = theta; } else { thetamax = theta; }
// //           theta = thetamin + randu() * (thetamax - thetamin);
// //           if ( (count) % 1000 == 0 ) {
// //             Rprintf("ESS for psi: %d iterations...\n", count);
// //           }
// //           if (count > 20000) stop("ESS is not converging...\n");
// //         }
// //       }
// //       psi = newpsi;
// //       loglam = newloglam;
// //       lam = newlam;
// //       
// //       postPsit = as<mat>(postPsi[t]);
// //       postPsit.row(s) = trans(psi.col(t));
// //       postPsi[t] = postPsit;
// //     }
// //     if(s == 0){ Rprintf("sampled psi\n"); }
// //     
// //     
// //     // update kappa_psi
// //     kappa_psi = cppf_gibbs_kappa_psit(psi, invRpsi, 2, 2);
// //     postKappa(s,1) = kappa_psi;
// //     if(s == 0){ Rprintf("sampled kappa_psi\n"); }
// //     
// //     
// //     if ( (s+1) % 10 == 0 ) {
// //       Rprintf("Generated %d samples...\n", s+1);
// //     }
// //   }
// //   
// //   
// //   return Rcpp::List::create(Rcpp::Named("postAlpha0tilde") = postAlpha0tilde,
// //                             Rcpp::Named("postAlpha") = postAlpha,
// //                             Rcpp::Named("postEta") = postEta,
// //                             Rcpp::Named("postLam") = postLam,
// //                             Rcpp::Named("postKappa") = postKappa,
// //                             Rcpp::Named("postSig2") = postSig2,
// //                             Rcpp::Named("postMu0") = postMu0,
// //                             Rcpp::Named("postDelta") = postDelta,
// //                             Rcpp::Named("postBeta") = postBeta,
// //                             Rcpp::Named("postPsi") = postPsi,
// //                             Rcpp::Named("Accprob") = accprob,
// //                             Rcpp::Named("sigma2") = sigma2,
// //                             Rcpp::Named("adapIter") = adapIter,
// //                             Rcpp::Named("COValpha") = COValpha,
// //                             Rcpp::Named("COVmu0") = COVmu0,
// //                             Rcpp::Named("COVbeta") = COVbeta);
// //   
// // }
// // 
// // 
// // 
// // // [[Rcpp::export]]
// // List cppf_fitMod_pref_b0t_psi_ja_fixC(
// //     int niter, mat distmat,
// //     List Xzoop, vec logYobl, vec logYsur, mat Spref,
// //     List Xwhale, mat Sdist, mat Wpam, vec dur_dist,
// //     vec f_dist_obs, List f_dist, List p_pam,
// //     mat indobl, mat indsur, vec indpref, vec indtwhale, mat inddist, vec indmdist, vec indmwhale,
// //     double alpha0tilde, vec alpha, mat eta, 
// //     double mu0, vec delta, double lam0, double lam1,
// //     double tau2, double kappa_eta, double phi_eta, double kappa_delta, double phi_delta,
// //     double sig2_obl, double sig2_sur,
// //     vec beta, vec psi, vec pii, vec c,
// //     double kappa_psi, double phi_psi,
// //     double cellarea,
// //     vec sigma2, mat COValpha, mat COVbeta, double COVmu0,
// //     bool updateCOV, int adaptInterval, double adaptFactorExponent, vec adapIter){
// //   
// //   double negativeInf = -std::numeric_limits<float>::infinity();;
// //   
// //   indobl = indobl-1, indsur = indsur-1, indpref = indpref-1, indtwhale = indtwhale-1;
// //   inddist = inddist-1, indmdist = indmdist-1, indmwhale = indmwhale-1;
// //   
// //   int nG = distmat.n_rows, n_al = alpha.size(), n_be = beta.size(), ntzoop = Xzoop.size(), ntwhale = Xwhale.size(), n_sigma2 = sigma2.size();
// //   
// //   mat postAlpha0tilde = zeros(niter, 1), postAlpha = zeros(niter, n_al), postLam = zeros(niter, 2);
// //   mat postDelta = zeros(niter, delta.size()), postMu0(niter, 1);
// //   mat postSig2 = zeros(niter, 1), postKappa = zeros(niter, 3);
// //   mat postBeta = zeros(niter, n_be);
// //   mat postEtat = zeros(niter, nG), postPsi = zeros(niter, nG);
// //   List postEta(ntzoop);
// //   for(int i = 0; i < ntzoop; i ++){ postEta[i] = zeros(niter, nG); }
// //   Rprintf("set up matrices for posterior samples\n");
// //   
// //   mat accprob = zeros(niter, n_sigma2);
// //   vec rhat = zeros(n_sigma2), gamma1 = zeros(n_sigma2), gamma2 = zeros(n_sigma2), dummyaccprob, postj;
// //   double c0 = 1, c1 = adaptFactorExponent, ropt = 0.234, logprob;
// //   mat cholCOValpha = trans( chol( sigma2[0] * ( COValpha + 0.0000000001 * diagmat(ones(n_al)) ) ) );
// //   mat cholCOVbeta = trans( chol( sigma2[1] * ( COVbeta + 0.0000000001 * diagmat(ones(n_be)) ) ) );
// //   double cholCOVmu0 = sqrt( sigma2[2] * COVmu0 );
// //   Rprintf("set up adaptive proposals\n");
// //   
// //   double newmu0;
// //   vec newalpha, newbeta, newlogZobl, newlogZsur, newloglam_pref, newlam_pref;
// //   mat newlogZ, newloglam, newlam;
// //   vec dummyvec;
// //   
// //   int count;
// //   vec xvals;
// //   double llprev, llthred, llnew, thetamin, thetamax, theta;
// //   vec prioretat, priorpsi, newpsi, priordelta, newdelta;
// //   mat neweta;
// //   bool accept;
// //   
// //   mat Reta = exp(-distmat / phi_eta), Rpsi = exp(-distmat / phi_psi), Rdelta = exp(-distmat / phi_delta);
// //   mat invReta = inv(Reta), invRpsi= inv(Rpsi), invRdelta = inv(Rdelta);
// //   mat cholReta = chol(Reta), cholRpsi = chol(Rpsi), cholRdelta = chol(Rdelta);
// //   rowvec oneinvReta = trans(ones(nG)) * invReta, oneinvRpsi = trans(ones(nG)) * invRpsi, oneinvRdelta = trans(ones(nG)) * invRdelta;
// //   double oneinvRetaone = sum( oneinvReta * ones(nG) ), oneinvRpsione = sum( oneinvRpsi * ones(nG) ), oneinvRdeltaone = sum( oneinvRdelta * ones(nG) );
// //   Rprintf("set up GP covariance matrices\n");
// //   
// //   mat logZ = cppf_logZ_a0t_etat_ja_pref(Xzoop, alpha, eta, delta);
// //   vec logZobl = cppf_logZobl(logZ, indobl);
// //   vec logZsur = cppf_logZsur(logZ, indsur);
// //   mat loglam = cppf_loglam_b0t_psi(Xwhale, logZ, beta, psi, indtwhale);
// //   mat lam = exp(loglam);
// //   vec loglam_pref = cppf_loglam_pref(mu0, delta);
// //   vec lam_pref = exp(loglam_pref);
// //   Rprintf("set up initial values\n");
// //   
// //   
// //   // Start MCMC
// //   for(int s = 0; s < niter; s++) {
// //     
// //     // update alpha
// //     if( updateCOV ){
// //       if( (s+1 >= adaptInterval) && (s+1 - (adaptInterval * trunc((s+1) / adaptInterval)) == 0) ){
// //         dummyaccprob = accprob.col(0);
// //         rhat[0] = sum( dummyaccprob.rows(s+1-adaptInterval, s-1) ) / (adaptInterval-1);
// //         gamma1[0] = 1 / pow(adapIter[0], c1);
// //         gamma2[0] = c0 * gamma1[0];
// //         sigma2[0] = exp( log(sigma2[0]) + gamma2[0] * (rhat[0] - ropt) );
// //         
// //         COValpha = COValpha + gamma1[0] * ( cov( postAlpha.rows(s+1-adaptInterval, s-1) ) - COValpha );
// //         cholCOValpha = trans( chol( sigma2[0] * ( COValpha + 0.0000000001 * diagmat(ones(n_al)) ) ) );
// //         
// //         adapIter[0] = adapIter[0] + 1;
// //       }
// //     }
// //     
// //     newalpha = alpha + cholCOValpha * randn(n_al);
// //     
// //     newlogZ = cppf_logZ_a0t_etat_ja_pref(Xzoop, newalpha, eta, delta);
// //     newlogZobl = cppf_logZobl(newlogZ, indobl);
// //     newlogZsur = cppf_logZsur(newlogZ, indsur);
// //     newloglam = cppf_loglam_b0t_psi(Xwhale, newlogZ, beta, psi, indtwhale);
// //     newlam = exp(newloglam);
// //     
// //     logprob = cppf_loglik_obl(logYobl, newlogZobl, sig2_obl) -
// //       cppf_loglik_obl(logYobl, logZobl, sig2_obl) +
// //       cppf_loglik_sur(logYsur, newlogZsur, lam0, lam1, sig2_sur) -
// //       cppf_loglik_sur(logYsur, logZsur, lam0, lam1, sig2_sur) +
// //       MVN_logh(newalpha, zeros(n_al), diagmat(ones(n_al))/100) -
// //       MVN_logh(alpha, zeros(n_al), diagmat(ones(n_al))/100);
// //     
// //     if (log(randu()) < logprob) {
// //       
// //       alpha = newalpha;
// //       logZ = newlogZ;
// //       logZobl = newlogZobl;
// //       logZsur = newlogZsur;
// //       loglam = newloglam;
// //       lam = newlam;
// //       
// //       accprob(s, 0) = 1;
// //     }
// //     postAlpha.row(s) = trans(alpha);
// //     if(s == 0){ Rprintf("sampled alpha\n"); }
// //     
// //     
// //     // Gibbs update for alpha0tilde
// //     alpha0tilde = cppf_gibbs_alpha0tilde(alpha.rows(0, ntzoop-1), tau2);
// //     postAlpha0tilde[s] = alpha0tilde;
// //     if(s == 0){ Rprintf("sampled alpha0tilde\n"); }
// //     
// //     
// //     
// //     // update eta via Elliptical slice sampling
// //     for(int t = 0; t < ntzoop; t ++){
// //       neweta = eta;
// //       
// //       llprev = cppf_loglik_obl(logYobl, logZobl, sig2_obl) +
// //         cppf_loglik_sur(logYsur, logZsur, lam0, lam1, sig2_sur);
// //       prioretat = sqrt(kappa_eta) * trans(cholReta) * randn(nG);
// //       
// //       thetamin = 0;
// //       thetamax = 2 * M_PI;
// //       theta = thetamin + randu() * (thetamax - thetamin);
// //       thetamin = theta - 2 * M_PI;
// //       thetamax = theta;
// //       
// //       llthred = llprev + log(randu());
// //       accept = false;
// //       count = 0;
// //       
// //       while(accept == false){
// //         count = count + 1;
// //         
// //         neweta.col(t) = eta.col(t) * cos(theta) + prioretat * sin(theta);
// //         neweta.col(t) = neweta.col(t) - mean(neweta.col(t));
// //         
// //         newlogZ = cppf_logZ_a0t_etat_ja_pref(Xzoop, alpha, neweta, delta);
// //         newlogZobl = cppf_logZobl(newlogZ, indobl);
// //         newlogZsur = cppf_logZsur(newlogZ, indsur);
// //         newloglam = cppf_loglam_b0t_psi(Xwhale, newlogZ, beta, psi, indtwhale);
// //         newlam = exp(newloglam);
// //         
// //         llnew = cppf_loglik_obl(logYobl, newlogZobl, sig2_obl) +
// //           cppf_loglik_sur(logYsur, newlogZsur, lam0, lam1, sig2_sur);
// //         
// //         if(llnew > llthred){
// //           llprev = llnew;
// //           accept = true;
// //         } else {
// //           if(theta < 0){ thetamin = theta; } else { thetamax = theta; }
// //           theta = thetamin + randu() * (thetamax - thetamin);
// //           if ( (count) % 1000 == 0 ) {
// //             Rprintf("ESS for eta: %d iterations...\n", count);
// //           }
// //           if (count > 20000) stop("ESS is not converging...\n");
// //         }
// //       }
// //       eta = neweta;
// //       logZ = newlogZ;
// //       logZobl = newlogZobl;
// //       logZsur = newlogZsur;
// //       loglam = newloglam;
// //       lam = newlam;
// //       
// //       postEtat = as<mat>(postEta[t]);
// //       postEtat.row(s) = trans(eta.col(t));
// //       postEta[t] = postEtat;
// //     }
// //     if(s == 0){ Rprintf("sampled eta\n"); }
// //     
// //     
// //     // update kappa_eta
// //     kappa_eta = cppf_gibbs_kappa_eta(eta, invReta, 2, 2);
// //     postKappa(s,0) = kappa_eta;
// //     if(s == 0){ Rprintf("sampled kappa_eta\n"); }
// //     
// //     
// //     // update lam0
// //     lam0 = cppf_gibbs_lam0(logYsur, logZsur, lam1, sig2_sur, 100);
// //     postLam(s, 0) = lam0;
// //     if(s == 0){ Rprintf("sampled lam0\n"); }
// //     
// //     
// //     // update lam1
// //     lam1 = cppf_gibbs_lam1(logYsur, logZsur, lam0, sig2_sur, 100);
// //     postLam(s, 1) = lam1;
// //     if(s == 0){ Rprintf("sampled lam1\n"); }
// //     
// //     
// //     // update sig2_sur
// //     sig2_sur = cppf_gibbs_sig2_sur(logYsur, logZsur, lam0, lam1, 2, 2);
// //     postSig2(s, 0) = sig2_sur;
// //     if(s == 0){ Rprintf("sampled sig2_sur\n"); }
// //     
// //     
// //     // update mu0
// //     if( updateCOV ){
// //       if( (s+1 >= adaptInterval) && (s+1 - (adaptInterval * trunc((s+1) / adaptInterval)) == 0) ){
// //         dummyaccprob = accprob.col(2);
// //         rhat[2] = sum( dummyaccprob.rows(s+1-adaptInterval, s-1) ) / (adaptInterval-1);
// //         gamma1[2] = 1 / pow(adapIter[2], c1);
// //         gamma2[2] = c0 * gamma1[2];
// //         sigma2[2] = exp( log(sigma2[2]) + gamma2[2] * (rhat[2] - ropt) );
// //         
// //         postj = postMu0.col(0);
// //         COVmu0 = COVmu0 + gamma1[2] * ( var( postj.rows(s+1-adaptInterval, s-1)) - COVmu0 );
// //         cholCOVmu0 = sqrt( sigma2[2] * COVmu0 );
// //         
// //         adapIter[2] = adapIter[2] + 1;
// //       }
// //     }
// //     
// //     newmu0 = mu0 + cholCOVmu0 * randn();
// //     
// //     newloglam_pref = cppf_loglam_pref(newmu0, delta);
// //     newlam_pref = exp(newloglam_pref);
// //     
// //     logprob = cppf_loglik_pref(newloglam_pref, newlam_pref, indpref, cellarea) - 
// //       cppf_loglik_pref(loglam_pref, lam_pref, indpref, cellarea) +
// //       normal_logh(newmu0, 0, 100) - normal_logh(mu0, 0, 100);
// //     
// //     if (log(randu()) < logprob) {
// //       
// //       mu0 = newmu0;
// //       loglam_pref = newloglam_pref;
// //       lam_pref = newlam_pref;
// //       
// //       accprob(s, 2) = 1;
// //     }
// //     postMu0(s,0) = mu0;
// //     if(s == 0){ Rprintf("sampled mu0\n"); }
// //     
// //     
// //     // update delta
// //     llprev = cppf_loglik_obl(logYobl, logZobl, sig2_obl) +
// //       cppf_loglik_sur(logYsur, logZsur, lam0, lam1, sig2_sur) +
// //       cppf_loglik_pref(loglam_pref, lam_pref, indpref, cellarea);
// //     priordelta = sqrt(kappa_delta) * trans(cholRdelta) * randn(nG);
// //     
// //     thetamin = 0;
// //     thetamax = 2 * M_PI;
// //     theta = thetamin + randu() * (thetamax - thetamin);
// //     thetamin = theta - 2 * M_PI;
// //     thetamax = theta;
// //     
// //     llthred = llprev + log(randu());
// //     accept = false;
// //     count = 0;
// //     
// //     while(accept == false){
// //       count = count + 1;
// //       
// //       newdelta = delta * cos(theta) + priordelta * sin(theta);
// //       newdelta = newdelta - mean(newdelta);
// //       
// //       newlogZ = cppf_logZ_a0t_etat_ja_pref(Xzoop, alpha, eta, newdelta);
// //       newlogZobl = cppf_logZobl(newlogZ, indobl);
// //       newlogZsur = cppf_logZsur(newlogZ, indsur);
// //       newloglam = cppf_loglam_b0t_psi(Xwhale, newlogZ, beta, psi, indtwhale);
// //       newlam = exp(newloglam);
// //       newloglam_pref = cppf_loglam_pref(mu0, newdelta);
// //       newlam_pref = exp(newloglam_pref);
// //       
// //       llnew = cppf_loglik_obl(logYobl, newlogZobl, sig2_obl) +
// //         cppf_loglik_sur(logYsur, newlogZsur, lam0, lam1, sig2_sur) +
// //         cppf_loglik_pref(newloglam_pref, newlam_pref, indpref, cellarea);
// //       
// //       if(llnew > llthred){
// //         llprev = llnew;
// //         accept = true;
// //       } else {
// //         if(theta < 0){ thetamin = theta; } else { thetamax = theta; }
// //         theta = thetamin + randu() * (thetamax - thetamin);
// //         if ( (count) % 1000 == 0 ) {
// //           Rprintf("ESS for delta: %d iterations...\n", count);
// //         }
// //         if (count > 20000) stop("ESS is not converging...\n");
// //       }
// //     }
// //     delta = newdelta;
// //     logZ = newlogZ;
// //     logZobl = newlogZobl;
// //     logZsur = newlogZsur;
// //     loglam = newloglam;
// //     lam = newlam;
// //     loglam_pref = newloglam_pref;
// //     lam_pref = newlam_pref;
// //     
// //     postDelta.row(s) = trans(delta);
// //     if(s == 0){ Rprintf("sampled delta\n"); }
// //     
// //     
// //     // update kappa_delta
// //     kappa_delta = cppf_gibbs_kappa_delta(delta, invRdelta, 2, 2);
// //     postKappa(s, 2) = kappa_delta;
// //     if(s == 0){ Rprintf("sampled kappa_delta\n"); }
// //     
// //     
// //     // update beta
// //     if( updateCOV ){
// //       if( (s+1 >= adaptInterval) && (s+1 - (adaptInterval * trunc((s+1) / adaptInterval)) == 0) ){
// //         dummyaccprob = accprob.col(1);
// //         rhat[1] = sum( dummyaccprob.rows(s+1-adaptInterval, s-1) ) / (adaptInterval-1);
// //         gamma1[1] = 1 / pow(adapIter[1], c1);
// //         gamma2[1] = c0 * gamma1[1];
// //         sigma2[1] = exp( log(sigma2[1]) + gamma2[1] * (rhat[1] - ropt) );
// //         
// //         COVbeta = COVbeta + gamma1[1] * ( cov( postBeta.rows(s+1-adaptInterval, s-1)) - COVbeta );
// //         cholCOVbeta = trans( chol( sigma2[1] * ( COVbeta + 0.0000000001 * diagmat(ones(n_be)) ) ) );
// //         
// //         adapIter[1] = adapIter[1] + 1;
// //       }
// //     }
// //     
// //     newbeta = beta + cholCOVbeta * randn(n_be);
// //     
// //     newloglam = cppf_loglam_b0t_psi(Xwhale, logZ, newbeta, psi, indtwhale);
// //     newlam = exp(newloglam);
// //     
// //     logprob = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) -
// //       cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //       cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea) -
// //       cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea) +
// //       MVN_logh(newbeta, zeros(n_be), diagmat(ones(n_be))/100) -
// //       MVN_logh(beta, zeros(n_be), diagmat(ones(n_be))/100);
// //     
// //     if (log(randu()) < logprob) {
// //       
// //       beta = newbeta;
// //       loglam = newloglam;
// //       lam = newlam;
// //       
// //       accprob(s, 1) = 1;
// //     }
// //     
// //     postBeta.row(s) = trans(beta);
// //     if(s == 0){ Rprintf("sampled beta\n"); }
// //     
// //     
// //     
// //     // update psi via elliptical slice sampling
// //     llprev = cppf_loglik_dist(loglam, lam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //       cppf_loglik_pam_ct(Wpam, lam, dur_dist, c, p_pam, cellarea);
// //     priorpsi = sqrt(kappa_psi) * trans(cholRpsi) * randn(nG);
// //     
// //     thetamin = 0;
// //     thetamax = 2 * M_PI;
// //     theta = thetamin + randu() * (thetamax - thetamin);
// //     thetamin = theta - 2 * M_PI;
// //     thetamax = theta;
// //     
// //     llthred = llprev + log(randu());
// //     accept = false;
// //     count = 0;
// //     
// //     while(accept == false){
// //       count = count + 1;
// //       
// //       newpsi = psi * cos(theta) + priorpsi * sin(theta);
// //       newpsi = newpsi - mean(newpsi);
// //       
// //       newloglam = cppf_loglam_b0t_psi(Xwhale, logZ, beta, newpsi, indtwhale);
// //       newlam = exp(newloglam);
// //       
// //       llnew = cppf_loglik_dist(newloglam, newlam, pii, f_dist_obs, f_dist, inddist, indmdist, indmwhale, cellarea) +
// //         cppf_loglik_pam_ct(Wpam, newlam, dur_dist, c, p_pam, cellarea);
// //       
// //       if(llnew > llthred){
// //         llprev = llnew;
// //         accept = true;
// //       } else {
// //         if(theta < 0){ thetamin = theta; } else { thetamax = theta; }
// //         theta = thetamin + randu() * (thetamax - thetamin);
// //         if ( (count) % 1000 == 0 ) {
// //           Rprintf("ESS for psi: %d iterations...\n", count);
// //         }
// //         if (count > 20000) stop("ESS is not converging...\n");
// //       }
// //     }
// //     psi = newpsi;
// //     loglam = newloglam;
// //     lam = newlam;
// //     
// //     postPsi.row(s) = trans(psi);
// //     if(s == 0){ Rprintf("sampled psi\n"); }
// //     
// //     
// //     // update kappa_psi
// //     kappa_psi = cppf_gibbs_kappa_psi(psi, invRpsi, 2, 2);
// //     postKappa(s,1) = kappa_psi;
// //     if(s == 0){ Rprintf("sampled kappa_psi\n"); }
// //     
// //     
// //     if ( (s+1) % 10 == 0 ) {
// //       Rprintf("Generated %d samples...\n", s+1);
// //     }
// //   }
// //   
// //   
// //   return Rcpp::List::create(Rcpp::Named("postAlpha0tilde") = postAlpha0tilde,
// //                             Rcpp::Named("postAlpha") = postAlpha,
// //                             Rcpp::Named("postEta") = postEta,
// //                             Rcpp::Named("postLam") = postLam,
// //                             Rcpp::Named("postKappa") = postKappa,
// //                             Rcpp::Named("postSig2") = postSig2,
// //                             Rcpp::Named("postMu0") = postMu0,
// //                             Rcpp::Named("postDelta") = postDelta,
// //                             Rcpp::Named("postBeta") = postBeta,
// //                             Rcpp::Named("postPsi") = postPsi,
// //                             Rcpp::Named("Accprob") = accprob,
// //                             Rcpp::Named("sigma2") = sigma2,
// //                             Rcpp::Named("adapIter") = adapIter,
// //                             Rcpp::Named("COValpha") = COValpha,
// //                             Rcpp::Named("COVmu0") = COVmu0,
// //                             Rcpp::Named("COVbeta") = COVbeta);
// //   
// // }


