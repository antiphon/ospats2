 #include <Rcpp.h>
#include "helpers.h"
using namespace Rcpp;

// first go, brute force

// [[Rcpp::export]]
NumericVector compute_Sh(NumericMatrix Z,
                         IntegerVector strat,
                         int nstrata){
  NumericVector Sh(nstrata);
  double M = Z.ncol();
  int N = Z.nrow();
  NumericVector sumy2(nstrata);
  NumericVector  sumy(nstrata);
  IntegerVector Nh(nstrata);
  int i, m, h;
  double Th;

  //Nh.fill(0);
  for(i = 0; i < N; i++) Nh(strat(i))++;

  //
  for(h = 0; h < nstrata; h++) {
    Th = 0;
    for(m = 0; m < M; m++) {
      sumy.fill(0);
      sumy2.fill(0);
      for(i = 0; i < N; i++) {
        sumy (strat(i)) += Z(i,m);
        sumy2(strat(i)) += Z(i,m) * Z(i,m);
      }
      Th += sqrt(sumy2(h) - sumy(h) * sumy(h) / Nh(h));
    }
    Sh(h) = Th/M * sqrt(Nh(h));
  }
  return Sh;
}

// Sh of one h

// [[Rcpp::export]]
double compute_Sh_h(NumericMatrix Z,
                         IntegerVector strat,
                         int h){
  double Sh;
  double M = Z.ncol();
  int N = Z.nrow();
  double sumy2;
  double  sumy;
  int Nh;
  int i, m;
  double Th;

  Nh = 0; //Nh.fill(0);
  for(i = 0; i < N; i++) if(strat(i) == h) Nh++;

  //
  Th = 0;
  for(m = 0; m < M; m++) {
    sumy = 0;
    sumy2 = 0;
    for(i = 0; i < N; i++) if(strat(i) == h){
      sumy  += Z(i,m);
      sumy2 += Z(i,m) * Z(i,m);
    }
    Th += sqrt(sumy2 - sumy * sumy / Nh);
  }
  Sh = Th/M * sqrt(Nh);
  return Sh;
}

// Simultaneously recompute Sh of two strata
// [[Rcpp::export]]
NumericVector compute_Sh_changes(NumericMatrix Z,
                    IntegerVector strat,
                    int h_from, int h_to){
  NumericVector Sh(2);
  double M = Z.ncol();
  int N = Z.nrow();
  NumericVector sumy2(2);
  NumericVector sumy (2);
  IntegerVector h_fromto = IntegerVector::create(h_from, h_to);
  IntegerVector Nh(2);
  int i, h, m;
  double Th;

  //Nh.fill(0);
  for(i = 0; i < N; i++)
    if(strat(i) == h_from) Nh(0)++; else if(strat(i)==h_to) Nh(1)++;

  //
  for(h = 0; h < 2; h++){
    Th = 0;
    for(m = 0; m < M; m++) {
      sumy.fill(0);
      sumy2.fill(0);
      for(i = 0; i < N; i++) if(strat(i) == h_fromto(h)){
        sumy(h)  += Z(i,m);
        sumy2(h) += Z(i,m) * Z(i,m);
      }
      Th += sqrt(sumy2(h) - sumy(h) * sumy(h) / Nh(h) );
    }
    Sh(h) = Th/M * sqrt(Nh(h));
  }
  return Sh;
}


// [[Rcpp::export]]
List ospats2_MC_c(
  NumericMatrix Z,
  IntegerVector strat_init,
  NumericVector Sh_init,
  int niter_inner,
  int niter_outer,
  double temperature,
  double coolingrate,
  int verbose
) {


  int i, j, ii, si, si_new;
  int it_out,it_in;

  int transfers;

  int n_strata = Sh_init.size();
  NumericVector Sh(n_strata);
  NumericVector Sh_best(n_strata);
  double Obj = sum(Sh_init); // make sure Sh's area already weighted by sqrt(Nh)

  int n = strat_init.size();
  IntegerVector strat(n);
  IntegerVector strat_best(n);

  NumericVector Sh_new(n_strata);
  NumericVector Sh_changes(2);

  double Odelta, Sh_no_i, Sh_add_i, prob, Obarfinal, Obarbest = R_PosInf;

  bool verbose2 = verbose > 1;
  bool verbose1 = verbose > 0;

  // order of visits
  IntegerVector u(n);
  for(i = 0; i < n; i++) u[i] = i;

  for(it_out = 0; it_out < niter_outer; it_out++) {
    std::copy(strat_init.begin(), strat_init.end(), strat.begin());
    std::copy(Sh_init.begin(), Sh_init.end(), Sh.begin());

    for(it_in = 0; it_in < niter_inner; it_in++) {
      Rcpp::checkUserInterrupt();
      transfers = 0;
      randomShuffle(u); // randomise order of visiting units

      for(ii = 0; ii < n; ii++){
        i  = u[ii];
        si = strat[i];
        for(si_new = 0; si_new < n_strata; si_new++) {
          if(si == si_new) continue;
          // is this strata better for i?
          strat(i) = si_new;
          //Sh_new = compute_Sh(Z, strat, n_strata);
          //Odelta = Sh_new(si_new) - Sh(si_new) + Sh_new(si) - Sh(si);
          Sh_no_i =  compute_Sh_h(Z, strat, si);
          Sh_add_i = compute_Sh_h(Z, strat, si_new);
          Odelta = Sh_add_i - Sh(si_new) + Sh_no_i - Sh(si);
          //Sh_changes = compute_Sh_changes(Z, strat, si, si_new);
          //Odelta = Sh_changes(1) - Sh(si_new) + Sh_changes(0) - Sh(si);
          if(Odelta < (-Obj * 1e-10) ) Odelta = 0;
          prob = exp(-Odelta/temperature);
          if(Rcpp::runif(1)[0] < prob) {
            strat[i]    = si_new;
            //Sh(si)     = Sh_new(si);
            //Sh(si_new) = Sh_new(si_new);
            Sh(si)     = Sh_no_i;
            Sh(si_new) = Sh_add_i;
            //Sh(si)       = Sh_changes(0);
            //Sh(si_new)   = Sh_changes(1);
            Obj = sum(Sh);
            transfers++;
            break; // for si_new loop
          }
          else{ // put back
            strat(i) = si;
          }
        }
      } // over units

      // cooling
      temperature *= coolingrate;
      if(verbose2) Rprintf("   tras'rs [%i]\n", transfers);
      if(transfers == 0) break;
    }// innter
    Obj = sum(Sh);
    Obarfinal = Obj/(1.0*n);

    if(Obarfinal < Obarbest) {
      Obarbest = Obarfinal;
      std::copy(strat.begin(), strat.end(), strat_best.begin());
      std::copy(Sh.begin(), Sh.end(), Sh_best.begin());
    }
    if(verbose1) Rprintf("Run %4i: Objective function [%7.4f]\n", it_out, Obarfinal);
  } // outer


  return List::create(Named("Sh_best")   = Sh_best,
                      Named("strat_best") = strat_best);
}
