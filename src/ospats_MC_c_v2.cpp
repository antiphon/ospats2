 #include <Rcpp.h>
#include "helpers.h"
using namespace Rcpp;

// [[Rcpp::export]]
List ospats2_MC_v2_c(
  NumericMatrix Z,
  IntegerVector strat_init,
  NumericVector Sh_init,
  int niter_inner,
  int niter_outer,
  double temperature,
  double coolingrate,
  int verbose
) {


  int i, m, ii, si, si_new;
  int it_out,it_in;

  int transfers;

  int n_strata = Sh_init.size();
  NumericVector Sh(n_strata);
  NumericVector Sh_best(n_strata);
  double Obj = sum(Sh_init); // make sure Sh's area already weighted by sqrt(Nh)

  int n = strat_init.size();
  IntegerVector strat(n);
  IntegerVector strat_best(n);


  // sum and squaresum matrices.
  int M = Z.ncol();

  NumericMatrix Ysum (n_strata, M);
  NumericMatrix Y2sum(n_strata, M);
  IntegerVector Nh0(n_strata), Nh(n_strata);

  double Odelta, Sh_no_i, Sh_add_i, prob, Obarfinal, Obarbest = R_PosInf;

  bool verbose2 = verbose > 1;
  bool verbose1 = verbose > 0;

  double x_p, x_m;

  // initial counts
  for(i=0; i < n; i++) Nh0(strat_init(i))++;

  // order of visits
  IntegerVector u(n);
  for(i = 0; i < n; i++) u[i] = i;

  for(it_out = 0; it_out < niter_outer; it_out++) {
    // start from the sample
    std::copy(strat_init.begin(), strat_init.end(), strat.begin());
    std::copy(Sh_init.begin(), Sh_init.end(), Sh.begin());
    std::copy(Nh0.begin(), Nh0.end(), Nh.begin());
    for(m = 0; m < M; m++) {
      for(si = 0; si < n_strata; si++) {
        Ysum(si,m) =0;
        Y2sum(si,m)=0;
      }
      for(i = 0; i < n; i++){
        Ysum(  strat(i), m) += Z(i,m);
        Y2sum( strat(i), m) += pow(Z(i,m), 2);
      }
    }
    // done init
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
          // Updates:
          x_m = 0;
          x_p = 0;
          for(m = 0; m < M; m++) {
            x_p += sqrt( Y2sum(si_new, m) + pow(Z(i,m),2) - pow(Ysum(si_new, m)+Z(i,m), 2)/(Nh(si_new)+1)  );
            x_m += sqrt( Y2sum(si    , m) - pow(Z(i,m),2) - pow(Ysum(si    , m)-Z(i,m), 2)/(Nh(si)-1)  );
          }
          Sh_add_i = sqrt(Nh(si_new)+1) * x_p /M;
          Sh_no_i  = sqrt(Nh(si)-1) * x_m / M;
          //
          Odelta = Sh_add_i - Sh(si_new) + Sh_no_i - Sh(si);
          if(Odelta < (-Obj * 1e-10) ) Odelta = 0;
          prob = exp(-Odelta/temperature);
          if(Rcpp::runif(1)[0] < prob) {
            strat(i)     = si_new;
            // update Ysum and Y2sum
            for(m = 0; m < M; m++) {
              Ysum(si    , m)  -= Z(i,m);
              Ysum(si_new, m)  += Z(i,m);
              Y2sum(si    , m) -= pow(Z(i,m), 2);
              Y2sum(si_new, m) += pow(Z(i,m), 2);
            }
            Sh(si)     = Sh_no_i;
            Sh(si_new) = Sh_add_i;
            Nh(si    )--;
            Nh(si_new)++;
            Obj = sum(Sh);
            transfers++;
            break; // for si_new loop
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
