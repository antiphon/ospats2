#include <Rcpp.h>
#include "helpers.h"
using namespace Rcpp;

// [[Rcpp::export]]
List ospats2_c(
  NumericMatrix D2,
  IntegerVector strat_init,
  NumericVector OA2_init,
  int niter_inner,
  int niter_outer,
  double temperature,
  double coolingrate,
  int verbose
) {


  int i, j, ii, si, si_new;
  int it_out,it_in;

  int transfers;

  int n_strata = OA2_init.size();
  NumericVector OA2(n_strata);
  NumericVector OA2_best(n_strata);
  double Obj = sum(sqrt(OA2_init));

  int n = strat_init.size();
  IntegerVector strat(n);
  IntegerVector strat_best(n);

  double Odelta, OA2_no_i, OA2_add_i, prob, Obarfinal, Obarbest = R_PosInf;

  bool verbose2 = verbose > 1;
  bool verbose1 = verbose > 0;

  // order of visits
  IntegerVector u(n);
  for(i = 0; i < n; i++) u[i] = i;

  for(it_out = 0; it_out < niter_outer; it_out++) {
    Rcpp::checkUserInterrupt();
    std::copy(strat_init.begin(), strat_init.end(), strat.begin());
    std::copy(OA2_init.begin(), OA2_init.end(), OA2.begin());

    for(it_in = 0; it_in < niter_inner; it_in++) {
      transfers = 0;
      randomShuffle(u); // randomise order of visiting units

      for(ii = 0; ii < n; ii++){
        i  = u[ii];
        si = strat[i];
        for(si_new = 0; si_new < n_strata; si_new++) {
          if(si == si_new) continue;
          // is this strata better for i?
          OA2_no_i  = OA2[si];
          OA2_add_i = OA2[si_new];
          for(j=0; j<n; j++) {
            if(strat[j] == si )          OA2_no_i -= D2(i,j);
            else if(strat[j] == si_new) OA2_add_i += D2(i,j);
          }
          Odelta = sqrt(OA2_no_i)- sqrt(OA2[si]) + sqrt(OA2_add_i) - sqrt(OA2[si_new]);
          if(Odelta < (-Obj * 1e-10) ) Odelta = 0;
          prob = exp(-abs(Odelta)/temperature);
          if(Rcpp::runif(1)[0] < prob) {
            strat[i]    = si_new;
            OA2[si]     = OA2_no_i;
            OA2[si_new] = OA2_add_i;
            Obj = sum(sqrt(OA2));
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
    Obj = sum(sqrt(OA2));
    Obarfinal = Obj/(1.0*n);

    if(Obarfinal < Obarbest) {
      Obarbest = Obarfinal;
      std::copy(strat.begin(), strat.end(), strat_best.begin());
      std::copy(OA2.begin(), OA2.end(), OA2_best.begin());
    }
    if(verbose1) Rprintf("Run %4i: Objective function [%7.4f]\n", it_out, Obarfinal);
  } // outer


  return List::create(Named("OA2_best")   = OA2_best,
                      Named("strat_best") = strat_best);
}


