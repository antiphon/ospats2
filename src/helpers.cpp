#include "helpers.h"

inline int randWrapper(const int n) { return floor(unif_rand()*n); }

// In place
void randomShuffle(Rcpp::IntegerVector a) {
  int n = a.size();
  int j;

  // Fisher-Yates Shuffle Algorithm
  for (int i = 0; i < n - 1; i++) {
    j = i + randWrapper(n - i);
    std::swap(a[i], a[j]);
  }
}
