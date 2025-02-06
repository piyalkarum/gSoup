#include <Rcpp.h>
#include <progress.hpp>
using namespace Rcpp;

// Helper function to compute the mean of a range, ignoring NAs
double mean_na_rm(const NumericVector& x) {
  NumericVector y = na_omit(x);
  if (y.size() == 0) return NA_REAL;
  return mean(y);
}

// Function to process a single vector (row or column)
NumericVector process_vector(const NumericVector& x, int width) {
  int n = x.size();
  NumericVector result(n, NA_REAL);

  for (int j = 0; j < n; ++j) {
    NumericVector left_window, right_window;

    // Left window
    if (j + width < n) {
      left_window = x[Range(j, j + width)];
    } else {
      left_window = x[Range(j, n - 1)];
    }
    double left_mean = mean_na_rm(left_window);

    // Right window
    if (j - width >= 0) {
      right_window = x[Range(j - width, j)];
    } else {
      right_window = x[Range(0, j)];
    }
    double right_mean = mean_na_rm(right_window);

    // Combine means
    result[j] = mean_na_rm(NumericVector::create(left_mean, right_mean));
  }

  return result;
}

// [[Rcpp::export]]
NumericMatrix wind_intpol_cpp(NumericMatrix mt, int wind = 50, int width = 3) {
  int nrow = mt.nrow();
  int ncol = mt.ncol();

  // Initialize progress bar
  Progress pb(wind, true);

  for (int i = 0; i < wind; ++i) {
    if (Progress::check_abort()) return NumericMatrix(0);

    // Process rows
    for (int r = 0; r < nrow; ++r) {
      NumericVector row = mt(r, _);
      NumericVector processed_row = process_vector(row, width);
      mt(r, _) = processed_row;
    }

    // Process columns
    for (int c = 0; c < ncol; ++c) {
      NumericVector col = mt(_, c);
      NumericVector processed_col = process_vector(col, width);
      mt(_, c) = processed_col;
    }

    // Update progress bar
    pb.increment();
  }

  return mt;
}
