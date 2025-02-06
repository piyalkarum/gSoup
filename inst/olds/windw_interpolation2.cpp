#include <Rcpp.h>
#include <progress.hpp>
using namespace Rcpp;

// Helper function to calculate mean while ignoring NA values
double mean_na(const NumericVector& x) {
  NumericVector valid = na_omit(x);
  if (valid.size() == 0) return NA_REAL;
  return mean(valid);
}

// Helper function to calculate max while ignoring NA values
double max_na(const NumericVector& x) {
  NumericVector valid = na_omit(x);
  if (valid.size() == 0) return NA_REAL;
  return max(valid);
}

// [[Rcpp::export]]
NumericMatrix wind_intpol(NumericMatrix mat, int wind = 50, int width = 3) {
  int nrow = mat.nrow();
  int ncol = mat.ncol();

  Progress pb(wind, true); // Initialize progress bar

  for (int w = 0; w < wind; ++w) {
    if (Progress::check_abort()) return mat; // Allow user interruption

    // Row-wise extrapolation
    NumericMatrix ht(nrow, ncol);
    for (int i = 0; i < nrow; ++i) {
      NumericVector row = mat(i, _);
      NumericVector mps(row.size(), NA_REAL);
      NumericVector mpsf(row.size(), NA_REAL);

      for (int j = 0; j < row.size(); ++j) {
        // Left direction
        if (j + width < row.size()) {
          mps[j] = NumericVector::is_na(row[j]) ? max_na(row[Range(j, j + 1)]) : mean_na(row[Range(j, j + width)]);
        } else {
          mps[j] = NumericVector::is_na(row[j]) ? max_na(row[Range(j, row.size() - 1)]) : mean_na(row[Range(j, row.size() - 1)]);
        }

        // Right direction
        if (j - width >= 0) {
          mpsf[j] = NumericVector::is_na(row[j]) ? max_na(row[Range(j - width, j)]) : mean_na(row[Range(j - width, j)]);
        } else {
          mpsf[j] = NumericVector::is_na(row[j]) ? max_na(row[Range(0, j)]) : mean_na(row[Range(0, j)]);
        }
      }

      for (int j = 0; j < row.size(); ++j) {
        ht(i, j) = mean_na(NumericVector::create(mps[j], mpsf[j]));
      }
    }

    // Column-wise extrapolation
    NumericMatrix vt(nrow, ncol);
    for (int j = 0; j < ncol; ++j) {
      NumericVector col = mat(_, j);
      NumericVector mps(col.size(), NA_REAL);
      NumericVector mpsf(col.size(), NA_REAL);

      for (int i = 0; i < col.size(); ++i) {
        // Left direction
        if (i + width < col.size()) {
          mps[i] = NumericVector::is_na(col[i]) ? max_na(col[Range(i, i + 1)]) : mean_na(col[Range(i, i + width)]);
        } else {
          mps[i] = NumericVector::is_na(col[i]) ? max_na(col[Range(i, col.size() - 1)]) : mean_na(col[Range(i, col.size() - 1)]);
        }

        // Right direction
        if (i - width >= 0) {
          mpsf[i] = NumericVector::is_na(col[i]) ? max_na(col[Range(i - width, i)]) : mean_na(col[Range(i - width, i)]);
        } else {
          mpsf[i] = NumericVector::is_na(col[i]) ? max_na(col[Range(0, i)]) : mean_na(col[Range(0, i)]);
        }
      }

      for (int i = 0; i < col.size(); ++i) {
        vt(i, j) = mean_na(NumericVector::create(mps[i], mpsf[i]));
      }
    }

    // Combine row and column results
    for (int i = 0; i < nrow; ++i) {
      for (int j = 0; j < ncol; ++j) {
        mat(i, j) = mean_na(NumericVector::create(ht(i, j), vt(i, j)));
      }
    }

    pb.increment(); // Update progress bar
  }

  return mat;
}
