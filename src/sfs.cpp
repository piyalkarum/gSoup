#include <Rcpp.h>
#include <vector>
#include <string>
#include <algorithm>
#include <unordered_map>

using namespace Rcpp;

// [[Rcpp::export]]
IntegerVector calculate_sfs_cpp(CharacterMatrix alignment) {
  int num_seqs = alignment.nrow();
  int num_sites = alignment.ncol();

  // Initialize a vector to store the site frequency spectrum
  IntegerVector sfs(num_seqs - 1);

  // Loop over each site (column)
  for (int i = 0; i < num_sites; i++) {
    std::unordered_map<std::string, int> allele_counts;
    int valid_alleles = 0;

    // Count alleles at the site, ignoring gaps
    for (int j = 0; j < num_seqs; j++) {
      std::string allele = as<std::string>(alignment(j, i));
      if (allele != "-") {
        allele_counts[allele]++;
        valid_alleles++;
      }
    }

    // Skip site if fewer than 2 valid sequences or monomorphic
    if (valid_alleles < 2 || allele_counts.size() < 2) continue;

    // Collect allele counts and sort descending
    std::vector<int> counts;
    for (auto const& pair : allele_counts) {
      counts.push_back(pair.second);
    }
    std::sort(counts.begin(), counts.end(), std::greater<int>());

    // Sum derived allele counts (all except the most frequent)
    int derived_count = 0;
    for (size_t k = 1; k < counts.size(); k++) {
      derived_count += counts[k];
    }

    // Make sure derived_count fits into the SFS vector
    if (derived_count > 0 && derived_count < num_seqs) {
      sfs[derived_count - 1]++;
    }
  }

  return sfs;
}