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
  
  // Loop over each site
  for (int i = 0; i < num_sites; i++) {
    // Create a map to count alleles
    std::unordered_map<std::string, int> allele_counts;
    
    // Loop over each sequence at the current site
    for (int j = 0; j < num_seqs; j++) {
      std::string allele = as<std::string>(alignment(j, i));
      
      // Exclude gaps
      if (allele != "-") {
        allele_counts[allele]++;
      }
    }
    
    // Only proceed if there are at least two valid alleles
    if (allele_counts.size() >= 2) {
      // Convert map values to a vector to find the derived allele counts
      std::vector<int> counts;
      for (auto const& pair : allele_counts) {
        counts.push_back(pair.second);
      }
      
      // Sort counts to identify the derived alleles (exclude the most frequent one)
      std::sort(counts.begin(), counts.end(), std::greater<int>());
      
      // Sum the counts of derived alleles (everything except the first, most frequent one)
      int derived_count = 0;
      for (size_t k = 1; k < counts.size(); k++) {
        derived_count += counts[k];
      }
      
      // Ensure the derived count is within the range of the SFS
      if (derived_count > 0 && derived_count < num_seqs) {
        // Increment the corresponding bin in the SFS
        sfs[derived_count - 1]++; // -1 because R indexing starts at 1
      }
    }
  }
  
  return sfs;
}
