#include <Rcpp.h>
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>

using namespace Rcpp;

// Function to calculate pairwise differences between two sequences
std::pair<int, int> pairwise_diff(const CharacterVector& seq1, const CharacterVector& seq2, bool pairwise_deletion) {
  int diff_count = 0;
  int total_count = 0;
  
  for (int i = 0; i < seq1.size(); ++i) {
    if (pairwise_deletion) {
      if (!CharacterVector::is_na(seq1[i]) && !CharacterVector::is_na(seq2[i]) && seq1[i] != "-" && seq2[i] != "-") {
        total_count++;
        if (seq1[i] != seq2[i]) {
          diff_count++;
        }
      }
    } else {
      if (!CharacterVector::is_na(seq1[i]) && !CharacterVector::is_na(seq2[i])) {
        total_count++;
        if (seq1[i] != seq2[i]) {
          diff_count++;
        }
      }
    }
  }
  return std::make_pair(diff_count, total_count);
}

// [[Rcpp::export]]
List calc_nuc_div_theta(CharacterMatrix dna_matrix, bool pairwise_deletion = true) {
  int num_seqs = dna_matrix.nrow();
  int seq_length = dna_matrix.ncol();
  double total_diff = 0;
  double total_sites = 0;

  // Calculate pairwise differences for all sequences
  for (int i = 0; i < num_seqs - 1; ++i) {
    for (int j = i + 1; j < num_seqs; ++j) {
      std::pair<int, int> diffs = pairwise_diff(dna_matrix(i, _), dna_matrix(j, _), pairwise_deletion);
      total_diff += diffs.first;
      total_sites += diffs.second;
    }
  }

  // Calculate nucleotide diversity (pi)
  double nucleotide_diversity = total_diff / total_sites;

  // Calculate the number of segregating sites
  int segregating_sites = 0;
  for (int i = 0; i < seq_length; ++i) {
    std::set<std::string> unique_bases;
    for (int j = 0; j < num_seqs; ++j) {
      if (dna_matrix(j, i) != "-" && !CharacterVector::is_na(dna_matrix(j, i))) {
        unique_bases.insert(Rcpp::as<std::string>(dna_matrix(j, i)));
      }
    }
    if (unique_bases.size() > 1) {
      segregating_sites++;
    }
  }

  // Calculate Watterson's theta
  double a1 = 0;
  for (int i = 1; i < num_seqs; ++i) {
    a1 += 1.0 / i;
  }
  double sequence_length = seq_length;
  double watterson_theta = segregating_sites / (a1 * sequence_length);

  return List::create(
    Named("pi") = nucleotide_diversity,
    Named("theta") = watterson_theta,
    Named("tot_diff") = total_diff,
    Named("seg_sit") = segregating_sites
  );
}
