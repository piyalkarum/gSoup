#include <Rcpp.h>
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>
#include <set>

using namespace Rcpp;

// Function to calculate pairwise differences with different models
std::pair<double, int> pairwise_diff(const CharacterVector& seq1, const CharacterVector& seq2, bool pairwise_deletion, std::string model) {
  int diff_count = 0;
  int trans_count = 0;
  int total_count = 0;
  
  for (int i = 0; i < seq1.size(); ++i) {
    if (pairwise_deletion) {
      if (!CharacterVector::is_na(seq1[i]) && !CharacterVector::is_na(seq2[i]) && seq1[i] != "-" && seq2[i] != "-") {
        total_count++;
        if (seq1[i] != seq2[i]) {
          diff_count++;
          if ((seq1[i] == "A" && seq2[i] == "G") || (seq1[i] == "G" && seq2[i] == "A") ||
              (seq1[i] == "C" && seq2[i] == "T") || (seq1[i] == "T" && seq2[i] == "C")) {
            trans_count++;
          }
        }
      }
    } else {
      if (!CharacterVector::is_na(seq1[i]) && !CharacterVector::is_na(seq2[i])) {
        total_count++;
        if (seq1[i] != seq2[i]) {
          diff_count++;
          if ((seq1[i] == "A" && seq2[i] == "G") || (seq1[i] == "G" && seq2[i] == "A") ||
              (seq1[i] == "C" && seq2[i] == "T") || (seq1[i] == "T" && seq2[i] == "C")) {
            trans_count++;
          }
        }
      }
    }
  }
  
  double distance = 0.0;
  if (model == "raw") {
    distance = static_cast<double>(diff_count) / total_count;
  } else if (model == "JC69") {
    double p = static_cast<double>(diff_count) / total_count;
    distance = -0.75 * std::log(1 - (4.0 / 3.0) * p);
  } else if (model == "K80") {
    double p = static_cast<double>(trans_count) / total_count;
    double q = static_cast<double>(diff_count - trans_count) / total_count;
    distance = -0.5 * std::log(1 - 2 * p - q) - 0.25 * std::log(1 - 2 * q);
  }
  
  return std::make_pair(distance, total_count);
}

// [[Rcpp::export]]
List calc_nuc_div_theta(CharacterMatrix dna_matrix, bool pairwise_deletion = true, std::string model = "raw") {
  int num_seqs = dna_matrix.nrow();
  int seq_length = dna_matrix.ncol();
  NumericVector sample_avg_distances(num_seqs, 0.0);
  NumericVector counts(num_seqs, 0.0);
  double total_diff = 0;
  double total_sites = 0;

  // Calculate pairwise distances and accumulate averages per sample
  for (int i = 0; i < num_seqs - 1; ++i) {
    for (int j = i + 1; j < num_seqs; ++j) {
      std::pair<double, int> diffs = pairwise_diff(dna_matrix(i, _), dna_matrix(j, _), pairwise_deletion, model);
      sample_avg_distances[i] += diffs.first;
      sample_avg_distances[j] += diffs.first;
      counts[i]++;
      counts[j]++;
      total_diff += diffs.first * diffs.second;
      total_sites += diffs.second;
    }
  }

  // Finalize average distances per sample
  for (int i = 0; i < num_seqs; ++i) {
    if (counts[i] > 0) {
      sample_avg_distances[i] /= counts[i];
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
    Named("seg_sit") = segregating_sites,
    Named("average_per_sample") = sample_avg_distances
  );
}