#include <Rcpp.h>
#include <map>
#include <string>
#include <algorithm>
#include <set>

using namespace Rcpp;

// Genetic code map
std::map<std::string, std::string> create_genetic_code() {
  std::map<std::string, std::string> genetic_code = {
    {"TTT", "F"}, {"TTC", "F"}, {"TTA", "L"}, {"TTG", "L"},
    {"CTT", "L"}, {"CTC", "L"}, {"CTA", "L"}, {"CTG", "L"},
    {"ATT", "I"}, {"ATC", "I"}, {"ATA", "I"}, {"ATG", "M"},
    {"GTT", "V"}, {"GTC", "V"}, {"GTA", "V"}, {"GTG", "V"},
    {"TCT", "S"}, {"TCC", "S"}, {"TCA", "S"}, {"TCG", "S"},
    {"CCT", "P"}, {"CCC", "P"}, {"CCA", "P"}, {"CCG", "P"},
    {"ACT", "T"}, {"ACC", "T"}, {"ACA", "T"}, {"ACG", "T"},
    {"GCT", "A"}, {"GCC", "A"}, {"GCA", "A"}, {"GCG", "A"},
    {"TAT", "Y"}, {"TAC", "Y"}, {"TAA", "X"}, {"TAG", "X"},
    {"CAT", "H"}, {"CAC", "H"}, {"CAA", "Q"}, {"CAG", "Q"},
    {"AAT", "N"}, {"AAC", "N"}, {"AAA", "K"}, {"AAG", "K"},
    {"GAT", "D"}, {"GAC", "D"}, {"GAA", "E"}, {"GAG", "E"},
    {"TGT", "C"}, {"TGC", "C"}, {"TGA", "X"}, {"TGG", "W"},
    {"CGT", "R"}, {"CGC", "R"}, {"CGA", "R"}, {"CGG", "R"},
    {"AGT", "S"}, {"AGC", "S"}, {"AGA", "R"}, {"AGG", "R"},
    {"GGT", "G"}, {"GGC", "G"}, {"GGA", "G"}, {"GGG", "G"}
  };
  return genetic_code;
}

// Function to translate codons to amino acids
std::string codon_to_aa(const std::string& codon, const std::map<std::string, std::string>& genetic_code) {
  auto it = genetic_code.find(codon);
  if (it != genetic_code.end()) {
    return it->second;
  } else {
    return "X";  // Unknown codon
  }
}

// Function to calculate synonymous and nonsynonymous differences and total sites
// [[Rcpp::export]]
NumericVector calculate_pi(CharacterMatrix alignment) {
  int n = alignment.nrow();
  int len = alignment.ncol();
  double total_syn_sites = 0.0;
  double total_nonsyn_sites = 0.0;
  double syn_diffs = 0.0;
  double nonsyn_diffs = 0.0;
  int valid_pairs = 0;

  std::map<std::string, std::string> genetic_code = create_genetic_code();

  for (int j = 0; j < n; j++) {
    for (int k = j + 1; k < n; k++) {
      double seq_syn_sites = 0.0;
      double seq_nonsyn_sites = 0.0;
      double seq_syn_diffs = 0.0;
      double seq_nonsyn_diffs = 0.0;

      for (int i = 0; i < len; i += 3) {
        if (i + 2 >= len) {
          continue; // Skip incomplete codons
        }

        std::string codon1 = as<std::string>(alignment(j, i)) + as<std::string>(alignment(j, i + 1)) + as<std::string>(alignment(j, i + 2));
        std::string codon2 = as<std::string>(alignment(k, i)) + as<std::string>(alignment(k, i + 1)) + as<std::string>(alignment(k, i + 2));

        // Convert to uppercase
        std::transform(codon1.begin(), codon1.end(), codon1.begin(), ::toupper);
        std::transform(codon2.begin(), codon2.end(), codon2.begin(), ::toupper);

        
        std::string aa1 = codon_to_aa(codon1, genetic_code);
        std::string aa2 = codon_to_aa(codon2, genetic_code);

        
        // Calculate synonymous and nonsynonymous sites for this codon
        double codon_syn_sites = 0.0;
        double codon_nonsyn_sites = 0.0;

        // Iterate through all three positions of the codon
        for (int pos = 0; pos < 3; pos++) {
          std::string mutated_codon1 = codon1;
          std::string mutated_codon2 = codon2;

          for (char base : {'A', 'C', 'G', 'T'}) {
            if (base != codon1[pos]) {
              mutated_codon1[pos] = base;
              std::string mutated_aa1 = codon_to_aa(mutated_codon1, genetic_code);
              if (mutated_aa1 == aa1) {
                codon_syn_sites += 1.0 / 3.0;
              } else {
                codon_nonsyn_sites += 1.0 / 3.0;
              }
            }

            if (base != codon2[pos]) {
              mutated_codon2[pos] = base;
              std::string mutated_aa2 = codon_to_aa(mutated_codon2, genetic_code);
              if (mutated_aa2 == aa2) {
                codon_syn_sites += 1.0 / 3.0;
              } else {
                codon_nonsyn_sites += 1.0 / 3.0;
              }
            }
          }
        }

        // Add to the sequence-level totals
        seq_syn_sites += codon_syn_sites / 2.0 ;
        seq_nonsyn_sites += codon_nonsyn_sites / 2.0 ;

        // Count differences
        if (aa1 != aa2 || aa1 == "X" || aa2 == "X") {
          seq_nonsyn_diffs += 1;
        } else {
        if (codon1 != codon2) {
          seq_syn_diffs += 1;
        }
        }
      }

      if (seq_syn_sites > 0 || seq_nonsyn_sites > 0) {
        total_syn_sites += seq_syn_sites;
        total_nonsyn_sites += seq_nonsyn_sites;
        syn_diffs += seq_syn_diffs;
        nonsyn_diffs += seq_nonsyn_diffs;
        valid_pairs++;
      }
    }
  }

  if (valid_pairs == 0 || total_syn_sites == 0 || total_nonsyn_sites == 0) {
    return NumericVector::create(NA_REAL, NA_REAL);
  }

  double pi_syn = syn_diffs /3.0 / total_syn_sites;
  double pi_nonsyn = nonsyn_diffs / 3.0 / total_nonsyn_sites;

  return NumericVector::create(pi_syn, pi_nonsyn);
}
