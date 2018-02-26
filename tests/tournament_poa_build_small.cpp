#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>

#include "poa_alignment.hpp"
#include "sequence.hpp"
#include "tournament_poa_build.hpp"

using namespace poa_alignment;
// std::vector<int> SmithWatermanNoSimd(Sequence query,
//                                     std::vector<Sequence> database, int
//                                     match,
//                                     int mismatch, int gap_penalty);

int EditDistance(std::string a, std::string b) {
  std::vector<std::vector<int>> dp(a.size() + 1,
                                   std::vector<int>(b.size() + 1, 1e9));
  dp[0][0] = 0;
  for (int i = 0; i < static_cast<int>(a.size()); ++i) {
    for (int j = 0; j < static_cast<int>(b.size()); ++j) {
      dp[i + 1][j] = std::min(dp[i][j] + 1, dp[i + 1][j]);
      dp[i][j + 1] = std::min(dp[i][j] + 1, dp[i][j + 1]);
      dp[i + 1][j + 1] = std::min(dp[i][j] + (a[i] != b[j]), dp[i + 1][j + 1]);
    }
  }
  return dp[a.size()][b.size()];
}

std::string ReadFile(const char *path) {
  std::ifstream file;
  file.open(path);
  std::stringstream file_stream;
  file_stream << file.rdbuf();
  file.close();
  return file_stream.str();
}

struct GraphStats {
  std::string consensus;
  int num_nodes = 0;
  int num_edges = 0;
  int total_weight = 0;
};

GraphStats GetGraphStats(const Graph &g) {
  GraphStats gs;
  auto nodes = TopologicalSort(g.start_nodes);
  gs.num_nodes = nodes.size();
  for (auto *node : nodes) {
    gs.num_edges += node->edges.size();
    for (auto edge : node->edges) {
      gs.total_weight += edge.second;
    }
  }

  auto consensus = FindConcensus(g.start_nodes);
  for (auto x : consensus) {
    gs.consensus += x->letter;
  }
  return gs;
}

GraphStats AlignNormal(std::vector<Sequence> seqs, int match, int mismatch,
                       int gap_penalty) {
  NodeStorage storage;
  std::vector<int> weights;
  for (int j = 1; j < static_cast<int>(seqs[0].quality.size()); ++j) {
    weights.push_back(seqs[0].quality[j] + seqs[0].quality[j - 1]);
  }
  Graph g(&storage, seqs[0], weights);

  for (int i = 1; i < static_cast<int>(seqs.size()); ++i) {
    std::vector<int> weights;
    for (int j = 1; j < static_cast<int>(seqs[i].quality.size()); ++j) {
      weights.push_back(seqs[i].quality[j] + seqs[i].quality[j - 1]);
    }
    AlignSequenceToGraph(g, seqs[i], weights, match, mismatch, gap_penalty);
  }

  return GetGraphStats(g);
}

int main() {
  int match = 5;
  int mismatch = 4;
  int gap_penalty = 6;
  const char sequences_path[] = "data/test_500.fq";

  auto fastq_file = ReadFile(sequences_path);
  auto sequences = ParseFastq(fastq_file);

  std::string x1 =
      "AATGATGCGCTTTGTTGGCGCGGTGGCTTGATGCAGGGGCTAATCGACCTCTGGCAACCACTTTTCCATGAC"
      "AGGAGTTGAATATGGCATTCAGTAATCCCTTCGATGATCCGCAGGGAGCGTTTTACATATTGCGCAATGCGC"
      "AGGGGCAATTCAGTCTGTGGCCGCAACAATGCGTCTTACCGGCAGGCTGGGACATTGTGTGTCAGCCGCAGT"
      "CACAGGCGTCCTGCCAGCAGTGGCTGGAAGCCCACTGGCGTACTCTGACACCGACGAATTTTACCCAGTTGC"
      "AGGGAGGCACAATGAGCCAGCATTTACCTTTGGTCGCCGCACAGCCCGGCATCTGGATGGCAGAAAAACTGT"
      "CAGAATTACCCTCCGCCTGGAGCGTGGCGCATTACGTTGAGTTAACCGGAGAGGTTGATTCGCCATTACTGG"
      "CCCGCGCGGTGGTTGCCGGACTAGCGCAAGCAGATACGC";
  std::string x2 =
      "ACATGATGCGCTTTGTTGGCGCGGTGGCTTGATGCAGGGGCTAATCGACCTCTGGCAACCACTTTTCCATGA"
      "CAGGAGTTGAATATGGCATTCAGTAATCCCTTCGATGATCCGCGCAGGGAGCGTTTTACATATTGCGCAATG"
      "CGCAGGGGCAATTCAGTCTGTGGCCGCAACAATGCGTCTTACCGGCAGGCTGGGACATTGTGTGTCAGCCGC"
      "AGTCACAGGCGTCCTGCCAGCAGTGGCTGGAAGCCCACTGGCGTACTCTGACACCGACGAATTTTACCCAGT"
      "TGCAGGAGGCACAATGAGCCAGCATTTACCTTTGGTCGCCGCACAGCCCGGCATCTGGATGGCAGAAAAACT"
      "GTCAGAATTACCCTCCGCCTGGAGCGTGGCGCATTACGTTGAGTTAACCGGAGAGGTTGATTCGCCATTACT"
      "GGCCCGCGCGGTGGTTGCCGGACTAGCGCAAGCAGATACGC";
  std::cerr << "bitno " << EditDistance(x1, x2) << "\n";

  for (int i = 2; i <= 54; i += 4) {
    std::vector<Sequence> seqs(sequences.begin(), sequences.begin() + i);

    std::random_shuffle(seqs.begin(), seqs.end());
    auto r = AlignNormal(seqs, match, mismatch, gap_penalty).consensus;
    Sequence reference;
    reference.sequence =
        "AATGATGCGCTTTGTTGGCGCGGTGGCTTGATGCAGGGGCTAATCGACCTCTGGCAACCACTTTTCCATG"
        "ACAGGAGTTGAATATGGCATTCAGTAATCCCTTCGATGATCCGCAGGGAGCGTTTTACATATTGCGCAAT"
        "GCGCAGGGGCAATTCAGTCTGTGGCCGCAACAATGCGTCTTACCGGCAGGCTGGGACATTGTGTGTCAGC"
        "CGCAGTCACAGGCGTCCTGCCAGCAGTGGCTGGAAGCCCACTGGCGTACTCTGACACCGACGAATTTTAC"
        "CCAGTTGCAGGGAGGCACAATGAGCCAGCATTTACCTTTGGTCGCCGCACAGCCCGGCATCTGGATGGCA"
        "GAAAAACTGTCAGAATTACCCTCCGCCTGGAGCGTGGCGCATTACGTTGAGTTAACCGGAGAGGTTGATT"
        "CGCCATTACTGGCCCGCGCGGTGGTTGCCGGACTAGCGCAAGCAGATACGC";

    std::random_shuffle(seqs.begin(), seqs.end());

    clock_t start1 = clock();
    auto gs1 = AlignNormal(seqs, match, mismatch, gap_penalty);
    clock_t end1 = clock();
    Sequence c;
    c.sequence = gs1.consensus;
    std::cout << i << " & ";
    std::cout << gs1.num_nodes << " & " << gs1.num_edges << " & "
              << gs1.total_weight << " & "
              << EditDistance(c.sequence, reference.sequence) << " & "
              << (double)(end1 - start1) / CLOCKS_PER_SEC << " & ";

    NodeStorage storage;
    clock_t start2 = clock();
    auto gs2 = GetGraphStats(
        GraphTournamentBuildPoa(seqs, match, mismatch, gap_penalty, storage));
    clock_t end2 = clock();
    Sequence c2;
    c2.sequence = gs2.consensus;
    std::cout << gs2.num_nodes << " & " << gs2.num_edges << " & "
              << gs2.total_weight << " & "
              << EditDistance(c2.sequence, reference.sequence) << " & "
              << (double)(end2 - start2) / CLOCKS_PER_SEC << std::endl;
  }

  return 0;
}
