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
std::vector<int> SmithWatermanNoSimd(Sequence query,
                                     std::vector<Sequence> database,
                                     const ScoreMatrix &matrix, int q, int r);

std::string ReadFile(const char *path) {
  std::ifstream file;
  file.open(path);
  std::stringstream file_stream;
  file_stream << file.rdbuf();
  file.close();
  return file_stream.str();
}

struct GraphStats {
  std::string concensus;
  int num_nodes = 0;
  int num_edges = 0;
  int total_weight = 0;
};

GraphStats GetGraphStats(const Graph &g, const ScoreMatrix &matrix) {
  GraphStats gs;
  auto nodes = TopologicalSort(g.start_nodes);
  gs.num_nodes = nodes.size();
  for (auto *node : nodes) {
    gs.num_edges += node->edges.size();
    for (auto edge : node->edges) {
      gs.total_weight += edge.second;
    }
  }

  auto concensus = FindConcensus(g.start_nodes);
  for (auto *node : concensus) {
    gs.concensus += matrix.position_to_letter(node->letter);
  }
  return gs;
}

GraphStats AlignNormal(std::vector<Sequence> seqs, const ScoreMatrix &matrix,
                       int penalty) {
  NodeStorage storage;
  for (auto &s : seqs) {
    TranslateSequence(&s, matrix);
  }
  Graph g(&storage, seqs[0]);
  for (int i = 1; i < static_cast<int>(seqs.size()); ++i) {
    AlignSequenceToGraph(
        g, seqs[i],
        std::vector<int>(static_cast<int>(seqs[i].sequence.size()) - 1, 1),
        matrix, penalty);
    std::cerr << "jesam " << i << "\n";
  }

  return GetGraphStats(g, matrix);
}

int main() {
  const char matrix_path[] = "data/matrix/blosum45.mat";
  const char sequences_path[] = "data/sequences";
  int penalty = 1;

  ScoreMatrix matrix(ReadFile(matrix_path));
  std::ifstream f(sequences_path);
  std::string sequence;
  std::vector<Sequence> sequences;
  while (f >> sequence) {
    Sequence seq;
    seq.sequence = sequence;
    sequences.push_back(seq);
  }
  sequences.erase(sequences.begin());

  for (int i = 4; i <= 32; i += 4) {
    std::vector<Sequence> seqs(sequences.begin(), sequences.begin() + i);

    std::random_shuffle(seqs.begin(), seqs.begin() + i);
    auto r = AlignNormal(seqs, matrix, penalty).concensus;
    Sequence reference;
    reference.sequence = r;

    std::random_shuffle(seqs.begin(), seqs.begin() + i);

    clock_t start1 = clock();
    auto gs1 = AlignNormal(seqs, matrix, penalty);
    clock_t end1 = clock();
    Sequence c;
    c.sequence = gs1.concensus;
    std::cout << i << " & ";
    std::cout << gs1.num_nodes << " & " << gs1.num_edges << " & "
              << gs1.total_weight << " & "
              << SmithWatermanNoSimd(c, std::vector<Sequence>{reference},
                                     matrix, 1, 1)[0] << " & "
              << (double)(end1 - start1) / CLOCKS_PER_SEC << " & ";

    NodeStorage storage;
    clock_t start2 = clock();
    auto gs2 = GetGraphStats(
        GraphTournamentBuildPoa(seqs, matrix, penalty, storage), matrix);
    clock_t end2 = clock();
    Sequence c2;
    c2.sequence = gs2.concensus;
    std::cout << gs2.num_nodes << " & " << gs2.num_edges << " & "
              << gs2.total_weight << " & "
              << SmithWatermanNoSimd(c2, std::vector<Sequence>{reference},
                                     matrix, 1, 1)[0] << " & "
              << (double)(end2 - start2) / CLOCKS_PER_SEC << std::endl;
  }

  return 0;
}
