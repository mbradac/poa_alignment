#include "tournament_poa_build.hpp"
#include <cstdlib>

namespace poa_alignment {
Graph GraphTournamentBuildPoa(std::vector<Sequence> sequences, int match,
                              int mismatch, int gap_penalty,
                              NodeStorage &storage) {
  int n = sequences.size();
  if (n == 1) {
    std::vector<int> weights;
    for (int j = 1; j < static_cast<int>(sequences[0].quality.size()); ++j) {
      weights.push_back(sequences[0].quality[j] + sequences[0].quality[j - 1]);
    }
    return Graph(&storage, sequences[0], weights);
  }

  int mid = (n + 1) / 2;
  auto g1 = GraphTournamentBuildPoa(
      std::vector<Sequence>(sequences.begin(), sequences.begin() + mid), match,
      mismatch, gap_penalty, storage);
  auto g2 = GraphTournamentBuildPoa(
      std::vector<Sequence>(sequences.begin() + mid, sequences.end()), match,
      mismatch, gap_penalty, storage);
  AlignGraphToGraph(g1, g2, match, mismatch, gap_penalty);
  return g1;
}
}
