#include "tournament_poa_build.hpp"
#include <cstdlib>

namespace poa_alignment {
Graph GraphTournamentBuildPoa(std::vector<Sequence> sequences,
                              const ScoreMatrix &matrix, int penalty,
                              NodeStorage &storage) {
  int n = sequences.size();
  if (n == 1) {
    TranslateSequence(&sequences[0], matrix);
    return Graph(&storage, sequences[0]);
  }

  int mid = (n + 1) / 2;
  auto g1 = GraphTournamentBuildPoa(
      std::vector<Sequence>(sequences.begin(), sequences.begin() + mid), matrix,
      penalty, storage);
  auto g2 = GraphTournamentBuildPoa(
      std::vector<Sequence>(sequences.begin() + mid, sequences.end()), matrix,
      penalty, storage);
  if (rand() % 2) {
    AlignGraphToGraph(g1, g2, matrix, penalty);
    return g1;
  }
  AlignGraphToGraph(g2, g1, matrix, penalty);
  return g2;
}
}
