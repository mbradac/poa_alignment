#include <vector>

#include "poa_alignment.hpp"
#include "sequence.hpp"

namespace poa_alignment {

Graph GraphTournamentBuildPoa(std::vector<Sequence> sequences,
                              const ScoreMatrix &matrix, int penalty,
                              NodeStorage &storage);
}
