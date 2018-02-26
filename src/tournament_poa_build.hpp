#include <vector>

#include "poa_alignment.hpp"
#include "sequence.hpp"

namespace poa_alignment {

Graph GraphTournamentBuildPoa(std::vector<Sequence> sequences, int match,
                              int mismatch, int gap_penalty,
                              NodeStorage &storage);
}
