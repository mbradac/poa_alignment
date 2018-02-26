#ifndef POA_ALIGNMENT_SEQUENCE_H_
#define POA_ALIGNMENT_SEQUENCE_H_

#include <string>
#include <vector>

namespace poa_alignment {

struct Sequence {
  std::string sequence;
  std::vector<int> quality;
};

std::vector<Sequence> ParseFastq(const std::string &fastq);
}

#endif  // POA_ALIGNMENT_SEQUENCE_H_
