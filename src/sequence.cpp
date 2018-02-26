#include "sequence.hpp"

#include <iostream>
#include <cstring>
#include <sstream>
#include <cassert>

namespace poa_alignment {

std::vector<Sequence> ParseFastq(const std::string &fastq) {
  std::istringstream fastq_stream(fastq);
  std::string line;
  Sequence sequence;
  std::vector<Sequence> sequences;

  while (getline(fastq_stream, line)) {
    Sequence sequence;

    getline(fastq_stream, line);
    assert(line.size() != 0U);
    sequence.sequence = line;

    getline(fastq_stream, line);
    assert(line.size() != 0U);

    getline(fastq_stream, line);
    assert(line.size() != 0U);
    for (auto c : line) {
      sequence.quality.push_back(c - '!' + 1);
    }
    sequences.push_back(sequence);
  }
  return sequences;
}
}
