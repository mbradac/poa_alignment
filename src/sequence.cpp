#include "sequence.hpp"

#include <cstring>
#include <sstream>
#include <cassert>

namespace poa_alignment {

std::vector<Sequence> ParseFasta(const std::string &fasta) {
  std::istringstream fasta_stream(fasta);
  std::string line;
  Sequence sequence;
  bool is_first = true;
  std::vector<Sequence> sequences;

  while (getline(fasta_stream, line)) {
    assert(line.size() != 0U);
    if (line[0] == '>') {
      if (!is_first) {
        assert(sequence.sequence != "");
        sequences.push_back(sequence);
        sequence.identifier = sequence.description = sequence.sequence = "";
      }
      std::istringstream line_stream(line.substr(1));
      line_stream >> sequence.identifier;
      getline(line_stream >> std::ws, sequence.description);
    } else {
      assert(!is_first);
      sequence.sequence += line;
    }
    is_first = false;
  }
  assert(sequence.sequence != "");
  sequences.push_back(sequence);
  return sequences;
}

ScoreMatrix::ScoreMatrix(const std::string &score_matrix) {
  std::istringstream matrix_stream(score_matrix);
  std::string line;
  memset(matrix_, -1, sizeof matrix_);
  memset(position_of_letter_, -1, sizeof position_of_letter_);

  assert(getline(matrix_stream, line));
  std::string letter;
  std::istringstream line_stream(line);
  alphabet_size_ = 0;
  while (line_stream >> letter) {
    assert(letter.size() == 1U);
    position_of_letter_[static_cast<int>(letter[0])] = alphabet_size_++;
  }
  for (int i = 0; i < alphabet_size_; ++i) {
    assert(getline(matrix_stream, line));
    std::istringstream line_stream(line);
    for (int j = 0; j < alphabet_size_; ++j) {
      assert((line_stream >> matrix_[i][j]));
    }
  }
  ++alphabet_size_;
}

void TranslateSequence(Sequence *sequence, const ScoreMatrix &matrix,
                       int length_multiple) {
  for (int i = 0; i < (int)sequence->sequence.size(); ++i) {
    int c = matrix.get_position(sequence->sequence[i]);
    assert(c != -1);
    sequence->sequence[i] = c;
  }
  for (int i = static_cast<int>(sequence->sequence.size());
       i % length_multiple != 0; ++i) {
    sequence->sequence.push_back(matrix.dummy_character());
  }
}
}
