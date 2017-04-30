#include "sequence.hpp"

#include <cstring>
#include <sstream>

namespace SmithWatermanSIMD {

int ParseFasta(const std::string &fasta, std::vector<Sequence> *sequences) {
  std::istringstream fasta_stream(fasta);
  std::string line;
  Sequence sequence;
  bool is_first = true;

  while (getline(fasta_stream, line)) {
    if (line.size() == 0U) return -1;
    if (line[0] == '>') {
      if (!is_first) {
        if (sequence.sequence == "") return -2;
        sequences->push_back(sequence);
        sequence.identifier = sequence.description = sequence.sequence = "";
      }
      std::istringstream line_stream(line.substr(1));
      line_stream >> sequence.identifier;
      getline(line_stream >> std::ws, sequence.description);
    } else {
      if (is_first == true) return -3;
      sequence.sequence += line;
    }
    is_first = false;
  }
  if (sequence.sequence == "") return -2;
  sequences->push_back(sequence);
  return 0;
}

int ScoreMatrix::Init(const std::string &score_matrix) {
  std::istringstream matrix_stream(score_matrix);
  std::string line;
  memset(matrix_, -1, sizeof matrix_);
  memset(position_of_letter_, -1, sizeof position_of_letter_);

  if (!getline(matrix_stream, line)) return -1;
  std::string letter;
  std::istringstream line_stream(line);
  alphabet_size_ = 0;
  while (line_stream >> letter) {
    if (letter.size() != 1U) return -2;
    position_of_letter_[(int)letter[0]] = alphabet_size_++;
  }
  for (int i = 0; i < alphabet_size_; ++i) {
    if (!getline(matrix_stream, line)) return -3;
    std::istringstream line_stream(line);
    for (int j = 0; j < alphabet_size_; ++j) {
      if (!(line_stream >> matrix_[i][j])) return -3;
    }
  }
  ++alphabet_size_;
  return 0;
}

int TranslateSequence(Sequence *sequence, const ScoreMatrix &matrix,
                      int length_multiple) {
  for (int i = 0; i < (int)sequence->sequence.size(); ++i) {
    int c = matrix.get_position(sequence->sequence[i]);
    if (c == -1) return -1;
    sequence->sequence[i] = c;
  }
  for (int i = (int)sequence->sequence.size(); i % length_multiple != 0; ++i) {
    sequence->sequence.push_back(matrix.dummy_character());
  }
  return 0;
}
}
