#ifndef POA_ALIGNMENT_SEQUENCE_H_
#define POA_ALIGNMENT_SEQUENCE_H_

#include <string>
#include <vector>

namespace poa_alignment {

struct Sequence {
  std::string identifier;
  std::string description;
  std::string sequence;
};

std::vector<Sequence> ParseFasta(const std::string &fasta);

class ScoreMatrix {
 public:
  const static int kMatrixSize = 257;
  ScoreMatrix(const std::string &score_matrix);
  ScoreMatrix &operator=(const ScoreMatrix &) = delete;
  ScoreMatrix(const ScoreMatrix &) = delete;
  int get_position(int c) const { return position_of_letter_[c]; }
  int get_score(int i, int j) const { return matrix_[i][j]; }
  int alphabet_size() const { return alphabet_size_; }
  int dummy_character() const { return alphabet_size_ - 1; }
  const int *get_matrix_row(int i) const { return matrix_[i]; }
  char position_to_letter(int x) const {
    for (int i = 0; i < kMatrixSize; ++i) {
      if (position_of_letter_[i] == x) return i;
    }
    return -1;
  }

 private:
  int alphabet_size_;
  int matrix_[kMatrixSize][kMatrixSize];
  int position_of_letter_[kMatrixSize];
};

void TranslateSequence(Sequence *sequence, const ScoreMatrix &matrix,
                       int length_multiple = 1);
}

#endif  // POA_ALIGNMENT_SEQUENCE_H_
