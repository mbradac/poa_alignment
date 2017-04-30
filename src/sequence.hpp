#ifndef MISLAV_BRADAC_SEMINAR_SEQUENCE_H_
#define MISLAV_BRADAC_SEMINAR_SEQUENCE_H_

#include <string>
#include <vector>

namespace SmithWatermanSIMD {

struct Sequence {
  std::string identifier;
  std::string description;
  std::string sequence;
};

int ParseFasta(const std::string &fasta, std::vector<Sequence> *sequences);

class ScoreMatrix {
 public:
  const static int kMatrixSize = 257;
  ScoreMatrix() {}
  ScoreMatrix &operator=(const ScoreMatrix &) = delete;
  ScoreMatrix(const ScoreMatrix &) = delete;
  int Init(const std::string &score_matrix);
  inline int get_position(int c) const { return position_of_letter_[c]; }
  inline int get_score(int i, int j) const { return matrix_[i][j]; }
  inline int alphabet_size() const { return alphabet_size_; }
  inline int dummy_character() const { return alphabet_size_ - 1; }
  inline const int *get_matrix_row(int i) const { return matrix_[i]; }

 private:
  int alphabet_size_;
  int matrix_[kMatrixSize][kMatrixSize];
  int position_of_letter_[kMatrixSize];
};

int TranslateSequence(Sequence *sequence, const ScoreMatrix &matrix,
                      int length_multiple);
}

#endif  // MISLAV_BRADAC_SEMINAR_SEQUENCE_H_
