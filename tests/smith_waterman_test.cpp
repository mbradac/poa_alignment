#include <cassert>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <vector>
#include "sequence.hpp"

using namespace SmithWatermanSIMD;

std::vector<short> SmithWatermanNoSimd(Sequence query,
                                       std::vector<Sequence> database,
                                       const ScoreMatrix &matrix, int q, int r);

std::string ReadFile(const char *path) {
  std::ifstream file;
  file.open(path);
  std::stringstream file_stream;
  file_stream << file.rdbuf();
  file.close();
  return file_stream.str();
}

int main() {
  const char query_path[] = "data/query/P19930.fasta";
  const char database_path[] = "data/database/uniprot_sprot196.fasta";
  const char matrix_path[] = "data/matrix/blosum50.mat";
  const char results_path[] = "data/results/P19930_sprot196_blosum50_r1q3";
  int r = 1;
  int q = 3;

  std::vector<Sequence> query_vector;
  std::string query_string = ReadFile(query_path);
  assert(ParseFasta(query_string, &query_vector) == 0);
  assert(query_vector.size() == 1U);
  std::vector<Sequence> database;
  std::string database_string = ReadFile(database_path);
  assert(ParseFasta(database_string, &database) == 0);
  std::string matrix_string = ReadFile(matrix_path);
  ScoreMatrix matrix;
  assert(matrix.Init(matrix_string) == 0);

  std::vector<short> results =
      SmithWatermanNoSimd(query_vector[0], database, matrix, q, r);
  std::vector<short> real_results;
  std::ifstream results_file;
  results_file.open(results_path);
  short result;
  while (results_file >> result) {
    real_results.push_back(result);
  }
  results_file.close();
  if (results.size() != real_results.size()) {
    printf("WRONG!\n");
    return 1;
  }
  for (int i = 0; i < (int)results.size(); ++i) {
    if (results[i] != real_results[i]) {
      printf("WRONG!\n");
      return 1;
    }
  }
  printf("OK\n");
  return 0;
}
