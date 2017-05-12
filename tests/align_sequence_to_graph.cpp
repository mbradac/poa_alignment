#include <cassert>
#include <cstdio>
#include "sequence.hpp"
#include "poa_alignment.hpp"
#include <cassert>
#include <sstream>
#include <fstream>
#include <stack>
#include <unordered_set>

using namespace poa_alignment;

std::string ReadFile(const char *path) {
  std::ifstream file;
  file.open(path);
  std::stringstream file_stream;
  file_stream << file.rdbuf();
  file.close();
  return file_stream.str();
}

// TODO: write real tests
void Run(const std::string &s1, const std::string &s2) {
  printf("Run\n===\n");
  Sequence sequence1;
  sequence1.sequence = s1;
  Sequence sequence2;
  sequence2.sequence = s2;
  const char matrix_path[] = "data/matrix/blosum50.mat";
  Storage<Node> storage;
  ScoreMatrix matrix(ReadFile(matrix_path));
  auto graph = GraphFromSequence(sequence1, matrix, storage);
  AlignSequenceToGraph(graph, sequence2, matrix, 1, storage);
  std::stack<Node *> jobs;
  std::unordered_set<Node *> visited;
  jobs.push(graph.first);
  while (!jobs.empty()) {
    auto *node = jobs.top();
    printf("%c: ", matrix.position_to_letter(node->letter));
    jobs.pop();
    for (auto *next : node->edges) {
      printf("%c, ", matrix.position_to_letter(next->letter));
      if (visited.find(next) == visited.end()) {
        visited.insert(next);
        jobs.push(next);
      }
    }
    printf("\n");
  }
}

int main() {
  Run("DDD", "CCP");
  Run("PPTTT", "PPKKK");
  Run("PDM", "PCM");
  Run("DDP", "CCP");
  Run("PKMCVRPQKNETC", "THKMDVRNETDM");
  return 0;
}
