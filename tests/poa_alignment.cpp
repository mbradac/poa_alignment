#include <cassert>
#include <cstdio>
#include <cassert>
#include <sstream>
#include <fstream>
#include <stack>
#include <unordered_map>
#include <algorithm>

#include "sequence.hpp"
#include "poa_alignment.hpp"

namespace {

using namespace poa_alignment;

std::string ReadFile(const char *path) {
  std::ifstream file;
  file.open(path);
  std::stringstream file_stream;
  file_stream << file.rdbuf();
  file.close();
  return file_stream.str();
}

void Check(Node *n1, Node *n2, std::unordered_map<Node *, Node *> &visited) {
  assert(n1->letter == n2->letter);
  if (n1->edges.size() != n2->edges.size()) {
    assert(n1->edges.size() == n2->edges.size());
  }
  assert(n1->edges.size() == n2->edges.size());
  auto it = visited.find(n1);
  if (it != visited.end()) {
    assert(it->second == n2);
    return;
  }
  visited[n1] = n2;
  auto cmp =
      [](const std::pair<Node *, int> &n, const std::pair<Node *, int> &m) {
        return n.first->letter < m.first->letter;
      };
  std::sort(n1->edges.begin(), n1->edges.end(), cmp);
  std::sort(n2->edges.begin(), n2->edges.end(), cmp);
  for (int i = 0; i < static_cast<int>(n1->edges.size()); ++i) {
    assert(n1->edges[i].second == n2->edges[i].second);
    Check(n1->edges[i].first, n2->edges[i].first, visited);
  }
}

void Run(const std::string &s1, const std::string &s2,
         std::vector<Node *> expected_graph) {
  Sequence sequence1;
  sequence1.sequence = s1;
  Sequence sequence2;
  sequence2.sequence = s2;
  const char matrix_path[] = "data/matrix/blosum50.mat";
  ScoreMatrix matrix(ReadFile(matrix_path));
  auto graph = GraphFromSequence(sequence1, matrix);
  AlignSequenceToGraph(graph, sequence2, matrix, 1);
  std::unordered_map<Node *, Node *> visited;
  assert(graph.start_nodes.size() == expected_graph.size());
  auto cmp =
      [](const Node *n1, const Node *n2) { return n1->letter < n2->letter; };
  sort(graph.start_nodes.begin(), graph.start_nodes.end(), cmp);
  sort(expected_graph.begin(), expected_graph.end(), cmp);
  for (int i = 0; i < static_cast<int>(graph.start_nodes.size()); ++i) {
    Check(graph.start_nodes[i], expected_graph[i], visited);
  }
}
}

#define CREATE(VAR, LETTER)                                              \
  auto *VAR = new Node(LETTER == -1 ? -1 : matrix.get_position(LETTER)); \
  storage.emplace_back(VAR);

void TestAlignSequenceToGraph(const ScoreMatrix &matrix,
                              std::vector<std::unique_ptr<Node>> &storage) {
  {
    CREATE(d1, 'D');
    CREATE(d2, 'D');
    CREATE(d3, 'D');
    CREATE(c1, 'C');
    CREATE(c2, 'C');
    CREATE(p1, 'P');
    std::vector<Node *> start_nodes = {d1, c1};
    d1->edges = {{d2, 1}};
    d2->edges = {{d3, 1}};
    c1->edges = {{c2, 1}};
    c2->edges = {{p1, 1}};
    Run("DDD", "CCP", start_nodes);
  }
  {
    CREATE(p1, 'P');
    CREATE(p2, 'P');
    CREATE(t1, 'T');
    CREATE(t2, 'T');
    CREATE(t3, 'T');
    CREATE(k1, 'K');
    CREATE(k2, 'K');
    CREATE(k3, 'K');
    std::vector<Node *> start_nodes = {p1};
    p1->edges = {{p2, 2}};
    p2->edges = {{t1, 1}, {k1, 1}};
    t1->edges = {{t2, 1}};
    t2->edges = {{t3, 1}};
    k1->edges = {{k2, 1}};
    k2->edges = {{k3, 1}};
    Run("PPTTT", "PPKKK", start_nodes);
  }
  {
    CREATE(p1, 'P');
    CREATE(d1, 'D');
    CREATE(m1, 'M');
    CREATE(c1, 'C');
    std::vector<Node *> start_nodes = {p1};
    p1->edges = {{d1, 1}, {c1, 1}};
    d1->edges = {{m1, 1}};
    c1->edges = {{m1, 1}};
    Run("PDM", "PCM", start_nodes);
  }
  {
    CREATE(p1, 'P');
    CREATE(p2, 'P');
    CREATE(c1, 'C');
    std::vector<Node *> start_nodes = {p1};
    p1->edges = {{p2, 1}};
    p2->edges = {{c1, 2}};
    Run("PPC", "PC", start_nodes);
  }
  {
    CREATE(d1, 'D');
    CREATE(d2, 'D');
    CREATE(p1, 'P');
    CREATE(c1, 'C');
    CREATE(c2, 'C');
    std::vector<Node *> start_nodes = {d1, c1};
    d1->edges = {{d2, 1}};
    d2->edges = {{p1, 1}};
    c1->edges = {{c2, 1}};
    c2->edges = {{p1, 1}};
    Run("DDP", "CCP", start_nodes);
  }
  {
    CREATE(p1, 'P');
    CREATE(k1, 'K');
    CREATE(m1, 'M');
    CREATE(c1, 'C');
    CREATE(v1, 'V');
    CREATE(r1, 'R');
    CREATE(p2, 'P');
    CREATE(q1, 'Q');
    CREATE(k2, 'K');
    CREATE(n1, 'N');
    CREATE(e1, 'E');
    CREATE(t1, 'T');
    CREATE(c2, 'C');
    CREATE(t2, 'T');
    CREATE(h1, 'H');
    CREATE(d1, 'D');
    CREATE(d2, 'D');
    CREATE(m2, 'M');
    std::vector<Node *> start_nodes = {p1, t2};
    p1->edges = {{k1, 1}};
    k1->edges = {{m1, 2}};
    m1->edges = {{c1, 1}, {d1, 1}};
    c1->edges = {{v1, 1}};
    v1->edges = {{r1, 2}};
    r1->edges = {{p2, 1}, {n1, 1}};
    p2->edges = {{q1, 1}};
    q1->edges = {{k2, 1}};
    k2->edges = {{n1, 1}};
    n1->edges = {{e1, 2}};
    e1->edges = {{t1, 2}};
    t1->edges = {{c2, 1}, {d2, 1}};
    t2->edges = {{h1, 1}};
    h1->edges = {{k1, 1}};
    d1->edges = {{v1, 1}};
    d2->edges = {{m2, 1}};
    Run("PKMCVRPQKNETC", "THKMDVRNETDM", start_nodes);
  }
}

void TestFindConcensus(const ScoreMatrix &matrix,
                       std::vector<std::unique_ptr<Node>> &storage) {
  {
    CREATE(p1, 'P');
    CREATE(k1, 'K');
    CREATE(m1, 'M');
    CREATE(c1, 'C');
    CREATE(v1, 'V');
    CREATE(r1, 'R');
    CREATE(p2, 'P');
    CREATE(q1, 'Q');
    CREATE(k2, 'K');
    CREATE(n1, 'N');
    CREATE(e1, 'E');
    CREATE(t1, 'T');
    CREATE(c2, 'C');
    CREATE(t2, 'T');
    CREATE(h1, 'H');
    CREATE(d1, 'D');
    CREATE(d2, 'D');
    CREATE(m2, 'M');
    std::vector<Node *> start_nodes = {p1, t2};
    p1->edges = {{k1, 2}};
    k1->edges = {{m1, 2}};
    m1->edges = {{c1, 1}, {d1, 2}};
    c1->edges = {{v1, 1}};
    v1->edges = {{r1, 2}};
    r1->edges = {{p2, 1}, {n1, 2}};
    p2->edges = {{q1, 1}};
    q1->edges = {{k2, 1}};
    k2->edges = {{n1, 1}};
    n1->edges = {{e1, 2}};
    e1->edges = {{t1, 2}};
    t1->edges = {{c2, 3}, {d2, 1}};
    t2->edges = {{h1, 1}, {e1, 3}};
    h1->edges = {{k1, 1}};
    d1->edges = {{v1, 1}};
    d2->edges = {{m2, 1}};
    auto concensus_nodes = FindConcensus(start_nodes);
    std::string concensus_string;
    for (auto *node : concensus_nodes) {
      concensus_string += matrix.position_to_letter(node->letter);
    }
    assert(concensus_string == "PKMDVRPQKNETC");
  }
  {
    CREATE(p1, 'P');
    CREATE(k1, 'K');
    CREATE(m1, 'M');
    CREATE(c1, 'C');
    CREATE(v1, 'V');
    CREATE(r1, 'R');
    CREATE(p2, 'P');
    CREATE(q1, 'Q');
    CREATE(k2, 'K');
    CREATE(n1, 'N');
    CREATE(e1, 'E');
    CREATE(t1, 'T');
    CREATE(c2, 'C');
    CREATE(t2, 'T');
    CREATE(h1, 'H');
    CREATE(d1, 'D');
    CREATE(d2, 'D');
    CREATE(m2, 'M');
    std::vector<Node *> start_nodes = {p1, t2};
    p1->edges = {{k1, 2}};
    k1->edges = {{m1, 2}};
    m1->edges = {{c1, 1}, {d1, 2}};
    c1->edges = {{v1, 1}};
    v1->edges = {{r1, 2}};
    r1->edges = {{p2, 1}, {n1, 2}};
    p2->edges = {{q1, 1}};
    q1->edges = {{k2, 1}};
    k2->edges = {{n1, 1}};
    n1->edges = {{e1, 2}};
    e1->edges = {{t1, 2}};
    t1->edges = {{c2, 3}, {d2, 1}};
    t2->edges = {{h1, 1}};
    h1->edges = {{k1, 1}};
    d1->edges = {{v1, 1}};
    d2->edges = {{m2, 1}};
    auto concensus_nodes = FindConcensus(start_nodes);
    std::string concensus_string;
    for (auto *node : concensus_nodes) {
      concensus_string += matrix.position_to_letter(node->letter);
    }
    assert(concensus_string == "PKMDVRNETC");
  }
}

int main() {
  const char matrix_path[] = "data/matrix/blosum50.mat";
  std::vector<std::unique_ptr<Node>> storage;
  ScoreMatrix matrix(ReadFile(matrix_path));
  TestAlignSequenceToGraph(matrix, storage);
  printf("OK!\n");
  return 0;
}
