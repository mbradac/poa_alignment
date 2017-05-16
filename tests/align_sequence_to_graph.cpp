#include <cassert>
#include <cstdio>
#include "sequence.hpp"
#include "poa_alignment.hpp"
#include <cassert>
#include <sstream>
#include <fstream>
#include <stack>
#include <unordered_map>
#include <algorithm>

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
  assert(n1->edges.size() == n2->edges.size());
  assert(n1->letter == n2->letter);
  auto it = visited.find(n1);
  if (it != visited.end()) {
    assert(it->second == n2);
    return;
  }
  visited[n1] = n2;
  auto cmp = [](const Node *n, const Node *m) { return n->letter < m->letter; };
  std::sort(n1->edges.begin(), n1->edges.end(), cmp);
  std::sort(n2->edges.begin(), n2->edges.end(), cmp);
  for (int i = 0; i < static_cast<int>(n1->edges.size()); ++i) {
    Check(n1->edges[i], n2->edges[i], visited);
  }
}

void Run(const std::string &s1, const std::string &s2, Node *expected_graph) {
  Sequence sequence1;
  sequence1.sequence = s1;
  Sequence sequence2;
  sequence2.sequence = s2;
  const char matrix_path[] = "data/matrix/blosum50.mat";
  Storage<Node> storage;
  ScoreMatrix matrix(ReadFile(matrix_path));
  auto graph = GraphFromSequence(sequence1, matrix, storage);
  AlignSequenceToGraph(graph, sequence2, matrix, 1, storage);
  std::unordered_map<Node *, Node *> visited;
  Check(graph.first, expected_graph, visited);
}
}

#define CREATE(VAR, LETTER)     \
  auto *VAR = storage.Create(); \
  VAR->letter = LETTER == -1 ? -1 : matrix.get_position(LETTER);

int main() {
  Storage<Node> storage;
  const char matrix_path[] = "data/matrix/blosum50.mat";
  ScoreMatrix matrix(ReadFile(matrix_path));
  {
    CREATE(start, -1);
    CREATE(d1, 'D');
    CREATE(d2, 'D');
    CREATE(d3, 'D');
    CREATE(c1, 'C');
    CREATE(c2, 'C');
    CREATE(p1, 'P');
    CREATE(end, -1);
    start->edges = {d1, c1};
    d1->edges = {d2};
    d2->edges = {d3};
    d3->edges = {end};
    c1->edges = {c2};
    c2->edges = {p1};
    p1->edges = {end};
    Run("DDD", "CCP", start);
  }
  {
    CREATE(start, -1);
    CREATE(p1, 'P');
    CREATE(p2, 'P');
    CREATE(t1, 'T');
    CREATE(t2, 'T');
    CREATE(t3, 'T');
    CREATE(k1, 'K');
    CREATE(k2, 'K');
    CREATE(k3, 'K');
    CREATE(end, -1);
    start->edges = {p1, p1};
    p1->edges = {p2};
    p2->edges = {t1, k1};
    t1->edges = {t2};
    t2->edges = {t3};
    t3->edges = {end};
    k1->edges = {k2};
    k2->edges = {k3};
    k3->edges = {end};
    Run("PPTTT", "PPKKK", start);
  }
  {
    CREATE(start, -1);
    CREATE(p1, 'P');
    CREATE(d1, 'D');
    CREATE(m1, 'M');
    CREATE(c1, 'C');
    CREATE(end, -1);
    start->edges = {p1, p1};
    p1->edges = {d1, c1};
    d1->edges = {m1};
    m1->edges = {end};
    c1->edges = {m1};
    Run("PDM", "PCM", start);
  }
  {
    CREATE(start, -1);
    CREATE(p1, 'P');
    CREATE(p2, 'P');
    CREATE(c1, 'C');
    CREATE(end, -1);
    start->edges = {p1, p2};
    p1->edges = {p2};
    p2->edges = {c1};
    c1->edges = {end};
    Run("PPC", "PC", start);
  }
  {
    CREATE(start, -1);
    CREATE(d1, 'D');
    CREATE(d2, 'D');
    CREATE(p1, 'P');
    CREATE(c1, 'C');
    CREATE(c2, 'C');
    CREATE(end, -1);
    start->edges = {d1, c1};
    d1->edges = {d2};
    d2->edges = {p1};
    p1->edges = {end};
    c1->edges = {c2};
    c2->edges = {p1};
    Run("DDP", "CCP", start);
  }
  {
    CREATE(start, -1);
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
    CREATE(end, -1);
    start->edges = {p1, t2};
    p1->edges = {k1};
    k1->edges = {m1};
    m1->edges = {c1, d1};
    c1->edges = {v1};
    v1->edges = {r1};
    r1->edges = {p2, n1};
    p2->edges = {q1};
    q1->edges = {k2};
    k2->edges = {n1};
    n1->edges = {e1};
    e1->edges = {t1};
    t1->edges = {c2, d2};
    c2->edges = {end};
    t2->edges = {h1};
    h1->edges = {k1};
    d1->edges = {v1};
    d2->edges = {m2};
    m2->edges = {end};
    Run("PKMCVRPQKNETC", "THKMDVRNETDM", start);
  }
  printf("OK!\n");
  return 0;
}
