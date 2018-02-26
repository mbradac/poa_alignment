#include <cstdio>
#include <sstream>
#include <fstream>
#include <stack>
#include <unordered_map>
#include <algorithm>
#include <iostream>
#include <cstdlib>

#include "sequence.hpp"
#include "poa_alignment.hpp"

#define Assert(x)   \
  do {              \
    if (!(x)) {     \
      std::abort(); \
    }               \
  } while (0)

namespace {

using namespace poa_alignment;

void Check(Node *n1, Node *n2, std::unordered_map<Node *, Node *> &visited) {
  Assert(n1->letter == n2->letter);
  if (n1->edges.size() != n2->edges.size()) {
    Assert(n1->edges.size() == n2->edges.size());
  }
  Assert(n1->edges.size() == n2->edges.size());
  auto it = visited.find(n1);
  if (it != visited.end()) {
    Assert(it->second == n2);
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
    Assert(n1->edges[i].second == n2->edges[i].second);
    Check(n1->edges[i].first, n2->edges[i].first, visited);
  }
}

void Run(const std::string &s1, const std::string &s2,
         std::vector<Node *> expected_graph, std::vector<int> weights,
         std::vector<Node *> expected_path) {
  Sequence sequence1;
  sequence1.sequence = s1;
  Sequence sequence2;
  sequence2.sequence = s2;
  NodeStorage storage;
  Graph graph(&storage, sequence1);
  auto path = AlignSequenceToGraph(graph, sequence2, weights, 5, 4, 6);
  std::unordered_map<Node *, Node *> visited;
  Assert(graph.start_nodes.size() == expected_graph.size());
  auto cmp =
      [](const Node *n1, const Node *n2) { return n1->letter < n2->letter; };
  std::sort(graph.start_nodes.begin(), graph.start_nodes.end(), cmp);
  std::sort(expected_graph.begin(), expected_graph.end(), cmp);
  for (int i = 0; i < static_cast<int>(graph.start_nodes.size()); ++i) {
    Check(graph.start_nodes[i], expected_graph[i], visited);
  }
  Assert(path.size() == expected_path.size());
  for (int i = 0; i < static_cast<int>(expected_path.size()); ++i) {
    Assert(visited[path[i]] == expected_path[i]);
  }
}

void Run(const std::string &s1, const std::string &s2,
         std::vector<Node *> expected_graph,
         std::vector<Node *> expected_path) {
  Run(s1, s2, expected_graph, std::vector<int>(s2.size() - 1, 1),
      expected_path);
}

#define CREATE(VAR, LETTER)     \
  auto *VAR = new Node(LETTER); \
  storage.emplace_back(VAR);

void TestAlignSequenceToGraph(std::vector<std::unique_ptr<Node>> &storage) {
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
    Run("DDD", "CCP", start_nodes, {c1, c2, p1});
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
    Run("PPTTT", "PPKKK", start_nodes, {p1, p2, k1, k2, k3});
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
    Run("PDM", "PCM", start_nodes, {p1, c1, m1});
  }
  {
    CREATE(p1, 'P');
    CREATE(p2, 'P');
    CREATE(c1, 'C');
    std::vector<Node *> start_nodes = {p1};
    p1->edges = {{p2, 1}};
    p2->edges = {{c1, 2}};
    Run("PPC", "PC", start_nodes, {p2, c1});
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
    Run("DDP", "CCP", start_nodes, {c1, c2, p1});
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
    Run("PKMCVRPQKNETC", "THKMDVRNETDM", start_nodes,
        {t2, h1, k1, m1, d1, v1, r1, n1, e1, t1, d2, m2});
  }
  {
    std::vector<int> weights = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
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
    k1->edges = {{m1, 4}};
    m1->edges = {{c1, 1}, {d1, 4}};
    c1->edges = {{v1, 1}};
    v1->edges = {{r1, 7}};
    r1->edges = {{p2, 1}, {n1, 7}};
    p2->edges = {{q1, 1}};
    q1->edges = {{k2, 1}};
    k2->edges = {{n1, 1}};
    n1->edges = {{e1, 9}};
    e1->edges = {{t1, 10}};
    t1->edges = {{c2, 1}, {d2, 10}};
    t2->edges = {{h1, 1}};
    h1->edges = {{k1, 2}};
    d1->edges = {{v1, 5}};
    d2->edges = {{m2, 11}};
    Run("PKMCVRPQKNETC", "THKMDVRNETDM", start_nodes, weights,
        {t2, h1, k1, m1, d1, v1, r1, n1, e1, t1, d2, m2});
  }
  printf("TestAlignSequenceToGraph OK!\n");
}

void TestFindConcensus(std::vector<std::unique_ptr<Node>> &storage) {
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
      concensus_string += node->letter;
    }
    Assert(concensus_string == "PKMDVRPQKNETC");
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
      concensus_string += node->letter;
    }
    Assert(concensus_string == "PKMDVRNETC");
  }
  printf("TestFindConcensus OK!\n");
}

void TestGraphDisassemble(std::vector<std::unique_ptr<Node>> &storage) {
  auto RemovePath =
      [](const std::vector<Node *> &path, std::vector<Node *> &start_nodes) {
        if (!path.size()) return;

        std::stack<Node *> jobs;
        std::unordered_map<Node *, int> remaining_edges;
        for (auto *start_node : start_nodes) {
          jobs.push(start_node);
        }
        while (!jobs.empty()) {
          auto *node = jobs.top();
          jobs.pop();
          for (auto _next : node->edges) {
            auto *next = _next.first;
            bool visited = remaining_edges[next];
            remaining_edges[next] += 1;
            if (!visited) {
              jobs.push(next);
            }
          }
        }

        int n = path.size();
        for (int i = 0; i < n; ++i) {
          path[i]->marked = true;
        }
        for (int i = 0; i < n - 1; ++i) {
          auto it = std::find_if(
              path[i]->edges.begin(), path[i]->edges.end(),
              [&](std::pair<Node *, int> a) { return a.first == path[i + 1]; });
          Assert(it != path[i]->edges.end());
          if (!--remaining_edges[it->first]) {
            start_nodes.push_back(it->first);
          }
          path[i]->edges.erase(it);
        }
      };

  auto chain_to_sequence_weight = [](const std::vector<Node *> &chain) {
    Sequence sequence;
    std::vector<int> weights;
    int n = chain.size();
    for (int i = 0; i < n; ++i) {
      sequence.sequence += chain[i]->letter;
      if (i + 1 == n) continue;
      for (auto edge : chain[i]->edges) {
        if (edge.first == chain[i + 1]) {
          weights.push_back(edge.second);
        }
      }
    }
    Assert(weights.size() + 1 == sequence.sequence.size());
    return std::make_pair(sequence, weights);
  };
  {
    CREATE(n2, 'G');
    CREATE(n3, 'T');
    CREATE(n4, 'C');
    CREATE(n9, 'C');
    CREATE(n5, 'A');
    CREATE(n6, 'T');
    CREATE(n0, 'C');
    CREATE(n1, 'C');
    CREATE(n7, 'G');
    CREATE(n8, 'G');
    n0->edges = {{n1, 1}};
    n1->edges = {{n2, 3}};
    n2->edges = {{n3, 3}, {n9, 1}, {n7, 1}};
    n4->edges = {{n3, 2}};
    n9->edges = {{n4, 1}};
    n5->edges = {{n6, 3}};
    n6->edges = {{n2, 2}};
    n7->edges = {{n8, 1}};
    std::vector<Node *> start_nodes = {n0, n5};

    int j = 0;
    std::vector<std::string> expected_concensus = {"CCGT", "ATG", "GCCT",
                                                   "GGG"};
    std::vector<std::vector<int>> weights = {
        {1, 3, 3}, {3, 2}, {1, 1, 2}, {1, 1}};
    while (true) {
      auto concensus = FindConcensus(start_nodes);
      if (concensus.size() <= 1U) break;
      auto got = chain_to_sequence_weight(concensus);
      std::string concensus_string = got.first.sequence;
      Assert(concensus_string == expected_concensus[j]);
      Assert(got.second == weights[j]);
      RemovePath(concensus, start_nodes);
      ++j;
    }
    Assert(j == static_cast<int>(expected_concensus.size()));
    printf("TestGraphDisassemble OK!\n");
  }
}
#undef CREATE

void TestAlignGraphToGraph() {
  NodeStorage storage;

  Graph g1(&storage);
  auto n2 = storage.AddNode('G');
  auto n3 = storage.AddNode('T');
  auto n4 = storage.AddNode('C');
  auto n9 = storage.AddNode('C');
  auto n5 = storage.AddNode('A');
  auto n6 = storage.AddNode('T');
  auto n0 = storage.AddNode('C');
  auto n1 = storage.AddNode('C');
  auto n7 = storage.AddNode('G');
  auto n8 = storage.AddNode('G');
  n0->edges = {{n1, 1}};
  n1->edges = {{n2, 3}};
  n2->edges = {{n3, 3}, {n9, 1}, {n7, 1}};
  n4->edges = {{n3, 2}};
  n9->edges = {{n4, 1}};
  n5->edges = {{n6, 3}};
  n6->edges = {{n2, 2}};
  n7->edges = {{n8, 1}};
  g1.start_nodes = std::vector<Node *>{n0, n5};

  Graph g2(&storage);
  auto bn0 = storage.AddNode('C');
  auto bn1 = storage.AddNode('C');
  auto bn2 = storage.AddNode('G');
  auto bn3 = storage.AddNode('C');
  auto bn4 = storage.AddNode('T');
  auto bn5 = storage.AddNode('A');
  auto bn6 = storage.AddNode('T');
  auto bn7 = storage.AddNode('G');
  bn0->edges = {{bn1, 2}};
  bn1->edges = {{bn2, 2}, {bn5, 1}};
  bn2->edges = {{bn3, 3}, {bn7, 2}};
  bn3->edges = {{bn4, 3}};
  bn5->edges = {{bn6, 1}};
  bn6->edges = {{bn3, 1}};
  g2.start_nodes = std::vector<Node *>{bn0};

  Graph g3(&storage);
  auto cn0 = storage.AddNode('C');
  auto cn1 = storage.AddNode('C');
  auto cn2 = storage.AddNode('G');
  auto cn3 = storage.AddNode('C');
  auto cn4 = storage.AddNode('T');
  auto cn5 = storage.AddNode('A');
  auto cn6 = storage.AddNode('T');
  auto cn7 = storage.AddNode('G');
  auto cn8 = storage.AddNode('A');
  auto cn9 = storage.AddNode('T');
  auto cn10 = storage.AddNode('C');
  auto cn11 = storage.AddNode('G');
  cn0->edges = {{cn1, 3}};
  cn1->edges = {{cn2, 5}, {cn5, 1}};
  cn2->edges = {{cn3, 4}, {cn4, 3}, {cn7, 3}};
  cn3->edges = {{cn4, 3}, {cn10, 1}};
  cn5->edges = {{cn6, 1}};
  cn6->edges = {{cn3, 1}};
  cn8->edges = {{cn9, 3}};
  cn9->edges = {{cn2, 2}};
  cn10->edges = {{cn4, 2}};
  cn7->edges = {{cn11, 1}};
  g3.start_nodes = std::vector<Node *>{cn0, cn8};

  AlignGraphToGraph(g2, g1, 2, 2, 1);

  Assert(g2.start_nodes.size() == g3.start_nodes.size());
  auto cmp = [](Node *n, Node *m) { return n->letter < m->letter; };

  std::sort(g2.start_nodes.begin(), g2.start_nodes.end(), cmp);
  std::sort(g3.start_nodes.begin(), g3.start_nodes.end(), cmp);
  std::unordered_map<Node *, Node *> visited;
  for (int i = 0; i < static_cast<int>(g2.start_nodes.size()); ++i) {
    Check(g2.start_nodes[i], g3.start_nodes[i], visited);
  }
  printf("TestAlignGraphToGraph OK!\n");
}
}

int main() {
  std::vector<std::unique_ptr<Node>> storage;
  TestAlignSequenceToGraph(storage);
  TestFindConcensus(storage);
  TestGraphDisassemble(storage);
  TestAlignGraphToGraph();
  printf("OK!\n");
  return 0;
}
