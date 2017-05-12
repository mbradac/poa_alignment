#include "poa_alignment.hpp"

#include <vector>
#include <memory>
#include <stack>
#include <unordered_map>
#include <unordered_set>
#include <cassert>

namespace std {
template <>
struct hash<std::pair<int, poa_alignment::Node *>> {
  size_t operator()(const std::pair<int, poa_alignment::Node *> &p) const {
    return p.first * 1337 + reinterpret_cast<size_t>(p.second);
  }
};
}

namespace poa_alignment {

std::vector<Node *> TopologicalSort(Node *start_node) {
  // Find number of incoming edges for every node.
  std::stack<Node *> jobs;
  std::unordered_map<Node *, int> remaining_edges;
  jobs.push(start_node);
  while (!jobs.empty()) {
    auto *node = jobs.top();
    jobs.pop();
    for (auto *next : node->edges) {
      bool visited = remaining_edges[next];
      remaining_edges[next] += 1;
      if (!visited) {
        jobs.push(next);
      }
    }
  }

  // Find topological sorting of graph.
  jobs.push(start_node);
  std::vector<Node *> sorted;
  while (!jobs.empty()) {
    auto *node = jobs.top();
    jobs.pop();
    sorted.push_back(node);
    for (auto &next : node->edges) {
      if (!--remaining_edges[next]) {
        jobs.push(next);
      }
    }
  }
  return sorted;
}

std::pair<Node *, Node *> GraphFromSequence(Sequence sequence,
                                            const ScoreMatrix &matrix,
                                            Storage<Node> &storage) {
  TranslateSequence(&sequence, matrix);
  auto *start_node = storage.Create();
  start_node->letter = -1;
  auto last_node = start_node;
  auto add_node = [&](char letter) {
    auto node = storage.Create();
    node->letter = letter;
    last_node->edges.push_back(node);
    last_node = node;
  };
  for (auto x : sequence.sequence) {
    add_node(x);
  }
  add_node(-1);
  return {start_node, last_node};
}

void AlignSequenceToGraph(std::pair<Node *, Node *> graph, Sequence sequence,
                          const ScoreMatrix &matrix, int gap_penalty,
                          Storage<Node> &storage) {
  assert(gap_penalty >= 0);
  Node *const start_node = graph.first;
  Node *const end_node = graph.second;
  TranslateSequence(&sequence, matrix);
  std::vector<Node *> nodes = TopologicalSort(start_node);
  std::unordered_map<std::pair<int, Node *>, int> dp;
  std::unordered_map<std::pair<int, Node *>, std::pair<int, Node *>> prev;
  auto update = [&](int next_i, Node *next_node, int i, Node *node, int d) {
    auto new_score = dp[{i, node}] + d;
    auto &old_score = dp[{next_i, next_node}];
    if (old_score < new_score) {
      old_score = new_score;
      prev[{next_i, next_node}] = {i, node};
    }
  };

  // Smith waterman.
  for (int i = 0; i < static_cast<int>(sequence.sequence.size()); ++i) {
    for (auto *node : nodes) {
      update(i + 1, node, i, node, -gap_penalty);
      for (auto *next : node->edges) {
        update(i, next, i, node, -gap_penalty);
        update(i + 1, next, i, node,
               matrix.get_score(sequence.sequence[i], node->letter));
      }
    }
  }

  for (auto x : dp) {
    printf("%d %p: %d\n", x.first.first, x.first.second, x.second);
  }

  // Find best score.
  std::pair<int, Node *> state(0, nullptr);
  for (auto score : dp) {
    if (dp[state] < score.second) {
      state = score.first;
    }
  }

  // Create nodes for ending of the sequence.
  auto *node = state.second ? prev[state].second : start_node;
  for (int i = state.first; i < static_cast<int>(sequence.sequence.size());
       ++i) {
    auto new_node = storage.Create();
    new_node->letter = sequence.sequence[i];
    node->edges.push_back(new_node);
    node = new_node;
  }
  if (state.first != static_cast<int>(sequence.sequence.size())) {
    node->edges.push_back(end_node);
  }
  if (!state.second) return;

  // Create nodes for the matched part.
  auto *graph_node = state.second;
  auto *sequence_node = state.second;
  while (true) {
    auto prev_state = prev[state];
    printf("matched\n");
    printf("====== %d %p : %d\n", state.first, state.second, dp[state]);
    printf("====== %d %p : %d\n", prev_state.first, prev_state.second,
           dp[prev_state]);
    if (state.first != prev_state.first && state.second != prev_state.second) {
      printf("oba\n");
      if (sequence_node != graph_node) {
        prev_state.second->edges.push_back(sequence_node);
      }
      sequence_node = graph_node = prev_state.second;
    } else if (state.second == prev_state.second) {
      printf("graph isti %d\n", prev_state.first);
      auto *new_node = storage.Create();
      new_node->letter = sequence.sequence[prev_state.first];
      new_node->edges.push_back(sequence_node);
      sequence_node = new_node;
    } else {
      printf("sequence isti\n");
      graph_node = prev_state.second;
    }
    state = prev_state;
    if (dp[state] == 0) break;
  }

  if (state.second == start_node) return;
  node = start_node;
  for (int i = 0; i < state.first; ++i) {
    auto new_node = storage.Create();
    new_node->letter = sequence.sequence[i];
    node->edges.push_back(new_node);
    node = new_node;
  }
  node->edges.push_back(sequence_node);
}
}
