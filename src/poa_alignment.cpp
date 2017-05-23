#include "poa_alignment.hpp"

#include <vector>
#include <memory>
#include <stack>
#include <unordered_map>
#include <unordered_set>
#include <cassert>
#include <algorithm>

namespace std {
template <>
struct hash<std::pair<int, poa_alignment::Node *>> {
  size_t operator()(const std::pair<int, poa_alignment::Node *> &p) const {
    return p.first * 1337 + reinterpret_cast<size_t>(p.second);
  }
};
}

namespace poa_alignment {

std::vector<Node *> TopologicalSort(const std::vector<Node *> &start_nodes) {
  // Find number of incoming edges for every node.
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

  // Find topological sorting of graph.
  for (auto *start_node : start_nodes) {
    jobs.push(start_node);
  }
  std::vector<Node *> sorted;
  while (!jobs.empty()) {
    auto *node = jobs.top();
    jobs.pop();
    sorted.push_back(node);
    for (auto _next : node->edges) {
      auto *next = _next.first;
      if (!--remaining_edges[next]) {
        jobs.push(next);
      }
    }
  }
  return sorted;
}

Graph GraphFromSequence(Sequence sequence, const ScoreMatrix &matrix) {
  TranslateSequence(&sequence, matrix);
  assert(sequence.sequence.size());
  Graph graph;
  graph.InsertSequence(nullptr, sequence.sequence);
  return graph;
}

void AlignSequenceToGraph(Graph &graph, Sequence sequence,
                          const ScoreMatrix &matrix, int gap_penalty) {
  if (!sequence.sequence.size()) return;
  assert(gap_penalty >= 0);
  TranslateSequence(&sequence, matrix);
  std::vector<Node *> nodes = TopologicalSort(graph.start_nodes);
  using State = std::pair<int, Node *>;
  std::unordered_map<State, int> dp;
  std::unordered_map<State, State> prev;
  auto update = [&](int next_i, Node *next_node, int i, Node *node, int d) {
    auto new_score = dp[{i, node}] + d;
    auto &old_score = dp[{next_i, next_node}];
    if (old_score < new_score) {
      old_score = new_score;
      prev[{next_i, next_node}] = {i, node};
    }
  };

  // Smith waterman.
  for (Node *start_node : graph.start_nodes) {
    auto score = matrix.get_score(sequence.sequence[0], start_node->letter);
    dp[{0, start_node}] = std::max(0, score);
  }
  for (int i = 0; i < static_cast<int>(sequence.sequence.size()); ++i) {
    for (auto *node : nodes) {
      for (auto next : node->edges) {
        update(i, next.first, i - 1, next.first, -gap_penalty);
        update(i, next.first, i, node, -gap_penalty);
        update(i, next.first, i - 1, node,
               matrix.get_score(sequence.sequence[i], next.first->letter));
      }
    }
  }

  // Find best score.
  // TODO: Change this once global alignment is added.
  std::pair<int, Node *> state(0, nullptr);
  for (auto score : dp) {
    if (dp[state] < score.second) {
      state = score.first;
    }
  }

  // Find alignment.
  std::vector<State> alignment_states;
  while (dp[state]) {
    alignment_states.push_back(state);
    state = prev[state];
  }
  std::reverse(alignment_states.begin(), alignment_states.end());

  if (alignment_states.size() == 0U) {
    // No alignment, insert whole sequence.
    graph.InsertSequence(nullptr, sequence.sequence);
    return;
  }

  // Insert unmatched beggining of the sequence.
  auto *sequence_node = graph.InsertSequence(
      nullptr, sequence.sequence.substr(0, alignment_states[0].first));

  // Insert nodes for the matched part.
  State prev_state(alignment_states[0].first - 1, nullptr);
  for (auto state : alignment_states) {
    if (state.first != prev_state.first && state.second != prev_state.second) {
      if (!sequence_node) {
        // Can happen only as a first node in a match if first letter of the
        // sequnce is matched, i. e. there is no unmatched beginning of the
        // sequence.
        sequence_node = state.second;
      } else if (sequence_node == prev_state.second) {
        auto got = graph.AddEdge(sequence_node, sequence.sequence[state.first]);
        sequence_node = got.first;
      } else {
        auto *next_node = [&]() {
          int letter = sequence.sequence[state.first];
          if (letter == state.second->letter) {
            return state.second;
          }
          return graph.AddNode(letter);
        }();
        sequence_node->edges.emplace_back(next_node, 1);
        sequence_node = next_node;
      }
    } else if (state.first != prev_state.first) {
      auto got = graph.AddEdge(sequence_node, sequence.sequence[state.first]);
      assert(got.second);
      sequence_node = got.first;
    }
    prev_state = state;
  }

  graph.InsertSequence(sequence_node,
                       sequence.sequence.substr(prev_state.first + 1));
}
}
