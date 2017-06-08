#include "poa_alignment.hpp"

#include <vector>
#include <memory>
#include <stack>
#include <unordered_map>
#include <unordered_set>
#include <cassert>
#include <algorithm>
#include <iostream>

namespace std {
template <>
struct hash<std::pair<int, poa_alignment::Node *>> {
  size_t operator()(const std::pair<int, poa_alignment::Node *> &p) const {
    return p.first * 1337 + reinterpret_cast<size_t>(p.second);
  }
};
}

enum class AlignmentType { LOCAL, FIXED, LEFT_FIXED, RIGHT_FIXED };

namespace poa_alignment {
using State = std::pair<int, Node *>;

void Graph::RemovePath(const std::vector<Node *> &path) {
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
    assert(it != path[i]->edges.end());
    if (!--remaining_edges[it->first]) {
      start_nodes.push_back(it->first);
    }
    path[i]->edges.erase(it);
  }
}

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

Graph::Graph(NodeStorage *storage, Sequence sequence,
             const std::vector<int> &weights)
    : storage_(*storage) {
  assert(sequence.sequence.size());
  InsertSequence(nullptr, sequence.sequence, weights);
}

auto CalculateAlignment(const std::vector<Node *> nodes,
                        const std::vector<Node *> start_nodes,
                        Sequence sequence, const ScoreMatrix &matrix,
                        int gap_penalty, AlignmentType type) {
  assert(gap_penalty >= 0);
  std::unordered_map<State, int> dp;
  std::unordered_map<State, State> prev;
  auto update = [&](int next_i, Node *next_node, int i, Node *node, int d) {
    if (type == AlignmentType::FIXED || type == AlignmentType::LEFT_FIXED) {
      dp.insert({{i, node}, -1e9});
      dp.insert({{next_i, next_node}, -1e9});
    }
    auto new_score = dp[{i, node}] + d;
    auto &old_score = dp[{next_i, next_node}];
    if (old_score < new_score) {
      old_score = new_score;
      prev[{next_i, next_node}] = {i, node};
    }
  };

  // Calculate alignment.
  for (Node *start_node : start_nodes) {
    auto score = matrix.get_score(sequence.sequence[0], start_node->letter);
    if (type == AlignmentType::FIXED || type == AlignmentType::LEFT_FIXED) {
      dp.insert({{0, start_node}, -1e9});
    }
    dp[{0, start_node}] = std::max(dp[{0, start_node}], score);
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
  auto state = [&]() {
    if (type == AlignmentType::FIXED || type == AlignmentType::RIGHT_FIXED) {
      return State(static_cast<int>(sequence.sequence.size()) - 1,
                   nodes.back());
    }
    State state(0, nullptr);
    for (auto score : dp) {
      if (dp[state] < score.second) {
        state = score.first;
      }
    }
    return state;
  }();

  // Find alignment.
  std::vector<State> alignment_states;
  while (
      ((type == AlignmentType::LOCAL || type == AlignmentType::RIGHT_FIXED) &&
       dp[state]) ||
      ((type == AlignmentType::FIXED || type == AlignmentType::LEFT_FIXED) &&
       state.second)) {
    alignment_states.push_back(state);
    state = prev[state];
  }
  std::reverse(alignment_states.begin(), alignment_states.end());

  return alignment_states;
}

// TODO: This API is terrible, implement GraphView class to pass subgraph to
// alignment function.
std::vector<Node *> AlignSequenceToGraph(Graph &graph, Node *start_node,
                                         Node *end_node, Sequence sequence,
                                         std::vector<int> weights,
                                         const ScoreMatrix &matrix,
                                         int gap_penalty) {
  if (!sequence.sequence.size()) return {};
  assert(weights.size() + 1 == sequence.sequence.size());
  auto nodes = [&]() {
    auto nodes = TopologicalSort(graph.start_nodes);
    auto begin_it = start_node
                        ? std::find(nodes.begin(), nodes.end(), start_node)
                        : nodes.begin();
    auto end_it =
        end_node ? std::next(std::find(nodes.begin(), nodes.end(), end_node))
                 : nodes.end();
    return std::vector<Node *>(begin_it, end_it);
  }();
  AlignmentType type = [&]() {
    if (start_node && end_node) return AlignmentType::FIXED;
    if (start_node) return AlignmentType::LEFT_FIXED;
    if (end_node) return AlignmentType::RIGHT_FIXED;
    return AlignmentType::LOCAL;
  }();
  auto start_nodes =
      start_node ? std::vector<Node *>{start_node} : graph.start_nodes;

  auto alignment_states = CalculateAlignment(nodes, start_nodes, sequence,
                                             matrix, gap_penalty, type);

  if (alignment_states.size() == 0U) {
    // No alignment, insert whole sequence.
    return graph.InsertSequence(nullptr, sequence.sequence, weights);
  }

  std::vector<Node *> sequence_nodes;
  // Insert unmatched beggining of the sequence.
  auto *sequence_node = [&]() -> Node *{
    if (!alignment_states[0].first) return nullptr;
    sequence_nodes = graph.InsertSequence(
        nullptr, sequence.sequence.substr(0, alignment_states[0].first),
        std::vector<int>(weights.begin(),
                         weights.begin() + alignment_states[0].first - 1));
    return sequence_nodes.back();
  }();

  // Insert nodes for the matched part.
  State prev_state(alignment_states[0].first - 1, nullptr);
  for (auto state : alignment_states) {
    if (state.first != prev_state.first && state.second != prev_state.second) {
      if (!sequence_node) {
        // Can happen only as a first node in a match if first letter of the
        // sequnce is matched, i. e. there is no unmatched beginning of the
        // sequence.
        sequence_node = state.second;
        sequence_nodes.push_back(sequence_node);
      } else if (sequence_node == prev_state.second) {
        auto got = graph.AddEdge(sequence_node, sequence.sequence[state.first],
                                 weights[state.first - 1]);
        sequence_node = got.first;
        sequence_nodes.push_back(sequence_node);
      } else {
        auto *next_node = [&]() {
          int letter = sequence.sequence[state.first];
          if (letter == state.second->letter) {
            return state.second;
          }
          return graph.AddNode(letter);
        }();
        auto found = false;
        for (auto &edge : sequence_node->edges) {
          if (edge.first == next_node) {
            edge.second += weights[state.first - 1];
            found = true;
          }
        }
        if (!found) {
          sequence_node->edges.emplace_back(next_node,
                                            weights[state.first - 1]);
        }
        sequence_node = next_node;
        sequence_nodes.push_back(sequence_node);
      }
    } else if (state.first != prev_state.first) {
      auto got = graph.AddEdge(sequence_node, sequence.sequence[state.first],
                               weights[state.first - 1]);
      assert(got.second);
      sequence_node = got.first;
      sequence_nodes.push_back(sequence_node);
    }
    prev_state = state;
  }

  // Insert end of the sequence.
  auto ending = graph.InsertSequence(
      sequence_node, sequence.sequence.substr(prev_state.first + 1),
      std::vector<int>(weights.begin() + prev_state.first, weights.end()));
  sequence_nodes.insert(sequence_nodes.end(), ending.begin(), ending.end());
  return sequence_nodes;
}

std::vector<Node *> AlignSequenceToGraph(Graph &graph, Sequence sequence,
                                         std::vector<int> weights,
                                         const ScoreMatrix &matrix,
                                         int gap_penalty) {
  return AlignSequenceToGraph(graph, nullptr, nullptr, sequence, weights,
                              matrix, gap_penalty);
}

std::vector<Node *> FindConcensus(const std::vector<Node *> &start_nodes) {
  std::vector<Node *> nodes = TopologicalSort(start_nodes);
  if (nodes.size() == 0U) return std::vector<Node *>();
  std::unordered_map<Node *, Node *> predecessors;
  std::unordered_map<Node *, int> scores;

  auto find_max_score_node = [&](const std::vector<Node *> &nodes,
                                 Node *node_on_path) {
    assert(!node_on_path || !node_on_path->marked);
    if (node_on_path) {
      for (auto *node : nodes) {
        predecessors[node] = nullptr;
      }
      for (auto score : scores) {
        if (score.first != node_on_path) {
          scores[score.first] = -1;
        }
      }
      for (auto edge : node_on_path->edges) {
        predecessors[edge.first] = node_on_path;
        scores[edge.first] = edge.second;
      }
    }

    Node *max_score_node = nullptr;
    for (auto *node : nodes) {
      if (scores[node] < 0) continue;
      if (predecessors[node] && !predecessors[node]->marked) {
        scores[node] += scores[predecessors[node]];
      }
      if (scores[node] > scores[max_score_node]) {
        max_score_node = node;
      }
      for (auto edge : node->edges) {
        auto node_score = node->marked ? 0 : scores[node];
        auto pred_score =
            predecessors[edge.first] && !predecessors[edge.first]->marked
                ? scores[predecessors[edge.first]]
                : 0;
        if (scores[edge.first] < edge.second ||
            (scores[edge.first] == edge.second && pred_score < node_score)) {
          scores[edge.first] = edge.second;
          predecessors[edge.first] = node;
        }
      }
    }
    return max_score_node;
  };

  auto *max_score_node = find_max_score_node(nodes, nullptr);
  if (!max_score_node) return {};
  while (!max_score_node->marked && max_score_node->edges.size()) {
    assert(max_score_node);
    nodes.erase(nodes.begin(),
                std::find(nodes.begin(), nodes.end(), max_score_node) + 1);
    max_score_node = find_max_score_node(nodes, max_score_node);
  }

  std::vector<Node *> concensus_nodes;
  for (auto *node = max_score_node; node; node = predecessors[node]) {
    concensus_nodes.push_back(node);
  }
  std::reverse(concensus_nodes.begin(), concensus_nodes.end());
  return concensus_nodes;
}

void AlignGraphToGraph(Graph &graph1, Graph &graph2, const ScoreMatrix &matrix,
                       int gap_penalty) {
  // TODO: this is also used in the test at the moment either expose it or
  // rewrite test
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
    assert(weights.size() + 1 == sequence.sequence.size());
    return std::make_pair(sequence, weights);
  };

  std::unordered_map<Node *, Node *> graph_nodes_mapping;

  while (true) {
    auto concensus = FindConcensus(graph2.start_nodes);
    if (concensus.size() <= 1U) break;
    auto got = chain_to_sequence_weight(concensus);
    graph2.RemovePath(concensus);
    auto start_node_it = graph_nodes_mapping.find(concensus[0]);
    Node *start_node = start_node_it == graph_nodes_mapping.end()
                           ? nullptr
                           : start_node_it->second;
    auto end_node_it = graph_nodes_mapping.find(concensus.back());
    Node *end_node = end_node_it == graph_nodes_mapping.end()
                         ? nullptr
                         : end_node_it->second;
    // TODO: send real subgraph here, don't fake it by sending whole graph and
    // start and end node
    auto path = AlignSequenceToGraph(graph1, start_node, end_node, got.first,
                                     got.second, matrix, gap_penalty);
    assert(concensus.size() == path.size());
    for (int i = 0; i < static_cast<int>(path.size()); ++i) {
      graph_nodes_mapping[concensus[i]] = path[i];
    }
  }
}
}
