#ifndef POA_ALIGNMENT_POA_ALIGNMENT_H_
#define POA_ALIGNMENT_POA_ALIGNMENT_H_

#include <vector>
#include <memory>
#include <cassert>

#include "sequence.hpp"

namespace poa_alignment {

struct Node {
  Node(char letter) : letter(letter) {}
  char letter = -1;
  // Next node and weight.
  std::vector<std::pair<Node *, int>> edges;
};

class NodeStorage {
 public:
  Node *AddNode(int letter) {
    storage_.emplace_back(new Node(letter));
    return storage_.back().get();
  }

 private:
  std::vector<std::unique_ptr<Node>> storage_;
};

class Graph {
 public:
  std::vector<Node *> start_nodes;

  // Storage must outlive Graph.
  Graph(NodeStorage *storage, Sequence sequence)
      : Graph(storage, sequence, [&]() {
          assert(sequence.sequence.size());
          return std::vector<int>(
              static_cast<int>(sequence.sequence.size()) - 1, 1);
        }()) {}
  Graph(NodeStorage *storage, Sequence sequence,
        const std::vector<int> &weights);

  Node *AddNode(int letter) { return storage_.AddNode(letter); }

  Node *AddStartNode(int letter) {
    Node *node = AddNode(letter);
    start_nodes.push_back(node);
    return node;
  }

  std::pair<Node *, bool> AddEdge(Node *from, int letter, int weight) {
    for (auto &edge : from->edges) {
      if (edge.first->letter == letter) {
        edge.second += weight;
        return {edge.first, false};
      }
    }
    auto *new_node = storage_.AddNode(letter);
    from->edges.emplace_back(new_node, weight);
    return {new_node, true};
  }

  void RemovePath(const std::vector<Node *> &_path);

  std::vector<Node *> InsertSequence(Node *prev, std::string sequence,
                                     const std::vector<int> &weights) {
    std::vector<Node *> sequence_nodes;
    if (sequence.size() == 0U) return {};
    if (!prev) {
      prev = AddStartNode(sequence[0]);
      sequence_nodes.push_back(prev);
      assert(weights.size() + 1 == sequence.size());
      sequence = sequence.substr(1);
    } else {
      assert(weights.size() == sequence.size());
    }
    int n = sequence.size();
    for (int i = 0; i < n; ++i) {
      prev = AddEdge(prev, sequence[i], weights[i]).first;
      sequence_nodes.push_back(prev);
    }
    return sequence_nodes;
  }

 private:
  NodeStorage &storage_;
};

std::vector<Node *> TopologicalSort(const std::vector<Node *> &start_nodes);

// void AlignGraphToGraph(Graph &graph1, Graph &graph2);

std::vector<Node *> AlignSequenceToGraph(Graph &graph, Sequence sequence,
                                         std::vector<int> weights,
                                         const ScoreMatrix &matrix,
                                         int gap_penalty);

std::vector<Node *> FindConcensus(const std::vector<Node *> &start_nodes);
}
#endif  // POA_ALIGNMENT_POA_ALIGNMENT_H_
