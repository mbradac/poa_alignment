#ifndef POA_ALIGNMENT_POA_ALIGNMENT_H_
#define POA_ALIGNMENT_POA_ALIGNMENT_H_

#include <vector>
#include <memory>

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
  // Storage must outlive Graph.
  Graph(NodeStorage *storage, Sequence sequence, const ScoreMatrix &matrix);
  std::vector<Node *> start_nodes;

  Node *AddNode(int letter) { return storage_.AddNode(letter); }

  Node *AddStartNode(int letter) {
    Node *node = AddNode(letter);
    start_nodes.push_back(node);
    return node;
  }

  std::pair<Node *, bool> AddEdge(Node *from, int letter) {
    for (auto &edge : from->edges) {
      if (edge.first->letter == letter) {
        edge.second += 1;
        return {edge.first, false};
      }
    }
    auto *new_node = storage_.AddNode(letter);
    from->edges.emplace_back(new_node, 1);
    return {new_node, true};
  }

  Node *InsertSequence(Node *prev, std::string sequence) {
    if (sequence.size() == 0U) return nullptr;
    if (!prev) {
      prev = AddStartNode(sequence[0]);
      sequence = sequence.substr(1);
    }
    for (char c : sequence) {
      prev = AddEdge(prev, c).first;
    }
    return prev;
  }

 private:
  NodeStorage &storage_;
};

std::vector<Node *> TopologicalSort(const std::vector<Node *> &start_nodes);

Graph GraphFromSequence(Sequence sequence, const ScoreMatrix &matrix);

void AlignSequenceToGraph(Graph &graph, Sequence sequence,
                          const ScoreMatrix &matrix, int gap_penalty);

std::vector<Node *> FindConcensus(const std::vector<Node *> &start_nodes);
}
#endif  // POA_ALIGNMENT_POA_ALIGNMENT_H_
