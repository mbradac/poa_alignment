#ifndef POA_ALIGNMENT_POA_ALIGNMENT_H_
#define POA_ALIGNMENT_POA_ALIGNMENT_H_

#include <vector>
#include <memory>

#include "sequence.hpp"

namespace poa_alignment {

template <typename T>
class Storage {
 public:
  T *Create() {
    storage_.emplace_back(new T());
    return storage_.back().get();
  }

 private:
  std::vector<std::unique_ptr<T>> storage_;
};

class Node {
 public:
  char letter = -1;
  std::vector<Node *> edges;

 private:
  Node() = default;
  friend class Storage<Node>;
};

std::vector<Node *> TopologicalSort(Node *start_node);

std::pair<Node *, Node *> GraphFromSequence(Sequence sequence,
                                            const ScoreMatrix &matrix,
                                            Storage<Node> &storage);

void AlignSequenceToGraph(std::pair<Node *, Node *> start_node,
                          Sequence sequence, const ScoreMatrix &matrix,
                          int gap_penalty, Storage<Node> &storage);
}
#endif  // POA_ALIGNMENT_POA_ALIGNMENT_H_
