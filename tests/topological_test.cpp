#include <cassert>
#include <memory>
#include <vector>
#include "poa_alignment.hpp"

using namespace poa_alignment;
struct Storage {
  Node *Create(int letter) {
    storage.emplace_back(new Node(letter));
    return storage.back().get();
  }
  std::vector<std::unique_ptr<Node>> storage;
};

int main() {
  Storage storage;
  auto node_a = storage.Create('a');
  auto node_b = storage.Create('b');
  node_b->edges.emplace_back(node_a, 0);
  auto node_c = storage.Create('c');
  node_c->edges.emplace_back(node_b, 0);
  node_c->edges.emplace_back(node_a, 0);
  auto node_d = storage.Create('d');
  node_d->edges.emplace_back(node_b, 0);
  node_d->edges.emplace_back(node_a, 0);
  node_d->edges.emplace_back(node_c, 0);
  auto node_e = storage.Create('e');
  node_e->edges.emplace_back(node_d, 0);
  node_e->edges.emplace_back(node_a, 0);
  node_e->edges.emplace_back(node_c, 0);
  node_e->edges.emplace_back(node_b, 0);
  auto sorted = TopologicalSort(std::vector<Node *>{node_e});
  assert(sorted.size() == 5U);
  assert(sorted[0]->letter == 'e');
  assert(sorted[1]->letter == 'd');
  assert(sorted[2]->letter == 'c');
  assert(sorted[3]->letter == 'b');
  assert(sorted[4]->letter == 'a');
  return 0;
}
