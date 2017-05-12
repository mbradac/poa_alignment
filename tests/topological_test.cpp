#include <cassert>
#include <memory>
#include "poa_alignment.hpp"

using namespace poa_alignment;

int main() {
  Storage<Node> storage;
  auto node_a = storage.Create();
  node_a->letter = 'a';
  auto node_b = storage.Create();
  node_b->letter = 'b';
  node_b->edges.push_back(node_a);
  auto node_c = storage.Create();
  node_c->letter = 'c';
  node_c->edges.push_back(node_b);
  node_c->edges.push_back(node_a);
  auto node_d = storage.Create();
  node_d->letter = 'd';
  node_d->edges.push_back(node_b);
  node_d->edges.push_back(node_a);
  node_d->edges.push_back(node_c);
  auto node_e = storage.Create();
  node_e->letter = 'e';
  node_e->edges.push_back(node_d);
  node_e->edges.push_back(node_a);
  node_e->edges.push_back(node_c);
  node_e->edges.push_back(node_b);
  auto sorted = TopologicalSort(node_e);
  assert(sorted.size() == 5U);
  assert(sorted[0]->letter == 'e');
  assert(sorted[1]->letter == 'd');
  assert(sorted[2]->letter == 'c');
  assert(sorted[3]->letter == 'b');
  assert(sorted[4]->letter == 'a');
  return 0;
}
