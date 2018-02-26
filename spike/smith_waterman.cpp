#include <algorithm>
#include <string>
#include <vector>
#include "sequence.hpp"

using namespace poa_alignment;

std::vector<int> SmithWatermanNoSimd(Sequence query,
                                     std::vector<Sequence> database, int match,
                                     int mismatch, int gap_penalty) {
  int query_length = query.sequence.size();
  std::vector<int> results;

  auto get_score = [&](char a, char b) { return a == b ? match : -mismatch; };

  for (Sequence &sequence : database) {
    std::vector<int> h(query_length + 1), f(query_length + 1);
    std::vector<int> prev_h(query_length + 1);
    results.push_back(0);
    for (int i = 0; i < static_cast<int>(sequence.sequence.size()); ++i) {
      int e = 0;
      for (int j = 0; j < query_length; ++j) {
        f[j + 1] =
            std::max(prev_h[j + 1] - gap_penalty, f[j + 1] - gap_penalty);
        e = std::max(h[j] - gap_penalty, e - gap_penalty);
        h[j + 1] = std::max(
            prev_h[j] + get_score(sequence.sequence[i], query.sequence[j]),
            static_cast<int>(std::max(e, std::max(f[j + 1], 0))));
        results.back() = std::max(results.back(), h[j + 1]);
      }
      prev_h = h;
    }
  }

  return results;
}
