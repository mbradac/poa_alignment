#include <algorithm>
#include <string>
#include <vector>
#include "sequence.hpp"

using namespace poa_alignment;

typedef long long llint;
std::vector<int> SmithWatermanNoSimd(Sequence query,
                                     std::vector<Sequence> database,
                                     const ScoreMatrix &matrix, int q, int r) {
  TranslateSequence(&query, matrix, 1);
  int query_length = query.sequence.size();
  std::vector<int> results;

  for (Sequence &sequence : database) {
    TranslateSequence(&sequence, matrix, 1);
    std::vector<int> h(query_length + 1), f(query_length + 1);
    std::vector<int> prev_h(query_length + 1);
    results.push_back(0);
    for (int i = 0; i < static_cast<int>(sequence.sequence.size()); ++i) {
      int e = 0;
      for (int j = 0; j < query_length; ++j) {
        f[j + 1] = std::max(prev_h[j + 1] - q, f[j + 1] - r);
        e = std::max(h[j] - q, e - r);
        h[j + 1] = std::max(
            prev_h[j] +
                matrix.get_score(sequence.sequence[i], query.sequence[j]),
            static_cast<int>(std::max(e, std::max(f[j + 1], 0))));
        results.back() = std::max(results.back(), h[j + 1]);
      }
      prev_h = h;
    }
  }

  return results;
}
