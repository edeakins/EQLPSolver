#ifndef UTIL_OCMATCHING_H
#define UTIL_OCMATCHING_H

#include <vector>

/* Maximum cardinality bipartite matching via Kuhn's augmenting-path
 * algorithm. Used by orbital crossover to match nonbasic residual (r)
 * variables with degenerate variables they can be swapped against, so that
 * degenerate crossover pivots can be performed combinatorially rather than by
 * running the simplex (see matching.pdf).
 *
 * adj[u] lists the right-node (residual) indices adjacent to left node u.
 * nRight is the number of right nodes (residuals).
 * On return, matchRight[v] holds the left node matched to right node v, or -1.
 * The function returns the number of matched right nodes.
 */
class OCBipartiteMatching {
 public:
  static int solve(const std::vector<std::vector<int> >& adj, int nRight,
                   std::vector<int>& matchRight) {
    const int nLeft = static_cast<int>(adj.size());
    matchRight.assign(nRight, -1);
    std::vector<char> used;
    int matched = 0;
    for (int u = 0; u < nLeft; ++u) {
      used.assign(nRight, 0);
      if (tryAugment(u, adj, used, matchRight)) ++matched;
    }
    return matched;
  }

 private:
  static bool tryAugment(int u, const std::vector<std::vector<int> >& adj,
                         std::vector<char>& used, std::vector<int>& matchRight) {
    for (size_t e = 0; e < adj[u].size(); ++e) {
      const int v = adj[u][e];
      if (used[v]) continue;
      used[v] = 1;
      if (matchRight[v] == -1 ||
          tryAugment(matchRight[v], adj, used, matchRight)) {
        matchRight[v] = u;
        return true;
      }
    }
    return false;
  }
};

#endif
