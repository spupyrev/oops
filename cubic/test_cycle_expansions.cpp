#include "graph_algorithms.h"
#include "logging.h"

#include <algorithm>
#include <array>
#include <set>
#include <string>
#include <utility>
#include <vector>

using Crossing = std::pair<int, int>;

/// Return whether all vertices incident with an edge lie in one component.
bool connectedSupport(const int n, const std::vector<EdgeTy>& edges) {
  std::vector<bool> active(n, false);
  for (const auto& [first, second] : edges) {
    active[first] = true;
    active[second] = true;
  }
  std::vector<bool> reached(n, false);
  const auto startIt = std::find(active.begin(), active.end(), true);
  if (startIt == active.end())
    return false;
  const int start = static_cast<int>(startIt - active.begin());
  reached[start] = true;
  bool changed = true;
  while (changed) {
    changed = false;
    for (const auto& [first, second] : edges) {
      if (reached[first] == reached[second])
        continue;
      reached[first] = reached[second] = true;
      changed = true;
    }
  }
  for (int vertex = 0; vertex < n; vertex++)
    if (active[vertex] && !reached[vertex])
      return false;
  return true;
}

/// One local drawing listed in an expansion lemma.
struct LocalDrawing {
  std::string boundaryOrder;
  std::vector<Crossing> crossings;
};

/// Verify a listed local drawing by planarity after replacing each crossing.
void checkLocalDrawing(const int k, const LocalDrawing& drawing) {
  CHECK(static_cast<int>(drawing.boundaryOrder.size()) == k);
  std::vector<EdgeTy> localEdges;
  auto makeEdge = [](const int u, const int v) {
    return EdgeTy(std::min(u, v), std::max(u, v));
  };
  for (int i = 0; i < k; i++)
    localEdges.push_back(makeEdge(i, (i + 1) % k));
  for (int i = 0; i < k; i++)
    localEdges.push_back(makeEdge(i, k + i));

  std::vector<bool> crossedEdges(2 * k, false);
  std::vector<EdgeTy> planarized;
  struct CrossingWheel {
    std::array<int, 4> rim;
    int apex;
  };
  std::vector<CrossingWheel> crossingWheels;
  int nextVertex = 2 * k + 1;
  for (const auto& [first, second] : drawing.crossings) {
    CHECK(0 <= first && first < 2 * k && 0 <= second && second < 2 * k);
    CHECK(first != second && !crossedEdges[first] && !crossedEdges[second]);
    CHECK(localEdges[first].first != localEdges[second].first &&
          localEdges[first].first != localEdges[second].second &&
          localEdges[first].second != localEdges[second].first &&
          localEdges[first].second != localEdges[second].second);
    crossedEdges[first] = crossedEdges[second] = true;
    const std::array<int, 4> rim = {
        nextVertex, nextVertex + 1, nextVertex + 2, nextVertex + 3};
    const int crossingApex = nextVertex + 4;
    nextVertex += 5;
    crossingWheels.push_back({rim, crossingApex});
    for (int i = 0; i < 4; i++)
      planarized.push_back(makeEdge(rim[i], rim[(i + 1) % 4]));
    for (int vertex : rim)
      planarized.push_back(makeEdge(crossingApex, vertex));
    const std::array<int, 4> endpoints = {
        localEdges[first].first, localEdges[second].first,
        localEdges[first].second, localEdges[second].second};
    for (int i = 0; i < 4; i++)
      planarized.push_back(makeEdge(endpoints[i], rim[i]));
  }

  // Add the uncrossed part of the local drawing and its prescribed boundary cycle.
  for (int i = 0; i < k; i++) {
    const int current = k + (drawing.boundaryOrder[i] - '0');
    const int next = k + (drawing.boundaryOrder[(i + 1) % k] - '0');
    planarized.push_back(makeEdge(current, next));
  }
  for (int edge = 0; edge < 2 * k; edge++) {
    const auto [u, v] = localEdges[edge];
    if (!crossedEdges[edge])
      planarized.push_back(makeEdge(u, v));
  }

  // For each crossing wheel separately, remove its apex and rim edges.  The
  // rest of the augmented drawing (including the other crossing wheels) must
  // remain one connected bridge attached to all four rim vertices.  Hence it
  // lies opposite the apex in every planar embedding.  Deleting the apex
  // leaves an empty face in which the two alternating edge continuations can
  // be crossed.
  for (const CrossingWheel& wheel : crossingWheels) {
    std::vector<EdgeTy> exterior;
    for (const EdgeTy& edge : planarized) {
      const bool firstOnRim =
          std::find(wheel.rim.begin(), wheel.rim.end(), edge.first) !=
          wheel.rim.end();
      const bool secondOnRim =
          std::find(wheel.rim.begin(), wheel.rim.end(), edge.second) !=
          wheel.rim.end();
      if (edge.first == wheel.apex || edge.second == wheel.apex ||
          (firstOnRim && secondOnRim))
        continue;
      exterior.push_back(edge);
    }
    CHECK(connectedSupport(nextVertex, exterior),
          "disconnected exterior around a crossing for boundary %s",
          drawing.boundaryOrder.c_str());
  }

  // The same bridge argument fixes the prescribed outer boundary.  Remove
  // the artificial boundary-cycle edges: the planarized local drawing remains
  // connected and is incident with every boundary point, while the outer apex will
  // form the only bridge on the other side of that cycle.
  std::vector<EdgeTy> interior;
  for (const EdgeTy& edge : planarized) {
    const bool boundaryEdge =
        k <= edge.first && edge.first < 2 * k &&
        k <= edge.second && edge.second < 2 * k;
    if (!boundaryEdge)
      interior.push_back(edge);
  }
  CHECK(connectedSupport(nextVertex, interior),
        "disconnected interior for boundary %s", drawing.boundaryOrder.c_str());
  for (int boundaryVertex = k; boundaryVertex < 2 * k; boundaryVertex++)
    CHECK(std::any_of(interior.begin(), interior.end(),
                      [&](const EdgeTy& edge) {
                        return edge.first == boundaryVertex || edge.second == boundaryVertex;
                      }),
          "unused boundary point %d for boundary %s", boundaryVertex - k,
          drawing.boundaryOrder.c_str());

  const int apex = 2 * k;
  for (int i = 0; i < k; i++) {
    const int current = k + (drawing.boundaryOrder[i] - '0');
    planarized.push_back(makeEdge(apex, current));
  }
  CHECK(isPlanar(nextVertex, planarized, 0),
        "nonplanar local drawing for boundary %s", drawing.boundaryOrder.c_str());
}

int main() {
  // Edge indices 0..k-1 are cycle edges c_i; k..2k-1 are attachment edges s_i.
  // Boundary orders around a regular neighborhood of the path used to
  // replace a 5-cycle.  Only the middle attachment s2 may cross locally.
  const std::vector<LocalDrawing> fiveCycleExpansions = {
      {"01234", {}}, {"01243", {{2, 4}}},
      {"01342", {{4, 7}}}, {"01432", {{1, 4}}}};

  // Boundary orders around the disk of the uncrossed edge introduced by
  // splitting a 4-cycle.
  // The target consecutive cycle edges c0,c1,c2 never cross; the only
  // possible crossing uses c3 and s1, which is uncrossed outside the disk.
  const std::vector<LocalDrawing> adjacentEdgeFourCycleExpansions = {
      {"0123", {}}, {"0132", {{3, 5}}},
      {"0231", {{3, 5}}}, {"0321", {}}};

  // In the 4-cycle expansion, c0 and c2 are created inside the two endpoint
  // disks.  The two nontrivial boundary orders use the crossing c1 x c3.
  const std::vector<LocalDrawing> fourCycleExpansions = {
      {"0123", {}}, {"0132", {{1, 3}}},
      {"0231", {{1, 3}}}, {"0321", {}}};

  // For a 6-cycle, s0 and s3 receive no local crossing.  The attachments
  // s1, s2, s4, and s5 are uncrossed outside the disk and may cross locally.
  const std::vector<LocalDrawing> sixCycleExpansions = {
      {"012345", {}}, {"012354", {{3, 5}}},
      {"012435", {{2, 4}}}, {"012453", {{2, 5}, {3, 11}}},
      {"012534", {{2, 11}}}, {"012543", {{2, 5}}},
      {"013245", {{1, 3}}}, {"013254", {{1, 3}, {5, 10}}},
      {"013425", {{4, 8}}}, {"013452", {{5, 8}}},
      {"013524", {{1, 3}, {5, 10}, {8, 11}}},
      {"013542", {{3, 11}, {5, 8}}}, {"014235", {{1, 10}}},
      {"014253", {{1, 10}, {2, 11}, {3, 5}}},
      {"014325", {{1, 4}}}, {"014352", {{1, 4}, {8, 11}}},
      {"014523", {{1, 3}, {5, 7}, {10, 11}}},
      {"014532", {{1, 5}, {3, 11}}}, {"015234", {{1, 11}}},
      {"015243", {{1, 5}, {4, 8}}},
      {"015324", {{1, 10}, {5, 7}}}, {"015342", {{1, 5}, {2, 4}}},
      {"015423", {{1, 3}, {5, 7}}}, {"015432", {{1, 5}}},
      {"021345", {{0, 2}}}, {"021354", {{0, 2}, {3, 5}}},
      {"021435", {{0, 8}, {2, 4}}},
      {"021453", {{0, 8}, {2, 5}, {3, 11}}},
      {"021534", {{0, 8}, {2, 11}}},
      {"021543", {{0, 8}, {2, 5}}}, {"023145", {{3, 7}}},
      {"023154", {{3, 7}, {5, 10}}}, {"023415", {{4, 7}}},
      {"023514", {{1, 5}, {3, 11}, {7, 10}}},
      {"024135", {{0, 8}, {2, 4}, {7, 10}}},
      {"024153", {{0, 8}, {2, 5}, {3, 11}, {7, 10}}},
      {"024315", {{0, 11}, {4, 8}}},
      {"024513", {{0, 2}, {3, 7}, {5, 8}, {10, 11}}},
      {"025134", {{0, 4}, {2, 7}, {8, 11}}},
      {"025143", {{0, 2}, {4, 7}, {5, 8}}},
      {"025314", {{0, 4}, {3, 7}, {8, 11}}},
      {"025413", {{0, 2}, {3, 7}, {5, 8}}},
      {"031245", {{0, 3}, {2, 7}}},
      {"031254", {{0, 3}, {2, 7}, {5, 10}}},
      {"031425", {{0, 2}, {3, 7}, {4, 8}}},
      {"031524", {{0, 3}, {2, 7}, {5, 10}, {8, 11}}},
      {"032145", {{0, 3}}}, {"032154", {{0, 3}, {5, 10}}},
      {"032415", {{0, 3}, {7, 10}}},
      {"032514", {{0, 3}, {5, 10}, {7, 11}}},
      {"034125", {{0, 4}, {2, 10}, {7, 8}}},
      {"034215", {{0, 4}, {2, 10}}},
      {"035124", {{0, 10}, {2, 11}, {3, 5}, {7, 8}}},
      {"035214", {{0, 10}, {2, 11}, {3, 5}}},
      {"041235", {{0, 10}}}, {"041325", {{0, 4}, {3, 7}}},
      {"042135", {{0, 10}, {2, 7}}},
      {"042315", {{0, 4}, {1, 3}}},
      {"043125", {{0, 4}, {2, 7}}}, {"043215", {{0, 4}}}};

  {
    constexpr int k = 6;
    std::set<std::string> expectedOrders;
    std::string permutation = "012345";
    do {
      std::string reverse = permutation;
      std::reverse(reverse.begin() + 1, reverse.end());
      expectedOrders.insert(std::min(permutation, reverse));
    } while (std::next_permutation(permutation.begin() + 1, permutation.end()));

    std::set<std::string> listedOrders;
    for (const LocalDrawing& drawing : sixCycleExpansions) {
      CHECK(listedOrders.insert(drawing.boundaryOrder).second,
            "duplicate 6-cycle boundary order");
      for (const auto& [first, second] : drawing.crossings) {
        for (int edge : {first, second}) {
          CHECK(edge < k || edge == 7 || edge == 8 || edge == 10 || edge == 11,
                "6-cycle local crossing uses attachment s0 or s3");
        }
      }
      checkLocalDrawing(k, drawing);
    }
    CHECK(listedOrders == expectedOrders, "6-cycle table does not cover all boundary orders");
  }

  const std::set<std::string> fiveCycleOrders = {
      "01234", "01243", "01342", "01432"};
  std::set<std::string> listedFiveCycleOrders;
  for (const LocalDrawing& drawing : fiveCycleExpansions) {
    CHECK(listedFiveCycleOrders.insert(drawing.boundaryOrder).second,
          "duplicate 5-cycle expansion boundary order");
    for (const auto& [first, second] : drawing.crossings)
      for (int edge : {first, second})
        CHECK(edge < 5 || edge == 7, "5-cycle local crossing uses an attachment other than s2");
    checkLocalDrawing(5, drawing);
  }
  CHECK(listedFiveCycleOrders == fiveCycleOrders,
        "5-cycle expansion table does not cover the path rotations");

  auto checkFourCycleExpansions = [](const std::vector<LocalDrawing>& drawings,
                                    const std::set<int>& allowedCrossedEdges) {
    const std::set<std::string> fourCycleBoundaryOrders = {
        "0123", "0132", "0231", "0321"};
    std::set<std::string> listedOrders;
    for (const LocalDrawing& drawing : drawings) {
      CHECK(listedOrders.insert(drawing.boundaryOrder).second,
            "duplicate 4-cycle boundary order");
      for (const auto& [first, second] : drawing.crossings) {
        for (int edge : {first, second})
          CHECK(allowedCrossedEdges.count(edge),
                "4-cycle local crossing uses an edge not allowed by the expansion lemma");
      }
      checkLocalDrawing(4, drawing);
    }
    CHECK(listedOrders == fourCycleBoundaryOrders,
          "4-cycle table does not cover all boundary orders");
  };
  checkFourCycleExpansions(fourCycleExpansions, {1, 3});
  checkFourCycleExpansions(adjacentEdgeFourCycleExpansions, {3, 5});
  return 0;
}
