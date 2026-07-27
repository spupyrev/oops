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

  // Boundary orders around a vertex obtained by contracting a 5-cycle.  The
  // first table may cross only the two consecutive attachments s0,s1.
  const std::vector<LocalDrawing> fiveCycleContractions = {
      {"01234", {}}, {"01243", {{2, 4}}}, {"01324", {{1, 3}}},
      {"01342", {{1, 5}, {2, 4}}}, {"01423", {{1, 3}, {4, 6}}},
      {"01432", {{1, 4}}}, {"02134", {{0, 2}}}, {"02143", {{2, 5}}},
      {"02314", {{3, 6}}}, {"02413", {{1, 4}, {2, 5}, {3, 6}}},
      {"03124", {{0, 2}, {3, 5}}}, {"03214", {{0, 3}}}};

  // These drawings preserve the star {c4,c0,s0}; all four other attachments
  // may receive one local crossing after all five were uncrossed outside.
  const std::vector<LocalDrawing> fiveCycleStarContractions = {
      {"01234", {}}, {"01243", {{2, 9}}}, {"01324", {{1, 3}}},
      {"01342", {{1, 3}, {7, 9}}}, {"01423", {{1, 9}}},
      {"01432", {{1, 9}, {3, 7}}}, {"02134", {{2, 6}}},
      {"02143", {{2, 6}, {8, 9}}}, {"02314", {{3, 6}}},
      {"02413", {{1, 3}, {6, 8}, {7, 9}}},
      {"03124", {{1, 3}, {6, 8}}}, {"03214", {{1, 8}, {3, 6}}}};

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

  auto checkFiveCycleContractions = [](const std::vector<LocalDrawing>& drawings,
                                       const std::set<int>& allowedCrossedAttachments,
                                       const std::set<int>& uncrossedTarget) {
    std::set<std::string> expectedOrders;
    std::string permutation = "01234";
    do {
      std::string reverse = permutation;
      std::reverse(reverse.begin() + 1, reverse.end());
      expectedOrders.insert(std::min(permutation, reverse));
    } while (std::next_permutation(permutation.begin() + 1, permutation.end()));

    std::set<std::string> listedOrders;
    for (const LocalDrawing& drawing : drawings) {
      CHECK(listedOrders.insert(drawing.boundaryOrder).second,
            "duplicate contracted 5-cycle boundary order");
      for (const auto& [first, second] : drawing.crossings) {
        for (int edge : {first, second}) {
          CHECK(edge < 5 || allowedCrossedAttachments.count(edge),
                "contracted 5-cycle drawing crosses an unavailable attachment");
          CHECK(!uncrossedTarget.count(edge),
                "contracted 5-cycle drawing crosses a target edge");
        }
      }
      checkLocalDrawing(5, drawing);
    }
    CHECK(listedOrders == expectedOrders,
          "contracted 5-cycle table does not cover every boundary order");
  };
  checkFiveCycleContractions(fiveCycleContractions, {5, 6}, {});
  checkFiveCycleContractions(
      fiveCycleStarContractions, {6, 7, 8, 9}, {0, 4, 5});

  // Check the weaker rotation-independent condition used for ordinary
  // 1-flexibility.  Generate every dihedral relabeling of the ordinary
  // table.  For any three available attachments, some relabeling must avoid
  // every unavailable attachment.  If an attachment is prescribed, it and
  // any three of the other four are available, and the prescribed attachment
  // must not cross locally.
  auto normalizeOrder = [](std::string order) {
    std::rotate(
        order.begin(), std::find(order.begin(), order.end(), '0'),
        order.end());
    std::string reverse = order;
    std::reverse(reverse.begin() + 1, reverse.end());
    return std::min(order, reverse);
  };
  auto mapIndex = [](const int index, const int shift, const bool reflect) {
    int result = shift + (reflect ? -index : index);
    result %= 5;
    return result < 0 ? result + 5 : result;
  };
  std::vector<LocalDrawing> relabeledFiveCycleContractions;
  for (const LocalDrawing& drawing : fiveCycleContractions) {
    for (int shift = 0; shift < 5; shift++) {
      for (bool reflect : {false, true}) {
        LocalDrawing relabeled;
        for (char label : drawing.boundaryOrder)
          relabeled.boundaryOrder.push_back(
              static_cast<char>('0' + mapIndex(label - '0', shift, reflect)));
        relabeled.boundaryOrder = normalizeOrder(relabeled.boundaryOrder);
        auto mapEdge = [&](const int edge) {
          if (edge >= 5)
            return 5 + mapIndex(edge - 5, shift, reflect);
          return reflect ? mapIndex(edge + 1, shift, true)
                         : mapIndex(edge, shift, false);
        };
        for (const auto& [first, second] : drawing.crossings)
          relabeled.crossings.emplace_back(
              mapEdge(first), mapEdge(second));
        relabeledFiveCycleContractions.push_back(std::move(relabeled));
      }
    }
  }

  auto hasUsableExpansion = [&](
                                const std::string& order,
                                const int availableAttachments,
                                const int preservedEdge) {
    return std::any_of(
        relabeledFiveCycleContractions.begin(),
        relabeledFiveCycleContractions.end(),
        [&](const LocalDrawing& drawing) {
          if (drawing.boundaryOrder != order)
            return false;
          for (const auto& [first, second] : drawing.crossings) {
            for (int edge : {first, second}) {
              if (edge == preservedEdge)
                return false;
              if (edge >= 5 &&
                  ((availableAttachments >> (edge - 5)) & 1) == 0)
                return false;
            }
          }
          return true;
        });
  };
  std::set<std::string> fiveCycleBoundaryOrders;
  for (const LocalDrawing& drawing : fiveCycleContractions)
    fiveCycleBoundaryOrders.insert(drawing.boundaryOrder);
  for (const std::string& order : fiveCycleBoundaryOrders) {
    for (int available = 0; available < (1 << 5); available++) {
      if (__builtin_popcount(available) != 3)
        continue;
      CHECK(hasUsableExpansion(order, available, -1),
            "three attachments do not cover boundary %s", order.c_str());
    }
    for (int attachment = 0; attachment < 5; attachment++) {
      for (int omitted = 0; omitted < 5; omitted++) {
        if (omitted == attachment)
          continue;
        const int available = ((1 << 5) - 1) ^ (1 << omitted);
        CHECK(hasUsableExpansion(order, available, 5 + attachment),
              "attachment %d is not preserved for boundary %s",
              attachment, order.c_str());
      }
    }
  }

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
