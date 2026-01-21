#pragma once

#include "common.h"
#include "logging.h"

// #include <unordered_set>
#include "satsuma/tsl.hpp"

struct CrossingPair {
  int32_t e1, e2;
  int32_t p1, p2;
  int64_t key;

  explicit CrossingPair(int m, int x1, int x2, int y1, int y2)
      : e1(x1), e2(x2), p1(y1), p2(y2) {
    if (e1 > e2) { std::swap(e1, e2); }
    if (p1 > p2) { std::swap(p1, p2); }
    // Normalize ordering so that (e1,e2) < (p1,p2) lexicographically
    CHECK(0 <= e1 && e1 < e2 && e2 < m);
    CHECK(0 <= p1 && p1 < p2 && p2 < m);
    if (std::make_pair(e1, e2) >= std::make_pair(p1, p2)) {
      std::swap(e1, p1);
      std::swap(e2, p2);
    }
    CHECK(std::make_pair(e1, e2) < std::make_pair(p1, p2));

    // Construct the key
    key = (int64_t)e1 * (int64_t)m * (int64_t)m * (int64_t)m +
          (int64_t)e2 * (int64_t)m * (int64_t)m +
          (int64_t)p1 * (int64_t)m +
          (int64_t)p2;
  }
};

struct CrossingTriple {
  int32_t a1, a2;
  int32_t b1, b2;
  int32_t c1, c2;
  int64_t key;

  explicit CrossingTriple(int m,
                          int x1, int x2,
                          int y1, int y2,
                          int z1, int z2)
      : a1(x1), a2(x2), b1(y1), b2(y2), c1(z1), c2(z2) {
    // normalize each pair
    CHECK(0 <= a1 && a1 < a2 && a2 < m);
    CHECK(0 <= b1 && b1 < b2 && b2 < m);
    CHECK(0 <= c1 && c1 < c2 && c2 < m);

    // normalize order of the three pairs lexicographically: (a1,a2) < (b1,b2) < (c1,c2)
    auto P = [](int u1, int u2) { return std::make_pair(u1, u2); };
    if (P(a1, a2) > P(b1, b2)) { std::swap(a1, b1); std::swap(a2, b2); }
    if (P(b1, b2) > P(c1, c2)) { std::swap(b1, c1); std::swap(b2, c2); }
    if (P(a1, a2) > P(b1, b2)) { std::swap(a1, b1); std::swap(a2, b2); }
    CHECK(P(a1, a2) < P(b1, b2) && P(b1, b2) < P(c1, c2));

    // pack 6 edge-indices base-m
    key = (((((int64_t)a1 * m + a2) * m + b1) * m + b2) * m + c1) * m + c2;
  }
};

class ForbiddenTuples {
public:
  ForbiddenTuples() = default;

  void insert(const CrossingPair& pair) {
    const auto [it, inserted] = forbiddenCrossingPairs.insert(pair.key);
    CHECK(inserted);
    clauses2.push_back(pair);
  }

  bool contains(const CrossingPair& pair) const {
    return forbiddenCrossingPairs.find(pair.key) != forbiddenCrossingPairs.end();
  }

  void insert(const CrossingTriple& triple) {
    const auto [it, inserted] = forbiddenCrossingTriples.insert(triple.key);
    CHECK(inserted);
    clauses3.push_back(triple);
  }

  bool contains(const CrossingTriple& triple) const {
    return forbiddenCrossingTriples.find(triple.key) != forbiddenCrossingTriples.end();
  }

  const std::vector<CrossingPair>& pairs() const {
    return clauses2;
  }

  const std::vector<CrossingTriple>& triples() const {
    return clauses3;
  }

private:
  // std::unordered_set<int64_t> forbiddenCrossingPairs;
  // std::unordered_set<int64_t> forbiddenCrossingTriples;
  tsl::robin_set<int64_t> forbiddenCrossingPairs;
  tsl::robin_set<int64_t> forbiddenCrossingTriples;

  std::vector<CrossingPair> clauses2;
  std::vector<CrossingTriple> clauses3;
};
