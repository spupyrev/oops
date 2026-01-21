#pragma once

#include <algorithm>
#include <cstring>
#include <iomanip>
#include <memory>
#include <queue>

namespace satsuma {

int sat_to_graph(int l) { return 2 * (abs(l) - 1) + (l < 0); }

int graph_to_sat(int vertex) {
  const bool is_neg = vertex % 2;
  int variable = floor(vertex / 2) + 1;
  return variable * (is_neg ? -1 : 1);
}

int graph_negate(int vertex) { return sat_to_graph(-graph_to_sat(vertex)); }

void terminate_with_error(std::string error_msg) {
  std::cerr << "c \nc " << error_msg << std::endl;
  exit(1);
}

static int hash32shift(int key) {
  key = ~key + (key << 15); // key = (key << 15) - key - 1;
  key = key ^ (key >> 12);
  key = key + (key << 2);
  key = key ^ (key >> 4);
  key = key * 2057; // key = (key + (key << 3)) + (key << 11);
  key = key ^ (key >> 16);
  return key;
}

// Hash function
struct any_hash {
  long operator()(const std::vector<int> &myVector) const {
    long answer = myVector.size();

    for (int i : myVector) {
      answer ^= hash32shift(i) + 0x9e3779b9 + (answer << 6) + (answer >> 2);
    }
    return answer;
  }
};

struct triple_hash {
  inline std::size_t operator()(const std::tuple<int, int, int> &v) const {
    return 48 * std::get<0>(v) + 24 * hash32shift(std::get<1>(v)) + hash32shift(std::get<2>(v));
  }
};

struct pair_hash {
  std::size_t operator()(const std::pair<int, int> &v) const { return v.first * 31 + hash32shift(v.second); }
};
} // namespace satsuma
