#pragma once

#include <algorithm>
#include <chrono>
#include <numeric>
#include <string>
#include <vector>
#include <iostream>
#include <sstream>

constexpr double EPS = 1e-8;
constexpr double PI = 3.14159265358979323846;

template<class T>
T Abs(const T& t) {
  if (t > 0)
    return t;
  return -t;
}

template<class T>
T Sgn(const T& t) {
  if (t > 0)
    return 1;
  if (t < 0)
    return -1;
  return 0;
}

template<class T>
T Sqr2(const T& t) {
  return ((t) * (t));
}

template <typename T>
std::string to_string(const T& n) {
  std::ostringstream ss;
  ss << n;
  return ss.str();
}

template <typename T>
std::string to_string(const std::vector<T>& vec, const std::string& separator) {
  std::string desc = "";

  for (auto p : vec) {
    if (desc.length() > 0)
      desc += separator;

    desc += to_string(p);
  }

  return desc;
}

template <typename T>
std::string to_string(const std::vector<T>& vec) {
  return to_string(vec, " ");
}

inline int to_int(const std::string& s) {
  int n;
  std::istringstream(s) >> n;
  return n;
}

inline double to_double(const std::string& s) {
  double n;
  std::istringstream(s) >> n;
  return n;
}

inline bool ends_with(const std::string& s, const std::string& postfix) {
  if (postfix.size() > s.size())
    return false;
  return std::equal(postfix.rbegin(), postfix.rend(), s.rbegin());
}

inline bool starts_with(const std::string& s, const std::string& prefix) {
  if (prefix.size() > s.size())
    return false;
  return std::equal(prefix.begin(), prefix.end(), s.begin());
}

std::string ms_to_str(uint64_t duration_ms);
std::string ms_to_str(std::chrono::time_point<std::chrono::steady_clock> startTime, 
                      std::chrono::time_point<std::chrono::steady_clock> endTime);

std::vector<std::string> SplitNotNull(const std::string& s, const std::string& c);
std::vector<int> SplitNotNullInt(const std::string& s, const std::string& c);

struct Rand {
  static size_t setSeed();
  static size_t setSeed(size_t seed);

  static double nextDouble();
  static bool nextBool();
  static bool check(double probability);
  static int next();
  static int next(int upperExclusive);
  static int next(int lowerInclusive, int upperExclusive);
  template<typename RandomIt>
  static void shuffle(RandomIt first, RandomIt last);
  static std::vector<int> permutation(size_t n);
};

template<typename RandomIt>
void Rand::shuffle(RandomIt first, RandomIt last) {
  typename std::iterator_traits<RandomIt>::difference_type i, n;
  n = last - first;
  for (i = n-1; i > 0; --i) {
    std::swap(first[i], first[Rand::next(i + 1)]);
  }
}

std::vector<int> identity(size_t n);
std::vector<int> identity(size_t start, size_t n);
std::vector<int> reverse(const std::vector<int>& p);
std::vector<int> inverse(const std::vector<int>& p);

template <typename T>
bool is_sorted(const std::vector<T>& vec) {
  return std::is_sorted(vec.begin(), vec.end());
}

template <typename T>
void remove_value(std::vector<T>& vec, const T& element) {
  vec.erase(std::remove(vec.begin(), vec.end(), element), vec.end());
}

inline bool contains(const std::vector<int>& vec, int element) {
  return std::find(vec.begin(), vec.end(), element) != vec.end();
}

inline bool contains(const std::vector<std::string>& vec, const char* element) {
  return std::find(vec.begin(), vec.end(), std::string(element)) != vec.end();
}

template <typename T>
void sort_unique(std::vector<T>& vec) {
  std::sort(vec.begin(), vec.end());
  vec.erase(std::unique(vec.begin(), vec.end()), vec.end());  
}

template <typename T>
bool equal_unsorted(const std::vector<T>& vec1, const std::vector<T>& vec2) {
  for (size_t i = 0; i < vec1.size(); i++) {
    if (std::find(vec2.begin(), vec2.end(), vec1[i]) == vec2.end())
      return false;
  }
  for (size_t i = 0; i < vec2.size(); i++) {
    if (std::find(vec1.begin(), vec1.end(), vec2[i]) == vec1.end())
      return false;
  }
  return true;
}

template <typename T>
double average(const std::vector<T>& vec) {
  if (vec.empty())
    return 0;
  double sum = std::accumulate(vec.begin(), vec.end(), 0.0);
  return sum / double(vec.size());
}

template <typename T>
double median(const std::vector<T>& vec) {
  const size_t size = vec.size();

  if (size == 0)
    return 0.0;
  else if (size == 1)
    return double(vec[0]);
  else if (size == 2)
    return double(vec[0] + vec[1]) / 2.0;

  std::vector<T> tmp(vec.begin(), vec.end());
  std::sort(tmp.begin(), tmp.end());

  if (size % 2 == 0)
    return double(tmp[size / 2 - 1] + tmp[size / 2]) / 2.0;
  else
    return double(tmp[size / 2]);
}

template <typename T>
T Max(const std::vector<T>& v) {
  return *std::max_element(v.begin(), v.end());
}

template <typename T>
T Min(const std::vector<T>& v) {
  return *std::min_element(v.begin(), v.end());
}

double confidence_interval(const std::vector<uint64_t>& vec);
double percentile(const std::vector<double>& vec, int value);

bool cross(const int l1, const int r1, const int l2, const int r2);
bool nest(const int l1, const int r1, const int l2, const int r2);

int Compare(const double numberA, const double numberB);
bool Equal(const double a, const double b, const double eps=EPS);
bool GreaterOrEqual(double numberA, double numberB);
bool Greater(double numberA, double numberB);
bool LessOrEqual(double numberA, double numberB);
bool Less(double numberA, double numberB);
