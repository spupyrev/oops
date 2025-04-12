#pragma once

#include <vector>

std::pair<int, std::vector<std::vector<int>>> breakSymmetries(int verbose, int numVars,
                                                              const std::vector<std::vector<int>> &clauses);
