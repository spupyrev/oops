#include "satsuma.h"
#include "cnf.h"
#include "cnf2wl.h"

#include <string>

void create_cnf2wl(cnf2wl& formula, int numVars, const std::vector<std::vector<int>> &clauses) {
  formula.reserve(numVars, clauses.size());
  std::vector<int> construct_clause;
  for (const auto& clause : clauses) {
    construct_clause.clear();
    construct_clause.insert(construct_clause.end(), clause.begin(), clause.end());
    formula.add_clause(construct_clause);
  }
}

std::pair<int, std::vector<std::vector<int>>>
applySatsumaSymmetry(int verbose, int numVars, std::vector<std::vector<int>> &clauses) {
  cnf2wl formula;
  create_cnf2wl(formula, numVars, clauses);
  clauses.clear();

  if (verbose >= 2) {
    std::clog << "c input cnf: #variables " << formula.n_variables() << " #clauses " << formula.n_clauses() 
              << " #arr " << formula.n_len() << "\n";
  }

  // call the main algorithm
  satsuma::preprocessor satsuma_preprocessor;
  satsuma_preprocessor.set_verbose(verbose);
  int newNumVars;
  std::vector<std::vector<int>> newClauses;
  satsuma_preprocessor.preprocess(formula, newNumVars, newClauses);

  return std::make_pair(newNumVars, newClauses);
}
