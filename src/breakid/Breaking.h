#pragma once

#include "Global.h"

namespace BreakID {

class Permutation;
class Specification;

class Breaker {
private:
  std::unordered_set<std::shared_ptr<Clause>, UVecHash, UvecEqual> clauses;
  std::unordered_set<std::shared_ptr<PBConstraint>, UVecHash, UvecEqual> pbConstrs;
  std::shared_ptr<Specification> originalTheory;
  uint nbExtraVars = 0;
  uint nbBinClauses = 0;
  uint nbRowClauses = 0;
  uint nbRegClauses = 0;

  void addBinary(uint l1, uint l2);
  void addTernary(uint l1, uint l2, uint l3);
  void addQuaternary(uint l1, uint l2, uint l3, uint l4);
  void add(std::shared_ptr<Clause> cl);
  void add(std::shared_ptr<Permutation> perm, std::vector<uint> &order, bool limitExtraConstrs);
  void addPB(std::shared_ptr<Permutation> perm, std::vector<uint> &order, bool limitExtraConstrs);
  void addShatter(std::shared_ptr<Permutation> perm, std::vector<uint> &order, bool limitExtraConstrs);

public:
  Breaker(std::shared_ptr<Specification> origTheo);

  ~Breaker() {};

  // Prints the current breaker
  void print();
  void result(int &newNumVars, std::vector<std::vector<int>> &newClauses);

  void addBinClause(uint l1, uint l2);
  void addRegSym(std::shared_ptr<Permutation> perm, std::vector<uint> &order);
  void addRowSym(std::shared_ptr<Permutation> perm, std::vector<uint> &order);

  uint getAuxiliaryNbVars();
  uint getTotalNbVars();
  uint getAddedNbClauses();
  uint getTotalNbClauses();

  uint getNbBinClauses();
  uint getNbRowClauses();
  uint getNbRegClauses();

  uint getTseitinVar();
};

} // namespace BreakID
