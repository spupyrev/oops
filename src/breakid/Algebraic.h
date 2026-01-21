#pragma once

#include "Global.h"

namespace BreakID {

class Specification;
class Breaker;

class Permutation : public std::enable_shared_from_this<Permutation> {
private:
  std::unordered_map<uint, uint> perm;
  std::vector<uint> cycleReprs; // smallest lit in each cycle
  uint maxCycleSize;
  size_t hash;

public:
  std::vector<uint> domain;
  std::vector<uint> posDomain;
  std::vector<uint> image;

  void addFromTo(uint from, uint to);
  void addCycle(std::vector<uint> &cyc);

  Permutation();
  Permutation(std::vector<std::pair<uint, uint>> &tuples);
  // Permutation constructed from swapping two rows.
  Permutation(std::vector<uint> &row1, std::vector<uint> &row2);

  ~Permutation() {};

  uint getImage(uint from);
  // return value is true iff the image is different from the original
  bool getImage(std::vector<uint> &orig, std::vector<uint> &img);
  bool getImage(std::map<uint, uint> &orig, std::map<uint, uint> &img);
  void getCycle(uint lit, std::vector<uint> &orb);
  bool isInvolution();
  bool permutes(uint lit);
  uint supportSize();
  bool isIdentity();

  void print(std::ostream &out);

  bool formsMatrixWith(std::shared_ptr<Permutation> other);
  std::pair<std::shared_ptr<Permutation>, std::shared_ptr<Permutation>> getLargest(std::shared_ptr<Permutation> other);
  void getSharedLiterals(std::shared_ptr<Permutation> other, std::vector<uint> &shared);
  std::vector<uint> &getCycleReprs();
  uint getMaxCycleSize();
  uint getNbCycles();

  bool equals(std::shared_ptr<Permutation> other);
};

class Matrix {
private:
  std::vector<std::vector<uint> *> rows; // TODO: refactor this as 1 continuous vector
  std::unordered_map<uint, uint> rowco;
  std::unordered_map<uint, uint> colco;

public:
  Matrix();
  ~Matrix();
  void print(std::ostream &out);

  void add(std::vector<uint> *row);
  uint nbColumns();
  uint nbRows();
  void tryToAddNewRow(std::shared_ptr<Permutation> p, uint rowIndex, Specification *theory);
  std::vector<uint> *getRow(uint rowindex);
  bool permutes(uint x);
  uint getLit(uint row, uint column);

  uint getRowNb(uint x);
  uint getColumnNb(uint x);

  std::shared_ptr<Permutation> testMembership(const std::shared_ptr<Permutation> p);
  std::shared_ptr<Permutation> getProductWithRowsWap(const std::shared_ptr<Permutation> p, uint r1,
                                                     uint r2); // return p*swap(r1,r2)
};

class Group {
private:
  std::vector<std::shared_ptr<Permutation>> permutations;
  std::vector<std::shared_ptr<Matrix>> matrices;
  std::unordered_set<uint> support;

  void cleanPermutations(std::shared_ptr<Matrix> matrix); // remove permutations implied by the matrix

public:
  // NOTE: if a group has a shared pointer to a theory, and a theory a shared pointer to a group, none of the memory
  // pointed to by these pointers will ever be freed :(
  Specification *theory; // non-owning pointer

  Group() {};

  ~Group() {};

  void add(std::shared_ptr<Permutation> p);
  void checkColumnInterchangeability(std::shared_ptr<Matrix> m);

  void print(std::ostream &out);

  std::shared_ptr<Matrix> getInitialMatrix();

  void addMatrices();
  void addMatrix(std::shared_ptr<Matrix> m); // cnf-parameter, otherwise we have to store a pointer to the cnf here :(
  uint getNbMatrices();
  uint getNbRowSwaps();

  std::shared_ptr<Matrix> getMatrix(uint idx);

  void getDisjointGenerators(std::vector<std::shared_ptr<Group>> &subgroups);
  uint getSize();

  bool permutes(uint lit);
  uint getSupportSize();

  void getOrderAndAddBinaryClausesTo(Breaker &brkr,
                                     std::vector<uint> &out_order); // returns a vector containing a lit for literals
                                                                    // relevant to construct sym breaking clauses
  void addBinaryClausesTo(Breaker &brkr, std::vector<uint> &out_order, const std::unordered_set<uint> &excludedLits);
  void addBreakingClausesTo(Breaker &brkr);

  void maximallyExtend(std::shared_ptr<Matrix> matrix, uint indexOfFirstNewRow);
};

} // namespace BreakID
