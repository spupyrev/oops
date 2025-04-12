#pragma once

#include "Global.h"

namespace BreakID {

class Graph;
class Permutation;
class Group;
class Matrix;
class Breaker;

class Specification {
  friend class Breaker;

protected:
  std::shared_ptr<Graph> graph;
  std::shared_ptr<Group> group;
  std::vector<std::string> originalSpec;

public:
  Specification();
  virtual ~Specification();

  virtual void print(std::ostream &out) = 0;
  virtual uint getSize() = 0;

  std::shared_ptr<Graph> getGraph();
  std::shared_ptr<Group> getGroup();

  virtual void cleanUp();

  virtual void setSubTheory(std::shared_ptr<Group> subgroup) = 0;

  virtual bool isSymmetry(Permutation &prm) = 0;
};

class CNF : public Specification {
  friend class Breaker;

private:
  std::unordered_set<std::shared_ptr<Clause>, UVecHash, UvecEqual>
      clauses; // must be an unordered_set, since we need to be able to test whether a clause exists to detect
               // symmetries

  void readCNF(int numVars_, const std::vector<std::vector<int>> &clauses_);

public:
  CNF(int numVars_, const std::vector<std::vector<int>> &clauses_);
  CNF(std::vector<std::shared_ptr<Clause>> &clss, std::shared_ptr<Group> grp);
  ~CNF();

  void print(std::ostream &out);
  uint getSize();

  void setSubTheory(std::shared_ptr<Group> subgroup);

  bool isSymmetry(Permutation &prm);
};

class LogicProgram : public Specification {
  friend class Breaker;

private:
  std::unordered_set<std::shared_ptr<Rule>, UVecHash, UvecEqual>
      rules; // must be an unordered_set, since we need to be able to test whether a rule exists to detect symmetries
  void readLogicProgram(std::istream &input);

public:
  LogicProgram(std::istream &input);
  LogicProgram(std::vector<std::shared_ptr<Rule>> &rls, std::shared_ptr<Group> grp);
  ~LogicProgram();

  void print(std::ostream &out);

  uint getSize();

  void setSubTheory(std::shared_ptr<Group> subgroup);

  bool isSymmetry(Permutation &prm);
};

class PB : public Specification {
  friend class Breaker;

private:
  std::unordered_set<std::shared_ptr<PBConstraint>, UVecHash, UvecEqual>
      constraints; // must be an unordered_set, since we need to be able to test whether a rule exists to detect
                   // symmetries
  void readPB(std::istream &input);
  // Takes as input a vector of variables (not lits) and a vector of (possibly neg) coefficients and a weight.
  // Transforms into a pbconstraint and adds.
  void processInequality(std::vector<int> lits, std::vector<int> coefs, int w);

public:
  PB(std::istream &input);
  PB(std::vector<std::shared_ptr<PBConstraint>> &rls, std::shared_ptr<Group> grp);

  void print(std::ostream &out);

  uint getSize();

  void setSubTheory(std::shared_ptr<Group> subgroup);

  bool isSymmetry(Permutation &prm);
};

} // namespace BreakID
