#include <fstream>
#include <memory>
#include <queue>

#include "Algebraic.h"
#include "Breaking.h"
#include "Global.h"
#include "Theory.h"

namespace BreakID {

Breaker::Breaker(std::shared_ptr<Specification> origTheo) : originalTheory(origTheo) {}

void Breaker::print() {
  // CNF CASE:
  std::cout << "c number of breaking clauses added: " << getAddedNbClauses() << "\n";
  std::cout << "c max original variable: " << nVars << "\n";
  std::cout << "c auxiliary variables: " << getAuxiliaryNbVars() << "\n";
  // original
  std::cout << "p cnf " << getTotalNbVars() << " " << getTotalNbClauses() << "\n";
  // originalTheory->print(std::cout);
  // new
  for (auto c : clauses) {
    c->print(std::cout);
  }
}

void Breaker::result(int &newNumVars, std::vector<std::vector<int>> &newClauses) {
  newNumVars = getTotalNbVars();
  newClauses.reserve(getAddedNbClauses());
  for (auto c : clauses) {
    std::vector<int> vec;
    c->to_vector(vec);
    newClauses.push_back(vec);
  }
}

void Breaker::add(std::shared_ptr<Clause> cl) { clauses.insert(cl); }

void Breaker::addBinary(uint l1, uint l2) {
  std::shared_ptr<Clause> toAdd(new Clause());
  toAdd->lits.push_back(l1);
  toAdd->lits.push_back(l2);
  add(toAdd);
}

void Breaker::addTernary(uint l1, uint l2, uint l3) {
  std::shared_ptr<Clause> toAdd(new Clause());
  toAdd->lits.push_back(l1);
  toAdd->lits.push_back(l2);
  toAdd->lits.push_back(l3);
  add(toAdd);
}

void Breaker::addQuaternary(uint l1, uint l2, uint l3, uint l4) {
  std::shared_ptr<Clause> toAdd(new Clause());
  toAdd->lits.push_back(l1);
  toAdd->lits.push_back(l2);
  toAdd->lits.push_back(l3);
  toAdd->lits.push_back(l4);
  add(toAdd);
}

void Breaker::addBinClause(uint l1, uint l2) {
  ++nbBinClauses;
  addBinary(l1, l2);
}

void Breaker::addRegSym(std::shared_ptr<Permutation> perm, std::vector<uint> &order) {
  uint current = getTotalNbClauses();
  if (useShatterTranslation) {
    addShatter(perm, order, true);
  } else if (nbPBConstraintsToGroup > 0) {
    addPB(perm, order, true);
  } else {
    add(perm, order, true);
  }
  nbRegClauses += getTotalNbClauses() - current;
}

void Breaker::addRowSym(std::shared_ptr<Permutation> perm, std::vector<uint> &order) {
  uint current = getTotalNbClauses();
  if (useShatterTranslation) {
    addShatter(perm, order, false);
  } else if (nbPBConstraintsToGroup > 0) {
    addPB(perm, order, true);
  } else {
    add(perm, order, false);
  }
  nbRowClauses += getTotalNbClauses() - current;
}

std::vector<uint> getVarsToBreakOn(std::shared_ptr<Permutation> perm, std::vector<uint> &order, bool limitExtraConstrs,
                                   int symBreakingFormLength) {
  std::unordered_set<uint> allowedLits; // which are not the last lit in their cycle, unless they map to their negation
  for (uint i = order.size(); i > 0; --i) {
    uint lit = order.at(i - 1);
    if (allowedLits.count(lit) == 0) { // we have a last lit of a cycle
      uint sym = perm->getImage(lit);
      while (sym != lit) { // add the other lits of the cycle and the negated cycle
        allowedLits.insert(sym);
        allowedLits.insert(neg(sym));
        sym = perm->getImage(sym);
      }
    }
  }

  std::vector<uint> result;
  result.reserve(order.size());
  for (auto l : order) {
    uint sym = perm->getImage(l);
    if (l != sym && allowedLits.count(l) > 0) {
      result.push_back(l);
    }
    if (limitExtraConstrs && (int)result.size() >= symBreakingFormLength) {
      break;
    }
    if (sym == neg(l)) {
      break;
    }
  }
  return result;
}

void Breaker::add(std::shared_ptr<Permutation> perm, std::vector<uint> &order, bool limitExtraConstrs) {
  std::vector<uint> varsToBreakOn = getVarsToBreakOn(perm, order, limitExtraConstrs, symBreakingFormLength);

  int nrExtraConstrs = 0;
  uint prevLit = 0;
  uint prevSym = 0;
  uint prevTst = 0; // previous tseitin
  for (auto l : varsToBreakOn) {
    uint sym = perm->getImage(l);
    uint tst = 0;
    if (nrExtraConstrs == 0) {
      // adding clause for l => sym :
      // ~l | sym
      addBinary(neg(l), sym);
    } else if (nrExtraConstrs == 1) {
      // adding clauses for (prevSym => prevLit) => tst and tst => (l => sym)
      tst = getTseitinVar();
      // prevSym | tst
      addBinary(prevSym, tst);
      // ~prevLit | tst
      addBinary(neg(prevLit), tst);
      // ~tst | ~l | sym
      addTernary(neg(tst), neg(l), sym);
      if (useFullTranslation) {
        // adding clauses for tst => (prevSym => prevLit)
        // ~tst | ~prevSym | prevLit
        addTernary(neg(tst), neg(prevSym), prevLit);
      }
    } else {
      // adding clauses for (prevSym => prevLit) & prevTst => tst and tst => (l => sym)
      tst = getTseitinVar();
      // prevSym | ~prevTst | tst
      addTernary(prevSym, neg(prevTst), tst);
      // ~prevLit | ~prevTst | tst
      addTernary(neg(prevLit), neg(prevTst), tst);
      // ~tst | ~l | sym
      addTernary(neg(tst), neg(l), sym);
      if (useFullTranslation) {
        // adding clauses for tst => prevTst and tst => (prevSym => prevLit)
        // ~tst | prevTst
        addBinary(neg(tst), prevTst);
        // ~tst | ~prevSym | prevLit
        addTernary(neg(tst), neg(prevSym), prevLit);
      }
    }
    ++nrExtraConstrs;

    prevLit = l;
    prevSym = sym;
    prevTst = tst;
  }
}

void Breaker::addShatter(std::shared_ptr<Permutation> perm, std::vector<uint> &order, bool limitExtraConstrs) {
  std::vector<uint> varsToBreakOn = getVarsToBreakOn(perm, order, limitExtraConstrs, symBreakingFormLength);

  int nrExtraConstrs = 0;
  uint prevLit = 0;
  uint prevSym = 0;
  uint prevTst = 0; // previous tseitin
  for (auto l : varsToBreakOn) {
    uint sym = perm->getImage(l);
    uint tst = 0;
    if (nrExtraConstrs == 0) {
      // adding clause for l => sym :
      // ~l | sym
      addBinary(neg(l), sym);
    } else if (nrExtraConstrs == 1) {
      tst = getTseitinVar();
      // clause(-z, -x, p[x], 0);
      addTernary(neg(prevLit), neg(l), sym);
      // clause(-z, vars+1, 0);
      addBinary(neg(prevLit), tst);
      // clause(p[z], -x, p[x], 0);
      addTernary(prevSym, neg(l), sym);
      // clause(p[z], vars+1, 0);
      addBinary(prevSym, tst);
    } else {
      tst = getTseitinVar();
      // clause(-vars, -z, -x, p[x], 0);
      addQuaternary(neg(prevTst), neg(prevLit), neg(l), sym);
      // clause(-vars, -z, vars+1, 0);
      addTernary(neg(prevTst), neg(prevLit), tst);
      // clause(-vars, p[z], -x, p[x], 0);
      addQuaternary(neg(prevTst), prevSym, neg(l), sym);
      // clause(-vars, p[z], vars+1, 0);
      addTernary(neg(prevTst), prevSym, tst);
    }
    ++nrExtraConstrs;

    prevLit = l;
    prevSym = sym;
    prevTst = tst;
  }
}

/*
Assume some literal-mapping symmetry s, and 3 ordered variables x<y<z. Bart's initial PB sbf is:

-4 x 4 s(x) -2 y 2 s(y) -1 z 1 s(z) >= 0

(though it's partial as there might be variables beyond x,y,z).

We can split this constraint by introducing auxiliary variables o, p, q, r. Then we get

1 o >= 1
-2 o -1 x 1 s(x) 1 p >= -1
-2 p -1 y 1 s(y) 1 q >= -1
-2 q -1 z 1 s(z) 1 r >= -1
-1 r >= -1 (trivial)

Eliminating o, p, q, r by canceling addition implies the above sbf. And eliminating any single auxiliary variable
multiplies some set of coefficients by 2.

Note that the core constraint

-2 o -1 x 1 s(x) 1 p >= -1

is logically equivalent to the sbf encoding with 3 3-clauses used in the CNF context.
 */

std::shared_ptr<PBConstraint> getPBConstraint(const std::vector<uint> &lits, uint firstTseitin, uint lastTseitin) {
  int weight = 0;
  uint coeff = 1 << (lits.size() / 2);
  if (firstTseitin != 0) {
    weight -= coeff;
  }
  if (lastTseitin != 0) {
    weight += 1;
  }
  std::shared_ptr<PBConstraint> pbc = std::make_shared<PBConstraint>(weight);
  if (firstTseitin != 0) {
    pbc->addTerm(firstTseitin, -coeff);
  }
  if (lastTseitin != 0) {
    pbc->addTerm(lastTseitin, 1);
  }
  for (uint i = 0; i < lits.size(); i += 2) {
    coeff = coeff >> 1; // divide coeff by 2
    pbc->addTerm(lits[i], -coeff);
    pbc->addTerm(lits[i + 1], coeff);
  }
  return pbc;
}

void Breaker::addPB(std::shared_ptr<Permutation> perm, std::vector<uint> &order, bool limitExtraConstrs) {
  std::vector<uint> varsToBreakOn = getVarsToBreakOn(perm, order, limitExtraConstrs, symBreakingFormLength);

  // Create sets of lits that make up each constraint
  std::vector<std::vector<uint>> constrLits; // lits for each pb constraint
  for (auto l : varsToBreakOn) {
    if (constrLits.size() == 0 || (int)constrLits.back().size() >= 2 * nbPBConstraintsToGroup) {
      constrLits.push_back({});
    }
    constrLits.back().push_back(l);
    constrLits.back().push_back(perm->getImage(l));
  }

  // For each set of lits, create the corresponding PB constraint
  // We optimize the number of tseitin variables
  if (constrLits.size() == 0) {
    return;
  }
  if (constrLits.size() == 1) {
    pbConstrs.insert(getPBConstraint(constrLits.front(), 0, 0));
    return;
  }

  // else, we have at least 2 constraints, of which the first and last form a special tseitin case
  uint firstTseitin = 0;
  uint lastTseitin = getTseitinVar();
  pbConstrs.insert(getPBConstraint(constrLits.front(), firstTseitin, lastTseitin));

  for (int i = 1; i < (int)constrLits.size() - 1; ++i) {
    firstTseitin = lastTseitin;
    lastTseitin = getTseitinVar();
    pbConstrs.insert(getPBConstraint(constrLits[i], firstTseitin, lastTseitin));
  }

  firstTseitin = lastTseitin;
  lastTseitin = 0;
  pbConstrs.insert(getPBConstraint(constrLits.back(), firstTseitin, lastTseitin));
}

uint Breaker::getAuxiliaryNbVars() { return nbExtraVars; }

uint Breaker::getTotalNbVars() { return nVars + nbExtraVars; }

uint Breaker::getAddedNbClauses() { return clauses.size() + pbConstrs.size(); }

uint Breaker::getTotalNbClauses() { return originalTheory->getSize() + getAddedNbClauses(); }

uint Breaker::getNbBinClauses() { return nbBinClauses; }

uint Breaker::getNbRowClauses() { return nbRowClauses; }

uint Breaker::getNbRegClauses() { return nbRegClauses; }

uint Breaker::getTseitinVar() {
  ++nbExtraVars;
  return encode(getTotalNbVars());
}
} // namespace BreakID
