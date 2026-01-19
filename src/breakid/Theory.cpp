#include <algorithm>
#include <bitset>
#include <iostream>
#include <iterator>
#include <sstream>

#include "Algebraic.h"
#include "Breaking.h"
#include "Graph.h"
#include "Theory.h"

//=========================CNF==================================================

namespace BreakID {

using namespace std;

void CNF::readCNF(int numVars_, const std::vector<std::vector<int>> &clauses_) {
  nVars = numVars_;
  clauses.reserve(clauses_.size());

  std::set<uint> inclause = std::set<uint>();
  for (auto clause_ : clauses_) {
    for (int l : clause_) {
      assert(l != 0);
      inclause.insert(encode(l));
    }
    if (inclause.size() == 0) {
      gracefulError("Theory can not contain empty clause.");
    }
    bool isTautology = false;
    for (auto lit : inclause) {
      if (inclause.count(neg(lit)) > 0) {
        isTautology = true;
        break;
      }
    }
    if (not isTautology) {
      std::shared_ptr<Clause> cl(new Clause(inclause));
      clauses.insert(cl);
    }
    inclause.clear();
  }
}

CNF::CNF(int numVars_, const std::vector<std::vector<int>> &clauses_) {
  readCNF(numVars_, clauses_);

  if (verbosity >= 2) {
    std::clog << "*** Creating first graph..." << std::endl;
  }
  graph = make_shared<Graph>(clauses);
  if (verbosity > 2) {
    std::clog << "**** Number of nodes: " << graph->getNbNodes() << std::endl;
    std::clog << "**** Number of edges: " << graph->getNbEdges() << std::endl;
  }

  group = make_shared<Group>();
  if (verbosity >= 2) {
    std::clog << "*** Detecting symmetry group..." << std::endl;
  }
  std::vector<std::shared_ptr<Permutation>> symgens;
  graph->getSymmetryGenerators(symgens);
  for (auto symgen : symgens) {
    group->add(symgen);
  }
}

CNF::CNF(std::vector<std::shared_ptr<Clause>> &clss, std::shared_ptr<Group> grp) {
  clauses.insert(clss.cbegin(), clss.cend());
  graph = make_shared<Graph>(clauses);
  group = grp;
  for (uint l = 0; l < 2 * nVars; ++l) {
    if (not grp->permutes(l)) {
      graph->setUniqueColor(l);
    }
  }
  for (uint m = 0; m < grp->getNbMatrices(); ++m) {
    auto mat = grp->getMatrix(m);
    for (uint r = 0; r < mat->nbRows() - 1; ++r) {
      getGraph()->setUniqueColor(*(mat->getRow(r)));
    }
  }
}

CNF::~CNF() {}

uint CNF::getSize() { return clauses.size(); }

void CNF::setSubTheory(std::shared_ptr<Group> subgroup) {
  // TODO: what is this method supposed to do: keep all clauses that are not mapped to themselves? Is it simply made
  // approximative on purpose or by accident?
  std::vector<std::shared_ptr<Clause>> subclauses;
  for (auto cl : clauses) {
    for (auto lit : cl->lits) {
      if (subgroup->permutes(lit)) {
        subclauses.push_back(cl);
        break;
      }
    }
  }
  subgroup->theory = new CNF(subclauses, subgroup);
}

bool CNF::isSymmetry(Permutation &prm) {
  for (auto cl : clauses) {
    std::shared_ptr<Clause> symmetrical(new Clause());
    if (!prm.getImage(cl->lits, symmetrical->lits)) {
      continue;
    }
    std::sort(symmetrical->lits.begin(), symmetrical->lits.end());
    if (clauses.count(symmetrical) == 0) {
      return false;
    }
  }
  return true;
}

/******************
 * SPECIFICATION
 *
 */

Specification::~Specification() {}

Specification::Specification() {}

std::shared_ptr<Graph> Specification::getGraph() { return graph; }

std::shared_ptr<Group> Specification::getGroup() { return group; }

void Specification::cleanUp() {
  graph.reset();
  group.reset();
}

/******************
 * LOGIC PROGRAM
 */

void checkVarExists(int lit) {
  if ((uint)abs(lit) > nVars) {
    nVars = abs(lit);
  }
}
void LogicProgram::readLogicProgram(std::istream &input) {
  string line;
  vector<uint> headlits;
  vector<uint> bodylits;
  vector<int> weights;

  // PARSE RULES:
  if (verbosity > 5) {
    std::clog << "***** Reading RULES: " << std::endl;
  }
  while (getline(input, line)) {
    originalSpec.push_back(line);
    if (line.size() == 0 || line.front() == '%') {
      continue;
    }
    if (line.front() == '0') {
      break;
    }
    istringstream iss(line);

    int bound = 0;
    int rule_type;
    iss >> rule_type;
    if (rule_type != 1 && rule_type != 2 && rule_type != 3 && rule_type != 5 && rule_type != 6 && rule_type != 8) {
      std::cerr << "UNSUPPORTED RULE: " << line << std::endl;
      gracefulError("Unsupported rule type: currently only supporting basic(1), constraint(2), choice(3), weight(5), "
                    "minimize(6) and disjunctive(8) rules. Use lp2normal to use this tool with a richer language");
      // TODO: Optimisation rules: fixed literals.
    }

    // Parse head of rule
    if (rule_type == 1 || rule_type == 2 || rule_type == 5) {
      int head;
      iss >> head;
      checkVarExists(head);
      headlits.push_back(encode(head));
    } else if (rule_type == 3 || rule_type == 8) {
      int nbheads;
      iss >> nbheads;
      int lit;
      for (int i = 0; i < nbheads; i++) {
        iss >> lit;
        checkVarExists(lit);
        headlits.push_back(encode(lit));
      }
    } else if (rule_type == 6) {
      int lit;
      iss >> lit;
      assert(lit == 0);
    }

    if (rule_type == 5) {
      iss >> bound;
    }

    // Parse body of rule
    int nblits;
    int nbneglits;
    iss >> nblits;
    iss >> nbneglits;

    if (rule_type == 2) {
      iss >> bound;
    }

    int i = 0;
    int lit;
    for (; i < nbneglits; i++) {
      iss >> lit;
      checkVarExists(lit);
      bodylits.push_back(encode(-lit));
    }
    for (; i < nblits; i++) {
      iss >> lit;
      checkVarExists(lit);
      bodylits.push_back(encode(lit));
    }

    // parse weights
    if (rule_type == 5 || rule_type == 6) {
      int weight;
      for (i = 0; i < nblits; i++) {
        iss >> weight;
        weights.push_back(weight);
      }
    }

    std::shared_ptr<Rule> r(new Rule(rule_type, headlits, bodylits, bound, weights));
    rules.insert(r);
    headlits.clear();
    bodylits.clear();
    weights.clear();
  }

  // PARSE SYMBOL TABLE:
  if (verbosity > 5) {
    std::clog << "***** Reading SYMBOL TABLE: " << std::endl;
  }
  while (getline(input, line)) {
    originalSpec.push_back(line);
    if (line.size() == 0 || line.front() == '%') {
      continue;
    }
    if (line.front() == '0') {
      break;
    }
    // Currently not storing this information
  }

  // PARSE CONSTRAINTS:
  if (verbosity > 5) {
    std::clog << "***** Reading CONSTRAINTS: " << std::endl;
  }
  bool negative = false;
  while (getline(input, line)) {
    originalSpec.push_back(line);
    if (line.size() == 0 || line.front() == '%') {
      continue;
    }
    if (line.front() == '0') {
      if (negative) {
        break;
      } else {
        continue;
      }
    }
    if (line.front() == 'B') {
      if (line.at(1) == '+') {
        negative = false;
        continue;
      }
      if (line.at(1) == '-') {
        negative = true;
        continue;
      }
    }
    int lit;
    istringstream iss(line);
    iss >> lit;
    if (negative) {
      // Note the negation caused by putting this in the body of a rule: negative constraint is positive body.
      bodylits.push_back(encode(lit));
    } else {
      bodylits.push_back(encode(-lit));
    }

    // Constraints are basic rules with empty head
    std::shared_ptr<Rule> r(new Rule(1, headlits, bodylits, 0, weights));
    rules.insert(r);
    headlits.clear();
    bodylits.clear();
  }
}

LogicProgram::LogicProgram(std::istream &input) {
  readLogicProgram(input);
  if (verbosity > 6) {
    std::clog << "Parsed logic program:\n";
    // print(std::clog);
  }
  if (verbosity >= 2) {
    std::clog << "*** Creating first graph..." << std::endl;
  }
  graph = make_shared<Graph>(rules);
  if (verbosity > 2) {
    std::clog << "**** Number of nodes: " << graph->getNbNodes() << std::endl;
    std::clog << "**** Number of edges: " << graph->getNbEdges() << std::endl;
  }

  group = make_shared<Group>();
  if (verbosity >= 2) {
    std::clog << "*** Detecting symmetry group..." << std::endl;
  }
  std::vector<std::shared_ptr<Permutation>> symgens;
  graph->getSymmetryGenerators(symgens);
  for (auto symgen : symgens) {
    group->add(symgen);
  }
}

LogicProgram::LogicProgram(std::vector<std::shared_ptr<Rule>> &rls, std::shared_ptr<Group> grp) {
  rules.insert(rls.cbegin(), rls.cend());
  graph = make_shared<Graph>(rules);
  group = grp;
  for (uint l = 0; l < 2 * nVars; ++l) {
    if (not grp->permutes(l)) {
      graph->setUniqueColor(l);
    }
  }
  for (uint m = 0; m < grp->getNbMatrices(); ++m) {
    auto mat = grp->getMatrix(m);
    for (uint r = 0; r < mat->nbRows() - 1; ++r) {
      getGraph()->setUniqueColor(*(mat->getRow(r)));
    }
  }
}

LogicProgram::~LogicProgram() {}

uint LogicProgram::getSize() { return rules.size(); }

void LogicProgram::setSubTheory(std::shared_ptr<Group> subgroup) {
  // TODO: what is this method supposed to do: keep all clauses that are not mapped to themselves? Is it simply made
  // approximative on purpose or by accident?
  std::vector<std::shared_ptr<Rule>> subrules;
  for (auto r : rules) {
    bool found = false;
    for (auto lit : r->headLits) {
      if (found) {
        break;
      }
      if (subgroup->permutes(lit)) {
        subrules.push_back(r);
        found = true;
      }
    }
    for (auto lit : r->bodyLits) {
      if (found) {
        break;
      }
      if (subgroup->permutes(lit)) {
        subrules.push_back(r);
        found = true;
      }
    }
  }
  subgroup->theory = new LogicProgram(subrules, subgroup);
}

bool LogicProgram::isSymmetry(Permutation &prm) {
  for (auto r : rules) {
    std::shared_ptr<Rule> symmetrical(new Rule());
    if (!(prm.getImage(r->headLits, symmetrical->headLits) || prm.getImage(r->bodyLits, symmetrical->bodyLits))) {
      continue;
    }
    symmetrical->weights = r->weights;
    symmetrical->ruleType = r->ruleType;
    symmetrical->bound = r->bound;
    std::sort(symmetrical->headLits.begin(), symmetrical->headLits.end());
    if (r->ruleType != 5 && r->ruleType != 6) {
      std::sort(symmetrical->bodyLits.begin(), symmetrical->bodyLits.end());
    }
    if (rules.count(symmetrical) == 0) {
      return false;
    }
  }
  return true;
};

/********************PB*/

void PB::processInequality(vector<int> lits, vector<int> coefs, int w) {
  auto constr = make_shared<PBConstraint>(w);
  for (size_t i = 0; i < lits.size(); i++) {
    constr->addTerm(encode(lits[i]), coefs[i]);
  }

  if (constr->getWeight() <= 0) {
    return;
  } // already satisfied.
  constraints.insert(constr);
}

void PB::readPB(std::istream &input) {
  for (string line; getline(input, line);) {
    if (line.empty())
      continue;
    else if (line[0] == '*') {
      if (line.substr(0, 13) == "* #variable= ") {
        istringstream is(line.substr(13));
        is >> nVars;
      }
    } else {
      string symbol;
      if (line.find(">=") != string::npos)
        symbol = ">=";
      else
        symbol = "=";
      assert(line.find(symbol) != string::npos);
      istringstream is(line.substr(0, line.find(symbol)));
      vector<int> lits;
      vector<int> coefs;
      int coef;
      string var;
      while (is >> coef >> var) {
        assert(coef != 0 && abs(coef) <= (int)1e9);
        int x = atoi(var.substr(1).c_str());
        lits.push_back(x);
        coefs.push_back(coef);
      }
      int w = atoi(line.substr(line.find("=") + 1).c_str());
      processInequality(lits, coefs, w);
      // Handle equality case with two constraints
      if (line.find(" = ") != string::npos) {
        for (int &coef : coefs)
          coef = -coef;
        w *= -1;
        processInequality(lits, coefs, w);
      }
    }
  }
}

PB::PB(std::istream &input) {
  readPB(input);
  if (verbosity > 6) {
    std::clog << "Parsed PB theory:\n";
    // print(std::clog);
  }
  if (verbosity >= 2) {
    std::clog << "*** Creating first graph..." << std::endl;
  }
  graph = make_shared<Graph>(constraints);
  if (verbosity > 2) {
    std::clog << "**** Number of nodes: " << graph->getNbNodes() << std::endl;
    std::clog << "**** Number of edges: " << graph->getNbEdges() << std::endl;
  }

  group = make_shared<Group>();
  if (verbosity >= 2) {
    std::clog << "*** Detecting symmetry group..." << std::endl;
  }
  std::vector<std::shared_ptr<Permutation>> symgens;
  graph->getSymmetryGenerators(symgens);
  for (auto symgen : symgens) {
    group->add(symgen);
  }
}

PB::PB(std::vector<std::shared_ptr<PBConstraint>> &constr, std::shared_ptr<Group> grp) {
  constraints.insert(constr.cbegin(), constr.cend());
  graph = make_shared<Graph>(constraints);
  group = grp;
  for (uint l = 0; l < 2 * nVars; ++l) {
    if (not grp->permutes(l)) {
      graph->setUniqueColor(l);
    }
  }
  for (uint m = 0; m < grp->getNbMatrices(); ++m) {
    auto mat = grp->getMatrix(m);
    for (uint r = 0; r < mat->nbRows() - 1; ++r) {
      getGraph()->setUniqueColor(*(mat->getRow(r)));
    }
  }
}

uint PB::getSize() { return constraints.size(); }

void PB::setSubTheory(std::shared_ptr<Group> subgroup) {
  std::vector<std::shared_ptr<PBConstraint>> subconsts;
  for (auto r : constraints) {
    for (auto lc : r->getTerms()) {
      if (subgroup->permutes(lc.first)) {
        subconsts.push_back(r);
        break;
      }
    }
  }
  subgroup->theory = new PB(subconsts, subgroup);
}

bool PB::isSymmetry(Permutation &prm) {
  for (auto c : constraints) {
    std::shared_ptr<PBConstraint> symmetrical = std::make_shared<PBConstraint>(c->getWeight());
    if (!prm.getImage(c->terms, symmetrical->terms)) {
      continue;
    }

    if (constraints.count(symmetrical) == 0) {
      return false;
    }
  }
  return true;
}
} // namespace BreakID
