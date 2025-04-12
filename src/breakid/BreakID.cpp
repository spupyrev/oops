#include <cstring>
#include <iterator>
#include <sstream>
#include <stdlib.h>

#include "Algebraic.h"
#include "Breaking.h"
#include "Global.h"
#include "Graph.h"
#include "Theory.h"

using namespace std;
using namespace BreakID;

namespace options {
// option strings:
string nointch = "-no-row";
string nobinary = "-no-bin";
string formlength = "-s";
string verbosity = "-v";
string timelim = "-t";
string nosmall = "-no-small";
string norelaxed = "-no-relaxed";
string onlybreakers = "-print-only-breakers";
string generatorfile = "-store-sym";
} // namespace options

// std::clog << "\nOptions:\n";
// std::clog << options::help << "\n  ";
// std::clog << "Display this help message instead of running BreakID.\n";
// std::clog << options::instancefile << "\n";
// std::clog << "Read instance from a file instead of input stream.\n";
// std::clog << options::nointch << "\n  ";
// std::clog << "Disable detection and breaking of row interchangeability.\n";
// std::clog << options::nobinary << "\n  ";
// std::clog << "Disable construction of additional binary symmetry breaking clauses based on stabilizer subgroups.\n";
// std::clog << options::nosmall << "\n  ";
// std::clog << "Disable compact symmetry breaking encoding, use Shatter's encoding instead.\n";
// std::clog << options::norelaxed << "\n  ";
// std::clog << "Disable relaxing constraints on auxiliary encoding variables, use longer encoding instead.\n";
// std::clog << options::formlength << " <default: " << symBreakingFormLength << ">\n  ";
// std::clog << "Limit the size of the constructed symmetry breaking formula's, measured as the number of auxiliary
// variables introduced. <-1> means no symmetry breaking.\n"; std::clog << options::timelim << " <default: " << timeLim
// << ">\n  "; std::clog << "Upper limit on time spent by Saucy detecting symmetry measured in seconds.\n"; std::clog <<
// options::verbosity << " <default: " << verbosity << ">\n  "; std::clog << "Verbosity of the output. <0> means no
// output other than the CNF augmented with symmetry breaking clauses.\n"; std::clog << options::fixedvars << "
// <default: none>\n  "; std::clog << "File with a list of variables that should be fixed, separated by whitespace.\n";
// std::clog << options::onlybreakers << "\n  ";
// std::clog << "Do not print original theory, only the symmetry breaking clauses.\n";
// std::clog << options::generatorfile << "\n  ";
// std::clog << "Store the detected symmetry generators in the given file.\n";
// std::clog << options::symmetryinput << " <default: none>\n  ";
// std::clog << "Pass a file with symmetry generators or row-interchangeable matrices to use as additional symmetry
// information. Same format as BreakID's output by "
//           << options::generatorfile << ".\n";
// std::clog << options::aspinput << "\n  ";
// std::clog << "Parse input in the smodels-lparse intermediate format instead of DIMACS.\n";
// std::clog << options::pbinput << " <default:unset>\n  ";
// std::clog
//         << "When this variable is set, we assume input in the OPB format. The value of the integer indicates the
//         number of variables grouped together in one PB lex-leader constraint. In particular: 1 means almost clausal
//         encoding and infinity means one PB constraint with exponentially sized coefficients. Special case is 0, which
//         takes the PB encoding of symmetry breaking CNF clauses.\n";

std::pair<int, std::vector<std::vector<int>>> breakSymmetries(int verbose, int numVars,
                                                              const std::vector<std::vector<int>> &clauses) {
  verbosity = verbose;
  nVars = 0;
  fixedLits.clear();
  startTime = 0;

  std::shared_ptr<Specification> theory;
  theory = std::make_shared<CNF>(numVars, clauses);

  if (verbosity >= 3) {
    theory->getGraph()->print();
  }

  if (verbosity >= 2) {
    std::clog << "**** symmetry generators detected: " << theory->getGroup()->getSize() << std::endl;
    if (verbosity > 2) {
      theory->getGroup()->print(std::clog);
    }
  }

  if (verbosity >= 2) {
    std::clog << "*** Detecting subgroups..." << std::endl;
  }

  vector<std::shared_ptr<Group>> subgroups;
  theory->getGroup()->getDisjointGenerators(subgroups);
  if (verbosity >= 2) {
    std::clog << "**** subgroups detected: " << subgroups.size() << std::endl;
  }

  if (verbosity > 2) {
    for (auto grp : subgroups) {
      std::clog << "group size: " << grp->getSize() << " support: " << grp->getSupportSize() << std::endl;
      if (verbosity > 2) {
        grp->print(std::clog);
      }
    }
  }

  theory->cleanUp(); // improve some memory overhead

  uint totalNbMatrices = 0;
  uint totalNbRowSwaps = 0;

  Breaker brkr(theory);
  for (auto grp : subgroups) {
    if (grp->getSize() > 1 && useMatrixDetection) {
      if (verbosity >= 2) {
        std::clog << "*** Detecting row interchangeability..." << std::endl;
      }
      theory->setSubTheory(grp);
      grp->addMatrices();
      totalNbMatrices += grp->getNbMatrices();
      totalNbRowSwaps += grp->getNbRowSwaps();
    }
    if (symBreakingFormLength > -1) {
      if (verbosity >= 2) {
        std::clog << "*** Constructing symmetry breaking formula..." << std::endl;
      }
      grp->addBreakingClausesTo(brkr);
    } // else no symmetry breaking formulas are needed :)
    grp.reset();
  }

  if (verbosity >= 2) {
    std::clog << "**** matrices detected: " << totalNbMatrices << std::endl;
    std::clog << "**** row swaps detected: " << totalNbRowSwaps << std::endl;
    std::clog << "**** extra binary symmetry breaking clauses added: " << brkr.getNbBinClauses() << "\n";
    std::clog << "**** regular symmetry breaking clauses added: " << brkr.getNbRegClauses() << "\n";
    std::clog << "**** row interchangeability breaking clauses added: " << brkr.getNbRowClauses() << "\n";
    std::clog << "**** total symmetry breaking clauses added: " << brkr.getAddedNbClauses() << "\n";
    std::clog << "**** auxiliary variables introduced: " << brkr.getAuxiliaryNbVars() << "\n";
  }

  int newNumVars;
  std::vector<std::vector<int>> newClauses;

  brkr.result(newNumVars, newClauses);

  return std::make_pair(newNumVars, newClauses);
}
