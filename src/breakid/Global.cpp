#include "Global.h"

#include <limits.h>

namespace BreakID {
using namespace std;

unsigned int nVars = 0;
std::vector<uint> fixedLits;
time_t startTime;

// OPTIONS:
bool useMatrixDetection = true;
bool useBinaryClauses = true;
bool useShatterTranslation = false;
bool useFullTranslation = false;
int symBreakingFormLength = 50;
int nbPBConstraintsToGroup = -1;
unsigned int verbosity = 2;
int timeLim = INT_MAX;

int timeLeft() {
  time_t now;
  time(&now);
  return timeLim - difftime(now, startTime);
}

bool timeLimitPassed() { return timeLeft() <= 0; }

void gracefulError(string str) {
  std::cerr << str << "\nExiting..." << endl;
  exit(1);
}

// Prints a clause c a rule falsevar <- (not c)
// afterwards, falsevar (a new variable) will be added to the "false" constraints.
void Clause::printAsRule(std::ostream &ostr, unsigned int falsevar) {
  ostr << "1 " << falsevar << " ";
  std::set<unsigned int> posBodyLits;
  std::set<unsigned int> negBodyLits;
  for (auto lit : lits) {
    auto decoded = decode(lit);
    if (decoded > 0) {
      posBodyLits.insert(decoded);
    } else {
      negBodyLits.insert(-decoded);
    }
  }
  ostr << (posBodyLits.size() + negBodyLits.size()) << " " << negBodyLits.size() << " ";
  for (auto decodedlit : negBodyLits) {
    ostr << decodedlit << " ";
  }
  for (auto decodedlit : posBodyLits) {
    ostr << decodedlit << " ";
  }
  ostr << std::endl;
}
} // namespace BreakID
