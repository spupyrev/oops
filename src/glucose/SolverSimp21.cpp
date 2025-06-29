#include "SolverSimp21.h"
#include "Sort.h"

#include <algorithm>
#include <cmath>
#include <signal.h>
#include <unistd.h>

using namespace std;
using namespace Simp21;

/// Options:
static const char *_cat = "CORE";

static DoubleOption opt_step_size(_cat, "step-size", "Initial step size", 0.40, DoubleRange(0, false, 1, false));
static DoubleOption opt_step_size_dec(_cat, "step-size-dec", "Step size decrement", 0.000001,
                                      DoubleRange(0, false, 1, false));
static DoubleOption opt_min_step_size(_cat, "min-step-size", "Minimal step size", 0.06,
                                      DoubleRange(0, false, 1, false));
static DoubleOption opt_var_decay(_cat, "var-decay", "The variable activity decay factor", 0.80,
                                  DoubleRange(0, false, 1, false));
static DoubleOption opt_clause_decay(_cat, "cla-decay", "The clause activity decay factor", 0.999,
                                     DoubleRange(0, false, 1, false));
static DoubleOption
    opt_random_var_freq(_cat, "rnd-freq",
                        "The frequency with which the decision heuristic tries to choose a random variable", 0,
                        DoubleRange(0, true, 1, true));
static DoubleOption opt_random_seed(_cat, "rnd-seed", "Used by the random variable selection", 91648253,
                                    DoubleRange(0, false, HUGE_VAL, false));
static IntOption opt_ccmin_mode(_cat, "ccmin-mode", "Controls conflict clause minimization (0=none, 1=basic, 2=deep)",
                                2, IntRange(0, 2));
static IntOption opt_phase_saving(_cat, "phase-saving",
                                  "Controls the level of phase saving (0=none, 1=limited, 2=full)", 2, IntRange(0, 2));
static BoolOption opt_rnd_init_act(_cat, "rnd-init", "Randomize the initial activity", false);
static IntOption opt_restart_first(_cat, "rfirst", "The base restart interval", 100, IntRange(1, INT32_MAX));
static DoubleOption opt_restart_inc(_cat, "rinc", "Restart interval increase factor", 2,
                                    DoubleRange(1, false, HUGE_VAL, false));
static DoubleOption opt_garbage_frac(_cat, "gc-frac",
                                     "The fraction of wasted memory allowed before a garbage collection is triggered",
                                     0.20, DoubleRange(0, false, HUGE_VAL, false));
static IntOption opt_chrono(_cat, "chrono", "Controls if to perform chrono backtrack", 100, IntRange(-1, INT32_MAX));
static IntOption opt_conf_to_chrono(_cat, "confl-to-chrono", "Controls number of conflicts to perform chrono backtrack",
                                    4000, IntRange(-1, INT32_MAX));

static IntOption opt_max_lbd_dup("DUP-LEARNTS", "lbd-limit",
                                 "specifies the maximum lbd of learnts to be screened for duplicates.", 12,
                                 IntRange(0, INT32_MAX));
static IntOption opt_min_dupl_app("DUP-LEARNTS", "min-dup-app",
                                  "specifies the minimum number of learnts to be included into db.", 3,
                                  IntRange(2, INT32_MAX));
static IntOption opt_dupl_db_init_size("DUP-LEARNTS", "dupdb-init", "specifies the initial maximal duplicates DB size.",
                                       500000, IntRange(1, INT32_MAX));

static IntOption opt_VSIDS_props_limit(
    "DUP-LEARNTS", "VSIDS-lim",
    "specifies the number of propagations after which the solver switches between LRB and VSIDS(in millions).", 30,
    IntRange(1, INT32_MAX));

// VSIDS_props_limit

//=================================================================================================
// Constructor/Destructor:

Solver::Solver()
    :
      // Parameters (user settable):
      //
      verbosity(0), step_size(opt_step_size), step_size_dec(opt_step_size_dec),
      min_step_size(opt_min_step_size), timer(5000), var_decay(opt_var_decay), clause_decay(opt_clause_decay),
      random_var_freq(opt_random_var_freq), random_seed(opt_random_seed), VSIDS(false), ccmin_mode(opt_ccmin_mode),
      phase_saving(opt_phase_saving), rnd_pol(false), rnd_init_act(opt_rnd_init_act), garbage_frac(opt_garbage_frac),
      restart_first(opt_restart_first), restart_inc(opt_restart_inc),

      // Parameters (the rest):
      //
      learntsize_factor((double)1 / (double)3), learntsize_inc(1.1),

      // Parameters (experimental):
      //
      learntsize_adjust_start_confl(100), learntsize_adjust_inc(1.5),

      // Statistics: 
      //
      solves(0), starts(0), decisions(0), rnd_decisions(0), propagations(0), conflicts(0), conflicts_VSIDS(0),
      dec_vars(0), clauses_literals(0), learnts_literals(0), max_literals(0), tot_literals(0), chrono_backtrack(0),
      non_chrono_backtrack(0),

      ok(true), cla_inc(1), var_inc(1), watches_bin(WatcherDeleted(ca)), watches(WatcherDeleted(ca)), qhead(0),
      simpDB_assigns(-1), simpDB_props(0), order_heap_CHB(VarOrderLt(activity_CHB)),
      order_heap_VSIDS(VarOrderLt(activity_VSIDS)), order_heap_distance(VarOrderLt(activity_distance)), remove_satisfied(true),

      core_lbd_cut(2), global_lbd_sum(0), lbd_queue(50), next_T2_reduce(10000), next_L_reduce(15000),
      confl_to_chrono(opt_conf_to_chrono), chrono(opt_chrono),

      counter(0),

      // Resource constraints:
      //
      conflict_budget(-1), propagation_budget(-1), asynch_interrupt(false),

      // simplfiy
      nbSimplifyAll(0), s_propagations(0),

      // simplifyAll adjust occasion
      curSimplify(1), nbconfbeforesimplify(1000), incSimplify(1000) {
  min_number_of_learnts_copies = opt_min_dupl_app;
  max_lbd_dup = opt_max_lbd_dup;
  dupl_db_init_size = opt_dupl_db_init_size;
  VSIDS_props_limit = opt_VSIDS_props_limit * 1000000;
  DISTANCE = true;
  cur_solustion_used = false;
  freeze_restart_num = 0;
  best_unsat_num = INT_MAX;
  max_trail = 0;
  timeout_ms = -1;
  start_time = std::chrono::system_clock::now();
}

Solver::~Solver() {}

// simplify All
CRef Solver::simplePropagate() {
  CRef confl = CRef_Undef;
  int num_props = 0;
  watches.cleanAll();
  watches_bin.cleanAll();

  while (qhead < trail.size()) {
    Lit p = trail[qhead++]; // 'p' is enqueued fact to propagate.
    vec<Watcher> &ws = watches[p];
    Watcher *i, *j, *end;
    num_props++;
    // First, propagate binary clauses
    vec<Watcher> &wbin = watches_bin[p];

    for (int k = 0; k < wbin.size(); k++) {
      Lit imp = wbin[k].blocker;

      if (value(imp) == l_False) {
        return wbin[k].cref;
      }

      if (value(imp) == l_Undef) {
        simpleUncheckEnqueue(imp, wbin[k].cref);
      }
    }

    for (i = j = (Watcher *)ws, end = i + ws.size(); i != end;) {
      // Try to avoid inspecting the clause:
      const Lit blocker = i->blocker;

      if (value(blocker) == l_True) {
        *j++ = *i++;
        continue;
      }

      // Make sure the false literal is data[1]:
      const CRef cr = i->cref;
      Clause &c = ca[cr];
      const Lit false_lit = ~p;

      if (c[0] == false_lit) {
        c[0] = c[1], c[1] = false_lit;
      }

      assert(c[1] == false_lit);
      //  i++;
      // If 0th watch is true, then clause is already satisfied.
      // However, 0th watch is not the blocker, make it blocker using a new watcher w
      // why not simply do i->blocker=first in this case?
      Lit first = c[0];

      if (first != blocker && value(first) == l_True) {
        i->blocker = first;
        *j++ = *i++;
        continue;
      } else {
        // ----------------- DEFAULT  MODE (NOT INCREMENTAL)
        for (int k = 2; k < c.size(); k++) {
          if (value(c[k]) != l_False) {
            // watcher i is abandonned using i++, because cr watches now ~c[k] instead of p
            // the blocker is first in the watcher. However,
            // the blocker in the corresponding watcher in ~first is not c[1]
            Watcher w(cr, first);
            i++;
            c[1] = c[k];
            c[k] = false_lit;
            watches[~c[1]].push(w);
            goto NextClause;
          }
        }
      }

      // Did not find watch -- clause is unit under assignment:
      i->blocker = first;
      *j++ = *i++;

      if (value(first) == l_False) {
        confl = cr;
        qhead = trail.size();

        // Copy the remaining watches:
        while (i < end) {
          *j++ = *i++;
        }
      } else {
        simpleUncheckEnqueue(first, cr);
      }

    NextClause:;
    }

    ws.shrink(i - j);
  }

  s_propagations += num_props;
  return confl;
}

void Solver::simpleUncheckEnqueue(Lit p, CRef from) {
  assert(value(p) == l_Undef);
  assigns[var(p)] = lbool(!sign(p)); // this makes a lbool object whose value is sign(p)
  vardata[var(p)].reason = from;
  trail.push_(p);
}

void Solver::cancelUntilTrailRecord() {
  for (int c = trail.size() - 1; c >= trailRecord; c--) {
    Var x = var(trail[c]);
    assigns[x] = l_Undef;
  }

  qhead = trailRecord;
  trail.shrink(trail.size() - trailRecord);
}

void Solver::litsEnqueue(int cutP, Clause &c) {
  for (int i = cutP; i < c.size(); i++) {
    simpleUncheckEnqueue(~c[i]);
  }
}

bool Solver::removed(CRef cr) { return ca[cr].mark() == 1; }

void Solver::simpleAnalyze(CRef confl, vec<Lit> &out_learnt, vec<CRef> &reason_clause, bool True_confl) {
  int pathC = 0;
  Lit p = lit_Undef;
  int index = trail.size() - 1;

  do {
    if (confl != CRef_Undef) {
      reason_clause.push(confl);
      Clause &c = ca[confl];

      // Special case for binary clauses
      // The first one has to be SAT
      if (p != lit_Undef && c.size() == 2 && value(c[0]) == l_False) {
        assert(value(c[1]) == l_True);
        Lit tmp = c[0];
        c[0] = c[1], c[1] = tmp;
      }

      // if True_confl==true, then choose p begin with the 1th index of c;
      for (int j = (p == lit_Undef && True_confl == false) ? 0 : 1; j < c.size(); j++) {
        Lit q = c[j];

        if (!seen[var(q)]) {
          seen[var(q)] = 1;
          pathC++;
        }
      }
    } else if (confl == CRef_Undef) {
      out_learnt.push(~p);
    }

    // if not break, while() will come to the index of trail blow 0, and fatal error occur;
    if (pathC == 0) {
      break;
    }

    // Select next clause to look at:
    while (!seen[var(trail[index--])])
      ;

    // if the reason cr from the 0-level assigned var, we must break avoid move forth further;
    // but attention that maybe seen[x]=1 and never be clear. However makes no matter;
    if (trailRecord > index + 1) {
      break;
    }

    p = trail[index + 1];
    confl = reason(var(p));
    seen[var(p)] = 0;
    pathC--;
  } while (pathC >= 0);
}

void Solver::simplifyLearnt(Clause &c) {
  original_length_record += c.size();
  trailRecord = trail.size(); // record the start pointer
  vec<Lit> falseLit;
  falseLit.clear();
  bool True_confl = false;
  int i, j;
  CRef confl;

  for (i = 0, j = 0; i < c.size(); i++) {
    if (value(c[i]) == l_Undef) {
      simpleUncheckEnqueue(~c[i]);
      c[j++] = c[i];
      confl = simplePropagate();

      if (confl != CRef_Undef) {
        break;
      }
    } else {
      if (value(c[i]) == l_True) {
        c[j++] = c[i];
        True_confl = true;
        confl = reason(var(c[i]));
        break;
      } else {
        falseLit.push(c[i]);
      }
    }
  }

  c.shrink(c.size() - j);

  if (confl != CRef_Undef || True_confl == true) {
    simp_learnt_clause.clear();
    simp_reason_clause.clear();

    if (True_confl == true) {
      simp_learnt_clause.push(c.last());
    }

    simpleAnalyze(confl, simp_learnt_clause, simp_reason_clause, True_confl);

    if (simp_learnt_clause.size() < c.size()) {
      for (i = 0; i < simp_learnt_clause.size(); i++) {
        c[i] = simp_learnt_clause[i];
      }

      c.shrink(c.size() - i);
    }
  }

  cancelUntilTrailRecord();
  simplified_length_record += c.size();
}

bool Solver::simplifyLearnt_core() {
  int ci, cj, li, lj;
  bool sat, false_lit;
  unsigned int nblevels;

  for (ci = 0, cj = 0; ci < learnts.size(); ci++) {
    CRef cr = learnts[ci];
    Clause &c = ca[cr];

    if (c.mark() == REMOVED) {
      continue;
    }

    if (c.done() || c.mark() == LONG) {
      learnts[cj++] = learnts[ci];
    } else {
      int saved_size = c.size();
      sat = false_lit = false;

      for (int i = 0; i < saved_size; i++) {
        if (value(c[i]) == l_True) {
          sat = true;
          break;
        } else if (value(c[i]) == l_False) {
          false_lit = true;
        }
      }

      if (sat) {
        removeClause(cr);
      } else {
        detachClause(cr, true);

        if (false_lit) {
          for (li = lj = 0; li < c.size(); li++) {
            if (value(c[li]) != l_False) {
              c[lj++] = c[li];
            }
          }

          c.shrink(li - lj);
        }

        simplifyLearnt(c);

        if (c.size() == 1) {
          // when unit clause occur, enqueue and propagate
          uncheckedEnqueue(c[0]);

          if (propagate() != CRef_Undef) {
            ok = false;
            return false;
          }

          // delete the clause memory in logic
          c.mark(REMOVED);
          ca.free(cr);
        } else {
          attachClause(cr);
          learnts[cj++] = learnts[ci];
          nblevels = computeLBD(c);

          if ((signed)nblevels < c.lbd()) {
            c.set_lbd(nblevels);
          }

          if (c.lbd() <= core_lbd_cut) {
            c.mark(SMALL);
          }

          c.done(true);
        }
      }
    }
  }

  learnts.shrink(ci - cj);
  return true;
}

int Solver::is_duplicate(std::vector<uint32_t> &c) {
  auto time_point_0 = std::chrono::high_resolution_clock::now();
  dupl_db_size++;
  int res = 0;
  int sz = c.size();
  std::vector<uint32_t> tmp(c);
  sort(tmp.begin(), tmp.end());
  uint64_t hash = 0;

  for (int i = 0; i < sz; i++) {
    hash ^= tmp[i] + 0x9e3779b9 + (hash << 6) + (hash >> 2);
  }

  int32_t head = tmp[0];
  auto it0 = ht.find(head);

  if (it0 != ht.end()) {
    auto it1 = ht[head].find(sz);

    if (it1 != ht[head].end()) {
      auto it2 = ht[head][sz].find(hash);

      if (it2 != ht[head][sz].end()) {
        it2->second++;
        res = it2->second;
      } else {
        ht[head][sz][hash] = 1;
      }
    } else {
      ht[head][sz][hash] = 1;
    }
  } else {
    ht[head][sz][hash] = 1;
  }

  auto time_point_1 = std::chrono::high_resolution_clock::now();
  duptime += std::chrono::duration_cast<std::chrono::microseconds>(time_point_1 - time_point_0);
  return res;
}

bool Solver::simplifyAll() {
  simplified_length_record = original_length_record = 0;

  if (!ok || propagate() != CRef_Undef) {
    return ok = false;
  }

  if (!simplifyLearnt_core()) {
    return ok = false;
  }

  checkGarbage();
  return true;
}
//=================================================================================================
// Minor methods:

// Creates a new SAT variable in the solver. If 'decision' is cleared, variable will not be
// used as a decision variable (NOTE! This has effects on the meaning of a SATISFIABLE result).
//
Var Solver::newVar(bool sign, bool dvar) {
  int v = nVars();
  watches_bin.init(mkLit(v, false));
  watches_bin.init(mkLit(v, true));
  watches.init(mkLit(v, false));
  watches.init(mkLit(v, true));
  assigns.push(l_Undef);
  vardata.push(mkVarData(CRef_Undef, 0));
  activity_CHB.push(0);
  activity_VSIDS.push(rnd_init_act ? drand(random_seed) * 0.00001 : 0);
  picked.push(0);
  conflicted.push(0);
  almost_conflicted.push(0);
  canceled.push(0);
  seen.push(0);
  seen2.push(0);
  polarity.push(sign);
  decision.push();
  trail.capacity(v + 1);
  setDecisionVar(v, dvar);
  activity_distance.push(0);
  return v;
}

bool Solver::addClause_(vec<Lit> &ps) {
  assert(decisionLevel() == 0);

  if (!ok) {
    return false;
  }

  // Check if clause is satisfied and remove false/duplicate literals:
  sort(ps);
  Lit p;
  int i, j;

  for (i = j = 0, p = lit_Undef; i < ps.size(); i++)
    if (value(ps[i]) == l_True || ps[i] == ~p) {
      return true;
    } else if (value(ps[i]) != l_False && ps[i] != p) {
      ps[j++] = p = ps[i];
    }

  ps.shrink(i - j);

  if (ps.size() == 0) {
    ok = false;
    return ok;
  } else if (ps.size() == 1) {
    uncheckedEnqueue(ps[0]);
    return ok = (propagate() == CRef_Undef);
  } else {
    CRef cr = ca.alloc(ps, false);
    clauses.push(cr);
    attachClause(cr);
  }

  return true;
}

void Solver::attachClause(CRef cr) {
  const Clause &c = ca[cr];
  assert(c.size() > 1);
  OccLists<Lit, vec<Watcher>, WatcherDeleted> &ws = c.size() == 2 ? watches_bin : watches;
  ws[~c[0]].push(Watcher(cr, c[1]));
  ws[~c[1]].push(Watcher(cr, c[0]));

  if (c.learnt()) {
    learnts_literals += c.size();
  } else {
    clauses_literals += c.size();
  }
}

void Solver::detachClause(CRef cr, bool strict) {
  const Clause &c = ca[cr];
  assert(c.size() > 1);
  OccLists<Lit, vec<Watcher>, WatcherDeleted> &ws = c.size() == 2 ? watches_bin : watches;

  if (strict) {
    remove(ws[~c[0]], Watcher(cr, c[1]));
    remove(ws[~c[1]], Watcher(cr, c[0]));
  } else {
    // Lazy detaching: (NOTE! Must clean all watcher lists before garbage collecting this clause)
    ws.smudge(~c[0]);
    ws.smudge(~c[1]);
  }

  if (c.learnt()) {
    learnts_literals -= c.size();
  } else {
    clauses_literals -= c.size();
  }
}

void Solver::removeClause(CRef cr) {
  Clause &c = ca[cr];

  detachClause(cr);

  // Don't leave pointers to free'd memory!
  if (locked(c)) {
    Lit implied = c.size() != 2 ? c[0] : (value(c[0]) == l_True ? c[0] : c[1]);
    vardata[var(implied)].reason = CRef_Undef;
  }

  c.mark(1);
  ca.free(cr);
}

bool Solver::satisfied(const Clause &c) const {
  for (int i = 0; i < c.size(); i++)
    if (value(c[i]) == l_True) {
      return true;
    }

  return false;
}

// Revert to the state at given level (keeping all assignment at 'level' but not beyond).
//
void Solver::cancelUntil(int bLevel) {
  if (decisionLevel() > bLevel) {
    add_tmp.clear();

    for (int c = trail.size() - 1; c >= trail_lim[bLevel]; c--) {
      Var x = var(trail[c]);

      if (level(x) <= bLevel) {
        add_tmp.push(trail[c]);
      } else {
        if (!VSIDS) {
          uint32_t age = conflicts - picked[x];

          if (age > 0) {
            double adjusted_reward = ((double)(conflicted[x] + almost_conflicted[x])) / ((double)age);
            double old_activity = activity_CHB[x];
            activity_CHB[x] = step_size * adjusted_reward + ((1 - step_size) * old_activity);

            if (order_heap_CHB.inHeap(x)) {
              if (activity_CHB[x] > old_activity) {
                order_heap_CHB.decrease(x);
              } else {
                order_heap_CHB.increase(x);
              }
            }
          }

          canceled[x] = conflicts;
        }

        assigns[x] = l_Undef;
        if (phase_saving > 1 || ((phase_saving == 1) && c > trail_lim.last())) {
          polarity[x] = sign(trail[c]);
        }

        insertVarOrder(x);
      }
    }

    qhead = trail_lim[bLevel];
    trail.shrink(trail.size() - trail_lim[bLevel]);
    trail_lim.shrink(trail_lim.size() - bLevel);

    for (int nLitId = add_tmp.size() - 1; nLitId >= 0; --nLitId) {
      trail.push_(add_tmp[nLitId]);
    }

    add_tmp.clear();
  }
}

//=================================================================================================
// Major methods:

Lit Solver::pickBranchLit() {
  Var next = var_Undef;
  Heap<VarOrderLt> &order_heap = DISTANCE ? order_heap_distance : ((!VSIDS) ? order_heap_CHB : order_heap_VSIDS);

  while (next == var_Undef || value(next) != l_Undef || !decision[next])
    if (order_heap.empty()) {
      return lit_Undef;
    } else {
      if (!VSIDS) {
        Var v = order_heap_CHB[0];
        uint32_t age = conflicts - canceled[v];

        while (age > 0) {
          double decay = pow(0.95, age);
          activity_CHB[v] *= decay;

          if (order_heap_CHB.inHeap(v)) {
            order_heap_CHB.increase(v);
          }

          canceled[v] = conflicts;
          v = order_heap_CHB[0];
          age = conflicts - canceled[v];
        }
      }

      next = order_heap.removeMin();
    }

  return mkLit(next, polarity[next]);
}

inline Solver::ConflictData Solver::FindConflictLevel(CRef cind) {
  ConflictData data;
  Clause &conflCls = ca[cind];
  data.nHighestLevel = level(var(conflCls[0]));

  if (data.nHighestLevel == decisionLevel() && level(var(conflCls[1])) == decisionLevel()) {
    return data;
  }

  int highestId = 0;
  data.bOnlyOneLitFromHighest = true;

  // find the largest decision level in the clause
  for (int nLitId = 1; nLitId < conflCls.size(); ++nLitId) {
    int nLevel = level(var(conflCls[nLitId]));

    if (nLevel > data.nHighestLevel) {
      highestId = nLitId;
      data.nHighestLevel = nLevel;
      data.bOnlyOneLitFromHighest = true;
    } else if (nLevel == data.nHighestLevel && data.bOnlyOneLitFromHighest == true) {
      data.bOnlyOneLitFromHighest = false;
    }
  }

  if (highestId != 0) {
    std::swap(conflCls[0], conflCls[highestId]);

    if (highestId > 1) {
      OccLists<Lit, vec<Watcher>, WatcherDeleted> &ws = conflCls.size() == 2 ? watches_bin : watches;
      remove(ws[~conflCls[highestId]], Watcher(cind, conflCls[1]));
      ws[~conflCls[0]].push(Watcher(cind, conflCls[1]));
    }
  }

  return data;
}

/*_________________________________________________________________________________________________
|
|  analyze : (confl : Clause*) (out_learnt : vec<Lit>&) (out_btlevel : int&)  ->  [void]
|
|  Description:
|    Analyze conflict and produce a reason clause.
|
|    Pre-conditions:
|      * 'out_learnt' is assumed to be cleared.
|      * Current decision level must be greater than root level.
|
|    Post-conditions:
|      * 'out_learnt[0]' is the asserting literal at level 'out_btlevel'.
|      * If out_learnt.size() > 1 then 'out_learnt[1]' has the greatest decision level of the
|        rest of literals. There may be others from the same level though.
|
|________________________________________________________________________________________________@*/
void Solver::analyze(CRef confl, vec<Lit> &out_learnt, int &out_btlevel, int &out_lbd) {
  int pathC = 0;
  Lit p = lit_Undef;
  // Generate conflict clause:
  //
  out_learnt.push(); // (leave room for the asserting literal)
  int index = trail.size() - 1;
  int nDecisionLevel = level(var(ca[confl][0]));
  assert(nDecisionLevel == level(var(ca[confl][0])));

  do {
    assert(confl != CRef_Undef); // (otherwise should be UIP)
    Clause &c = ca[confl];

    // For binary clauses, we don't rearrange literals in propagate(), so check and make sure the first is an implied
    // lit.
    if (p != lit_Undef && c.size() == 2 && value(c[0]) == l_False) {
      assert(value(c[1]) == l_True);
      Lit tmp = c[0];
      c[0] = c[1], c[1] = tmp;
    }

    // Update LBD if improved.
    if (c.learnt() && c.mark() != SMALL) {
      int lbd = computeLBD(c);

      if (lbd < c.lbd()) {
        if (c.lbd() <= 30) {
          c.removable(false); // Protect once from reduction.
        }

        c.set_lbd(lbd);

        if (lbd <= core_lbd_cut) {
          c.mark(SMALL);
        } else if (lbd <= 6 && c.mark() == LONG) {
          c.mark(MIDSZ);
        }
      }

      if (c.mark() == MIDSZ) {
        c.touched() = conflicts;
      } else if (c.mark() == LONG) {
        claBumpActivity(c);
      }
    }

    for (int j = (p == lit_Undef) ? 0 : 1; j < c.size(); j++) {
      Lit q = c[j];

      if (!seen[var(q)] && level(var(q)) > 0) {
        if (VSIDS) {
          varBumpActivity(var(q), .5);
          add_tmp.push(q);
        } else {
          conflicted[var(q)]++;
        }

        seen[var(q)] = 1;

        if (level(var(q)) >= nDecisionLevel) {
          pathC++;
        } else {
          out_learnt.push(q);
        }
      }
    }

    // Select next clause to look at:
    do {
      while (!seen[var(trail[index--])])
        ;

      p = trail[index + 1];
    } while (level(var(p)) < nDecisionLevel);

    confl = reason(var(p));
    seen[var(p)] = 0;
    pathC--;
  } while (pathC > 0);

  out_learnt[0] = ~p;
  // Simplify conflict clause:
  //
  int i, j;
  out_learnt.copyTo(analyze_toclear);

  if (ccmin_mode == 2) {
    uint32_t abstract_level = 0;

    for (i = 1; i < out_learnt.size(); i++) {
      abstract_level |= abstractLevel(var(out_learnt[i])); // (maintain an abstraction of levels involved in conflict)
    }

    for (i = j = 1; i < out_learnt.size(); i++)
      if (reason(var(out_learnt[i])) == CRef_Undef || !litRedundant(out_learnt[i], abstract_level)) {
        out_learnt[j++] = out_learnt[i];
      }
  } else if (ccmin_mode == 1) {
    for (i = j = 1; i < out_learnt.size(); i++) {
      Var x = var(out_learnt[i]);

      if (reason(x) == CRef_Undef) {
        out_learnt[j++] = out_learnt[i];
      } else {
        Clause &c = ca[reason(var(out_learnt[i]))];

        for (int k = c.size() == 2 ? 0 : 1; k < c.size(); k++)
          if (!seen[var(c[k])] && level(var(c[k])) > 0) {
            out_learnt[j++] = out_learnt[i];
            break;
          }
      }
    }
  } else {
    i = j = out_learnt.size();
  }

  max_literals += out_learnt.size();
  out_learnt.shrink(i - j);
  tot_literals += out_learnt.size();
  out_lbd = computeLBD(out_learnt);

  if (out_lbd <= 6 && out_learnt.size() <= 30) // Try further minimization?
    if (binResMinimize(out_learnt)) {
      out_lbd = computeLBD(out_learnt); // Recompute LBD if minimized.
    }

  // Find correct backtrack level:
  //
  if (out_learnt.size() == 1) {
    out_btlevel = 0;
  } else {
    int max_i = 1;

    // Find the first literal assigned at the next-highest level:
    for (int i = 2; i < out_learnt.size(); i++)
      if (level(var(out_learnt[i])) > level(var(out_learnt[max_i]))) {
        max_i = i;
      }

    // Swap-in this literal at index 1:
    Lit p = out_learnt[max_i];
    out_learnt[max_i] = out_learnt[1];
    out_learnt[1] = p;
    out_btlevel = level(var(p));
  }

  if (VSIDS) {
    for (int i = 0; i < add_tmp.size(); i++) {
      Var v = var(add_tmp[i]);

      if (level(v) >= out_btlevel - 1) {
        varBumpActivity(v, 1);
      }
    }

    add_tmp.clear();
  } else {
    seen[var(p)] = true;

    for (int i = out_learnt.size() - 1; i >= 0; i--) {
      Var v = var(out_learnt[i]);
      CRef rea = reason(v);

      if (rea != CRef_Undef) {
        const Clause &reaC = ca[rea];

        for (int i = 0; i < reaC.size(); i++) {
          Lit l = reaC[i];

          if (!seen[var(l)]) {
            seen[var(l)] = true;
            almost_conflicted[var(l)]++;
            analyze_toclear.push(l);
          }
        }
      }
    }
  }

  for (int j = 0; j < analyze_toclear.size(); j++) {
    seen[var(analyze_toclear[j])] = 0; // ('seen[]' is now cleared)
  }
}

// Try further learnt clause minimization by means of binary clause resolution.
bool Solver::binResMinimize(vec<Lit> &out_learnt) {
  // Preparation: remember which false variables we have in 'out_learnt'.
  counter++;

  for (int i = 1; i < out_learnt.size(); i++) {
    seen2[var(out_learnt[i])] = counter;
  }

  // Get the list of binary clauses containing 'out_learnt[0]'.
  const vec<Watcher> &ws = watches_bin[~out_learnt[0]];
  int to_remove = 0;

  for (int i = 0; i < ws.size(); i++) {
    Lit the_other = ws[i].blocker;

    // Does 'the_other' appear negatively in 'out_learnt'?
    if (seen2[var(the_other)] == counter && value(the_other) == l_True) {
      to_remove++;
      seen2[var(the_other)] = counter - 1; // Remember to remove this variable.
    }
  }

  // Shrink.
  if (to_remove > 0) {
    int last = out_learnt.size() - 1;

    for (int i = 1; i < out_learnt.size() - to_remove; i++)
      if (seen2[var(out_learnt[i])] != counter) {
        out_learnt[i--] = out_learnt[last--];
      }

    out_learnt.shrink(to_remove);
  }

  return to_remove != 0;
}

// Check if 'p' can be removed. 'abstract_levels' is used to abort early if the algorithm is
// visiting literals at levels that cannot be removed later.
bool Solver::litRedundant(Lit p, uint32_t abstract_levels) {
  analyze_stack.clear();
  analyze_stack.push(p);
  int top = analyze_toclear.size();

  while (analyze_stack.size() > 0) {
    assert(reason(var(analyze_stack.last())) != CRef_Undef);
    Clause &c = ca[reason(var(analyze_stack.last()))];
    analyze_stack.pop();

    // Special handling for binary clauses like in 'analyze()'.
    if (c.size() == 2 && value(c[0]) == l_False) {
      assert(value(c[1]) == l_True);
      Lit tmp = c[0];
      c[0] = c[1], c[1] = tmp;
    }

    for (int i = 1; i < c.size(); i++) {
      Lit p = c[i];

      if (!seen[var(p)] && level(var(p)) > 0) {
        if (reason(var(p)) != CRef_Undef && (abstractLevel(var(p)) & abstract_levels) != 0) {
          seen[var(p)] = 1;
          analyze_stack.push(p);
          analyze_toclear.push(p);
        } else {
          for (int j = top; j < analyze_toclear.size(); j++) {
            seen[var(analyze_toclear[j])] = 0;
          }

          analyze_toclear.shrink(analyze_toclear.size() - top);
          return false;
        }
      }
    }
  }

  return true;
}

/*_________________________________________________________________________________________________
|
|  analyzeFinal : (p : Lit)  ->  [void]
|
|  Description:
|    Specialized analysis procedure to express the final conflict in terms of assumptions.
|    Calculates the (possibly empty) set of assumptions that led to the assignment of 'p', and
|    stores the result in 'out_conflict'.
|________________________________________________________________________________________________@*/
void Solver::analyzeFinal(Lit p, vec<Lit> &out_conflict) {
  out_conflict.clear();
  out_conflict.push(p);

  if (decisionLevel() == 0) {
    return;
  }

  seen[var(p)] = 1;

  for (int i = trail.size() - 1; i >= trail_lim[0]; i--) {
    Var x = var(trail[i]);

    if (seen[x]) {
      if (reason(x) == CRef_Undef) {
        assert(level(x) > 0);
        out_conflict.push(~trail[i]);
      } else {
        Clause &c = ca[reason(x)];

        for (int j = c.size() == 2 ? 0 : 1; j < c.size(); j++)
          if (level(var(c[j])) > 0) {
            seen[var(c[j])] = 1;
          }
      }

      seen[x] = 0;
    }
  }

  seen[var(p)] = 0;
}

void Solver::uncheckedEnqueue(Lit p, int level, CRef from) {
  assert(value(p) == l_Undef);
  Var x = var(p);

  if (!VSIDS) {
    picked[x] = conflicts;
    conflicted[x] = 0;
    almost_conflicted[x] = 0;
    uint32_t age = conflicts - canceled[var(p)];

    if (age > 0) {
      double decay = pow(0.95, age);
      activity_CHB[var(p)] *= decay;

      if (order_heap_CHB.inHeap(var(p))) {
        order_heap_CHB.increase(var(p));
      }
    }
  }

  assigns[x] = lbool(!sign(p));
  vardata[x] = mkVarData(from, level);
  trail.push_(p);
}

/*_________________________________________________________________________________________________
|
|  propagate : [void]  ->  [Clause*]
|
|  Description:
|    Propagates all enqueued facts. If a conflict arises, the conflicting clause is returned,
|    otherwise CRef_Undef.
|
|    Post-conditions:
|      * the propagation queue is empty, even if there was a conflict.
|________________________________________________________________________________________________@*/
CRef Solver::propagate() {
  CRef confl = CRef_Undef;
  int num_props = 0;
  watches.cleanAll();
  watches_bin.cleanAll();

  while (qhead < trail.size()) {
    const Lit p = trail[qhead++]; // 'p' is enqueued fact to propagate.
    const Lit false_lit = ~p;
    const int currLevel = level(var(p));
    Watcher *i, *j, *end;
    num_props++;

    // Propagate binary clauses first
    const vec<Watcher> &ws_bin = watches_bin[p]; 
    for (int k = 0; k < ws_bin.size(); k++) {
      const Lit the_other = ws_bin[k].blocker;

      if (value(the_other) == l_False) {
        return ws_bin[k].cref;
      } else if (value(the_other) == l_Undef) {
        uncheckedEnqueue(the_other, currLevel, ws_bin[k].cref);
      }
    }

    vec<Watcher> &ws = watches[p];
    for (i = j = (Watcher *)ws, end = i + ws.size(); i != end;) {
      const Lit blocker = i->blocker;

      // Try to avoid inspecting the clause:
      if (value(blocker) == l_True) {
        *j++ = *i++;
        continue;
      }

      // Make sure the false literal is data[1]:
      const CRef cr = i->cref;
      Clause &c = ca[cr];

      if (c[0] == false_lit) {
        c[0] = c[1]; 
        c[1] = false_lit;
      }

      i++;
      // If 0th watch is true, then clause is already satisfied.
      const Lit first = c[0];
      Watcher w(cr, first);

      if (first != blocker && value(first) == l_True) {
        *j++ = w;
        continue;
      }

      // Look for new watch:
      const int c_size = c.size();
      for (int k = 2; k < c_size; k++) {
        if (value(c[k]) != l_False) {
          c[1] = c[k];
          c[k] = false_lit;
          watches[~c[1]].push(w);
          goto NextClause;
        }
      }

      // Did not find watch -- clause is unit under assignment:
      *j++ = w;

      if (value(first) == l_False) {
        confl = cr;
        qhead = trail.size();

        // Copy the remaining watches:
        while (i < end) {
          *j++ = *i++;
        }
      } else {
        if (currLevel == decisionLevel()) {
          uncheckedEnqueue(first, currLevel, cr);
        } else {
          int nMaxLevel = currLevel;
          int nMaxInd = 1;

          // pass over all the literals in the clause and find the one with the biggest level
          for (int nInd = 2; nInd < c_size; nInd++) {
            int nLevel = level(var(c[nInd]));

            if (nLevel > nMaxLevel) {
              nMaxLevel = nLevel;
              nMaxInd = nInd;
            }
          }

          if (nMaxInd != 1) {
            std::swap(c[1], c[nMaxInd]);
            j--;
            watches[~c[1]].push(w);
          }

          uncheckedEnqueue(first, nMaxLevel, cr);
        }
      }

    NextClause:;
    }

    ws.shrink(i - j);
  }

  propagations += num_props;
  simpDB_props -= num_props;
  return confl;
}

/*_________________________________________________________________________________________________
|
|  reduceDB : ()  ->  [void]
|
|  Description:
|    Remove half of the learnt clauses, minus the clauses locked by the current assignment. Locked
|    clauses are clauses that are reason to some assignment. Binary clauses are never removed.
|________________________________________________________________________________________________@*/
struct reduceDB_lt {
  ClauseAllocator &ca;
  reduceDB_lt(ClauseAllocator &ca_) : ca(ca_) {}
  bool operator()(CRef x, CRef y) const { return ca[x].activity() < ca[y].activity(); }
};

struct reduceDB_lt2 {
  ClauseAllocator &ca;
  reduceDB_lt2(ClauseAllocator &ca_) : ca(ca_) {}
  bool operator()(CRef x, CRef y) const {
    return ca[y].lbd() < ca[x].lbd();
  }
};

void Solver::reduceDB() {
  int i, j;
  vec<CRef> learntTmp;

  for (i = j = 0; i < learnts.size(); i++) {
    CRef cr = learnts[i];

    if (ca[cr].mark() == LONG) {
      learntTmp.push(cr);
    } else {
      learnts[j++] = cr;
    }
  }

  sort(learntTmp, reduceDB_lt(ca));

  if (learntTmp.size() > 100 && conflicts < 600000 && nFreeVars() > 20000) {
    int b = learntTmp.size() / 16;
    sort((CRef *)learntTmp + 7 * b, b + b / 2, reduceDB_lt2(ca));
  }

  int limit = learntTmp.size() / 2;

  for (i = 0; i < learntTmp.size(); i++) {
    Clause &c = ca[learntTmp[i]];

    if (c.removable() && !locked(c) && i < limit) {
      removeClause(learntTmp[i]);
    } else {
      if (!c.removable()) {
        limit++;
      }

      c.removable(true);
      learnts[j++] = learntTmp[i];
    }
  }

  learnts.shrink(learnts.size() - j);
  checkGarbage();
}

void Solver::reduceDB_Tier2() {
  int i, j;

  for (i = j = 0; i < learnts.size(); i++) {
    Clause &c = ca[learnts[i]];

    if (c.mark() == MIDSZ)
      if (!locked(c) && c.touched() + 30000 < conflicts) {
        c.mark(LONG);
        c.activity() = 0;
        claBumpActivity(c);
      }
  }
}

void Solver::removeSatisfied(vec<CRef> &cs) {
  int i, j;

  for (i = j = 0; i < cs.size(); i++) {
    Clause &c = ca[cs[i]];

    if (satisfied(c)) {
      removeClause(cs[i]);
    } else {
      cs[j++] = cs[i];
    }
  }

  cs.shrink(i - j);
}

void Solver::rebuildOrderHeap() {
  vec<Var> vs;

  for (Var v = 0; v < nVars(); v++)
    if (decision[v] && value(v) == l_Undef) {
      vs.push(v);
    }

  order_heap_CHB.build(vs);
  order_heap_VSIDS.build(vs);
  order_heap_distance.build(vs);
}

/*_________________________________________________________________________________________________
|
|  simplify : [void]  ->  [bool]
|
|  Description:
|    Simplify the clause database according to the current top-level assigment. Currently, the only
|    thing done here is the removal of satisfied clauses, but more things can be put here.
|________________________________________________________________________________________________@*/
bool Solver::simplify() {
  assert(decisionLevel() == 0);

  if (!ok || propagate() != CRef_Undef) {
    return ok = false;
  }

  if (nAssigns() == simpDB_assigns || (simpDB_props > 0)) {
    return true;
  }

  removeSatisfied(learnts);

  if (remove_satisfied) { // Can be turned off.
    removeSatisfied(clauses);
  }

  checkGarbage();
  rebuildOrderHeap();
  simpDB_assigns = nAssigns();
  simpDB_props = clauses_literals + learnts_literals; // (shouldn't depend on stats really, but it will do for now)
  return true;
}

void Solver::rand_rephase() {
  int var_nums = nVars();

  if (rand() % 100 < 50) {
    for (int i = 0; i < var_nums; ++i) {
      polarity[i] = !cur_solution[i];
    }

    return;
  }

  int pick_rand = rand() % 1000;

  if ((pick_rand -= 100) < 0) {
    for (int i = 0; i < var_nums; ++i) {
      polarity[i] = !best_solution[i];
    }
  } else if ((pick_rand -= 300) < 0) {
    for (int i = 0; i < var_nums; ++i) {
      polarity[i] = !cur_solution[i];
    }

    cur_solustion_used = true;
  }
  // top_trail 200
  else if ((pick_rand -= 300) < 0) {
    for (int i = 0; i < var_nums; ++i) {
      polarity[i] = !top_trail_soln[i];
    }
  }
  // reverse
  else if ((pick_rand -= 50) < 0) {
    for (int i = 0; i < var_nums; ++i) {
      polarity[i] = !polarity[i];
    }
  } else if ((pick_rand -= 25) < 0) {
    for (int i = 0; i < var_nums; ++i) {
      polarity[i] = best_solution[i];
    }
  } else if ((pick_rand -= 25) < 0) {
    for (int i = 0; i < var_nums; ++i) {
      polarity[i] = top_trail_soln[i];
    }
  }
  // 150
  else if ((pick_rand -= 140) < 0) {
    for (int i = 0; i < var_nums; ++i) {
      polarity[i] = rand() % 2 == 0 ? 1 : 0;
    }
  } else if ((pick_rand -= 5) < 0) {
    for (int i = 0; i < var_nums; ++i) {
      polarity[i] = 1;
    }
  } else if ((pick_rand -= 5) < 0) {
    for (int i = 0; i < var_nums; ++i) {
      polarity[i] = 0;
    }
  }
  // 50
  else {
    // do nothing
  }
}
/*_________________________________________________________________________________________________
|
|  search : (nof_conflicts : int) (params : const SearchParams&)  ->  [lbool]
|
|  Description:
|    Search for a model the specified number of conflicts.
|
|  Output:
|    'l_True' if a partial assigment that is consistent with respect to the clauseset is found. If
|    all variables are decision variables, this means that the clause set is satisfiable. 'l_False'
|    if the clause set is unsatisfiable. 'l_Undef' if the bound on number of conflicts is reached.
|________________________________________________________________________________________________@*/
lbool Solver::search(int &nof_conflicts) {
  assert(ok);
  int backtrack_level;
  int lbd;
  vec<Lit> learnt_clause;
  bool cached = false;
  starts++;
  freeze_restart_num--;
  bool can_propagate_solution = true;
  unsigned int state_starts_delta = 1500;

  if (starts > state_starts_delta) {
    rand_rephase();
  }

  // simplify
  //
  if ((signed)conflicts >= curSimplify * nbconfbeforesimplify) {
    nbSimplifyAll++;

    if (!simplifyAll()) {
      return l_False;
    }

    curSimplify = (conflicts / nbconfbeforesimplify) + 1;
    nbconfbeforesimplify += incSimplify;
  }

  int chrono_lim = 10000;

  if (VSIDS && nof_conflicts > 1000000) {
    nof_conflicts = 1000000;
  }

  for (;;) {
    CRef confl = propagate();

    if (confl != CRef_Undef) {
      // CONFLICT
      if (VSIDS) {
        if (--timer == 0 && var_decay < 0.95) {
          timer = 5000, var_decay += 0.01;
        }
      } else if (step_size > min_step_size) {
        step_size -= step_size_dec;
      }

      conflicts++;
      nof_conflicts--;

      if (nof_conflicts < -100000) {
        cancelUntil(0);
        return l_Undef;
      }

      ConflictData data = FindConflictLevel(confl);

      if (data.nHighestLevel == 0) {
        return l_False;
      }

      if (data.bOnlyOneLitFromHighest) {
        chrono_lim--;
        cancelUntil(data.nHighestLevel - 1);
        continue;
      }

      learnt_clause.clear();
      DISTANCE = 0;
      analyze(confl, learnt_clause, backtrack_level, lbd);

      // check chrono backtrack condition
      if (chrono_lim > 0 && (confl_to_chrono < 0 || confl_to_chrono <= (int)conflicts) && chrono > -1 &&
          (decisionLevel() - backtrack_level) >= chrono) {
        chrono_lim--;
        ++chrono_backtrack;
        cancelUntil(data.nHighestLevel - 1);
      } else { // default behavior
        ++non_chrono_backtrack;
        cancelUntil(backtrack_level);
      }

      lbd--;
      lbd_sum += lbd;

      if (VSIDS) {
        cached = false;
        conflicts_VSIDS++;
        lbd_queue.push(lbd);
        global_lbd_sum += (lbd > 50 ? 50 : lbd);
        // global_lbd_sum += lbd;
      }

      if (learnt_clause.size() == 1) {
        uncheckedEnqueue(learnt_clause[0]);
      } else {
        CRef cr = ca.alloc(learnt_clause, true);
        ca[cr].set_lbd(lbd);
        // duplicate learnts
        int id = 0;

        if (lbd <= (signed)max_lbd_dup) {
          std::vector<uint32_t> tmp;

          for (int i = 0; i < learnt_clause.size(); i++) {
            tmp.push_back(learnt_clause[i].x);
          }

          id = is_duplicate(tmp);

          if ((unsigned)id == min_number_of_learnts_copies + 1) {
            duplicates_added_conflicts++;
          }

          if ((unsigned)id == min_number_of_learnts_copies) {
            duplicates_added_tier2++;
          }
        }

        // duplicate learnts
        learnts.push(cr);

        if ((lbd <= core_lbd_cut) || (id == (signed)min_number_of_learnts_copies + 1)) {
          ca[cr].mark(SMALL);
        } else if ((lbd <= 6) || (id == (signed)min_number_of_learnts_copies)) {
          ca[cr].mark(MIDSZ);
          ca[cr].touched() = conflicts;
        } else {
          ca[cr].mark(LONG);
          claBumpActivity(ca[cr]);
        }

        attachClause(cr);
        uncheckedEnqueue(learnt_clause[0], backtrack_level, cr);
      }

      if (VSIDS) {
        varDecayActivity();
      }

      claDecayActivity();

      // the top_trail_soln should be update after each conflict
      if (trail.size() > max_trail) {
        max_trail = trail.size();

        for (int i = 0; i < nVars(); ++i) {
          lbool value_i = value(i);

          if (value_i == l_Undef) {
            top_trail_soln[i] = !polarity[i];
          } else {
            top_trail_soln[i] = value_i == l_True ? 1 : 0;
          }
        }
      }
    } else {
      // NO CONFLICT
      if (starts > state_starts_delta) {
        if (can_propagate_solution && freeze_restart_num < 1 && cur_solustion_used &&
            (trail.size() > (int)(0.4 * nVars()) || trail.size() > (int)(0.9 * max_trail))) {
          can_propagate_solution = false;
          cur_solustion_used = false;
          freeze_restart_num = 300; //  int     restarts_gap        = 300;
          propagate_solution();
        }
      }

      bool restart = false;

      if (!VSIDS) {
        restart = nof_conflicts <= 0;
      } else if (!cached) {
        restart = lbd_queue.full() && (lbd_queue.avg() * 0.8 > global_lbd_sum / conflicts_VSIDS);
        cached = true;
      }

      if (restart) {
        lbd_queue.clear();
        cached = false;
        // Reached bound on number of conflicts:
        cancelUntil(0);
        return l_Undef;
      }

      // Simplify the set of problem clauses:
      if (decisionLevel() == 0 && !simplify()) {
        return l_False;
      }

      if (conflicts >= next_T2_reduce) {
        next_T2_reduce = conflicts + 10000;
        reduceDB_Tier2();
      }

      if (conflicts >= next_L_reduce) {
        next_L_reduce = conflicts + 15000;
        reduceDB();
      }

      Lit next = lit_Undef;
      {
        // New variable decision:
        decisions++;
        next = pickBranchLit();

        if (next == lit_Undef)
        // Model found:
        {
          return l_True;
        }
      }
      // Increase decision level and enqueue 'next'
      newDecisionLevel();
      uncheckedEnqueue(next, decisionLevel());
    }
  }
}

/*
  Finite subsequences of the Luby-sequence:

  0: 1
  1: 1 1 2
  2: 1 1 2 1 1 2 4
  3: 1 1 2 1 1 2 4 1 1 2 1 1 2 4 8
  ...


 */
static double luby(double y, int x) {
  // Find the finite subsequence that contains index 'x', and the
  // size of that subsequence:
  int size, seq;

  for (size = 1, seq = 0; size < x + 1; seq++, size = 2 * size + 1)
    ;

  while (size - 1 != x) {
    size = (size - 1) >> 1;
    seq--;
    x = x % size;
  }

  return pow(y, seq);
}

static bool switch_mode = false;

uint32_t Solver::reduceduplicates() {
  uint32_t removed_duplicates = 0;
  std::vector<std::vector<uint64_t>> tmp;

  for (auto &outer_mp : ht) {                // variables
    for (auto &inner_mp : outer_mp.second) { // sizes
      for (auto &in_in_mp : inner_mp.second) {
        if (in_in_mp.second >= 2) {
          // min_number_of_learnts_copies
          tmp.push_back({(uint64_t)outer_mp.first, inner_mp.first, in_in_mp.first, in_in_mp.second});
        }
      }
    }
  }

  removed_duplicates = dupl_db_size - tmp.size();
  ht.clear();

  for (size_t i = 0; i < tmp.size(); i++) {
    ht[tmp[i][0]][tmp[i][1]][tmp[i][2]] = tmp[i][3];
  }

  return removed_duplicates;
}

// NOTE: assumptions passed in member-variable 'assumptions'.
lbool Solver::solve_() {
  model.clear();
  conflict.clear();

  if (!ok) {
    return l_False;
  }

  lbd_sum = 0;
  solves++;
  max_learnts = nClauses() * learntsize_factor;
  learntsize_adjust_confl = learntsize_adjust_start_confl;
  learntsize_adjust_cnt = (int)learntsize_adjust_confl;
  lbool status = l_Undef;
  cur_solution.resize(nVars());
  best_solution.resize(nVars());
  top_trail_soln.resize(nVars());

  add_tmp.clear();
  VSIDS = true;
  int init = 10000;

  while (status == l_Undef && init > 0) {
    status = search(init);
  }

  VSIDS = false;
  duplicates_added_conflicts = 0;
  duplicates_added_minimization = 0;
  duplicates_added_tier2 = 0;
  dupl_db_size = 0;
  size_t dupl_db_size_limit = dupl_db_init_size;
  // Search:
  int curr_restarts = 0;
  uint64_t curr_props = 0;
  uint32_t removed_duplicates = 0;
  //  int smallbd=0;
  int switch_distance = 500; // starts

  while (status == l_Undef && !isTimeout()) {
    if (dupl_db_size >= dupl_db_size_limit) {
      if (starts < 500) {
        removed_duplicates = reduceduplicates();
        dupl_db_size_limit *= 1.1;
        dupl_db_size -= removed_duplicates;
      }
    }

    if (propagations - curr_props > VSIDS_props_limit) {
      curr_props = propagations;
      switch_mode = true;
      VSIDS_props_limit = VSIDS_props_limit + VSIDS_props_limit / 10;
    }

    if (VSIDS) {
      int weighted = INT32_MAX;
      status = search(weighted);
    } else {
      int nof_conflicts = luby(restart_inc, curr_restarts) * restart_first;
      curr_restarts++;
      status = search(nof_conflicts);
    }

    switch_distance--;

    if (switch_distance < 0) {
      switch_distance = 500;

      if (VSIDS) {
        VSIDS = false;
      } else {
        VSIDS = true;
        picked.clear();
        conflicted.clear();
        almost_conflicted.clear();
        canceled.clear();
      }
    }
  }

  if (status == l_True) {
    // Extend & copy model:
    model.growTo(nVars());

    for (int i = 0; i < nVars(); i++) {
      model[i] = value(i);
    }
  } else if (status == l_False && conflict.size() == 0) {
    ok = false;
  }

  cancelUntil(0);
  return status;
}

//=================================================================================================
// Garbage Collection methods:

void Solver::relocAll(ClauseAllocator &to) {
  // All watchers:
  //
  watches.cleanAll();
  watches_bin.cleanAll();

  for (int v = 0; v < nVars(); v++)
    for (int s = 0; s < 2; s++) {
      Lit p = mkLit(v, s);
      vec<Watcher> &ws = watches[p];

      for (int j = 0; j < ws.size(); j++) {
        ca.reloc(ws[j].cref, to);
      }

      vec<Watcher> &ws_bin = watches_bin[p];

      for (int j = 0; j < ws_bin.size(); j++) {
        ca.reloc(ws_bin[j].cref, to);
      }
    }

  // All reasons:
  for (int i = 0; i < trail.size(); i++) {
    Var v = var(trail[i]);

    if (reason(v) != CRef_Undef && (ca[reason(v)].reloced() || locked(ca[reason(v)]))) {
      ca.reloc(vardata[v].reason, to);
    }
  }

  // All learnt:
  for (int i = 0; i < learnts.size(); i++) {
    ca.reloc(learnts[i], to);
  }

  // All original:
  int i, j;

  for (i = j = 0; i < clauses.size(); i++)
    if (ca[clauses[i]].mark() != 1) {
      ca.reloc(clauses[i], to);
      clauses[j++] = clauses[i];
    }

  clauses.shrink(i - j);
}

void Solver::garbageCollect() {
  // Initialize the next region to a size corresponding to the estimated utilization degree. This
  // is not precise but should avoid some unnecessary reallocations for the new region:
  ClauseAllocator to(ca.size() - ca.wasted());
  relocAll(to);

  to.moveTo(ca);
}

void Solver::propagate_solution() {
  int var_nums = nVars();
  int t_sz = trail.size();
  int idx = 0;
  int viewList_sz = t_sz;
  vector<Lit> viewList(var_nums + 2);

  for (int i = 0; i < t_sz; ++i) {
    viewList[i] = trail[i];
  }

  int undef_nums = 0;
  vector<int> undef_vars(var_nums - t_sz + 2);
  vector<int> idx_undef_vars(var_nums + 1, -1); // undef_vars' idx is not -1

  for (int i = 0; i < var_nums; ++i)
    if (value(i) == l_Undef) {
      idx_undef_vars[i] = undef_nums;
      undef_vars[undef_nums++] = i;
    } else {
      cur_solution[i] = (value(i) == l_True) ? 1 : 0;
    }

  while (undef_nums > 0) {
    while (idx < viewList_sz) {
      Lit p = viewList[idx++];
      vec<Watcher> &ws_bin = watches_bin[p];
      int ws_bin_sz = ws_bin.size();

      for (int k = 0; k < ws_bin_sz; k++) {
        Lit the_other = ws_bin[k].blocker;
        Var the_other_var = var(the_other);

        if (idx_undef_vars[the_other_var] > -1) {
          // no conflict and can decide.
          cur_solution[the_other_var] = sign(the_other) ? 0 : 1;
          viewList[viewList_sz++] = the_other;
          int end_var = undef_vars[--undef_nums];
          int idx_end_var = idx_undef_vars[the_other_var];
          undef_vars[idx_end_var] = end_var;
          idx_undef_vars[end_var] = idx_end_var;
          idx_undef_vars[the_other_var] = -1;
        }
      }

      if (undef_nums == 0) {
        break;
      }

      vec<Watcher> &ws = watches[p];
      Watcher *i, *j, *end;

      for (i = j = (Watcher *)ws, end = i + ws.size(); i != end;) {
        // Make sure the false literal is data[1]:
        CRef cr = i->cref;
        Clause &c = ca[cr];
        Lit false_lit = ~p;

        if (c[0] == false_lit) {
          c[0] = c[1], c[1] = false_lit;
        }

        i++;
        // If 0th watch is true, then clause is already satisfied.
        Lit first = c[0];
        Var first_var = var(first);
        Watcher w(cr, first);

        if (idx_undef_vars[first_var] == -1 && cur_solution[first_var] == (!sign(first))) {
          *j++ = w;
          continue;
        }

        int c_sz = c.size();

        for (int k = 2; k < c_sz; ++k) {
          Lit tmp_lit = c[k];
          Var tmp_var = var(tmp_lit);

          if (idx_undef_vars[tmp_var] == -1 && cur_solution[tmp_var] == sign(tmp_lit)) {
          } else {
            c[1] = c[k];
            c[k] = false_lit;
            watches[~c[1]].push(w);
            // next clause
            goto check_next_clause;
          }
        }

        *j++ = w;

        if (idx_undef_vars[first_var] == -1 && cur_solution[first_var] == sign(first)) {
          continue;
        } else {
          // unit can assign
          cur_solution[first_var] = sign(first) ? 0 : 1;
          viewList[viewList_sz++] = first;
          int end_var = undef_vars[--undef_nums];
          int idx_end_var = idx_undef_vars[first_var];
          undef_vars[idx_end_var] = end_var;
          idx_undef_vars[end_var] = idx_end_var;
          idx_undef_vars[first_var] = -1;
        }

      check_next_clause:;
      }

      ws.shrink(i - j);
    }

    if (undef_nums == 0) {
      break;
    }

    int choosevar_idx = rand() % undef_nums;
    Var choosevar = undef_vars[choosevar_idx];
    Lit choose = mkLit(choosevar, polarity[choosevar]);
    cur_solution[choosevar] = sign(choose) ? 0 : 1;
    viewList[viewList_sz++] = choose;
    int end_var = undef_vars[--undef_nums];
    int idx_end_var = idx_undef_vars[choosevar];
    undef_vars[idx_end_var] = end_var;
    idx_undef_vars[end_var] = idx_end_var;
    idx_undef_vars[choosevar] = -1;
  }

  int unsat_num = 0;

  for (int i = 0; i < clauses.size(); i++) {
    Clause &c = ca[clauses[i]];

    for (int j = 0; j < c.size(); j++) {
      Lit p = c[j];

      if (sign(p)) { // negative
        if (!cur_solution[var(p)]) {
          goto nextc;
        }
      } else if (cur_solution[var(p)]) {
        goto nextc;
      }
    }

    unsat_num++;
  nextc:;
  }

  if (unsat_num <= best_unsat_num) {
    for (int i = 0; i < var_nums; ++i) {
      best_solution[i] = cur_solution[i];
    }

    best_unsat_num = unsat_num;
  }
}
