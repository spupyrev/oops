#!/usr/bin/env python3
"""oops_split.py — deterministic non-1-planarity proof via cross1 case-split.

Strategy: pure cross1 lex-leader case-split under Aut(G). Deterministic case
analysis, not a heuristic search.

Sound for UNSAT direction. Termination is bounded a priori by the recursion
tree, which is determined entirely by Aut(G): each A-child fixes one more edge,
so a branch bottoms out once every edge is fixed. -total-timeout is the
wallclock backstop for hard/rigid instances the method cannot close.

Proof outline (sound by drawing-level Aut argument + encoder faithfulness):

  ROOT: G {nic|ic|1-planar} ∧ -Cno-transitive  =? UNSAT.
    1. Anchor: Case B (∀e: cross1(e)=true) UNSAT (1s) ⇒ some edge uncrossed.
    2. WLOG cross1(rep_i)=false for some edge orbit O_i  (exhaustive split
       over edge orbits at level 0).
    3. Recurse: at every node with "fixed-false" edges N, compute orbits of
       (E \\ N \\ P) under Stab(N), where P is "fixed-true". For each orbit O_i,
       lex-leader case: cumulative-positives of earlier orbits + cross1(rep_i)
       fixed false. Plus B-case (all remaining = true).
    4. Each leaf closes via CDCL, or the branch bottoms out with all edges fixed.

Requirements:
  - The `oops` binary must exist at the repository root (built via `make`).
  - networkx (for automorphism computation): `pip install networkx`.
  - Inputs are .g6/.s6 (one graph per line) or .cfg multi-graph files;
    `-part` selects one graph by index.

Per-node CDCL logs are written under the system temp dir in
oops_split_<part>_<pid>/.

Usage:
  ./oops_split.py -i=graphs.g6 -part=0
  ./oops_split.py -i=graphs.cfg -part=0 -node-timeout=60 -total-timeout=3600
"""
import argparse
import os
import subprocess
import sys
import tempfile
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent  # this script lives in <repo>/scripts/
try:
    import networkx as nx
except ImportError:
    sys.exit("error: networkx is required (pip install networkx)")

OOPS = str(REPO / "oops")


def parse_graph(path, part):
    """Read the part-th graph from the input file. Returns (n_vertices, edges).

    Dispatches on extension: .g6/.s6 (one graph per line, via networkx) or a
    .cfg multi-graph file. Vertices are labelled 0..n-1, matching oops, so the
    endpoint-labelled -fix-cross1 units the driver emits are unambiguous.
    """
    ext = Path(path).suffix.lower()
    if ext in (".g6", ".s6"):
        with open(path, "rb") as f:
            lines = [ln.strip() for ln in f if ln.strip()]
        if part >= len(lines):
            raise IndexError(f"part={part} but file has {len(lines)} graphs")
        reader = nx.from_graph6_bytes if ext == ".g6" else nx.from_sparse6_bytes
        G = reader(lines[part])
        edges = sorted((min(u, v), max(u, v)) for u, v in G.edges() if u != v)
        return G.number_of_nodes(), edges
    return parse_cfg_part(path, part)


def parse_cfg_part(cfg_path, part):
    """Read the part-th graph from a cfg file. Returns (n_vertices, edges)."""
    with open(cfg_path) as f:
        lines = f.readlines()
    blocks = []
    cur = None
    for ln in lines:
        s = ln.strip()
        if s.startswith("graph ") or s.startswith("digraph "):
            cur = []
            blocks.append(cur)
        elif s == "}" and cur is not None:
            cur = None
        elif cur is not None:
            cur.append(s)
    if part >= len(blocks):
        raise IndexError(f"part={part} but cfg has {len(blocks)} graphs")
    block = blocks[part]
    name_to_id = {}
    edges = []
    for s in block:
        if not s or s == "{" or ";" not in s:
            continue
        body = s.rstrip(";").strip()
        if "--" not in body and "->" not in body:
            tok = body.split()[0]
            if tok in ("node", "edge"):
                continue
            if tok not in name_to_id:
                name_to_id[tok] = len(name_to_id)
        else:
            sep = "--" if "--" in body else "->"
            u, v = (t.strip() for t in body.split(sep, 1))
            assert u in name_to_id, f"unknown node {u!r}"
            assert v in name_to_id, f"unknown node {v!r}"
            iu, iv = sorted((name_to_id[u], name_to_id[v]))
            if iu != iv:
                edges.append((iu, iv))
    return len(name_to_id), edges

# Global semaphore: caps the number of concurrent oops invocations. Initialized
# in main() once we know -threads.
_oops_sem = None

# Global wallclock deadline (epoch seconds). 0 means no limit. Set in main()
# from -total-timeout; consulted at every decision point in `prove`.
_deadline = 0


def _budget_exceeded():
    return _deadline > 0 and time.time() >= _deadline


# Global abort signal. The method only ever proves UNSAT (every leaf UNSAT).
# The moment any node reaches a terminal non-UNSAT verdict (SAT / CRASH), the
# overall answer is decided and no further search can change it, so we stop
# spawning and splitting. Sound for the UNSAT direction: aborting never
# converts a genuine result into a false UNSAT, since declaring UNSAT still
# requires every leaf to have returned UNSAT.
#
# `_abort_verdict` records the deciding verdict. This MUST be consulted before
# concluding UNSAT: an early break out of an as_completed loop can leave the
# deciding SAT/CRASH child unrecorded in a `results` dict, so a caller that
# aggregates only `results` could wrongly fall through to UNSAT. SAT is
# definitive (a drawing exists) so it wins over CRASH if both occur.
_abort = threading.Event()
_abort_verdict = None
_abort_lock = threading.Lock()


def _signal_abort(verdict):
    global _abort_verdict
    with _abort_lock:
        if _abort_verdict is None or verdict == "SAT":
            _abort_verdict = verdict
    _abort.set()


# ---------------------------------------------------------------------------
# Aut and orbits
# ---------------------------------------------------------------------------

class AutTimeout(Exception):
    """Raised when compute_aut exceeds its time budget."""


def compute_aut(n_v, edges, timeout=120):
    """Return list of automorphisms (each a tuple sigma[v] of length n_v).

    networkx's VF2 isomorphism enumeration is worst-case exponential and can
    hang on larger / less-symmetric graphs. Guard it with a SIGALRM timeout so
    a pathological graph raises AutTimeout instead of stalling indefinitely.
    Must run in the main thread (signal.alarm is main-thread only); compute_aut
    is only called from main() before any worker threads start.
    """
    import signal

    def _timeout(signum, frame):
        raise AutTimeout(f"compute_aut exceeded {timeout}s")

    old = signal.signal(signal.SIGALRM, _timeout)
    signal.alarm(int(timeout))
    try:
        G = nx.Graph()
        G.add_nodes_from(range(n_v))
        G.add_edges_from(edges)
        gm = nx.algorithms.isomorphism.GraphMatcher(G, G)
        return [tuple(it[v] for v in range(n_v))
                for it in gm.isomorphisms_iter()]
    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, old)


def stab(auts, edges, fixed_edge_idxs):
    """Subgroup that fixes each edge in `fixed_edge_idxs` setwise."""
    fixed_sets = [frozenset(edges[i]) for i in fixed_edge_idxs]
    return [s for s in auts
            if all(frozenset({s[v] for v in fs}) == fs for fs in fixed_sets)]


def edge_orbits(auts, edges, exclude):
    """Orbits of (edges \\ exclude) under `auts`. Returns list of sorted lists,
    sorted lex by representative."""
    edge_to_idx = {(min(u, v), max(u, v)): i for i, (u, v) in enumerate(edges)}
    excluded = set(exclude)
    remaining = set(range(len(edges))) - excluded
    orbs = []
    while remaining:
        seed = min(remaining)
        orb = {seed}
        frontier = {seed}
        while frontier:
            nxt = set()
            for e in frontier:
                u, v = edges[e]
                for sigma in auts:
                    a, b = sigma[u], sigma[v]
                    eim = edge_to_idx[(min(a, b), max(a, b))]
                    if eim not in orb:
                        nxt.add(eim)
                        orb.add(eim)
            frontier = nxt
        orbs.append(sorted(orb))
        remaining -= orb
    orbs.sort(key=lambda O: O[0])
    return orbs


# ---------------------------------------------------------------------------
# CDCL invocation
# ---------------------------------------------------------------------------

def run_cdcl(part, cfg, edges, fix_neg, fix_pos, budget, log_path,
             encoding="plain", no_transitive=True):
    """Invoke oops with the given cross1 fixings. Returns verdict.

    Caps concurrent oops processes via the global semaphore.
    `edges` is the edge list (index -> (u, v)); fixings are edge indices and
    are passed to oops as endpoint-labelled units 'u:v+/-', so we never depend
    on oops's internal edge ordering.
    `encoding` ∈ {"nic", "ic", "plain"} → adds -nic, -ic, or nothing
    (plain 1-planarity).
    `no_transitive` → adds -Cno-transitive (sound for UNSAT only; on SAT
    branches the encoder asserts and we report CRASH).
    """
    def tok(e, s):
        u, v = edges[e]
        return f"{u}:{v}{s}"
    fix_str = ";".join([tok(e, "+") for e in sorted(fix_pos)] +
                       [tok(e, "-") for e in sorted(fix_neg)])
    cmd = [
        OOPS, f"-i={cfg}", "-unsat=1", f"-part={part}",
        f"-timeout={budget}", "-verbose=1", "-colors=false",
    ]
    if encoding != "plain":
        cmd.append(f"-{encoding}")
    if no_transitive:
        cmd.append("-Cno-transitive")
    if fix_str:
        cmd.append(f"-fix-cross1={fix_str}")
    with _oops_sem:
        # Re-check the deadline AFTER acquiring the slot. The semaphore wait
        # can be long (many subproblems queued); we don't want to launch a
        # fresh CDCL run after the budget has already been exceeded.
        if _budget_exceeded():
            return ("TIMEOUT", 0.0)
        if _abort.is_set():
            return ("ABORTED", 0.0)
        # `t0` is set after acquiring the semaphore so `elapsed` reflects only
        # the CDCL run, not time spent queued waiting for a free slot.
        t0 = time.time()
        with open(log_path, "wb") as f:
            # oops honors its own -timeout=budget. The Python-side backstop
            # (budget*2+60) guards against oops ignoring its own deadline.
            # We poll instead of blocking so that a global abort (some other
            # node reached a terminal verdict) can kill this oops promptly
            # rather than waiting out its full node-timeout.
            proc = subprocess.Popen(cmd, stdout=f, stderr=subprocess.STDOUT)
            hard_deadline = t0 + budget * 2 + 60
            while True:
                try:
                    proc.wait(timeout=1)
                    break
                except subprocess.TimeoutExpired:
                    if _abort.is_set():
                        proc.kill()
                        proc.wait()
                        return ("ABORTED", time.time() - t0)
                    if time.time() >= hard_deadline:
                        proc.kill()
                        proc.wait()
                        raise
        elapsed = time.time() - t0
    log = log_path.read_text(errors="replace")
    if "non-1-planar = 1" in log:
        return ("UNSAT", elapsed)
    if "1-planar = 1" in log and "non-1-planar = 0" in log:
        return ("SAT", elapsed)
    # A planar graph is trivially 1-planar, i.e. satisfiable. oops takes an
    # early-exit and prints '#planar = 1' without running the SAT encoding.
    if "planar = 1" in log and "1-planar = 0" in log and "non-1-planar = 0" in log:
        return ("SAT", elapsed)
    if "unknown = 1" in log:
        return ("TIMEOUT", elapsed)
    # Assertion failure (e.g., -Cno-transitive on a SAT branch produces
    # 'result.order[index].size() <= 1 failed'). Don't infer SAT or UNSAT —
    # the encoding is unsound for that branch; the caller should retry without
    # the offending flag.
    if "assertion " in log and "failed" in log:
        return ("CRASH", elapsed)
    # No verdict line and no detected crash: most likely Python watchdog
    # killed oops before it could write a verdict. Treat as TIMEOUT.
    return ("TIMEOUT", elapsed)


# ---------------------------------------------------------------------------
# Recursion
# ---------------------------------------------------------------------------

class Node:
    """A subproblem: F ∧ (cross1(e)=false for e in fix_neg) ∧
       (cross1(e)=true for e in fix_pos)."""
    __slots__ = ("fix_neg", "fix_pos", "depth", "label")

    def __init__(self, fix_neg, fix_pos, depth, label):
        self.fix_neg = list(fix_neg)
        self.fix_pos = list(fix_pos)
        self.depth = depth
        self.label = label


def lex_leader_children(node, orbits):
    """Generate B-case + A_0..A_{k-1} children (lex-leader case-split).

    B-case: all remaining edges crossed. A_i: representative of orbit i is
    uncrossed, and every earlier orbit's edges are fixed crossed. The
    accumulated positives make the children partition the space and hand each
    child a more constrained (faster-to-close) formula.
    """
    children = []
    # B-case: all remaining edges = positive.
    all_remaining = sorted({e for O in orbits for e in O})
    children.append(Node(
        fix_neg=node.fix_neg,
        fix_pos=node.fix_pos + all_remaining,
        depth=node.depth + 1,
        label=f"{node.label}.B",
    ))
    # A_i cases.
    cumulative_pos = list(node.fix_pos)
    for i, O in enumerate(orbits):
        rep = O[0]
        children.append(Node(
            fix_neg=node.fix_neg + [rep],
            fix_pos=list(cumulative_pos),
            depth=node.depth + 1,
            label=f"{node.label}.A{i}",
        ))
        cumulative_pos.extend(O)
    return children


def prove(node, *, part, cfg, edges, auts, node_timeout,
          log_dir, encoding="plain", no_transitive=True, indent=""):
    """Recursively prove `node` is UNSAT.
    Returns "UNSAT" | "SAT" | "CRASH" | "INCONCLUSIVE" | "ABORTED"."""
    if _abort.is_set():
        return "ABORTED"
    if _budget_exceeded():
        print(f"{indent}[{time.strftime('%H:%M:%S')}] {node.label}: "
              f"BUDGET — total-timeout exceeded before launch", flush=True)
        return "INCONCLUSIVE"
    log_path = log_dir / f"{node.label}.log"
    print(f"{indent}[{time.strftime('%H:%M:%S')}] start {node.label} "
          f"(d={node.depth}, |neg|={len(node.fix_neg)}, |pos|={len(node.fix_pos)})",
          flush=True)

    # Direct CDCL attempt.
    verdict, rt = run_cdcl(part, cfg, edges, node.fix_neg, node.fix_pos,
                           node_timeout, log_path, encoding=encoding,
                           no_transitive=no_transitive)
    print(f"{indent}[{time.strftime('%H:%M:%S')}] {node.label}: "
          f"{verdict} in {rt:.1f}s", flush=True)
    if verdict == "UNSAT":
        return verdict
    # Terminal non-UNSAT verdicts decide the whole run: no amount of further
    # splitting can turn the answer back into UNSAT. Signal a global abort so
    # sibling/other subtrees stop instead of fanning out further.
    if verdict in ("SAT", "CRASH"):
        _signal_abort(verdict)
        return verdict
    # This node's oops was killed because another node already decided the run.
    if verdict == "ABORTED":
        return "ABORTED"
    # TIMEOUT: split. The recursion is naturally bounded — each A-child adds one
    # edge to fix_neg, so once every edge is fixed there are no orbits left
    # (handled below). -total-timeout is the wallclock backstop for hard/rigid
    # instances the method cannot close.
    if _abort.is_set():
        return "ABORTED"
    if _budget_exceeded():
        print(f"{indent}  BUDGET — total-timeout exceeded before split", flush=True)
        return "INCONCLUSIVE"

    # Compute orbits of the remaining edges to split on, under the
    # sub-automorphism group that fixes fix_neg.
    sub_auts = stab(auts, edges, node.fix_neg)
    orbits = edge_orbits(sub_auts, edges, exclude=node.fix_neg + node.fix_pos)
    if not orbits:
        # No more edges to split; CDCL just couldn't close in budget.
        # Caller may want to retry with bigger budget.
        print(f"{indent}  no orbits left; CDCL INCONCLUSIVE", flush=True)
        return "INCONCLUSIVE"

    children = lex_leader_children(node, orbits)
    print(f"{indent}  |Stab|={len(sub_auts)}; orbits={[len(O) for O in orbits]}; "
          f"{len(children)} children", flush=True)

    # Run children in parallel. We don't cap workers here — the global
    # _oops_sem limits actual oops concurrency. The "extra" workers are just
    # idle Python threads waiting for the semaphore, which is cheap.
    results = {}
    with ThreadPoolExecutor(max_workers=len(children)) as ex:
        future_to_child = {
            ex.submit(prove, c, part=part, cfg=cfg, edges=edges, auts=auts,
                      node_timeout=node_timeout, log_dir=log_dir,
                      encoding=encoding, no_transitive=no_transitive,
                      indent=indent + "  "): c
            for c in children
        }
        for fut in as_completed(future_to_child):
            c = future_to_child[fut]
            results[c.label] = fut.result()
            # A terminal non-UNSAT verdict anywhere aborts the whole run; stop
            # waiting on the remaining siblings (already-running oops calls will
            # finish quickly and their children short-circuit on the flag).
            if _abort.is_set():
                break

    vals = results.values()
    if any(v == "SAT" for v in vals):
        return "SAT"
    if any(v == "CRASH" for v in vals):
        return "CRASH"
    if _abort.is_set() or any(v == "ABORTED" for v in vals):
        return "ABORTED"
    if any(v == "INCONCLUSIVE" for v in vals):
        return "INCONCLUSIVE"
    return "UNSAT"


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    # Input selection, single-dash to mirror oops (-i=..., -part=...).
    ap.add_argument("-i", dest="cfg", default=str(REPO / "benchmarks/wiki_nonplanar.cfg"),
                    help="input graph (.g6/.s6/.cfg); relative paths resolve "
                    "against the repo root")
    ap.add_argument("-part", type=int, default=0,
                    help="index of the graph to select from the input file")
    ap.add_argument("-node-timeout", type=int, default=300,
                    help="per-node CDCL budget in seconds")
    ap.add_argument("-total-timeout", type=int, default=0,
                    help="overall wallclock budget in seconds (0 = no limit). "
                    "When exceeded, no new CDCL invocations are started; "
                    "in-flight ones finish and the run reports INCONCLUSIVE.")
    ap.add_argument("-threads", type=int,
                    default=int(os.environ.get("THREADS", os.cpu_count() or 1)),
                    help="max concurrent oops invocations (default: CPU count)")
    # Encoding, mirroring oops: plain 1-planarity by default; -ic / -nic
    # select the IC / NIC variants.
    ap.add_argument("-ic", dest="encoding", action="store_const", const="ic",
                    default="plain", help="enforce IC-1-planarity")
    ap.add_argument("-nic", dest="encoding", action="store_const", const="nic",
                    help="enforce NIC-1-planarity")
    # Direction, mirroring oops: -unsat is the default (adds -Cno-transitive,
    # fast and sound for the UNSAT/non-1-planarity direction). -sat drops
    # -Cno-transitive so satisfiable branches don't assert (slower).
    ap.add_argument("-unsat", dest="no_transitive", action="store_true",
                    default=True, help="UNSAT direction (default): "
                    "adds -Cno-transitive, fast, non-1-planarity proofs only")
    ap.add_argument("-sat", dest="no_transitive", action="store_false",
                    help="omit -Cno-transitive (slower, but sound for "
                    "satisfiable branches too)")
    args = ap.parse_args()

    # Ctrl-C handling. The main thread otherwise blocks in ThreadPoolExecutor
    # shutdown (joining workers) and cannot raise KeyboardInterrupt until the
    # pool drains — so a plain Ctrl-C appears dead. Instead, route SIGINT into
    # the existing abort path: workers poll _abort every 1s in run_cdcl, kill
    # their oops child, and return, so the pool drains within ~1s and the run
    # stops gracefully. A second Ctrl-C restores the default handler for a hard
    # exit if the first is not enough.
    import signal

    def _on_sigint(signum, frame):
        print("\n[Interrupt] Ctrl-C — aborting; workers stopping...", flush=True)
        _abort.set()
        signal.signal(signal.SIGINT, signal.SIG_DFL)
    signal.signal(signal.SIGINT, _on_sigint)

    # The recursion nests one blocked thread per open subproblem; a large tree
    # can create tens of thousands of them. Shrink the per-thread stack from the
    # ~8MB default so they fit in address space (real oops work runs in
    # subprocesses, so the Python stacks stay shallow).
    threading.stack_size(512 * 1024)

    global _oops_sem, _deadline
    _oops_sem = threading.Semaphore(args.threads)
    print(f"# threads: {args.threads} (max concurrent oops invocations)", flush=True)

    t_total = time.time()
    if args.total_timeout > 0:
        _deadline = t_total + args.total_timeout
        print(f"# total-timeout: {args.total_timeout}s "
              f"(deadline at {time.strftime('%H:%M:%S', time.localtime(_deadline))})",
              flush=True)
    # Pre-flight: the oops binary must exist and actually run. Without this,
    # a missing binary or a runtime loader error (e.g. a libstdc++ version
    # mismatch) would produce no verdict line and be silently misread as a
    # per-node TIMEOUT, yielding a bogus INCONCLUSIVE instead of a clear error.
    if not Path(OOPS).exists():
        sys.exit(f"error: oops binary not found at {OOPS} (build it with `make`)")
    try:
        probe = subprocess.run([OOPS, "-help"], capture_output=True, timeout=30)
    except Exception as e:
        sys.exit(f"error: could not run {OOPS}: {e}")
    if probe.returncode != 0:
        msg = probe.stderr.decode(errors="replace").strip().splitlines()
        hint = ("\n  " + msg[0]) if msg else ""
        sys.exit(f"error: {OOPS} failed to run (exit {probe.returncode}).{hint}")

    n_v, edges = parse_graph(args.cfg, args.part)
    print(f"# graph: |V|={n_v}, |E|={len(edges)}", flush=True)

    # Pre-flight: Aut + edge orbits.
    t0 = time.time()
    try:
        auts = compute_aut(n_v, edges)
    except AutTimeout as e:
        print(f"\nAUT-TIMEOUT: {e}. The automorphism computation (networkx VF2)")
        print(f"  did not finish; the method cannot run on this graph.")
        sys.exit(5)
    print(f"# |Aut|={len(auts)} (computed in {time.time()-t0:.2f}s)", flush=True)
    full_orbits = edge_orbits(auts, edges, exclude=[])
    print(f"# edge orbits under Aut: {len(full_orbits)}; "
          f"sizes={[len(o) for o in full_orbits]}", flush=True)

    log_dir = Path(tempfile.gettempdir()) / f"oops_split_{args.part}_{os.getpid()}"
    log_dir.mkdir(parents=True, exist_ok=True)
    print(f"# logs: {log_dir}", flush=True)

    # Anchor: Case B = all edges crossed.
    # In a 1-planar drawing each edge is crossed at most once and each crossing
    # pairs exactly two edges, so "all edges crossed" needs a perfect matching
    # of E into crossing pairs — impossible when |E| is odd. So for odd |E|
    # Case B is UNSAT by parity, unconditionally (no encoder-faithfulness
    # assumption needed). For even |E| we still refute it via CDCL, always with
    # -Cno-transitive (Case B is a UNSAT proof, and -Cno-transitive is sound for
    # UNSAT).
    print(f"\n[Anchor] Case B: ∀e: cross1(e)=true ⇒ UNSAT", flush=True)
    if len(edges) % 2 == 1:
        print(f"  Case B: UNSAT by parity (|E|={len(edges)} is odd; no perfect "
              f"matching of edges into crossing pairs)", flush=True)
    else:
        verdict, rt = run_cdcl(args.part, args.cfg, edges, [], list(range(len(edges))),
                               args.node_timeout, log_dir / "anchor_B.log",
                               encoding=args.encoding,
                               no_transitive=True)
        print(f"  Case B: {verdict} in {rt:.1f}s", flush=True)
        if verdict != "UNSAT":
            print(f"\n  Case B not UNSAT — skewness might equal |E|.")
            print(f"  cross1(e)=false anchor not sound; aborting.")
            sys.exit(2)

    # Root: case-split over edge orbits under full Aut.
    # For edge-transitive G this is just one A-case (cross1(0)=false) plus B.
    # For non-edge-transitive G this enumerates one anchor per orbit.
    print(f"\n[Root] case-split on edge orbits under Aut", flush=True)
    root = Node(fix_neg=[], fix_pos=[], depth=0, label="R")
    children = lex_leader_children(root, full_orbits)
    # Skip the B-case at root (== anchor we just verified).
    children = [c for c in children if not c.label.endswith(".B")]
    print(f"  {len(children)} root children (one per edge orbit)", flush=True)

    t0 = time.time()
    results = {}
    with ThreadPoolExecutor(max_workers=len(children)) as ex:
        future_to_child = {
            ex.submit(prove, c, part=args.part, cfg=args.cfg, edges=edges,
                      auts=auts, node_timeout=args.node_timeout,
                      log_dir=log_dir, encoding=args.encoding,
                      no_transitive=args.no_transitive,
                      indent="  "): c
            for c in children
        }
        for fut in as_completed(future_to_child):
            c = future_to_child[fut]
            results[c.label] = fut.result()
            if _abort.is_set():
                break
    elapsed = time.time() - t0

    print(f"\n[Result] root children:", flush=True)
    for label, v in sorted(results.items()):
        print(f"  {label}: {v}", flush=True)

    print(f"\nTotal recursion wallclock: {elapsed/60:.1f} min", flush=True)
    print(f"Total wallclock (incl. anchor + Aut): {(time.time()-t_total)/60:.1f} min", flush=True)

    enc = args.encoding.upper()
    vals = results.values()
    # Consult the recorded deciding verdict FIRST. When the run aborted, an
    # early break can leave the deciding SAT/CRASH child absent from `results`,
    # so aggregating `results` alone could wrongly fall through to UNSAT. The
    # abort verdict is always recorded by the node that triggered it.
    if _abort_verdict == "SAT" or any(v == "SAT" for v in vals):
        print(f"\nSAT: graph (part={args.part}) IS {enc}-1-planar (a drawing exists).")
        print(f"  The case-split method proves NON-1-planarity; on a 1-planar")
        print(f"  graph it cannot succeed. Run oops directly to get the drawing:")
        print(f"    oops -i={args.cfg} -part={args.part} -o=out.svg")
        sys.exit(0)
    if _abort_verdict == "CRASH" or any(v == "CRASH" for v in vals):
        print(f"\nCRASH: oops asserted on a branch (expected when -Cno-transitive")
        print(f"  meets a satisfiable sub-formula). This strongly indicates the")
        print(f"  graph IS {enc}-1-planar, so the non-1-planarity proof cannot close.")
        print(f"  Confirm with a direct run:  oops -i={args.cfg} -part={args.part} -o=out.svg")
        sys.exit(4)
    # ABORTED can only reach here without a SAT/CRASH (both handled above) when
    # the -total-timeout fired, so it means the same as INCONCLUSIVE: no verdict.
    if _abort.is_set() or any(v in ("INCONCLUSIVE", "ABORTED") for v in vals):
        print(f"\nINCONCLUSIVE: a node timed out and could not be split further,")
        print(f"  or the -total-timeout budget was exhausted.")
        print(f"  Suggest: increase -node-timeout or -total-timeout. Note the")
        print(f"  method needs a rich Aut(G); it degenerates on rigid graphs (|Aut|=1).")
        sys.exit(2)
    print(f"\nUNSAT: graph (part={args.part}) is NOT {enc}-1-planar")
    sys.exit(1)


if __name__ == "__main__":
    main()
