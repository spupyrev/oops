#!/usr/bin/env python3

import itertools
import pathlib
import subprocess
import tempfile

def boundary_orders(size):
    for tail in itertools.permutations(range(1, size)):
        order = (0,) + tail
        if order <= (0,) + tuple(reversed(tail)):
            yield order


def augmented_patch(kind, order):
    if kind == "BALOM":
        order_n = 15
        edges = {(min(i, (i + 1) % 5), max(i, (i + 1) % 5))
                 for i in range(5)}
        edges.update(
            (min(first, second), max(first, second))
            for first, second in (
                (0, 5), (2, 6), (5, 7), (7, 6),
                (1, 8), (5, 9), (7, 10),
                (6, 11), (3, 12), (4, 13)))
        first_boundary, apex = 8, 14
    else:
        order_n = 17
        edges = {(min(i, (i + 1) % 5), max(i, (i + 1) % 5))
                 for i in range(5)}
        edges.update(
            (min(first, second), max(first, second))
            for first, second in (
                (0, 5), (2, 6), (5, 7), (7, 8), (8, 6),
                (1, 9), (5, 10), (7, 11), (8, 12),
                (6, 13), (3, 14), (4, 15)))
        first_boundary, apex = 9, 16

    fixed = set()
    for index in range(len(order)):
        first = first_boundary + order[index]
        second = first_boundary + order[(index + 1) % len(order)]
        edges.add((min(first, second), max(first, second)))
        fixed.add(tuple(sorted((first, second))))
        edges.add((first_boundary + index, apex))
        fixed.add((first_boundary + index, apex))
    return order_n, edges, sorted(fixed)


def graph6(order, edges):
    values = []
    value = 0
    bits = 0
    for second in range(1, order):
        for first in range(second):
            value = (value << 1) | ((first, second) in edges)
            bits += 1
            if bits == 6:
                values.append(value + 63)
                value = 0
                bits = 0
    if bits:
        values.append((value << (6 - bits)) + 63)
    return bytes([order + 63] + values) + b"\n"


def is_expandable(oops, kind, order):
    graph_order, edges, fixed = augmented_patch(kind, order)
    literals = ";".join(f"{first}:{second}-" for first, second in fixed)
    with tempfile.NamedTemporaryFile(suffix=".g6") as stream:
        stream.write(graph6(graph_order, edges))
        stream.flush()
        result = subprocess.run(
            (
                str(oops), f"-i={stream.name}", "-colors=0",
                "-unsat=1", f"-fix-cross1={literals}",
            ),
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        )
    if result.returncode != 0:
        raise RuntimeError(result.stdout)
    if "#unknown = 1" in result.stdout:
        raise RuntimeError(result.stdout)
    if "#planar = 1" in result.stdout or "#1-planar = 1" in result.stdout:
        return True
    if "#non-1-planar = 1" in result.stdout:
        return False
    raise RuntimeError(f"OOPS returned no verdict:\n{result.stdout}")


def main():
    oops = pathlib.Path("./oops").resolve()
    expected_bad = {
        "BALOM": set(),
        "BALOOM": {
            "0245136", "0246135", "0264153", "0315426", "0316425",
        },
    }
    for kind, size in (("BALOM", 6), ("BALOOM", 7)):
        actual_bad = {
            "".join(map(str, order))
            for order in boundary_orders(size)
            if not is_expandable(oops, kind, order)
        }
        if actual_bad != expected_bad[kind]:
            raise RuntimeError(
                f"{kind}: expected {sorted(expected_bad[kind])}, "
                f"found {sorted(actual_bad)}")
        print(
            f"{kind}: {len(list(boundary_orders(size))) - len(actual_bad)} "
            f"expandable boundary orders")


if __name__ == "__main__":
    main()
