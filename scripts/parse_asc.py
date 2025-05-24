#!/usr/bin/env python3

file_path = "/home/spupyrev/downloads/14_5_4.asc"

def extract_edges_from_block(lines):
    graph = {}
    for line in lines:
        if ':' not in line or not line.strip()[0].isdigit():
            continue
        node, neighbors = line.split(':')
        node = int(node.strip())
        neighbors = list(map(int, neighbors.strip().split()))
        graph[node] = neighbors

    # Extract unique edges
    edges = set()
    for node, neighbors in graph.items():
        for neighbor in neighbors:
            edge = tuple(sorted((node, neighbor)))
            edges.add(edge)
    return sorted(edges)


def parse_graphs_from_file(filename):
    graphs = []
    current_block = []

    with open(filename, 'r') as f:
        for line in f:
            stripped = line.strip()
            if stripped.startswith("Taillenweite:"):
                # Process the current graph block
                edges = extract_edges_from_block(current_block)
                #print(edges)
                graphs.append(edges)
                current_block = []
            else:
                current_block.append(line)

    return graphs


# --- Main ---
all_graphs = parse_graphs_from_file(file_path)

for i, edges in enumerate(all_graphs, start=1):
    print("graph k{} {{".format(i))
    mx = 0
    for u, v in edges:
        mx = max(mx, max(u, v))
    for i in range(mx):
        print(f"  {i};")
    for u, v in edges:
        assert 0 < u and u <= mx
        assert 0 < v and v <= mx
        print(f"  {u-1} -- {v-1};")
    print("}")
    print()
