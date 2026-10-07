reads = [
    "ATTCA",
    "ATTGA",
    "CATTG",
    "CTTAT",
    "GATTG",
    "TATTT",
    "TCATT",
    "TCTTA",
    "TGATT",
    "TTATT",
    "TTCAT",
    "TTCTT",
    "TTGAT",
]

k = 3
kmers = []

for item in reads:
    for kmer in range(len(item) - k + 1):
        kmers.append(item[kmer:kmer+k])

edges = []

for kmer in kmers:
    prefix = kmer[:-1]
    suffix = kmer[1:]
    edges.append((prefix, suffix))

with open("question2.dot", "w") as file:
    file.write("digraph debruin {\n")

    for prefix, suffix in edges:
        file.write(f'   "{prefix}" -> "{suffix}";\n')

    file.write("}\n")


unique_edges = list(set(edges))

graph = {}

for prefix, suffix in unique_edges:
    if prefix not in graph:
        graph[prefix] = []
    graph[prefix].append(suffix)

def find_walk(current_node, path, edge_counts):
    if all(count > 0 for count in edge_counts.values()):
        return path
    
    for next_node in graph.get(current_node, []):
        edge = (current_node, next_node)

        if edge_counts[edge]<4:
            edge_counts[edge] += 1
            result = find_walk(
                next_node,
                path + [next_node],
                edge_counts
            )

            if result is not None:
                return result

            edge_counts[edge] -= 1

    return None


edge_counts = {edge: 0 for edge in unique_edges}

path = find_walk("AT", ["AT"], edge_counts)

if path is not None:
    genome = path[0]

    for node in path[1:]:
        genome += node [-1]
    
    print("Possible genome:", genome)

else:
    print("No possible genome found")