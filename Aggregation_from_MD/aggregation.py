import numpy as np
from collections import defaultdict

def calculate_aggregates(pdb_frame):
    def extract_coordinates_and_resrange(lines):
        coords = []
        resids = set()
        for line in lines:
            if line.startswith(("ATOM", "HETATM")):
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    coords.append([x, y, z])
                    resid = int(line[22:26])
                    resids.add(resid)
                except ValueError:
                    continue
        if coords and resids:
            return np.array(coords), min(resids), max(resids)
        else:
            return np.array([]), None, None

    sugars = []
    res_ranges = []
    current = []

    for line in pdb_frame:
        if line.startswith(("ATOM", "HETATM")):
            current.append(line)
        elif line.startswith(("TER", "ENDMDL")):
            coords, rmin, rmax = extract_coordinates_and_resrange(current)
            if len(coords) > 0:
                sugars.append(coords)
                res_ranges.append((rmin, rmax))
            current = []

    if current:
        coords, rmin, rmax = extract_coordinates_and_resrange(current)
        if len(coords) > 0:
            sugars.append(coords)
            res_ranges.append((rmin, rmax))

    n = len(sugars)
    contact_graph = defaultdict(set)

    for i in range(n):
        for j in range(i + 1, n):
            a1 = sugars[i][:, None, :]
            a2 = sugars[j][None, :, :]
            delta = a1 - a2
            distances = np.linalg.norm(delta, axis=-1)
            if np.any(distances < 3.5):
                contact_graph[i].add(j)
                contact_graph[j].add(i)

    visited = set()
    aggregates = []

    for i in range(n):
        if i not in visited:
            stack = [i]
            group = set()
            while stack:
                node = stack.pop()
                if node not in visited:
                    visited.add(node)
                    group.add(node)
                    stack.extend(contact_graph[node])
            aggregates.append(group)

    true_aggregates = [group for group in aggregates if len(group) >= 2]

    for group in true_aggregates:
        descriptions = []
        for i in sorted(group):
            rmin, rmax = res_ranges[i]
            if rmin == rmax:
                descriptions.append(f"{rmin} (sugar {i})")
            else:
                descriptions.append(f"{rmin}–{rmax} (sugar {i})")
        print("Aggregate:", ", ".join(descriptions))

    if not true_aggregates:
        return 0, 0.0
    else:
        sizes = [len(g) for g in true_aggregates]
        return sum(sizes) / len(sizes), max(sizes)
