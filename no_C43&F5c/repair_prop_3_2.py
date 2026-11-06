"""Exact independent reconstruction of the objectives and degree assumptions.

The certificate data are treated only as encoded unlabelled/typed hypergraphs.
The formal claims being checked are specified in this module.
"""
from fractions import Fraction
from itertools import combinations
from math import comb

BETA = Fraction(4641, 10000)


def _edges(flag):
    return {frozenset(e) for e in flag[2][0][1]}


def exact_objective(cert):
    """Uncoloured edge density for 3.1; targeted rooted pattern density for 3.2."""
    m = len(cert['positives'])
    if m not in (0, 7):
        raise ValueError('Only Propositions 3.1 and 3.2 are supported')
    result = []
    for f in cert['base flags']:
        n = f[0]
        if n != 6:
            raise ValueError('Unexpected base flag size')
        edges = _edges(f)
        if m == 0:
            result.append(Fraction(len(edges), comb(n, 3)))
            continue
        count = 0
        vertices = frozenset(range(n))
        for r in combinations(range(n), 3):
            root = frozenset(r)
            complement = vertices - root
            if root not in edges or complement not in edges:
                continue
            high = sum(
                sum(frozenset((v, a, b)) in edges
                    for a, b in combinations(r, 2)) >= 2
                for v in complement
            )
            if high == 2:
                count += 1
        result.append(Fraction(count, comb(n, 3)))
    return result


def exact_positives(cert):
    """Exact rooted-pattern degree assumption (degree - 4641/10000).

    A factor flag has four vertices, one distinguished (root) and three
    unlabelled. Its isomorphism class is determined by the number of edges
    through the root and whether the opposite triple is an edge.
    """
    groups = cert['factor flags']
    if len(groups) != 1 or len(groups[0]) != 7:
        raise ValueError('Expected seven rooted four-vertex factor flags')
    classification = {}
    for j, f in enumerate(groups[0]):
        if f[0] != 4 or len(f[1]) != 1:
            raise ValueError('Incorrect factor-flag format')
        root = f[1][0]
        edges = _edges(f)
        private = frozenset(set(range(4)) - {root})
        key = (sum(root in e for e in edges), private in edges)
        if key in classification:
            raise ValueError('Repeated factor-flag class')
        classification[key] = j
    if len(classification) != 7:
        raise ValueError('Missing rooted class')
    n = len(cert['base flags'])
    rows = [[Fraction(0) for _ in range(n)] for _ in range(7)]
    for i, f in enumerate(cert['base flags']):
        if f[0] != 6:
            raise ValueError('Expected 6-vertex base flag')
        edges = _edges(f)
        for root in range(6):
            rest = [v for v in range(6) if v != root]
            for left in combinations(rest, 2):
                private = frozenset(v for v in rest if v not in left)
                root_count = sum(frozenset((root, a, b)) in edges
                                 for a, b in combinations(private, 2))
                key = (root_count, private in edges)
                j = classification.get(key)
                if j is None:
                    raise ValueError('Unexpected rooted factor class (forbidden K4?)')
                weight = (int(frozenset((root, *left)) in edges) - BETA) / 60
                rows[j][i] += weight
    return rows
