#!/usr/bin/env python3
"""Fast, independent, all-flag exact verification for Propositions 3.1 and 3.2.

This checks every constraint, objective (3.2), positivity vector (3.2),
multipliers, and positive semidefiniteness, using exact Python Fractions.
Unlike the general verifier, here flag-isomorphism classes are encoded by
small bit masks, so all 818 six-vertex constraints can be checked promptly.
No SageMath installation is required. Trust the input pickle.
"""
import argparse
from collections import Counter
from fractions import Fraction
from itertools import combinations, permutations
from math import factorial, comb
from pathlib import Path
import pickle

from repair_prop_3_2 import exact_objective, exact_positives

TRIPLES = list(combinations(range(6), 3))
TRI_INDEX = {a: i for i, a in enumerate(TRIPLES)}


def matrix_from_upper(triangle, n):
    result = [[Fraction(0) for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for j in range(i, n):
            # Formula in the original verifier.py.
            idx = i * (2 * n - i + 1) // 2 + (j - i)
            result[i][j] = result[j][i] = triangle[idx]
    if len(triangle) != n * (n + 1) // 2:
        raise ValueError('Malformed triangular matrix')
    return result


def is_psd_exact(M):
    """Exact LDL^T test, including zero pivots."""
    n = len(M)
    L = [[Fraction(0)] * n for _ in range(n)]
    D = [Fraction(0)] * n
    for j in range(n):
        L[j][j] = Fraction(1)
        d = M[j][j] - sum(L[j][k]**2 * D[k] for k in range(j))
        if d < 0:
            return False
        D[j] = d
        for i in range(j + 1, n):
            residual = M[i][j] - sum(L[i][k] * L[j][k] * D[k]
                                           for k in range(j))
            if d == 0:
                if residual != 0:
                    return False
            else:
                L[i][j] = residual / d
    return True


def local_edge_mask(vertices, edges):
    n = len(vertices)
    edge_set = {frozenset(edge) for edge in edges}
    return sum((1 << j) for j, trip in enumerate(combinations(range(n), 3))
               if frozenset(vertices[x] for x in trip) in edge_set)


def canonical_mask(mask, size, type_size):
    if size - type_size == 1:
        return mask
    if size - type_size != 2:
        raise ValueError('Expected one or two unlabelled vertices')
    tris = list(combinations(range(size), 3))
    triindex = {t: i for i, t in enumerate(tris)}
    swap = list(range(size))
    swap[-1], swap[-2] = swap[-2], swap[-1]
    other = 0
    for i, trip in enumerate(tris):
        if mask & (1 << i):
            other |= (1 << triindex[tuple(sorted(swap[x] for x in trip))])
    return min(mask, other)


def class_lookup(typed_flags):
    """Map each type's exact rooted flag to its canonical small bitmask."""
    all_lookups = []
    for key, flags in typed_flags.items():
        t = len(flags[0][1])
        n = flags[0][0]
        if (t, n) not in [(2, 4), (4, 5)]:
            raise ValueError(f'Unexpected flag type: ({t}, {n})')
        lookup = {}
        for j, flag in enumerate(flags):
            if flag[0] != n or len(flag[1]) != t:
                raise ValueError('Inconsistent group dimensions')
            roots = list(flag[1])
            unlabeled = sorted(set(range(n)) - set(roots))
            local = roots + unlabeled
            mask = canonical_mask(local_edge_mask(local, flag[2][0][1]), n, t)
            if mask in lookup:
                raise ValueError(f'Duplicate typed flags in group {key}: {mask}')
            lookup[mask] = j
        all_lookups.append((t, n, lookup))
    return all_lookups


def precompute_assignments(t, flag_size):
    """All ordered type embeddings and disjoint private vertex choices."""
    k = flag_size - t
    local_triples = list(combinations(range(flag_size), 3))
    assigns = []
    for roots in permutations(range(6), t):
        rest = sorted(set(range(6)) - set(roots))
        for private1 in combinations(rest, k):
            rem = tuple(v for v in rest if v not in private1)
            for private2 in combinations(rem, k):
                embeddings = (roots + private1, roots + private2)
                mappings = []
                for local in embeddings:
                    mappings.append(tuple(
                        TRI_INDEX[tuple(sorted(local[a] for a in triple))]
                        for triple in local_triples
                    ))
                assigns.append(tuple(mappings))
    expected = factorial(6) // factorial(6 - t) * comb(6 - t, k) * comb(6 - t - k, k)
    if len(assigns) != expected:
        raise ValueError('Incorrect denominator or count')
    return assigns, expected


def proj(mask, mask_map):
    val = 0
    for j, glob_idx in enumerate(mask_map):
        val |= ((mask >> glob_idx) & 1) << j
    return val


def exact_inner_sum(global_mask, matrices, lookups, assignments):
    total = Fraction(0)
    for X, (t, size, lookup), (actions, denom) in zip(matrices, lookups, assignments):
        counts = Counter()
        for left_map, right_map in actions:
            first = canonical_mask(proj(global_mask, left_map), size, t)
            second = canonical_mask(proj(global_mask, right_map), size, t)
            i = lookup.get(first)
            j = lookup.get(second)
            # If the induced root is not of this type, neither product flag
            # can match; this assignment contributes zero.
            if i is not None and j is not None:
                counts[i, j] += 1
        total += sum((c * X[i][j] for (i, j), c in counts.items()), Fraction(0)) / denom
    return total


def verify(cert, limit=None):
    n = len(cert['base flags'])
    if n != 818 or any(g[0] != 6 for g in cert['base flags']):
        raise ValueError('Expected the 818 six-vertex base flags')
    result = cert['result']
    if any(not isinstance(x, (int, Fraction)) or x < 0 for x in cert['e vector']):
        raise AssertionError('A positivity assumption multiplier is negative or inexact')
    if len(cert['e vector']) != len(cert['positives']):
        raise AssertionError('Assumption/multiplier dimension mismatch')
    if len(cert['slack vector']) != n or len(cert['target']) != n:
        raise AssertionError('Invalid constraint-vector dimensions')
    if any(s < 0 for s in cert['slack vector']):
        raise AssertionError('Negative slack in certificate')
    if len(cert['positives']) == 7:
        assert cert['positives'] == exact_positives(cert), 'Positivity vector does not match exact Prop 3.2'
        assert cert['target'] == exact_objective(cert), 'Objective does not match exact Prop 3.2'
        print('PASS: Every objective and exact positivity coefficient matches.')
    elif len(cert['positives']) != 0:
        raise ValueError('This fast validator supports Props 3.1 and 3.2 only')

    lookups = class_lookup(cert['typed flags'])
    if len(lookups) != len(cert['X matrices']):
        raise AssertionError('Matrix groups mismatch')
    matrices = []
    assignments = []
    for k, (t, fsize, lookup) in enumerate(lookups):
        X = matrix_from_upper(cert['X matrices'][k], len(lookup))
        if not is_psd_exact(X):
            raise AssertionError(f'Matrix {k + 1} not PSD over Q')
        matrices.append(X)
        assignments.append(precompute_assignments(t, fsize))
    print(f'PASS: {len(matrices)} matrices are exactly positive semidefinite.')

    if limit is None:
        indices = range(n)
    else:
        a, b = limit
        if not (1 <= a <= b <= n):
            raise ValueError(f'Invalid range {limit}')
        indices = range(a - 1, b)
    for i in indices:
        graph = cert['base flags'][i]
        global_mask = sum(1 << TRI_INDEX[tuple(edge)] for edge in graph[2][0][1])
        inner = exact_inner_sum(global_mask, matrices, lookups, assignments)
        positive = sum((m * p[i] for m, p in zip(cert['e vector'], cert['positives'])), Fraction(0))
        sign = 1 if cert['maximize'] else -1
        lhs = sign * (result - cert['target'][i])
        rhs = inner + positive + cert['slack vector'][i]
        if lhs != rhs:
            raise AssertionError(f'Constraint {i + 1} fails: LHS={lhs}, RHS={rhs}, delta={lhs-rhs}')
    print(f'PASS: All {len(indices)} exact rational identities hold.')
    print(f'PASS: All {len(cert["e vector"])} multipliers and all {n} slacks are nonnegative.')
    print('CERTIFICATE CHECKED SUCCESSFULLY')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('certificate', type=Path)
    parser.add_argument('--range', nargs=2, type=int, metavar=('START', 'END'))
    args = parser.parse_args()
    with args.certificate.open('rb') as f:
        cert = pickle.load(f)
    verify(cert, args.range)


if __name__ == '__main__':
    main()
