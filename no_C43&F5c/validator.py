import pickle
from math import comb, factorial
from fractions import Fraction
from itertools import permutations, combinations
import argparse


class flag:
    def __init__(self, size, edges=[], ftype=[], coloured=[]):
        self.size = size
        self.edges = edges
        self.ftype = ftype
        self.coloured = set(coloured)

    def subflag(self, vertices):
        new_edges = [edge for edge in self.edges
                     if all(v in vertices for v in edge)]
        new_coloured = self.coloured & set(vertices)
        relabel = {vertices[i]: i for i in range(len(vertices))}
        new_edges_1 = [[relabel[v] for v in edge] for edge in new_edges]
        new_ftype_1 = [relabel[v] for v in self.ftype]
        new_coloured_1 = [relabel[v] for v in new_coloured]
        return flag(
            size=len(vertices),
            edges=new_edges_1,
            ftype=new_ftype_1,
            coloured=new_coloured_1
        )


# ---------------------------- Auxiliary functions ----------------------------

def is_isomorphism(bijection, flag1, flag2):
    """Checks if a vertex bijection between two flags is an isomorphism.

    """
    # Check type and colour mapping
    if flag2.ftype != [bijection[v] for v in flag1.ftype]:
        return False
    if flag2.coloured != {bijection[v] for v in flag1.coloured}:
        return False
    # Check edge mapping
    mapped_edges = set()
    for edge in flag1.edges:
        mapped_edge = frozenset(bijection[v] for v in edge)
        mapped_edges.add(mapped_edge)
    if mapped_edges != {frozenset(edge) for edge in flag2.edges}:
        return False
    return True


def are_isomorphic(flag1, flag2):
    """Checks if two flags are isomorphic.

    """
    # Check number of edges and coloured vertices
    if len(flag1.edges) != len(flag2.edges):
        return False
    if len(flag1.coloured) != len(flag2.coloured):
        return False
    V = range(flag1.size)
    # Check degree sequence
    if {len([e for e in flag1.edges if v in e]) for v in V} != {len([e for e in flag2.edges if v in e]) for v in V}:
        return False
    # Check every permutation
    for permutation in permutations(V):
        bijection = {v: u for v, u in zip(V, permutation)}
        if is_isomorphism(bijection, flag1, flag2):
            return True
    return False


def change_base(quantum, base_flags):
    """Given a qunatum flag and some base flags of the same size, expresses the
    quantum flag as a linear combination of the base flags. Returns the coefficients
    of this vector.

    """
    coefficients = [Fraction(0) for _ in range(len(base_flags))]
    for coef, small in quantum:
        for idx, big in enumerate(base_flags):
            # Obtain the density of 'small' in 'big'
            copies = 0
            # Check every small-size subset of big for isomorphism
            for U in combinations(range(big.size), small.size):
                if are_isomorphic(big.subflag(U), small):
                    copies += 1
            density = Fraction(copies, comb(big.size, small.size))
            coefficients[idx] += coef * density
    return coefficients


def get_densities(quantum_flag, flags, G):
    """Given a typed quantum flag, a list of typed flags, and a host 0-flag,
    calculates the list of densities of the type-average of the product of the
    quantum flag with each of the flags in the host flag.

    """
    n = G.size
    n_flags = len(flags)
    n_qflag = len(quantum_flag)
    flags_size = flags[0].size
    quantum_size = quantum_flag[0][1].size
    type_size = len(flags[0].ftype)
    V = range(n)
    M = [[0 for _ in range(n_flags)] for _ in range(n_qflag)]
    for ftype in permutations(V, type_size):
        new_flag = flag(G.size, G.edges, list(ftype), G.coloured)
        unlabeled_vertices = [v for v in V if v not in ftype]
        for V1 in combinations(unlabeled_vertices, quantum_size - type_size):
            flag1 = new_flag.subflag(ftype + V1)
            for i, flag_a in enumerate(quantum_flag):
                if are_isomorphic(flag1, flag_a[1]):
                    remaining_vertices = [v for v in unlabeled_vertices if v not in V1]
                    for V2 in combinations(remaining_vertices, flags_size - type_size):
                        flag2 = new_flag.subflag(ftype + V2)
                        for j, flag_b in enumerate(flags):
                            if are_isomorphic(flag2, flag_b):
                                M[i][j] += 1
                                break
                    break
    denominator = comb(n, type_size) * factorial(type_size)
    denominator *= comb(n - type_size, quantum_size - type_size)
    denominator *= comb(n - quantum_size, flags_size - type_size)
    densities = [0 for _ in range(n_flags)]
    for i in range(n_flags):
        numerator = sum([quantum_flag[j][0] * M[j][i] for j in range(n_qflag)])
        densities[i] = Fraction(numerator, denominator)
    return densities


def average(typed_flag):
    """Given a typed flag, calculates the probability that the subjacent 0-flag,
    together with a uniformly at random type is isomorphic to the flag.

    """
    n, edges, coloured = typed_flag.size, typed_flag.edges, typed_flag.coloured
    type_size = len(typed_flag.ftype)
    numerator = 0
    for ftype in permutations(range(n), type_size):
        if are_isomorphic(typed_flag, flag(n, edges, list(ftype), coloured)):
            numerator += 1
    return Fraction(numerator, comb(n, type_size) * factorial(type_size))


def recover_from_certificate(file_name):
    """Given a certificate, recovers the stored objects.

    """
    with open(file_name, 'rb') as handle:
        certificate = pickle.load(handle)
        base_flags_0 = certificate['base flags']
        target = certificate['target']
        positives = certificate['positives']
        factor_flags_0 = certificate['factor flags']

    # Get list of base flags and typed flags in format
    base_flags = [flag(size=f[0],
                       edges=[list(e) for e in f[2][0][1]], ftype=[],
                       coloured=[x[0] for x in (f[2][2][1] if
                                                len(f[2]) > 1 else [])])
                  for f in base_flags_0]

    factor_flags = []
    for a in factor_flags_0:
        b = []
        for c in a:
            b += [flag(size=c[0],
                       edges=[list(d) for d in c[2][0][1]],
                       ftype=[d for d in c[1]],
                       coloured=[x[0] for x in (c[2][2][1] if
                                                len(c[2]) > 1 else [])])]
        factor_flags += [b]
    return base_flags, target, factor_flags, positives


# ---------------------------- Objectives and assumptions ----------------------------

def prop_3_1():
    """Proposition 3.1: maximize the ordinary 3-edge density."""
    return [(Fraction(1), flag(3, edges=[[0, 1, 2]]))], []


def prop_3_2(base_flags):
    """Proposition 3.2: exact six-vertex objective and degree assumption.

    The objective follows the Sage notebook's construction: choose a rooted
    edge, require its complementary three vertices also to span an edge,
    and require exactly two of these three vertices to have at least two
    edges meeting the root triple in two vertices. Average over all 20
    possible root triples. All arithmetic is rational.

    The degree positivity assumption is degree - (4641/10000). For a full
    rooted 3-flag basis, its edge coefficient is 1 - beta and its nonedge
    coefficient is -beta (not 1 and -beta).
    """
    beta = Fraction(4641, 10000)
    target = []
    V = frozenset(range(6))
    for G in base_flags:
        if G.size != 6:
            raise ValueError('Proposition 3.2 requires six-vertex base flags')
        edges = {frozenset(e) for e in G.edges}
        count = 0
        for root in combinations(range(6), 3):
            rest = V - set(root)
            if frozenset(root) not in edges or frozenset(rest) not in edges:
                continue
            high_root_degrees = sum(
                sum(frozenset((v, a, b)) in edges
                    for a, b in combinations(root, 2)) >= 2
                for v in rest
            )
            if high_root_degrees == 2:
                count += 1
        target.append(Fraction(count, 20))
    degree_minus_beta = [
        (1 - beta, flag(3, edges=[[0, 1, 2]], ftype=[0])),
        (-beta, flag(3, ftype=[0]))
    ]
    return target, [degree_minus_beta]




def prop_3_3():
    """Returns the objective function as a quantum flag and the positivity
    assumptions as a list of quantum flags of Proposition 3.3.

    """
    objective = [
        (1, flag(3, edges=[[0, 1, 2]])),
        (1, flag(3, edges=[[0, 1, 2]], coloured=[1, 2])),
        (-(1 - Fraction(1, 4000)), flag(3, coloured=[2]))
    ]

    # Define assumptions

    # Local optimality
    Cp0_minus_Bp0 = [
        (1, flag(3, edges=[[0, 1, 2]], ftype=[0], coloured=[2])),
        (-1, flag(3, edges=[[0, 1, 2]], ftype=[0]))
    ]
    Cp1_minus_Bp1 = [
        (1, flag(3, edges=[[0, 1, 2]], ftype=[2], coloured=[2])),
        (-1, flag(3, edges=[[0, 1, 2]], ftype=[2], coloured=[1, 2]))
    ]

    # Degree regularity
    degeq_00 = [
        (1, flag(4, ftype=[2, 3], edges=[[0, 1, 2]], coloured=[0, 1])),
        (1, flag(4, ftype=[1, 3], edges=[[0, 1, 2]], coloured=[0])),
        (1, flag(4, ftype=[1, 0], edges=[[1, 2, 3]], coloured=[])),
        (1, flag(4, ftype=[2, 3], edges=[[0, 1, 2], [0, 2, 3]], coloured=[0, 1])),
        (1, flag(4, ftype=[2, 1], edges=[[0, 1, 2], [0, 2, 3]], coloured=[0])),
        (1, flag(4, ftype=[2, 1], edges=[[0, 2, 3], [1, 2, 3]], coloured=[0])),
        (1, flag(4, ftype=[0, 1], edges=[[0, 1, 2], [0, 2, 3]], coloured=[])),
        (1, flag(4, ftype=[3, 2], edges=[[0, 1, 3], [0, 2, 3], [1, 2, 3]], coloured=[0, 1])),
        (1, flag(4, ftype=[3, 1], edges=[[0, 1, 3], [0, 2, 3], [1, 2, 3]], coloured=[0])),
        (1, flag(4, ftype=[3, 0], edges=[[0, 1, 3], [0, 2, 3], [1, 2, 3]], coloured=[])),
        (-1, flag(4, ftype=[3, 2], edges=[[0, 1, 2]], coloured=[0, 1])),
        (-1, flag(4, ftype=[3, 1], edges=[[0, 1, 2]], coloured=[0])),
        (-1, flag(4, ftype=[0, 1], edges=[[1, 2, 3]], coloured=[])),
        (-1, flag(4, ftype=[3, 2], edges=[[0, 1, 2], [0, 2, 3]], coloured=[0, 1])),
        (-1, flag(4, ftype=[1, 2], edges=[[0, 1, 2], [0, 2, 3]], coloured=[0])),
        (-1, flag(4, ftype=[1, 2], edges=[[0, 2, 3], [1, 2, 3]], coloured=[0])),
        (-1, flag(4, ftype=[1, 0], edges=[[0, 1, 2], [0, 2, 3]], coloured=[])),
        (-1, flag(4, ftype=[2, 3], edges=[[0, 1, 3], [0, 2, 3], [1, 2, 3]], coloured=[0, 1])),
        (-1, flag(4, ftype=[1, 3], edges=[[0, 1, 3], [0, 2, 3], [1, 2, 3]], coloured=[0])),
        (-1, flag(4, ftype=[0, 3], edges=[[0, 1, 3], [0, 2, 3], [1, 2, 3]], coloured=[])),
    ]

    minus_degeq_00 = [(-coef, flag) for coef, flag in degeq_00]

    degeq_01 = [
        (1, flag(4, ftype=[3, 0], edges=[[1, 2, 3]], coloured=[0, 1, 2])),
        (1, flag(4, ftype=[2, 0], edges=[[1, 2, 3]], coloured=[0, 1])),
        (1, flag(4, ftype=[1, 0], edges=[[1, 2, 3]], coloured=[0])),
        (1, flag(4, ftype=[3, 0], edges=[[0, 2, 3], [1, 2, 3]], coloured=[0, 1, 2])),
        (1, flag(4, ftype=[2, 1], edges=[[0, 1, 2], [0, 2, 3]], coloured=[0, 1])),
        (1, flag(4, ftype=[2, 0], edges=[[0, 2, 3], [1, 2, 3]], coloured=[0, 1])),
        (1, flag(4, ftype=[2, 0], edges=[[0, 2, 3], [1, 2, 3]], coloured=[0])),
        (1, flag(4, ftype=[3, 0], edges=[[0, 1, 3], [0, 2, 3], [1, 2, 3]], coloured=[0, 1, 2])),
        (1, flag(4, ftype=[3, 0], edges=[[0, 1, 3], [0, 2, 3], [1, 2, 3]], coloured=[0, 1])),
        (1, flag(4, ftype=[3, 0], edges=[[0, 1, 3], [0, 2, 3], [1, 2, 3]], coloured=[0])),
        (-1, flag(4, ftype=[3, 0], edges=[[0, 1, 2]], coloured=[0, 1, 2])),
        (-1, flag(4, ftype=[3, 0], edges=[[0, 1, 2]], coloured=[0, 1])),
        (-1, flag(4, ftype=[3, 0], edges=[[0, 1, 2]], coloured=[0])),
        (-1, flag(4, ftype=[3, 0], edges=[[0, 1, 2], [0, 2, 3]], coloured=[0, 1, 2])),
        (-1, flag(4, ftype=[3, 0], edges=[[0, 1, 2], [0, 2, 3]], coloured=[0, 1])),
        (-1, flag(4, ftype=[2, 0], edges=[[0, 1, 2], [0, 1, 3]], coloured=[0, 1])),
        (-1, flag(4, ftype=[1, 0], edges=[[0, 1, 2], [0, 2, 3]], coloured=[0])),
        (-1, flag(4, ftype=[3, 0], edges=[[0, 1, 2], [0, 1, 3], [0, 2, 3]], coloured=[0, 1, 2])),
        (-1, flag(4, ftype=[2, 0], edges=[[0, 1, 2], [0, 1, 3], [0, 2, 3]], coloured=[0, 1])),
        (-1, flag(4, ftype=[1, 0], edges=[[0, 1, 2], [0, 1, 3], [0, 2, 3]], coloured=[0]))
    ]

    minus_degeq_01 = [(-coef, flag) for coef, flag in degeq_01]

    degeq_11 = [
        (1, flag(4, ftype=[1, 0], edges=[[1, 2, 3]], coloured=[0, 1, 2, 3])),
        (1, flag(4, ftype=[1, 0], edges=[[1, 2, 3]], coloured=[0, 1, 2])),
        (1, flag(4, ftype=[1, 0], edges=[[1, 2, 3]], coloured=[0, 1])),
        (1, flag(4, ftype=[0, 1], edges=[[0, 1, 2], [0, 2, 3]], coloured=[0, 1, 2, 3])),
        (1, flag(4, ftype=[0, 1], edges=[[0, 1, 2], [0, 2, 3]], coloured=[0, 1, 2])),
        (1, flag(4, ftype=[2, 0], edges=[[0, 2, 3], [1, 2, 3]], coloured=[0, 1, 2])),
        (1, flag(4, ftype=[0, 1], edges=[[0, 1, 2], [0, 2, 3]], coloured=[0, 1])),
        (1, flag(4, ftype=[3, 0], edges=[[0, 1, 3], [0, 2, 3], [1, 2, 3]], coloured=[0, 1, 2, 3])),
        (1, flag(4, ftype=[0, 1], edges=[[0, 1, 2], [0, 1, 3], [0, 2, 3]], coloured=[0, 1, 2])),
        (1, flag(4, ftype=[0, 1], edges=[[0, 1, 2], [0, 1, 3], [0, 2, 3]], coloured=[0, 1])),
        (-1, flag(4, ftype=[0, 1], edges=[[1, 2, 3]], coloured=[0, 1, 2, 3])),
        (-1, flag(4, ftype=[0, 1], edges=[[1, 2, 3]], coloured=[0, 1, 2])),
        (-1, flag(4, ftype=[0, 1], edges=[[1, 2, 3]], coloured=[0, 1])),
        (-1, flag(4, ftype=[1, 0], edges=[[0, 1, 2], [0, 2, 3]], coloured=[0, 1, 2, 3])),
        (-1, flag(4, ftype=[1, 0], edges=[[0, 1, 2], [0, 2, 3]], coloured=[0, 1, 2])),
        (-1, flag(4, ftype=[0, 2], edges=[[0, 2, 3], [1, 2, 3]], coloured=[0, 1, 2])),
        (-1, flag(4, ftype=[1, 0], edges=[[0, 1, 2], [0, 2, 3]], coloured=[0, 1])),
        (-1, flag(4, ftype=[0, 3], edges=[[0, 1, 3], [0, 2, 3], [1, 2, 3]], coloured=[0, 1, 2, 3])),
        (-1, flag(4, ftype=[1, 0], edges=[[0, 1, 2], [0, 1, 3], [0, 2, 3]], coloured=[0, 1, 2])),
        (-1, flag(4, ftype=[1, 0], edges=[[0, 1, 2], [0, 1, 3], [0, 2, 3]], coloured=[0, 1]))
    ]

    minus_degeq_11 = [(-coef, flag) for coef, flag in degeq_11]

    # Best fit assumption (C - 43099222842404160803560878960795552015959/109451433676840840832138301242637359506224 >= 0)
    C_minus_2_5 = [
        (-Fraction(43099222842404160803560878960795552015959, 109451433676840840832138301242637359506224), flag(3, coloured=[0, 1, 2])),
        (-Fraction(43099222842404160803560878960795552015959, 109451433676840840832138301242637359506224), flag(3, coloured=[0, 1])),
        (-Fraction(43099222842404160803560878960795552015959, 109451433676840840832138301242637359506224), flag(3, coloured=[0])),
        (-Fraction(43099222842404160803560878960795552015959, 109451433676840840832138301242637359506224), flag(3, coloured=[])),
        (-Fraction(43099222842404160803560878960795552015959, 109451433676840840832138301242637359506224), flag(3, [[0, 1, 2]], coloured=[0, 1, 2])),
        (-Fraction(43099222842404160803560878960795552015959, 109451433676840840832138301242637359506224), flag(3, [[0, 1, 2]], coloured=[0, 1])),
        (1-Fraction(43099222842404160803560878960795552015959, 109451433676840840832138301242637359506224), flag(3, [[0, 1, 2]], coloured=[0])),
        (-Fraction(43099222842404160803560878960795552015959, 109451433676840840832138301242637359506224), flag(3, [[0, 1, 2]]))
    ]

    positivity_assumptions = [
        # Local optimality
        Cp0_minus_Bp0, Cp1_minus_Bp1,
        # Degree regularity
        degeq_00, minus_degeq_00,
        degeq_01, minus_degeq_01,
        degeq_11, minus_degeq_11,
        # Best fit assumption
        C_minus_2_5
    ]
    return objective, positivity_assumptions


# ---------------------------- Main ----------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Independently validate the objective and positivity assumptions in flag-algebra certificates.")
    parser.add_argument("certificate", help="trusted certificate pickle; can include a directory path")
    parser.add_argument("--range", nargs=2, type=int, metavar=("START", "END"),
                        dest="selection", help="1-based inclusive flag-index range (default: all flags)")
    args = parser.parse_args()
    name = __import__('os').path.basename(args.certificate)
    if name not in ('prop_3_1_k4f5c.pickle',
                    'prop_3_2_k4f5c.pickle',
                    'prop_3_3_k4f5c.pickle'):
        parser.error('Unknown proposition certificate name')

    base_flags, target, factor_flags, positives = recover_from_certificate(args.certificate)
    first, last = (0, len(base_flags))
    if args.selection:
        a, b = args.selection
        if not (1 <= a <= b <= len(base_flags)):
            parser.error(f'--range must satisfy 1 <= START <= END <= {len(base_flags)}')
        first, last = a - 1, b
    if name == 'prop_3_1_k4f5c.pickle':
        objective, positivity_assumptions = prop_3_1()
    elif name == 'prop_3_2_k4f5c.pickle':
        objective, positivity_assumptions = prop_3_2(base_flags)
    else:
        objective, positivity_assumptions = prop_3_3()

    if isinstance(objective[0], tuple):
        objective = change_base(objective, base_flags)
    if target[first:last] != objective[first:last]:
        for idx in range(first, last):
            if target[idx] != objective[idx]:
                raise SystemExit(f'FAILED: objective at base graph {idx + 1}: '
                                 f'{target[idx]} != {objective[idx]}')
    print(f'PASS: Objective matches for flags {first + 1}--{last}.')

    if len(factor_flags) != len(positivity_assumptions):
        raise SystemExit('FAILED: number of positivity assumptions does not match.')
    expected_rows = sum(len(group) for group in factor_flags)
    if len(positives) != expected_rows:
        raise SystemExit('FAILED: number of positivity coefficient rows does not match.')
    n_prev = 0
    for j, assumption in enumerate(positivity_assumptions):
        rows = [get_densities(assumption, factor_flags[j], G)
                for G in base_flags[first:last]]
        actual = list(map(list, zip(*rows)))
        expected = [col[first:last] for col in
                    positives[n_prev:n_prev + len(factor_flags[j])]]
        if actual != expected:
            for local_idx, (got, want) in enumerate(zip(actual, expected)):
                if got != want:
                    k = next(i for i, (a, b) in enumerate(zip(got, want)) if a != b)
                    raise SystemExit('FAILED: positivity assumption {}, factor {}, base graph {}: '
                                     'expected {}, certificate {}'.format(
                                         j + 1, local_idx + 1, first + k + 1,
                                         got[k], want[k]))
            raise SystemExit(f'FAILED: positivity assumption {j + 1}.')
        print(f'PASS: Positivity assumption {j + 1} matches ({len(actual)} rows).')
        n_prev += len(factor_flags[j])
    print('SUCCESS: Independently reconstructed objectives and positivity assumptions match.')


if __name__ == '__main__':
    main()
