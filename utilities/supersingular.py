"""
Helper functions for the supersingular elliptic curve computations in FESTA
"""

# Sage Imports
from sage.all import ZZ, factor, inverse_mod, CRT

# Local imports
from montgomery_isogenies.kummer_line import KummerLine

from utilities.order import has_order_D
from utilities.pairing import weil_pairing_pari
from utilities.discrete_log import BiDLP
from utilities.fast_roots import sqrt_Fp2, quadratic_roots

# =========================================== #
#   Extract coefficent from Montgomery curve  #
# =========================================== #


def montgomery_coefficient(E):
    a_inv = E.a_invariants()
    A = a_inv[1]
    if a_inv != (0, A, 0, 1, 0):
        raise ValueError(
            "Parent function assumes the curve E is in the Montgomery model"
        )
    return A


# =========================================== #
#   Isomorphisms between Montgomery curves    #
# =========================================== #


def montgomery_curve_isomorphism(EA, EB, T1, T2):
    """
    Given two isomorphic Montgomery curves EA and EB, it computes an
    isomorphism phi : EA -> EB, phi(T1) and phi(T2)
    """

    i = EA.base_ring().gen()
    assert i**2 == -1
    A = montgomery_coefficient(EA)
    B = montgomery_coefficient(EB)

    # r = 0
    if A == B:
        return EB(*T1.xy()), EB(*T2.xy())
    if A == -B:
        # u ^ 2 = -1
        T1_p = EB(-T1[0], -i * T1[1])
        T2_p = EB(-T2[0], -i * T2[1])
        return T1_p, T2_p

    r1, r2 = quadratic_roots(B, 1)
    A_inv = 1 / A

    u_sqrd = (3 * r1 + B) * A_inv
    if u_sqrd.is_square() and (r1 * (3 * r1 + 2 * B) + 1) == u_sqrd * u_sqrd:
        r = r1
        u = sqrt_Fp2(u_sqrd)
    else:
        u_sqrd = (3 * r2 + B) * A_inv
        r = r2
        u = sqrt_Fp2(u_sqrd)

    u_cube = u_sqrd * u
    T1_p = EB(u_sqrd * T1[0] + r, u_cube * T1[1])
    T2_p = EB(u_sqrd * T2[0] + r, u_cube * T2[1])

    return T1_p, T2_p


# =========================================== #
# Compute points of order D and Torsion Bases #
# =========================================== #


def generate_kummer_point(E, x_start=0):
    """
    Generate points on a curve E with x-coordinate
    i + x for x in Fp and i is the generator of Fp^2
    such that i^2 = -1.
    """
    F = E.base_ring()
    one = F.one()

    if x_start:
        x = F.gen() + x_start
    else:
        x = F.gen() + one

    K = KummerLine(E)
    A = montgomery_coefficient(E)

    # Try 1000 times then give up, just protection
    # for infinite loops
    for _ in range(1000):
        y2 = x * (x**2 + A * x + 1)
        if y2.is_square():
            yield K(x)
        x += one

    raise ValueError(
        "Generated 1000 points, something is probably going wrong somewhere."
    )


def clear_cofactor(P, k, even_power=None):
    # In the case where we remove a known
    # large even factor, it's faster to
    # do this first to save on addition in
    # xDBLADD
    if even_power:
        P = P.double_iter(even_power)
        k >>= even_power
    return k*P


def generate_point_order_D(E, D, x_start=0, even_power=None):
    """
    Input:  An elliptic curve E / Fp2
            An integer D dividing (p +1)
    Output: A point P of order D.
    """
    p = E.base().characteristic()
    k = (p + 1) // D

    Ps = generate_kummer_point(E, x_start=x_start)
    for G in Ps:
        P = clear_cofactor(G, k, even_power=even_power)

        # Case when we randomly picked
        # a point in the n-torsion
        if P.is_zero():
            continue

        # Check that P has order exactly D
        if has_order_D(P, D):
            yield P.curve_point()

    raise ValueError(f"Never found a point P of order D.")


def compute_point_order_D(E, D, x_start=0, even_power=None):
    """
    Wrapper function around a generator which returns the first
    point of order D
    """
    return generate_point_order_D(E, D, x_start=x_start, even_power=even_power).__next__()


def compute_linearly_independent_point_with_pairing(E, P, D, x_start=0, even_power=None):
    """
    Input:  An elliptic curve E / Fp2
            A point P ∈ E[D]
            An integer D dividing (p +1)
    Output: A point Q such that E[D] = <P, Q>
            The Weil pairing e(P,Q)
    """
    Qs = generate_point_order_D(E, D, x_start=x_start, even_power=even_power)
    for Q in Qs:
        # Make sure the point is linearly independent
        pair = weil_pairing_pari(P, Q, D)
        if has_order_D(pair, D, multiplicative=True):
            Q._order = ZZ(D)
            return Q, pair
    raise ValueError("Never found a point Q linearly independent to P")


def compute_linearly_independent_point(E, P, D, x_start=0, even_power=None):
    """
    Wrapper function around `compute_linearly_independent_point_with_pairing`
    which only returns a linearly independent point
    """
    Q, _ = compute_linearly_independent_point_with_pairing(E, P, D, x_start=x_start, even_power=even_power)
    return Q


def torsion_basis_with_pairing(E, D, even_power=None):
    """
    Generate basis of E(Fp^2)[D] of supersingular curve

    While computing E[D] = <P, Q> we naturally compute the
    Weil pairing e(P,Q), which we also return as in some cases
    the Weil pairing is then used when solving the BiDLP
    """
    p = E.base().characteristic()

    # Ensure D divides the curve's order
    if (p + 1) % D != 0:
        print(f"{factor(D) = }")
        print(f"{factor(p+1) = }")
        raise ValueError(f"D must divide the point's order")

    P = compute_point_order_D(E, D, even_power=even_power)
    Q, ePQ = compute_linearly_independent_point_with_pairing(E, P, D, x_start=P[0], even_power=even_power)

    return P, Q, ePQ


def torsion_basis(E, D, even_power=None):
    """
    Wrapper function around torsion_basis_with_pairing which only
    returns the torsion basis <P,Q> = E[D]
    """
    P, Q, _ = torsion_basis_with_pairing(E, D, even_power=even_power)
    return P, Q


# =========================================== #
#   Entangled torsion basis for fast E[2^k]   #
# =========================================== #


def precompute_elligator_tables(F):
    """
    Precomputes tables of quadratic non-residue or
    quadratic residue in Fp2. Used to compute entangled
    torsion bases following https://ia.cr/2017/1143
    """
    u = 2 * F.gen()

    T1 = dict()
    T2 = dict()
    # TODO: estimate how large r should be
    for r in range(1, 30):
        v = 1 / (1 + u * r**2)
        if v.is_square():
            T2[r] = v
        else:
            T1[r] = v
    return T1, T2


def entangled_torsion_basis(E, elligator_tables, cofactor):
    """
    Optimised algorithm following https://ia.cr/2017/1143
    which modifies the elligator method of hashing to points
    to find points P,Q of order k*2^b. Clearing the cofactor
    gives the torsion basis without checking the order or
    computing a Weil pairing.

    To do this, we need tables TQNR, TQR of pairs of values
    (r, v) where r is an integer and v = 1/(1 + ur^2) where
    v is either a quadratic non-residue or quadratic residue
    in Fp2 and u = 2i = 2*sqrt(-1).
    """
    F = E.base_ring()
    p = F.characteristic()
    p_sqrt = (p + 1) // 4

    i = F.gen()
    u = 2 * i
    u0 = 1 + i

    TQNR, TQR = elligator_tables

    # Pick the look up table depending on whether
    # A = a + ib is a QR or NQR
    A = E.a_invariants()[1]
    if (0, A, 0, 1, 0) != E.a_invariants():
        raise ValueError("The elliptic curve E must be in Montgomery form")
    if A.is_square():
        T = TQNR
    else:
        T = TQR

    # Look through the table to find point with
    # rational (x,y)
    y = None
    for r, v in T.items():
        x = -A * v

        t = x * (x**2 + A * x + 1)

        # Break when we find rational y: t = y^2
        c, d = t.list()
        z = c**2 + d**2
        s = z**p_sqrt
        if s**2 == z:
            y = sqrt_Fp2(t)
            break

    if y is None:
        raise ValueError("Never found a y-coordinate, increase the lookup table size")

    z = (c + s) // 2
    alpha = z**p_sqrt
    beta = d / (2 * alpha)

    if alpha**2 == z:
        y = F([alpha, beta])
    else:
        y = -F([beta, alpha])

    S1 = E([x, y])
    S2 = E([u * r**2 * x, u0 * r * y])

    return cofactor * S1, cofactor * S2


# =============================================== #
#   Compute canonical representation for kernel   #
# =============================================== #


def compute_canonical_kernel(imP, imQ, D, basis=None, ePQ=None):
    """
    Algorithm 6 of the FESTA paper

    Given a scalar multiple of the of the image points phi(P)
    and phi(Q) in E', such that:
        <P,Q> = E[D] and phi : E -> E'
    compute the canonical kernel representation
        ker(phi) = <P + [s]Q>
    """
    E = imP.curve()
    assert E == imQ.curve()

    if basis == None:
        P, Q, ePQ = torsion_basis_with_pairing(E, D)
    else:
        P, Q = basis
        if ePQ is None:
            ePQ = weil_pairing_pari(P, Q, D)

    # phi(P) = [a1]P + [b1]Q
    a1, b1 = BiDLP(imP, P, Q, D, ePQ=ePQ)

    # phi(Q) = [a2]P + [b2]Q
    a2, b2 = BiDLP(imQ, P, Q, D, ePQ=ePQ)

    t1_list = []
    t2_list = []
    di_list = []

    # CRT recovery of t1, t2
    for pi, ei in D.factor():
        di = pi**ei
        di_list.append(di)

        if a1 % pi == 0:
            t1, t2 = 0, inverse_mod(a2, di)
        else:
            t1, t2 = inverse_mod(a1, di), 0

        t1_list.append(t1)
        t2_list.append(t2)

    t1 = CRT(t1_list, di_list)
    t2 = CRT(t2_list, di_list)

    return (t1 * b1 + t2 * b2) % D


# =============================================== #
#  Ensure Basis <P,Q> of E[2^k] has (0,0) under Q #
# =============================================== #


def fix_torsion_basis_renes(P, Q, k):
    """
    Set the torsion basis P,Q such that
    2^(k-1)Q = (0,0) to ensure that (0,0)
    is never a kernel of a two isogeny
    """
    cofactor = 2 ** (k - 1)

    R = cofactor * P
    if R[0] == 0:
        return Q, P
    R = cofactor * Q
    if R[0] == 0:
        return P, Q
    return P, P + Q
