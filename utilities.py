# Sage imports
from sage.all import (
    ZZ,
    GF,
    Matrix,
    prod, 
    proof,
    cached_function,
    cached_method,
)

# import pari for fast dlog
import cypari2

# ============================== #
#  Utility functions for masking #
# ============================== #

def canonical_matrix(M):
    """
    We recover M up to a minus sign, we so
    ensure the sign of all matrices are 
    canonical
    """
    alpha = M[0][0]
    if alpha < -alpha:
        return -M
    return M    

def random_diag_matrix(R, alpha=None):
    r"""
    Compute a random, diagonal, unitary 
    invertible matrix M \in Z / (2^b) Z
    """
    if alpha is None:
        alpha = R.random_element()

    # Ensure alpha is even
    alpha = alpha + ZZ(alpha % 2) + 1
    M = Matrix(R, 2, 2, [alpha, 0, 0, ~alpha])

    return canonical_matrix(M)

def random_circulant_matrix(R, beta=None):
    r"""
    Compute a random, circulant, unitary 
    invertible matrix M \in Z / (2^b) Z
    """
    if beta is None:
        beta = R.random_element()
    b = 4*beta
    aa = R(b*b + 1)
    a = aa.sqrt() # TODO: fast sqrt mod 2^k?
    M = Matrix(R, 2, 2, [a, b, b, a])

    return canonical_matrix(M)

def random_matrix(R, ele=None, diag=True):
    """
    Compute a random masking matrix which
    is commutative, invertible and unitary
    in Z / (N Z)
    """
    if diag:
        return random_diag_matrix(R, alpha=ele)
    return random_circulant_matrix(R, beta=ele)

def mask_torsion_points(A, P, Q):
    """
    Evaluate a 2x2 matrix A on the vector of elliptic
    curve points [P, Q] by acting with scalar multiplication
    """
    # Extract out scalars
    a00, a01, a10, a11 = A.list()

    # Act with scalar multiplication
    R = a00*P + a01*Q
    S = a10*P + a11*Q

    return R, S

# ================================================== #
#  Code to check whether a group element has order D #
# ================================================== #

def batch_cofactor_mul_generic(G_list, pis, group_action, lower, upper):
    """
    Input:  A list of elements `G_list`, such that
                G is the first entry and the rest is empty
                in the sublist G_list[lower:upper]
            A list `pis` of primes p such that
                their product is D
            The `group_action` of the group
            Indices lower and upper
    Output: None

    NOTE: G_list is created in place
    """

    # check that indices are valid
    if lower > upper:
        raise ValueError(f"Wrong input to cofactor_multiples()")

    # last recursion step does not need any further splitting
    if upper - lower == 1:
        return

    # Split list in two parts,
    # multiply to get new start points for the two sublists,
    # and call the function recursively for both sublists.
    mid = lower + (upper - lower + 1) // 2
    cl, cu = 1, 1
    for i in range(lower, mid):
        cu = cu * pis[i]
    for i in range(mid, upper):
        cl = cl * pis[i]
    # cl = prod(pis[lower:mid])
    # cu = prod(pis[mid:upper])

    G_list[mid] = group_action(G_list[lower], cu)
    G_list[lower] = group_action(G_list[lower], cl)

    batch_cofactor_mul_generic(G_list, pis, group_action, lower, mid)
    batch_cofactor_mul_generic(G_list, pis, group_action, mid, upper)


@cached_function
def has_order_constants(D):
    """
    Helper function, finds constants to
    help with has_order_D
    """
    D = ZZ(D)
    pis = [p for p, _ in D.factor()]
    D_radical = prod(pis)
    Dtop = D // D_radical
    return Dtop, pis


def has_order_D(G, D, multiplicative=False):
    """
    Given an element G in a group, checks if the
    element has order exactly D. This is much faster
    than determining its order, and is enough for 
    many checks we need when computing the torsion
    basis.

    We allow both additive and multiplicative groups
    which means we can use this when computing the order
    of points and elements in Fp^k when checking the 
    multiplicative order of the Weil pairing output
    """
    # For the case when we work with elements of Fp^k
    if multiplicative:
        group_action = lambda a, k: a**k
        is_identity = lambda a: a == 1
        identity = 1
    # For the case when we work with elements of E / Fp^k
    else:
        group_action = lambda a, k: k * a
        is_identity = lambda a: a.is_zero()
        identity = G.curve()(0)

    if is_identity(G):
        return False

    D_top, pis = has_order_constants(D)

    # If G is the identity after clearing the top
    # factors, we can abort early
    Gtop = group_action(G, D_top)
    if is_identity(Gtop):
        return False

    G_list = [identity for _ in range(len(pis))]
    G_list[0] = Gtop

    # Lastly we have to determine whether removing any prime 
    # factors of the order gives the identity of the group
    if len(pis) > 1:
        batch_cofactor_mul_generic(G_list, pis, group_action, 0, len(pis))
        if not all([not is_identity(G) for G in G_list]):
            return False

    return True


# ===================================== #
#  Fast DLP solving using Weil pairing  #
# ===================================== #

# Make instance of Pari
pari = cypari2.Pari()

def discrete_log_pari(a, base, order):
    """
    Wrapper around pari discrete log. Works like a.log(b), 
    but allows us to use the optional argument order. This
    is important as we skip Pari attempting to factor the
    full order of Fp^k, which is slow.
    """
    x = pari.fflog(a, base, order)
    return ZZ(x)

def weil_pairing_pari(P, Q, D, check=False):
    """
    Wrapper around Pari's implementation of the Weil pairing
    Allows the check of whether P,Q are in E[D] to be optional
    """
    if check:
        nP, nQ = D*P, D*Q
        if nP.is_zero() or nQ.is_zero():
            raise ValueError("points must both be n-torsion")

    return pari.ellweilpairing(P.curve(), P, Q, D)

def tate_pairing_pari(P, Q, D):
    """
    Wrapper around Pari's implementation of the Tate pairing
    NOTE: this is not the reduced Tate pairing, so to make this
    match with SageMath you need

    P.tate_pairing(Q, D, k) == pari.elltatepairing(E, P, Q, D)**((p^k - 1) / D)
    """
    E = P.curve()
    return pari.elltatepairing(E, P, Q, D)

def _precompute_baby_steps(base, step, e):
    """
    Helper function to compute the baby steps for 
    pohlig_hellman_base and windowed_pohlig_hellman.
    """
    baby_steps = [base]
    for _ in range(e):
        base = base**step
        baby_steps.append(base)
    return baby_steps

def pohlig_hellman_base(a, base, e):
    """
    Solve the discrete log for a = base^x for 
    elements base,a of order 2^e using the 
    Pohlig-Hellman algorithm.
    """
    baby_steps = _precompute_baby_steps(base, 2, e)

    dlog = 0
    exp  = 2**(e-1)

    # Solve the discrete log mod 2, only 2 choices
    # for each digit!
    for i in range(e):
        if a**exp != 1:
            a /= baby_steps[i]
            dlog += (2**i)

        if a ==  1:
            break

        exp //= 2

    return dlog

def windowed_pohlig_hellman(a, base, e, window):
    """
    Solve the discrete log for a = base^x for 
    elements base,a of order 2^e using the 
    windowed Pohlig-Hellman algorithm following
    https://ia.cr/2016/963.

    Algorithm runs recursively, computing in windows
    l^wi for window=[w1, w2, w3, ...]. 

    Runs the base case when window = []
    """  
    # Base case when all windows have been used
    if not window:
        return pohlig_hellman_base(a, base, e)

    # Collect the next window
    w, window = window[0], window[1:]
    step = 2**w

    # When the window is not a divisor of e, we compute
    # e mod w and solve for both the lower e - e mod w
    # bits and then the e mod w bits at the end.
    e_div_w, e_rem = divmod(e, w)
    e_prime = e - e_rem

    # First force elements to have order e - e_rem
    a_prime = a**(2**e_rem)
    base_prime = base**(2**e_rem)

    # Compute base^(2^w*i) for i in (0, ..., e/w-1)
    baby_steps = _precompute_baby_steps(base_prime, step, e_div_w)
    
    # Initialise some pieces
    dlog = 0
    if e_prime:
        exp = 2**(e_prime - w)

        # Work in subgroup of size 2^w
        s = base_prime**(exp)

        # Windowed Pohlig-Hellman to recover dlog as
        # alpha = sum l^(i*w) * alpha_i
        for i in range(e_div_w):
            # Solve the dlog in 2^w
            ri = a_prime**exp
            alpha_i = windowed_pohlig_hellman(ri, s, w, window)

            # Update a value and dlog computation
            a_prime /= baby_steps[i]**(alpha_i)
            dlog  += alpha_i * step**i

            if a_prime == 1:
                break

            exp //= step

    # TODO: 
    # I don't know if this is a nice way to do
    # this last step... Works well enough but could
    # be improved I imagine...
    exp = 2**e_prime
    if e_rem:
        base_last = base**exp
        a_last = a / base**dlog
        dlog_last = pohlig_hellman_base(a_last, base_last, e_rem)
        dlog += exp * dlog_last

    return dlog

def BiDLP(R, P, Q, D, ePQ=None):
    """
    Given a basis P,Q of E[D] finds
    a,b such that R = [a]P + [b]Q.

    Uses the fact that 
        e([a]P + [b]Q, [c]P + [d]Q) = e(P,Q)^(ad-bc)

    Optional: include the pairing e(P,Q) which can be precomputed
    which is helpful when running multiple BiDLP problems with P,Q
    as input. This happens, for example, during compression.
    """
    # e(P,Q)
    if ePQ:
        pair_PQ = ePQ
    else:
        pair_PQ = weil_pairing_pari(P, Q, D)

    # Write R = aP + bQ for unknown a,b
    # e(R, Q) = e(P, Q)^a
    pair_a = weil_pairing_pari(R, Q, D)

    # e(R,-P) = e(P, Q)^b
    pair_b = weil_pairing_pari(R, -P, D)

    # Now solve the dlog in Fq
    a = discrete_log_pari(pair_a, pair_PQ, D)
    b = discrete_log_pari(pair_b, pair_PQ, D)

    return a, b

def BiDLP_power_two(R, P, Q, e, window, ePQ=None):
    r"""
    Same as the above, but uses optimisations using that
    D = 2^e.

    First, rather than use the Weil pairing, we can use
    the Tate pairing which is approx 2x faster. 

    Secondly, rather than solve the discrete log naively,
    we use an optimised windowed Pohlig-Hellman.

    NOTE: this could be optimised further, following the
    SIKE key-compression algorithms and for the Tate pairing
    we could reuse Miller loops to save even more time, but
    that seemed a little overkill for a SageMath PoC

    Finally, as the Tate pairing produces elements in \mu_n
    we also have fast inversion from conjugation, but SageMath
    has slow conjugation, so this doesn't help for now.
    """
    p = R.curve().base_ring().characteristic()
    D = 2**e
    exp = (p**2 - 1) // D

    # e(P,Q)
    if ePQ:
        pair_PQ = ePQ
    else:
        pair_PQ = tate_pairing_pari(P, Q, D)**exp

    # Write R = aP + bQ for unknown a,b
    # e(R, Q) = e(P, Q)^a
    pair_a = tate_pairing_pari(Q, -R, D)**exp

    # e(R,-P) = e(P, Q)^b
    pair_b = tate_pairing_pari(P, R, D)**exp

    # Now solve the dlog in Fq
    a = windowed_pohlig_hellman(pair_a, pair_PQ, e, window)
    b = windowed_pohlig_hellman(pair_b, pair_PQ, e, window)

    return a, b

def DLP_power_two(R, P, Q, e, window, ePQ=None, first=True):
    r"""
    This is the same as BiDLP but it only returns either a or b 
    depending on whether first is true or false.
    This is used in compression, where we only send 3 of the 4 
    scalars from BiDLP
    """
    p = R.curve().base_ring().characteristic()
    D = 2**e
    exp = (p**2 - 1) // D

    # e(P,Q)
    if ePQ:
        pair_PQ = ePQ
    else:
        pair_PQ = tate_pairing_pari(P, Q, D)**exp

    if first:
        pair_a = tate_pairing_pari(Q, -R, D)**exp
        x = windowed_pohlig_hellman(pair_a, pair_PQ, e, window)

    else:
        pair_b = tate_pairing_pari(P, R, D)**exp
        x = windowed_pohlig_hellman(pair_b, pair_PQ, e, window)

    return x

# ================================================ #
#     Compute optimised strategy for (2,2)-chain   #
# ================================================ #

def optimised_strategy(n):
    """
    Algorithm 60: https://sike.org/files/SIDH-spec.pdf
    Shown to be appropriate for (l,l)-chains in 
    https://ia.cr/2023/508
    
    Note: the costs we consider are:
       eval_c: the cost of one Richelot isogeny
       mul_x:  the cost of one divisor doubling
       
    We find that a single doubling has approximately a 
    cost of 0.175*evaluation which suggests the params
    
    eval_c = 1.000
    mul_c  = 0.175
    
    However, the first doubling has a slightly higher cost
    due to the parsing of the polynomials to extract coefficients
    and so repeated doubling is cheaper than a single doubling.
    
    The best case would be to work this into the strategy generation,
    but while these costs are so noisey from SageMath benchmarks we
    find there's little use in overoptimising here to save a few ms.
    
    Once a proper, constant-time implementation has been developed
    we can re-address these costs and find the true optimal strategy 
    of the (2,2)-chain for decryption.
    """

    eval_c = 1.000
    mul_c  = 0.175

    S = {1:[]}
    C = {1:0 }
    for i in range(2, n+1):
        b, cost = min(((b, C[i-b] + C[b] + b*mul_c + (i-b)*eval_c) for b in range(1,i)), key=lambda t: t[1])
        S[i] = [b] + S[i-b] + S[b]
        C[i] = cost

    return S[n]

# ====================================== #
#  Compute f^(-1) mod g for f,g in R[X]  #
# ====================================== #

def invert_mod_polynomial_quadratic(f, g):
    """
    Given polynomials f, g with deg(g) = 2
    and deg(f) < deg(g) compute f^(-1) mod g

    Inverse is found using linear algebra.
    Cost: 9M 1I
    """
    R = f.parent()

    f0, f1 = f[0], f[1]
    g0, g1, g2 = g[0], g[1], g[2]

    f0_g2 = f0*g2

    A =  f0_g2
    B = -f1*g0
    C =  f1*g2
    D =  f0_g2 - f1*g1

    inv_det = ((A*D - C*B)).inverse()
    inv_det *= g2

    h0 =  D * inv_det
    h1 = -C * inv_det

    return R([h0, h1])

def invert_mod_polynomial_quartic(f, g):
    """
    Given polynomials f, g with deg(g) = 4
    and deg(f) < deg(g) compute f^(-1) mod g

    Inverse is found using linear algebra.
    Cost: 48M 2I

    - 4M 1I (ensure g is monic)
    - 12M   (compute matrix coefficients)
    - 28M   (compute matrix determinant)
    - 4M 1I (solve for inverse)

    - 48M 2I (total)
    """
    R = f.parent()

    # Extract out the coefficients
    f0, f1, f2, f3 = f.list()
    g0, g1, g2, g3, g4 = g.monic().list()

    # First, make g monic
    # 1I 4M
    if not g4.is_one():
        inv_g4 = 1/g4
        g0 *= inv_g4
        g1 *= inv_g4
        g2 *= inv_g4
        g3 *= inv_g4

    # Compute the matrix coefficients
    # M = [[a1, a2, a3, a4]
    #      [b1, b2, b3, b4]
    #      [c1, c2, c3, c4]
    #      [d1, d2, d3, d4]

    # Linear shift and feedback
    # 12M
    a1, b1, c1, d1 = f0, f1, f2, f3

    a2 =    - d1*g0
    b2 = a1 - d1*g1
    c2 = b1 - d1*g2
    d2 = c1 - d1*g3

    a3 =    - d2*g0
    b3 = a2 - d2*g1
    c3 = b2 - d2*g2
    d3 = c2 - d2*g3

    a4 =    - d3*g0
    b4 = a3 - d3*g1
    c4 = b3 - d3*g2
    d4 = c3 - d3*g3


    # Compute the determinant of the matrix
    # 28M
    b1_c2 = b1*c2
    b1_c3 = b1*c3
    b1_c4 = b1*c4

    b2_c1 = b2*c1
    b2_c3 = b2*c3
    b2_c4 = b2*c4

    b3_c1 = b3*c1    
    b3_c2 = b3*c2
    b3_c4 = b3*c4
    
    b4_c1 = b4*c1  
    b4_c2 = b4*c2
    b4_c3 = b4*c3

    D1 = (b3_c4 - b4_c3)*d2 + (b4_c2 - b2_c4)*d3 + (b2_c3 - b3_c2)*d4
    D2 = (b4_c3 - b3_c4)*d1 - (b4_c1 - b1_c4)*d3 + (b3_c1 - b1_c3)*d4
    D3 = (b2_c4 - b4_c2)*d1 + (b4_c1 - b1_c4)*d2 + (b1_c2 - b2_c1)*d4
    D4 = (b3_c2 - b2_c3)*d1 - (b3_c1 - b1_c3)*d2 + (b2_c1 - b1_c2)*d3

    det  = a1*D1 + a2*D2 + a3*D3 + a4*D4

    # Invert the determinant
    # 1I
    det_inv = 1/det

    # Compute solution
    # 4M
    h0 = D1*det_inv
    h1 = D2*det_inv
    h2 = D3*det_inv
    h3 = D4*det_inv

    return R([h0,h1,h2,h3])

# ========================== #
#  Cast Sage types to bytes  #
# ========================== #

# TODO: I should just push this to Sage,
# or at least some kind of nice API of field 
# elements to bytes.

def integer_to_bytes(n, byte_len=None):
    """
    Represent an integer n as bytes in big-endian
    """
    n = int(n)
    if byte_len is None:
        byte_len = (n.bit_length() + 7) // 8
    return n.to_bytes(byte_len, 'big')

def bytes_to_integer(b):
    """
    Represent bytes as an integer in big-endian
    """
    return ZZ(int.from_bytes(b, 'big'))

# ============================================ #
#     Fast square root and quadratic roots     #
# ============================================ #

def sqrt_Fp2(a):
    """
    Efficiently computes the sqrt
    of an element in Fp2 using that
    we always have a prime p such that
    p ≡ 3 mod 4.
    """
    Fp2 = a.parent()
    p = Fp2.characteristic()
    i = Fp2.gen() # i = √-1

    a1 = a ** ((p - 3) // 4)
    x0 = a1 * a
    alpha = a1 * x0

    if alpha == -1:
        x = i * x0
    else:
        b = (1 + alpha) ** ((p - 1) // 2)
        x = b * x0

    return x

def quadratic_roots(b, c):
    """
    Computes roots to the quadratic polynomial

        f = x^2 + b * x + c

    Using the quadratic formula

    Just like in school!
    """
    d2 = b**2 - 4 * c
    d = sqrt_Fp2(d2)
    return ((-b + d) / 2, -(b + d) / 2)

# ========================== #
#     Speed up SageMath!     #
# ========================== #

def speed_up_sagemath():
    """
    First we set proof.all(False) for general speed ups which
    keeping everything correct (enough)
    
    Then, we apply a monkey patch to cache the vector_space which
    helps with the performance of polynomial computations
    """
    # Skips strong primality checks and other slow things
    proof.all(False)
    
    # Cache vector spaces to improve performance
    p = 2**127 - 1 # Arbitrary large prime
    to_patch = [GF(3), GF(3**2), GF(p), GF(p**2)]
    for x in to_patch:
        type(x).vector_space = cached_method(type(x).vector_space)

# ====================== #
#     Print Debugging    #
# ====================== #

def print_info(str, banner="="):
    """
    Print information with a banner to help
    with visibility during debug printing
    """
    print(banner * 80)
    print(f"{str}".center(80))
    print(banner * 80)
