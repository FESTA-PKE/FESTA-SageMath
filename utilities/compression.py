"""
This file is written to hold the compression helper functions for the 
main (de)compression functions of FESTA(+).
"""

# Sage imports
from sage.all import EllipticCurve, inverse_mod

# Local imports
from utilities.supersingular import entangled_torsion_basis
from utilities.pairing import tate_pairing_pari
from utilities.discrete_log import (
    BiDLP_power_two,
    DLP_power_two,
    windowed_pohlig_hellman,
)
from utilities.utils import (
    integer_to_bytes,
    bytes_to_integer,
)

# ======================================================== #
#     Helper Functions for compression and decompression   #
# ======================================================== #


def tuple_to_bytes(v, byte_len):
    """
    Takes a tuple of integers (a0, a1, a2, ...) and represents
    them as bytes
    """
    return b"".join(integer_to_bytes(vi, byte_len=byte_len) for vi in v)


def tuple_from_bytes(v_bytes, byte_len):
    """
    Takes bytes representing a tuple of integers and parses
    them to recover the integers (a0, a1, a2, ...)
    """
    v_bytes = [v_bytes[i : i + byte_len] for i in range(0, len(v_bytes), byte_len)]
    v = [bytes_to_integer(vi_bytes) for vi_bytes in v_bytes]
    return v


def field_element_to_bytes(x, byte_len):
    """
    Represent an element x = a + i*b in Fp^2 as bytes
    """
    # x = a + b*i
    a_int, b_int = x.list()
    # convert to bytes
    a_bytes = integer_to_bytes(a_int, byte_len=byte_len)
    b_bytes = integer_to_bytes(b_int, byte_len=byte_len)
    # concatenate bytes
    return a_bytes + b_bytes


def field_element_from_bytes(Fp2, x_bytes, byte_len):
    """
    Recover an element x = a + i*b in Fp^2 from bytes
    """
    a_bytes = x_bytes[:byte_len]
    b_bytes = x_bytes[byte_len:]
    a_int = bytes_to_integer(a_bytes)
    b_int = bytes_to_integer(b_bytes)
    return Fp2([a_int, b_int])


# ======================================================== #
#         Compress and Decompress Montgomery Curve         #
# ======================================================== #


def compress_curve(E, p_byte_len):
    """
    Given a Montgomery curve: y^2 = x(x^2 + Ax + 1) represent
    it by its Montgomery constant A in bytes
    """
    a_inv = E.a_invariants()
    A = a_inv[1]
    assert a_inv == (0, A, 0, 1, 0), "Curve must be in Montgomery form"

    return field_element_to_bytes(A, p_byte_len)


def decompress_curve(F, A_bytes, p_byte_len):
    """
    Given the Montgomery constant A represented as bytes
    recover the Montgomery curve: y^2 = x(x^2 + Ax + 1)
    """
    A = field_element_from_bytes(F, A_bytes, p_byte_len)
    return EllipticCurve(F, [0, A, 0, 1, 0])


# ======================================================== #
#          Compress and Decompress Torsion Basis           #
# ======================================================== #


def compress_two_torsion_basis(
    E, R, S, b, elligator_tables, cofactor, dlog_window, basis_byte_len
):
    """
    Given a torsion basis E[2^b] = <R, S> compress the points
    as bytes by first finding a canonical torsion basis
    E[2^b] = <P, Q> and finding integers R = [aR]P + [bR]Q and
    S = [aS]P + [bS]Q

    The tuple (aR, bR, aS, bS) can be converted to bytes.

    NOTE: Computes the E[2^b] torsion basis using the optimised method
    of https://ia.cr/2017/1143 and for the discrete log uses the Tate
    pairing and optimised windowed Pohlig-Hellman by using that we're
    working with elements which order a power of two.
    """
    P, Q = entangled_torsion_basis(E, elligator_tables, cofactor)

    # Write the points R,S in their canonical representation
    p = E.base_ring().characteristic()
    D = 2**b
    e = (p**2 - 1) // D

    ePQ = tate_pairing_pari(P, Q, D) ** e
    aR, bR = BiDLP_power_two(R, P, Q, b, dlog_window, ePQ=ePQ)

    # From the Weil pairing dlog we will have:
    #     x = (aR * bS - bR * aS) mod 2^b
    # When aR is invertible, we send aS and recover bS
    # otherwise we send bS and invert bR to recover aS.
    # This means we only have to solve one discrete log
    # for S (rather than two for BiDLP)
    return_a = aR % 2
    xS = DLP_power_two(S, P, Q, b, dlog_window, ePQ=ePQ, first=return_a)

    # Package BiDLP output into vector
    dlp_vec = (aR, bR, xS)
    return tuple_to_bytes(dlp_vec, basis_byte_len)


def recover_lost_scalar(R, S, P, Q, b, deg_phi, aR, bR, xS, window):
    """
    This is a helper function to recover the lost fourth
    scalar such that imR = [aR]P + [bR]Q and imS = [aS]P + [bS]Q
    given either the tuple (aR, bR, aS) or (aR, bR, bS).

    To do this, we use that given imR, imS = phi(R), phi(S)
    compatibility with the Weil pairing means:
       e(R, S)^deg(phi) = e(imR, imS)

    For the canonical torsion basis E[D] = <P, Q>
       imR = [aR]P + [bR]Q , imS = [aS]P + [bS]Q
    We have:
       e(imR, imS) = e(P, Q)^(aR*bS - bR*aS)

    Solving the discrete log for
       e(R, S)^deg(phi) = e(P, Q)^(aR*bS - bR*aS)

    Assuming aR is invertible we have that
    e(R, S) = e(P, Q)^x

       x = (aR*bS - bR*aS) / deg(phi)
       bS = (x*deg(phi) + bR*aS) / aR

    When bR is invertible, we instead recover aS

       aS = (aR*bS - x*deg(phi)) / bR
    """
    p = R.base_ring().characteristic()
    D = 2**b
    exp = (p**2 - 1) // D

    ePQ = tate_pairing_pari(P, Q, D) ** exp
    eRS = tate_pairing_pari(R, S, D) ** exp
    x = windowed_pohlig_hellman(eRS, ePQ, b, window)

    if aR % 2 != 0:
        aS = xS
        bS = (x * deg_phi + bR * aS) * inverse_mod(aR, D)
        bS = bS % D
    else:
        bS = xS
        aS = (aR * bS - x * deg_phi) * inverse_mod(bR, D)
        aS = aS % D

    return aS, bS


def decompress_two_torsion_basis(
    E,
    dlp_bytes,
    preimage_data,
    b,
    elligator_tables,
    cofactor,
    dlog_window,
    basis_byte_len,
):
    """
    Given the vector (aR, bR, aS, bS) represented as bytes,
    we recover points R = [aR]P + [bR]Q and S = [aS]P + [bS]Q
    where <P, Q> = E[2^b] is the canonical torsion basis

    NOTE: Computes the E[2^b] torsion basis using the optimised method
    of https://ia.cr/2017/1143
    """
    P, Q = entangled_torsion_basis(E, elligator_tables, cofactor)

    # Extract out the DLP from bytes, compression means we lose
    # either aS or bS
    aR, bR, xS = tuple_from_bytes(dlp_bytes, basis_byte_len)

    # recover (aS, bS) from (aR, bR, xS) and the Weil pairing
    R, S, deg_phi = preimage_data
    aS, bS = recover_lost_scalar(R, S, P, Q, b, deg_phi, aR, bR, xS, dlog_window)

    # Reconstruct points
    imR = aR * P + bR * Q
    imS = aS * P + bS * Q

    return imR, imS


# ======================================================== #
#         Compress Montgomery Curve and Torsion Basis      #
# ======================================================== #


def compress_curve_and_two_torsion_basis(
    E, R, S, b, elligator_tables, cofactor, dlog_window, p_byte_len, basis_byte_len
):
    """
    Helper function which given a curve E and a torsion basis E[2^b] = <R,S>
    applies the above compression techniques and outputs a single byte string
    """
    A_bytes = compress_curve(E, p_byte_len)
    RS_bytes = compress_two_torsion_basis(
        E, R, S, b, elligator_tables, cofactor, dlog_window, basis_byte_len
    )
    return A_bytes + RS_bytes


def decompress_curve_and_two_torsion_basis(
    F,
    ERS_bytes,
    preimage_data,
    b,
    elligator_tables,
    cofactor,
    dlog_window,
    p_byte_len,
    basis_byte_len,
):
    """
    Helper function which given the bytes corresponding to (A, (aR, bR, aS, bS))
    recovers the elliptic curve y^2 = x(x^2 + Ax + 1) and the torsion basis
    E[2^b] = <R, S> such that R = [aR]P + [bR]Q and S = [aS]P + [bS]Q
    where <P, Q> = E[2^b] is the canonical torsion basis
    """
    A_bytes = ERS_bytes[: 2 * p_byte_len]
    dlp_bytes = ERS_bytes[2 * p_byte_len :]

    E = decompress_curve(F, A_bytes, p_byte_len)
    R, S = decompress_two_torsion_basis(
        E,
        dlp_bytes,
        preimage_data,
        b,
        elligator_tables,
        cofactor,
        dlog_window,
        basis_byte_len,
    )

    return E, R, S
