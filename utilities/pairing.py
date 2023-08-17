# import pari for fast dlog
import cypari2

# Make instance of Pari
pari = cypari2.Pari()

# ===================================== #
#   Directly access pairings from Pari  #
# ===================================== #


def weil_pairing_pari(P, Q, D, check=False):
    """
    Wrapper around Pari's implementation of the Weil pairing
    Allows the check of whether P,Q are in E[D] to be optional
    """
    if check:
        nP, nQ = D * P, D * Q
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
