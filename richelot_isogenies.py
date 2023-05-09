"""
The file implements the (2^b,2^b)-isogeny between elliptic products and some helper
functions which FESTA call to compute and evaluate this isogeny.

The (2^b,2^b)-isogeny chain is computed in three steps:

1. A gluing isogeny which takes a pair of supersingular elliptic curves (E, E') 
   and pairs of points (P, Q), (P', Q') and outputs a Jacobian of a hyperelliptic 
   curve J(H) and a pair of divisors (D, D').

2. An isogeny between Jacobians of hyperelliptic curves with domain J(H) and kernel
   (D, D') and outputs the codomain J(H')

3. A splitting isogeny which takes J(H) and outputs a pair of supersingular elliptic
   curves.

The formula used follow https://ia.cr/2022/1283 for the compact explicit formula 
for the gluing isogeny and follow Benjamin Smith's thesis for the Richelot 
correspondence and Splitting isogeny. http://iml.univ-mrs.fr/~kohel/phd/thesis_smith.pdf

This code has been adapted from the file:

  https://github.com/jack4818/Castryck-Decru-SageMath/blob/main/richelot_aux.py

which was originally used for the SageMath implementation of the Castryck-Decru
attack, as well as the variant of the attack described by Oudompheng. See:

  https://ia.cr/2022/975
  https://www.normalesup.org/~oudomphe/textes/202208-castryck-decru-shortcut.pdf
"""

# Sage imports
from sage.all import (
    PolynomialRing,
    EllipticCurve,
    Matrix,
    vector,
)

# local imports
from divisor_arithmetic import affine_dbl_iter, affine_add
from supersingular import weil_pairing_pari
from utilities import quadratic_roots, invert_mod_polynomial_quadratic, invert_mod_polynomial_quartic

def FromProdToJac(P2, Q2, R2, S2):
    """
    Compact explicit formula for the gluing isogeny is derived in
    https://ia.cr/2022/1283
    """
    Fp2 = P2.curve().base()
    
    # NOTE: we use a generic implementation of the PolynomialRing
    # as this is faster than NTL in this scenario due to the slowness
    # of coefficient extraction and polynomial construction!
    Rx = PolynomialRing(Fp2, name="x", implementation="generic")
    x = Rx.gens()[0]

    a1, a2, a3 = P2[0], Q2[0], (P2 + Q2)[0]
    b1, b2, b3 = R2[0], S2[0], (R2 + S2)[0]

    # Compute coefficients
    M = Matrix(Fp2, [
        [a1*b1, a1, b1],
        [a2*b2, a2, b2],
        [a3*b3, a3, b3]])
    R, S, T = M.inverse() * vector(Fp2, [1,1,1])
    RD = R * M.determinant()
    da = (a1 - a2)*(a2 - a3)*(a3 - a1)
    db = (b1 - b2)*(b2 - b3)*(b3 - b1)

    R_inv = 1/R
    RD_inv = 1/RD
    s1, t1 = - da * RD_inv, db * RD_inv
    s2, t2 = - T * R_inv, -S * R_inv

    s1_inv = 1/s1
    a1_t = (a1 - s2) * s1_inv
    a2_t = (a2 - s2) * s1_inv
    a3_t = (a3 - s2) * s1_inv
    h = s1 * (x**2 - a1_t) * (x**2 - a2_t) * (x**2 - a3_t)

    # We need the image of (P_c, P) and (Q_c, Q) in J
    # The image of (P_c, P) is the image of P_c as a divisor on H
    # plus the image of P as a divisor on H.
    def isogeny(pair):
        # Argument members may be None to indicate the zero point.

        # The projection maps are:
        # H->C: (xC = s1/x²+s2, yC = s1 y)
        # so we compute Mumford coordinates of the divisor f^-1(P_c): a(x), y-b(x)
        Pc, P = pair
        if Pc:
            xPc, yPc = Pc.xy()
            uC = s1 * x**2 + s2 - xPc
            vC = Rx(yPc * s1_inv)

            DC = (uC, vC)

        # Same for E
        # H->E: (xE = t1 x² + t2, yE = t1 y/x^3)
        if P:
            xP, yP = P.xy()
            uE = (xP - t2) * x**2 - t1
            vE = yP * x**3 / t1
            vE = vE % uE

            DE = (uE, vE)

        if Pc and P:
            return affine_add(h, DC, DE)
        if Pc:
            return DC
        if P:
            return DE

    return h, isogeny

class RichelotCorr:
    """
    The Richelot correspondence between hyperelliptic
    curves y²=g1*g2*g3 and y²=h1*h2*h3=hnew(x)

    It is defined by equations:
        g1(x1) h1(x2) + g2(x1) h2(x2) = 0
    and y1 y2 = g1(x1) h1(x2) (x1 - x2)

    Given a divisor D in Mumford coordinates:
        U(x) = x^2 + u1 x + u0 = 0
        y = V(x) = v1 x + v0
    Let xa and xb be the symbolic roots of U.
    Let s, p by the sum and product (s=-u1, p=u0)

    Then on x-coordinates, the image of D is defined by equation:
        (g1(xa) h1(x) + g2(xa) h2(x))
      * (g1(xb) h1(x) + g2(xb) h2(x))
    which is a symmetric function of xa and xb.
    This is a non-reduced polynomial of degree 4.

    Write gred = g-U = g1*x + g0
    then gred(xa) gred(xb) = g1^2*p + g1*g0*s + g0^2
    and  g1red(xa) g2red(xb) + g1red(xb) g2red(xa)
       = 2 g11 g21 p + (g11*g20+g10*g21) s + 2 g10*g20

    On y-coordinates, the image of D is defined by equations:
           V(xa) y = Gred1(xa) h1(x) (xa - x)
        OR V(xb) y = Gred1(xb) h1(x) (xb - x)
    If we multiply:
    * y^2 has coefficient V(xa)V(xb)
    * y has coefficient h1(x) * (V(xa) Gred1(xb) (x-xb) + V(xb) Gred1(xa) (x-xa))
      (x-degree 3)
    * 1 has coefficient Gred1(xa) Gred1(xb) h1(x)^2 (x-xa)(x-xb)
                      = Gred1(xa) Gred1(xb) h1(x)^2 U(x)
      (x-degree 4)
    """
    def __init__(self, G1, G2, H1, H2, hnew):
        self.G1 = G1
        self.G2 = G2
        self.H1 = H1
        self.H11 = H1*H1
        self.H12 = H1*H2
        self.H22 = H2*H2
        self.hnew = hnew
        self.x = hnew.parent().gen()

    def map(self, D):
        "Computes (non-monic) Mumford coordinates for the image of D"
        U, V = D
        U = U.monic()
        V = V % U
        
        # Sum and product of (xa, xb)
        s, p = -U[1], U[0]
        
        # Compute X coordinates (non reduced, degree 4)
        g1red = self.G1 - U
        g2red = self.G2 - U
        g11, g10 = g1red[1], g1red[0]
        g21, g20 = g2red[1], g2red[0]

        # Precompute and reuse some multiplications
        tt = g11*p
        t0 = g11*tt + g11*g10*s + g10*g10
        t1 = g21*tt
        t2 = g10*g20

        # see above
        Px = t0 * self.H11 \
           + (t1 + t1 + (g11*g20 + g21*g10)*s + t2 + t2) * self.H12 \
           + (g21*(g21*p + g20*s) + g20*g20) * self.H22

        # Compute Y coordinates (non reduced, degree 3)
        v1, v0 = V[1], V[0]
        # coefficient of y^2 is V(xa)V(xb)
        Py2 = v1*v1*p + v1*v0*s + v0*v0
        
        # coefficient of y is h1(x) * (V(xa) Gred1(xb) (x-xb) + V(xb) Gred1(xa) (x-xa))
        # so we need to symmetrize:
        # V(xa) Gred1(xb) (x-xb)
        # = (v1*xa+v0)*(g11*xb+g10)*(x-xb)
        # = (v1*g11*p + v1*g10*xa + v0*g11*xb + v0*g10)*x
        # - xb*(v1*g11*p + v1*g10*xa + v0*g11*xb + v0*g10)
        # Symmetrizing xb^2 gives u1^2-2*u0

        # Precomputing some values we will reuse
        z0 = v1*g11*p
        z1 = v1*g10
        z2 = v0*g11
        z3 = v0*g10

        Py1 = (z0 + z0 + z3 + z3 + s*(z1 + z2))*self.x \
            - (p*(z1 + z1 - z2 - z2) + s*(z0 + s*z2 + z3))
        Py1 *= self.H1
        # coefficient of 1 is Gred1(xa) Gred1(xb) h1(x)^2 U(x)
        Py0 = self.H11 * U * t0

        # Now reduce the divisor, and compute Cantor reduction.
        # Py2 * y^2 + Py1 * y + Py0 = 0
        # y = - (Py2 * hnew + Py0) / Py1
        Py1inv = invert_mod_polynomial_quartic(Py1, Px)
        Py = (- Py1inv * (Py2 * self.hnew + Py0)) % Px

        Dx = ((self.hnew - Py*Py) // Px)
        Dy = (-Py) % Dx
        return (Dx, Dy)

def FromJacToJac(h, D11, D12, D21, D22):
    """
    Isogeny between J(H) -> J(H') with kernel ((D11, D12), (D21, D22))
    where D1, D2 are divisors on J(H) and Dij are the Mumford coordinates
    of the divisors Di
    """
    R,x = h.parent().objgen()

    D1 = (D11, D12)
    D2 = (D21, D22)

    G1, G2 = D1[0].monic(), D2[0].monic()
    G3, r3 = h.quo_rem(G1 * G2)
    assert r3 == 0

    delta = Matrix(G.padded_list(3) for G in (G1,G2,G3))
    # H1 = 1/det (G2[1]*G3[0] - G2[0]*G3[1])
    #        +2x (G2[2]*G3[0] - G3[2]*G2[0])
    #        +x^2(G2[1]*G3[2] - G3[1]*G2[2])
    # The coefficients correspond to the inverse matrix of delta.
    delta = delta.inverse()
    H1 = -delta[0][0]*x**2 + 2*delta[1][0]*x - delta[2][0]
    H2 = -delta[0][1]*x**2 + 2*delta[1][1]*x - delta[2][1]
    H3 = -delta[0][2]*x**2 + 2*delta[1][2]*x - delta[2][2]

    # New hyperelliptic curve H'
    hnew = H1*H2*H3

    # Class to compute the evaluation of the isogeny
    R = RichelotCorr(G1, G2, H1, H2, hnew)
    return hnew, R.map

def FromJacToProd(G1, G2, G3):
    """
    Construct the "split" isogeny from Jac(y^2 = G1*G2*G3)
    to a product of elliptic curves.

    This computation is the same as Benjamin Smith
    see 8.3 in http://iml.univ-mrs.fr/~kohel/phd/thesis_smith.pdf
    """
    h = G1*G2*G3
    x = h.parent().gen()

    M = Matrix(G.padded_list(3) for G in (G1,G2,G3))
    # Find homography
    u, v, w = M.right_kernel().gen()
    d = u/2
    ad, b = quadratic_roots(-v, w*d/2)
    a = ad/d

    # Apply transform G(x) -> G((a*x+b)/(x+d))*(x+d)^2
    # The coefficients of x^2 are M * (1, a, a^2)
    # The coefficients of 1 are M * (d^2, b*d, b^2)
    H11, H21, H31 = M * vector([1, a, a*a])
    H10, H20, H30 = M * vector([d*d, b*d, b*b])

    p1 = (H11*x+H10)*(H21*x+H20)*(H31*x+H30)
    p2 = (H11+H10*x)*(H21+H20*x)*(H31+H30*x)

    # We will need to map to actual elliptic curve
    p1norm = (x + H10*H21*H31)*(x + H20*H11*H31)*(x + H30*H11*H21)
    p2norm = (x + H11*H20*H30)*(x + H21*H10*H30)*(x + H31*H10*H20)
    E1 = EllipticCurve([0, p1norm[2], 0, p1norm[1], p1norm[0]])
    E2 = EllipticCurve([0, p2norm[2], 0, p2norm[1], p2norm[0]])

    def morphE1(x, y):
        # from y^2=p1 to y^2=p1norm
        return (H11*H21*H31*x, H11*H21*H31*y)
    
    def morphE2(x, y):
        # from y^2=p1 to y^2=p2norm
        return (H10*H20*H30*x, H10*H20*H30*y)
    
    # The morphisms are:
    # inverse homography:
    # H->H2: x, y => ((b-dx) / (x-a), y/(x-a)^3)
    # then H2->E1:(x,y) => (x^2,y)
    #   or H2->E2:(x,y) => (1/x^2,y/x^3)

    def isogeny(D):
        # To map a divisor, perform the change of coordinates
        # on Mumford coordinates
        U, V = D
        # apply homography
        # y = v1 x + v0 =>
        U_ = U[0] * (x+d)**2 + U[1]*(a*x+b)*(x+d) + U[2]*(a*x+b)**2
        V_ = V[0] * (x+d)**3 + V[1]*(a*x+b)*(x+d)**2
        V_ = V_ % U_
        v1, v0 = V_[1], V_[0]
        # prepare symmetric functions
        s = - U_[1] / U_[2]
        p = U_[0] / U_[2]
        # compute Mumford coordinates on E1
        # Points x1, x2 map to x1^2, x2^2
        U1 = x**2 - (s*s - 2*p)*x + p**2
        # y = v1 x + v0 becomes (y - v0)^2 = v1^2 x^2
        # so 2v0 y-v0^2 = p1 - v1^2 xH^2 = p1 - v1^2 xE1
        V1 = (p1 - v1**2 * x + v0**2) / (2*v0)
        # Reduce Mumford coordinates to get a E1 point
        V1 = V1 % U1
        U1red = (p1 - V1**2) // U1
        xP1 = -U1red[0] / U1red[1]
        yP1 = V1(xP1)

        # Same for E2
        # Points x1, x2 map to 1/x1^2, 1/x2^2
        U2 = x**2 - (s*s-2*p)/p**2*x + 1/p**2
        # yE = y1/x1^3, xE = 1/x1^2
        # means yE = y1 x1 xE^2
        # (yE - y1 x1 xE^2)(yE - y2 x2 xE^2) = 0
        # p2 - yE (x1 y1 + x2 y2) xE^2 + (x1 y1 x2 y2 xE^4) = 0
        V21 = x**2 * (v1 * (s*s-2*p) + v0*s)
        V20 = p2 + x**4 * (p*(v1**2*p + v1*v0*s + v0**2))
        # V21 * y = V20
        V21 = V21 % U2
        V21inv = invert_mod_polynomial_quadratic(V21, U2)
        V2 = (V21inv * V20) % U2

        # Reduce coordinates
        U2red = (p2 - V2**2) // U2
        xP2 = -U2red[0] / U2red[1]
        yP2 = V2(xP2)

        return E1(morphE1(xP1, yP1)), E2(morphE2(xP2, yP2))

    return isogeny, (E1, E2)

def _check_maximally_isotropic(P, Q, R, S, a):
    """
    Ensures that the kernel ((P, R), (Q, S)) is maximally
    isotropic by ensuring that:

    Isotropic:
    (P, Q) and (R, S) are valid torsion bases for E[2^a] and
    E'[2^a]

    Maximally isotropic:
    e(P, Q) * e(R, S) = 1

    Returns the scaled points 2^(a-1)(P, Q, R, S) on success
    for use in the gluing isogeny
    """
    # Scale to put points in E[2] and E'[2]
    k = 2**(a-1)
    P2 = k*P
    Q2 = k*Q
    R2 = k*R
    S2 = k*S

    two_torsion = (P2, Q2, R2, S2)

    # Ensure both bases are linearly independent
    if P2 == Q2 or R2 == S2:
        return False, None

    # Ensure all points have order dividing 2^a
    if not all((2*X).is_zero() for X in two_torsion):
        return False, None
        
    # Ensure all points have order exactly 2^a
    if any(X.is_zero() for X in two_torsion):
        return False, None

    # Ensure kernel is maximally isotropic
    ePQ = weil_pairing_pari(P, Q, 2*k) # 2k = 2^a
    eRS = weil_pairing_pari(R, S, 2*k)
    if ePQ*eRS != 1:
        return False, None

    return True, two_torsion

def split_richelot_chain(P, Q, R, S, a, strategy):
    r"""
    Given curves C, E and points (P, Q) \in E
                                 (R, S) \in E'
    
    Computes a (2^a, 2^a)-chain of isogenies with kernel:

    ker(Phi) = <(P, R), (Q, S)> : E x E' -> E'' x E''' 

    If the codomain splits, returns the chain and the codomain, 
    and None otherwise

    We expect this to fail only on malformed ciphertexts!
    """
    # We will output the final codomain as well as the isogeny
    # Phi : (E1, E2) -> J -> J -> ... -> J -> J -> (E3, E4)
    # Phi : Glue -> Richelot Chain -> Split
    richelot_chain = []

    # ========================= #
    #  Gluing step              
    #  (E, C) --> Jacobian      #
    # ========================= #

    check, (P2, Q2, R2, S2) = _check_maximally_isotropic(P, Q, R, S, a)
    if not check:
        return None, None

    # Compute the gluing map f : E x E' -> Jac(h)
    h, f = FromProdToJac(P2, Q2, R2, S2)

    # Evaluate the action of f on the points of
    # order 2^b
    D1 = f((P, R))
    D2 = f((Q, S))

    # Store the isogeny for later evaluations
    richelot_chain.append(f)

    # Bookkeeping variables
    strat_idx = 0
    level = [0]

    # We will precompute and push through elements
    # of the kernel, keep track of them with `ker`
    # and `kernel_elements`
    ker = (D1, D2)
    kernel_elements = [ker]

    # ======================================= #
    #  Middle Steps                           #
    #  Jacobian ---> Jacobian                 #
    # ======================================= #
    for i in range(1, a - 1):
        prev = sum(level)
        ker = kernel_elements[-1]
        while prev != (a - 1 - i):
            level.append(strategy[strat_idx])
            # Perform repeated doublings to compute
            # D_new = 2^strategy[strat_idx] D
            D1, D2 = ker
            D1 = affine_dbl_iter(h, D1[0], D1[1], strategy[strat_idx])
            D2 = affine_dbl_iter(h, D2[0], D2[1], strategy[strat_idx])
            ker = (D1, D2)

            # Update kernel elements and bookkeeping variables
            kernel_elements.append(ker)
            prev += strategy[strat_idx]
            strat_idx += 1

        # Compute the next step in the isogeny with the divisors D1, D2
        D1, D2 = ker
        h, f = FromJacToJac(h, D1[0], D1[1], D2[0], D2[1])
        
        # Update the chain of isogenies
        richelot_chain.append(f)

        # Remove elements from lists
        kernel_elements.pop()
        level.pop()

        # Push the kernel elements through the last step in the isogeny chain
        kernel_elements = [(f(D1), f(D2)) for D1, D2 in kernel_elements]

    # print(f"Middle steps took: {time.time() - t0:10f}")

    # Now we are left with a quadratic splitting: is it singular?
    D1, D2 = kernel_elements[-1]
    G1, G2 = D1[0], D2[0]
    G3, r3 = h.quo_rem(G1 * G2)
    assert r3 == 0

    delta = Matrix(G.padded_list(3) for G in (G1,G2,G3))
    if delta.determinant():
        # Determinant is non-zero, no splitting
        return None, None

    # ======================================= #
    #  Splitting Step                         #
    #  Jacobian ---> (E3, E4)                 #
    # ======================================= #
    f, h = FromJacToProd(G1, G2, G3)
    richelot_chain.append(f)
    return richelot_chain, h

def compute_richelot_chain(ker_Phi, b, strategy):
    """
    Helper function which takes as input a kernel for
    a (2^b,2^b)-isogeny and returns the isogeny which
    can be used to evaluate points.
    """
    # Unpack kernel generator
    glue_P1, glue_Q1, glue_P2, glue_Q2 = ker_Phi

    chain, _ = split_richelot_chain(
        glue_P1, glue_Q1, glue_P2, glue_Q2, b, strategy
    )
    if chain is None:
        raise ValueError("No splitting, ciphertext must be malformed")
    return chain

def evaluate_richelot_chain(Phi, X):
    """
    Evaluate the isogeny Phi : E1 x E2 -> E3 x E4 on a
    pair of points ((P1, P2), (Q1, Q2)) and obtain the
    points ((P3, P4), (Q3, Q4))
    """
    for f in Phi:
        X = f(X)
    return X

def compute_Li_from_richelot_chain(E0, Phi, L1_preimage, L2_preimage):
    """
    Given the (2^b, 2^b)-isogeny between elliptic products and pairs
    of points L1_preimage, L2_preimage compute L1 and L2 by 
    Phi(L1_preimage) and Phi(L2_preimage).

    Finally, use an isomorphism to map the target points back to the 
    starting curve E0. We don't know until after evaluating which
    point this will be, as there may be some twisting automorphism
    along the way.
    """
    L1_left, L1_right = evaluate_richelot_chain(Phi, L1_preimage)
    L2_left, L2_right = evaluate_richelot_chain(Phi, L2_preimage)

    # Pick the point with the right codomain to account for 
    # twisting
    if L1_left.curve().is_isomorphic(E0):
        iso_to_E0 = L2_left.curve().isomorphism_to(E0)
        L1 = iso_to_E0(L1_left)
        L2 = iso_to_E0(L2_left) 
    else:
        iso_to_E0 = L2_right.curve().isomorphism_to(E0)
        L1 = iso_to_E0(L1_right)
        L2 = iso_to_E0(L2_right) 

    return L1, L2