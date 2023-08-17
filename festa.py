# Python imports
from random import randint

# Sage Imports
from sage.all import EllipticCurve, GF, ceil, Matrix, ZZ, Zmod, inverse_mod

# pycryptodome imports
try:
    from Crypto.Hash import SHAKE128, SHAKE256
except ImportError as e:
    print(
        "Error importing SHAKE from pycryptodome. Have you tried installing requirements?"
    )
    print("sage -pip install -r requirements.txt")
    print(f"ImportError: {e}\n")
    exit()

# Local Imports
from utilities.supersingular import (
    torsion_basis,
    torsion_basis_with_pairing,
    precompute_elligator_tables,
    entangled_torsion_basis,
    compute_canonical_kernel,
    montgomery_curve_isomorphism,
)
from montgomery_isogenies.isogenies_x_only import (
    isogeny_from_scalar_x_only,
    evaluate_isogeny_x_only,
    random_isogeny_x_only,
)
from richelot_isogenies.richelot_isogenies import (
    compute_richelot_chain,
    compute_Li_from_richelot_chain,
)
from utilities.compression import (
    compress_curve_and_two_torsion_basis,
    decompress_curve_and_two_torsion_basis,
)
from utilities.masking import (
    random_matrix,
    mask_torsion_points,
    canonical_matrix,
)
from utilities.discrete_log import (
    BiDLP,
    DLP_power_two,
    BiDLP_power_two,
)
from utilities.strategy import optimised_strategy
from utilities.utils import (
    integer_to_bytes,
    bytes_to_integer,
)


class FESTA:
    """
    The FESTA PKE

    Exported functions:

    - keygen
    - export_public_key
    - encrypt
    - decrypt
    """

    def __init__(self, params, diag=True):
        """
        FESTA initialisation is used to set the main FESTA parameters
        along with auxiliary parameters which can be precomputed to save
        time during keygen, encryption and decryption.
        """
        # Extract parameters from parameter set
        self.p = params["p"]
        self.lambda_security = params["lambda_security"]

        # When we target higher security levels we use SHAKE256
        if self.lambda_security > 128:
            self.SHAKE = SHAKE256
        else:
            self.SHAKE = SHAKE128

        # (l^b, l^b)-Richelot chain
        self.l = params["l"]
        self.b = params["b"]
        self.l_power = self.l**self.b
        self.strategy = optimised_strategy(self.b - 1)
        self.Rlb = Zmod(self.l_power)

        # phi_i isogeny degrees
        self.d1 = params["d1"]
        self.dA = params["dA"]
        self.d2 = params["d2"]
        self.dA1 = params["dA1"]
        self.dA2 = params["dA2"]

        # Precompute scalars which are used for
        # clearing odd or even order from points
        d1_inverse = inverse_mod(self.d1, self.l_power)
        d2_inverse = inverse_mod(self.d2, self.l_power)
        self.clear_d1 = d1_inverse * self.d1
        self.clear_d2 = d2_inverse * self.d2
        self.clear_lb = inverse_mod(self.l_power, self.d2) * self.l_power

        # Extract mi, Ti to check params are valid
        m1, m2 = params["m1"], params["m2"]
        T1, T2 = params["T1"], params["T2"]

        N1 = m1**2 * self.d1 * self.dA1
        N2 = m2**2 * self.d2 * self.dA2

        # l^b = m1^2 T1 + m2^2 T2
        assert self.l_power == N1 + N2
        # T1*T2 = d1*dA*d2
        assert T1 * T2 == self.d1 * self.dA * self.d2

        # gi are scalars which are used for scaling the points
        # used in the (l,l)-chain kernel
        g1 = m2 * self.dA2 * self.d2
        g1_inv = inverse_mod(g1, self.l_power)
        self.g2 = m1 * self.d1 * g1_inv

        # ci are scalars used to recover si in inverse trapdoor
        self.c1 = -self.l_power * self.d1
        self.c2 = self.l_power * self.d2

        # Scalars used to recover the masking matrix B
        self.d1d2 = self.d1 * self.d2
        self.inv_scalar = inverse_mod(self.d1**2 * self.d2 * m1, self.l_power)

        # For the OAEP transform, we pad m with k 0-bits
        self.k = self.d1.nbits() - self.lambda_security - 1

        # Construct the FESTA base field and the starting
        # supersingular elliptic curve E0
        self.F = GF(self.p**2, name="z2", modulus=[1, 0, 1])
        self.E0 = EllipticCurve(self.F, [0, 6, 0, 1, 0])
        self.E0.set_order((self.p + 1) ** 2)

        # For fast 2^k torsion basis and dlogs
        self.cofactor = (self.p + 1) // self.l_power
        self.elligator_tables = precompute_elligator_tables(self.F)
        self.window = params["window"]

        # For magic square root in splitting
        self.N_constant = self.F(N1 + N2) / self.F(N1 - N2)

        # Pre-compute canonical torsion basis for E0[l^b] and E0[d1]
        # self.Pb, self.Qb = torsion_basis(self.E0, self.l_power)
        self.Pb, self.Qb = entangled_torsion_basis(
            self.E0, self.elligator_tables, self.cofactor
        )
        self.Pd1, self.Qd1 = torsion_basis(self.E0, self.d1, even_power=self.b)

        # Compute the byte length of integers for compression
        # and decompression
        self.p_byte_len = (self.p.nbits() + 7) // 8
        self.l_power_byte_len = (self.b + 7) // 8
        self.pk_bytes = 2 * self.p_byte_len + 3 * self.l_power_byte_len

        # Pick whether masking matrices are diagonal or circulant
        self.diag = diag

        # Internal values for key generation
        self.pk = None
        self.sk = None

    # ================================================ #
    #      Helper functions for the OAEP transform     #
    # ================================================ #

    def _extract_Zd(self, shake, d):
        """
        Extract bytes from a shake hash and reduce modulo d
        """
        # We sample lambda / 2 extra bits to remove sampling bias
        n_bits = ceil(self.lambda_security / 2 + d.nbits())

        # SHAKE API extracts bytes, not bits
        n_bytes = (n_bits + 7) // 8
        x_bytes = shake.read(n_bytes)

        # Reduce the result modulo d
        return ZZ(bytes_to_integer(x_bytes))

    def _extract_Mb(self, shake):
        """
        Extract bytes from a shake hash and reduce modulo 2^b
        Return a matrix in Mb
        """
        if self.l != 2:
            raise NotImplementedError("_extract_Mb is only implemented for l=2")

        # Extract out bytes from shake
        byte_len = (self.b + 7) // 8
        matrix_ele = bytes_to_integer(shake.read(byte_len))

        # Reduce to be exactly b-bits
        mask = self.b % 8
        if mask:
            matrix_ele = matrix_ele >> mask

        # Compute random matrix. Second element is generated canonically
        # given the first for both diagonal and circulant matrices.
        M = random_matrix(self.Rlb, ele=matrix_ele, diag=self.diag)

        return M

    def G(self, x, X):
        """
        Hash function: Z_d2 x M -> Z_d1
        Denoted G() in the paper
        """
        # Convert elements to bytes
        x_bytes = integer_to_bytes(x)

        # Convert matrix to bytes
        X_bytes = b"".join(integer_to_bytes(x) for x in X.list())

        # Initiate SHAKE128 with domain separator
        shake = SHAKE128.new(b"FESTA-G" + x_bytes + X_bytes)

        # Compute an element from SHAKE128
        return self._extract_Zd(shake, self.d1)

    def H(self, x):
        """
        Hash function: Z_d1 -> Z_d2 x M
        Denoted H() in the paper
        """
        # Convert Z_d1 to bytes for hashing
        x_bytes = integer_to_bytes(x)

        # Initiate SHAKE128 with domain separator
        shake = SHAKE128.new(b"FESTA-H" + x_bytes)

        # Compute elements from SHAKE128
        r = self._extract_Zd(shake, self.d2)
        R = self._extract_Mb(shake)

        return r, R

    # ================================================ #
    #    Compression and Decompression of PK and CT    #
    # ================================================ #

    def compress_public_key(self, EA, R, S):
        """
        The public key EA, R, S is compressed by:

        - Representing the elliptic curve EA by the Montgomery coefficient A
        - Representing points with integers a,b such that R = [a]P + [b]Q for
          P,Q the canonical basis of EA[l^b]
        """
        ERS_bytes = compress_curve_and_two_torsion_basis(
            EA,
            R,
            S,
            self.b,
            self.elligator_tables,
            self.cofactor,
            self.window,
            self.p_byte_len,
            self.l_power_byte_len,
        )
        return ERS_bytes

    def decompress_public_key(self, pk_bytes):
        """
        The public key is decompressed by parsing the bytes to get byte
        representations of A_bytes, for the Montgomery coefficient of EA
        and the tuple of integers (aR, bR) and (aS, bS) such that
        R = [a]P + [b]Q for the canonical basis of EA[l^b]
        """
        assert (
            len(pk_bytes) == self.pk_bytes
        ), "Length of compressed public key is invalid"
        preimage_data = (self.Pb, self.Qb, self.dA)
        EA, R, S = decompress_curve_and_two_torsion_basis(
            self.F,
            pk_bytes,
            preimage_data,
            self.b,
            self.elligator_tables,
            self.cofactor,
            self.window,
            self.p_byte_len,
            self.l_power_byte_len,
        )
        return EA, R, S

    def compress_ciphertext(self, ct):
        """
        Compression of the ciphertext is essentially two copies of the
        public key compression, where (Ei, Ri, Si) are compressed as above
        and the two byte strings are concatenated.
        """
        E1, R1, S1, E2, R2, S2 = ct

        # Compress (E1, R1, S1)
        ERS1_bytes = compress_curve_and_two_torsion_basis(
            E1,
            R1,
            S1,
            self.b,
            self.elligator_tables,
            self.cofactor,
            self.window,
            self.p_byte_len,
            self.l_power_byte_len,
        )

        # Compress (E2, R2, S2)
        ERS2_bytes = compress_curve_and_two_torsion_basis(
            E2,
            R2,
            S2,
            self.b,
            self.elligator_tables,
            self.cofactor,
            self.window,
            self.p_byte_len,
            self.l_power_byte_len,
        )

        return ERS1_bytes + ERS2_bytes

    def decompress_ciphertext(self, ct_bytes):
        """
        Decompression of the ciphertext is essentially two copies of the
        decompression of the public key, where (Ei, Ri, Si) are decompressed
        from the supplied byte string
        """
        if not self.pk:
            raise ValueError("A valid public key is required for compression")

        assert (
            len(ct_bytes) == 2 * self.pk_bytes
        ), "Length of compressed ciphertext is invalid, aborting"

        # Split bytes into (E1, R1, S1) and (E2, R2, S2)
        ERS1_bytes = ct_bytes[: self.pk_bytes]
        ERS2_bytes = ct_bytes[self.pk_bytes :]

        # Grab preimages of points which we will
        # decompress
        _, R_preimage, S_preimage = self.pk
        preimage_data_1 = (self.Pb, self.Qb, self.d1)
        preimage_data_2 = (R_preimage, S_preimage, self.d2)

        # Decompress (E1, R1, S1)
        E1, R1, S1 = decompress_curve_and_two_torsion_basis(
            self.F,
            ERS1_bytes,
            preimage_data_1,
            self.b,
            self.elligator_tables,
            self.cofactor,
            self.window,
            self.p_byte_len,
            self.l_power_byte_len,
        )

        # Decompress (E2, R2, S2)
        E2, R2, S2 = decompress_curve_and_two_torsion_basis(
            self.F,
            ERS2_bytes,
            preimage_data_2,
            self.b,
            self.elligator_tables,
            self.cofactor,
            self.window,
            self.p_byte_len,
            self.l_power_byte_len,
        )

        return (E1, R1, S1, E2, R2, S2)

    # ================================================ #
    #         FESTA Trapdoor Function Evaluation       #
    # ================================================ #

    def _trapdoor_eval(self, pk, s1, s2, B):
        """
        The FESTA trapdoor function following algorithm 5 in
        the paper.

        Input: the public key (EA, R, S)
               two integers s1, s2
               the masking matrix B

        Output: two elliptic curves E1 and E2 along with masked torsion bases
                <R1, S1> = E1[l^b] and <R2, S2> = E2[l^b]
        """
        # Extract data from pk
        EA, R, S = pk

        # Compute phi_1, the action of phi_1 on E0[2^b]
        # E1 = E0 / <P + [s1] Q>, E0[d1] = <P, Q>
        phi_1, E1 = isogeny_from_scalar_x_only(
            self.E0, self.d1, s1, basis=(self.Pd1, self.Qd1)
        )
        imPb, imQb = evaluate_isogeny_x_only(
            phi_1, self.Pb, self.Qb, self.l_power, self.d1
        )

        # Compute phi_2, the action of phi_2 on the public key points R, S
        # E2 = EA / <P + [s1] Q>, EA[d2] = <P, Q>
        phi_2, E2 = isogeny_from_scalar_x_only(EA, self.d2, s2)
        imR, imS = evaluate_isogeny_x_only(phi_2, R, S, self.l_power, self.d2)

        # Mask torsion points with the matrix B
        R1, S1 = mask_torsion_points(B, imPb, imQb)
        R2, S2 = mask_torsion_points(B, imR, imS)

        return (E1, R1, S1, E2, R2, S2)

    # ========================================== #
    #   Helper functions for trapdoor_inverse    #
    # ========================================== #

    def _recover_B_diag(self, R, P, Q):
        """
        Helper function to recover the elements of the
        masking matrix when it is known to be in diagonal
        form
        """
        # Recover the first coefficient of the matrix by
        # solving the discrete log
        alpha = DLP_power_two(R, P, Q, self.b, self.window)

        # Scale out by the constant
        # inv_scalar = (d1**2 * d2 * m1)^(-1) MOD l^b
        alpha = self.inv_scalar * alpha

        # We use unitary matricies where beta = 1/alpha
        # this saves one dlog
        beta = inverse_mod(alpha, self.l_power)

        # Construct diagonal B from alpha and beta
        return Matrix(self.Rlb, 2, 2, [alpha, 0, 0, beta])

    def _recover_B_circulant(self, R, P, Q):
        """
        Helper function to recover the elements of the
        masking matrix when it is known to be in circulant
        form
        """
        # Recover the first coefficient of the matrix by
        # solving the discrete log
        a, b = BiDLP_power_two(R, P, Q, self.b, self.window)
        # Scale out by the constant
        # inv_scalar = (d1**2 * d2 * m1)^(-1) MOD l^b
        a, b = self.inv_scalar * a, self.inv_scalar * b

        # Construct the circulant matrix B from a and b
        return Matrix(self.Rlb, 2, 2, [a, b, b, a])

    def _recover_B(self, phi_1_dual_R_scaled, phi_A1_Pb, phi_A1_Qb):
        """
        Helper function to recover the masking matrix B from the
        output of the (l,l)-chain during the inverse trapdoor function
        calculation
        """
        if self.diag:
            B = self._recover_B_diag(phi_1_dual_R_scaled, phi_A1_Pb, phi_A1_Qb)
        else:
            B = self._recover_B_circulant(phi_1_dual_R_scaled, phi_A1_Pb, phi_A1_Qb)

        return canonical_matrix(B)

    def _recover_si_and_B(self, L1, L2, im_basis_bd1_E0, im_basis_d2_EA):
        """
        Wrapper function which recovers s1, s2 and the masking matrix B given
        the output of the (l,l)-chain during the inverse trapdoor function
        calculation
        """
        # Scale point to get order 2^b, d1 and d2 points
        # on EA_prime, from the image of the (2,2)-isogeny
        imPd1d2 = self.l_power * L1
        imQd1d2 = self.l_power * L2

        imPb = self.d1d2 * L1  # order 2^b on EA'
        imPd1 = self.d2 * imPd1d2  # order d1 on EA'
        imQd1 = self.d2 * imQd1d2  # order d1 on EA'
        imPd2 = self.d1 * imPd1d2  # order d2 on EA'
        imQd2 = self.d1 * imQd1d2  # order d2 on EA'

        # Image of generators of E0[l^b + d1] on EA_prime
        imPd1b, imQd1b = im_basis_bd1_E0

        # Scale points to get order 2^b points on EA_prime
        # from the canonical torsion basis
        phi_A1_Pb = self.clear_d1 * imPd1b
        phi_A1_Qb = self.clear_d1 * imQd1b

        # Scale points to get order d1 points on EA_prime
        # from the canonical torsion basis
        phi_A1_Pd1 = self.l_power * imPd1b
        phi_A1_Qd1 = self.l_power * imQd1b

        # Image of generators of EA[d2] on EA_prime
        phi_A2_dual_Pd2, phi_A2_dual_Qd2 = im_basis_d2_EA

        s1 = compute_canonical_kernel(
            imPd1, imQd1, self.d1, basis=(phi_A1_Pd1, phi_A1_Qd1)
        )
        s2 = compute_canonical_kernel(
            imPd2, imQd2, self.d2, basis=(phi_A2_dual_Pd2, phi_A2_dual_Qd2)
        )

        B = self._recover_B(imPb, phi_A1_Pb, phi_A1_Qb)
        return s1, s2, B

    @staticmethod
    def recover_correct_image(L1, L2, EA_prime):
        L1_left, L1_right = L1.points()
        L2_left, L2_right = L2.points()

        # Pick the point with the right codomain to account for
        # twisting
        if L1_left.curve().j_invariant() == EA_prime.j_invariant():
            L1, L2 = montgomery_curve_isomorphism(
                L1_left.curve(), EA_prime, L1_left, L2_left
            )
        else:
            L1, L2 = montgomery_curve_isomorphism(
                L1_right.curve(), EA_prime, L1_right, L2_right
            )

        return L1, L2

    # =========================================== #
    #     FESTA Inverse of Trapdoor Function      #
    # =========================================== #

    def _trapdoor_inverse(self, c):
        """
        The FESTA inverted trapdoor function following algorithm 7 in
        the FESTA paper.  Given a tuple c and the secret key, sk, this
        function recovers the integers s1, s2 and the masking matrix B
        """
        if self.sk is None:
            raise ValueError(
                "Cannot invert without first running keygen. \
                             Question: what was c encrypted with?"
            )

        # Extract data from sk
        (
            A,  # Secret scaling matrix
            EA_prime,  # Middle curve E0 -> EA_prime -> EA
            im_basis_bd1_E0,  # phi_A1(E0[l^b] + E0[d1])
            im_basis_d2_EA,  # phi_A2_dual(EA[d2])
        ) = self.sk

        # Extract data from c
        E1, R1, S1, E2, R2, S2 = c

        # Unmask R2, S2
        R2_prime, S2_prime = mask_torsion_points(A.inverse(), R2, S2)

        # Compute the torsion basis needed for the evaluations
        Pd1_1, Qd1_1 = torsion_basis(E1, self.d1, even_power=self.b)
        Pd2_2, Qd2_2 = torsion_basis(E2, self.d2, even_power=self.b)

        # ker(Φ) = <(P1, P2), (Q1, Q2)>
        #        = <[m2 d2 dA]R1, [m1 d1]R2), [m2 d2 dA]S1, [m1 d1]S2>
        glue_P1 = R1
        glue_Q1 = S1
        glue_P2 = self.g2 * R2_prime
        glue_Q2 = self.g2 * S2_prime

        # Prepare points to map through Phi
        L1_1, L1_2 = Pd1_1 + R1, Pd2_2
        L2_1, L2_2 = Qd1_1 + S1, Qd2_2

        # Package the kernel generators of (l^b, l^b)-chain
        ker_Phi = (glue_P1, glue_Q1, glue_P2, glue_Q2)

        # Compute the (l^b, l^b)-chain
        Phi = compute_richelot_chain(ker_Phi, self.b, self.N_constant, self.strategy)

        # Compute L1 and L2 by evaluating Phi on the following points
        # L1 = Phi(L1_1, L1_2)
        # L1 = Phi(L2_1, L2_2)
        L1_preimage = (L1_1, L1_2)
        L2_preimage = (L2_1, L2_2)
        L1, L2 = compute_Li_from_richelot_chain(EA_prime, Phi, L1_preimage, L2_preimage)

        # Recover the secrets s1, s2 and B from the data of L1 and L2
        s1, s2, B = self._recover_si_and_B(L1, L2, im_basis_bd1_E0, im_basis_d2_EA)
        return s1, s2, B

    # =============================== #
    #          FESTA KeyGen           #
    # =============================== #

    def keygen(self):
        """
        Compute the key pair

        pk = (EA, R, S) in compressed form
        Where:
           EA is the image of the secret isogeny phi_A
           <R, S> = EA[l^b] is the masked torsion basis

        sk = (A, EA_prime, im_basis_bd1_E0, im_basis_d2_EA)
        Where:
            A is the secret masking matrix
            EA_prime is the intermediate supersingular curve

            im_basis_bd1_E0 is the image of E0[l^b + d1] on EA_prime
            im_basis_d2_EA is the image of EA[d2
            ] on EA_prime
        """
        # Random masking matrix
        A = random_matrix(self.Rlb, diag=self.diag)

        # Compute a random isogeny phi_A of degree dA1 * dA2
        # The isogeny phi_A is equal to phi_A2_tilde \circ phi_A1_tilde
        # Step one: compute the isogeny phi_A1 : E0 -> EA_prime of degree dA1
        # together with phi_A1(Pb) and phi_A1(Qb) and phi_A1(Pd1) and phi_A1(Qd1)
        phi_A1, EA_prime = random_isogeny_x_only(self.E0, self.dA1, even_power=self.b)
        imPd1b, imQd1b = evaluate_isogeny_x_only(
            phi_A1,
            self.Pb + self.Pd1,
            self.Qb + self.Qd1,
            self.l_power * self.d1,
            self.dA1,
        )

        # Remove the odd-order from the points to recover
        # phi_A1(Pb) and phi_A1(Qb)
        imPb = self.clear_d1 * imPd1b
        imQb = self.clear_d1 * imQd1b

        # Step two: compute the isogeny phi_A2 : EA_prime -> EA of degree dA2
        # together with the images of
        # phi_A2(imPb + PA_prime_d2) and phi_A2(imQb + QA_prime_d2)
        PA_prime_d2, QA_prime_d2 = torsion_basis(EA_prime, self.d2, even_power=self.b)
        phi_A2, EA = random_isogeny_x_only(EA_prime, self.dA2)
        imPbd2, imQbd2 = evaluate_isogeny_x_only(
            phi_A2,
            imPb + PA_prime_d2,
            imQb + QA_prime_d2,
            self.l_power * self.d2,
            self.dA2,
        )

        # Remove the odd-order from the points
        # to recover phi_A2 \circ phi_A1(Pb) and phi_A2 \circ phi_A1(Qb)
        imPb = self.clear_d2 * imPbd2
        imQb = self.clear_d2 * imQbd2

        # Remove the even order from the points
        # to recover phi_A2 on EA_prime[d2]
        imPAd2 = self.clear_lb * imPbd2
        imQAd2 = self.clear_lb * imQbd2

        # We now want the image of the canonical basis EA[d2] on EA_prime
        # phi_A2_dual(PA_d2) = [a1] PA_prime_d2 + [b1] QA_prime_d2
        # phi_A2_dual(QA_d2) = [a2] QA_prime_d2 + [b2] QA_prime_d2
        #
        # phi_A2(phi_A2_dual(PA_d2)) = [dA2] PA_d2 = [dA2] ([a1] phi_A2(PA_prime_d2) + [b1] phi_A2(QA_prime_d2))
        # So we can compute ai, bi from
        #
        # PA_d2 = [a1] imPAd2 + [b1] imQAd2

        # Compute canonical torsion basis on EA for the dlogs
        PA_d2, QA_d2, ePQA = torsion_basis_with_pairing(EA, self.d2, even_power=self.b)

        a1, b1 = BiDLP(PA_d2, imPAd2, imQAd2, self.d2, ePQ=ePQA)
        a2, b2 = BiDLP(QA_d2, imPAd2, imQAd2, self.d2, ePQ=ePQA)

        imPA_prime_d2 = a1 * PA_prime_d2 + b1 * QA_prime_d2
        imQA_prime_d2 = a2 * PA_prime_d2 + b2 * QA_prime_d2

        # Now we can recover the scalars from encryption by computing
        # canonical bases on EA_prime with the following images:
        im_basis_bd1_E0 = (imPd1b, imQd1b)
        im_basis_d2_EA = (imPA_prime_d2, imQA_prime_d2)

        # Mask torsion points using the matrix A
        R, S = mask_torsion_points(A, imPb, imQb)

        # Set key pair
        self.pk = (EA, R, S)
        self.pk_compressed = self.compress_public_key(EA, R, S)
        self.sk = (
            A,  # Secret scaling matrix
            EA_prime,  # Middle curve E0 -> EA_prime -> EA
            im_basis_bd1_E0,  # phi_A1(E0[l^b] + E0[d1])
            im_basis_d2_EA,  # phi_A2_dual(EA[d2])
        )

        return

    def export_public_key(self):
        """
        Helper function to act as a "transport" for the
        public key to a user for encryption
        """
        return self.pk_compressed

    # =============================== #
    #          FESTA Encrypt           #
    # =============================== #

    def encrypt(self, pk_compressed, m):
        """
        OAEP transform for the trapdoor function following
        Algorithm 3 of the FESTA paper.
        """
        if m > 2**self.lambda_security:
            raise ValueError("Message is too long")

        # Create random tokens for OAEP
        r = randint(0, self.d2 - 1)
        R = random_matrix(self.Rlb, diag=self.diag)

        # Pad out the bottom bits of m_prime
        m_prime = m << (self.k)
        assert m_prime < self.d1, "Padding was performed incorrectly..."

        # s is an element of Z_d1
        s = Zmod(self.d1)(m_prime + self.G(r, R))

        # Extract x and X from s
        x, X = self.H(s)

        # t is an element of Z_d2
        t = Zmod(self.d2)(x + r)
        T = X * R
        T = canonical_matrix(T)

        # Perform encryption with one way function
        pk = self.decompress_public_key(pk_compressed)
        c = self._trapdoor_eval(pk, s, t, T)
        c_bytes = self.compress_ciphertext(c)

        return c_bytes

    # =============================== #
    #          FESTA Decrypt           #
    # =============================== #

    def decrypt(self, c_bytes):
        """
        OAEP transform for the inverse trapdoor function following
        Algorithm 4 of the FESTA paper.
        """
        # Decompress the ciphertext
        c = self.decompress_ciphertext(c_bytes)

        # Recover the secret values with the inverse function
        s, t, T = self._trapdoor_inverse(c)

        # Recover m from various values
        x, X = self.H(s)
        r = Zmod(self.d2)(t - x)
        R = X.inverse() * T
        R = canonical_matrix(R)

        # Check m_prime has 0 as the k-lowest bits
        m_prime = Zmod(self.d1)(s - self.G(r, R))
        m_prime = int(m_prime)
        if m_prime & ((1 << self.k) - 1):
            print(f"Non zero bits in the padding")
            return False

        # Extract m by removing the bottom k bits
        m = m_prime >> self.k

        return m
