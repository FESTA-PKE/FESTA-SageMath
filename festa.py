# Python imports
from random import randint

# Sage Imports
from sage.all import EllipticCurve, GF, ceil, Matrix, ZZ, Zmod, inverse_mod

# pycryptodome imports
try:
    from Crypto.Hash import SHAKE128, SHAKE256
except ImportError as e:
    print("Error importing SHAKE from pycryptodome. Have you tried installing requirements?")
    print("sage -pip install -r requirements.txt")
    print(f"ImportError: {e}\n")
    exit()

# Local Imports
from supersingular import (
    torsion_basis,
    precompute_elligator_tables,
    entangled_torsion_basis,
    compute_canonical_kernel,
    montgomery_to_weierstrass_model
)
from isogenies_x_only import (
    isogeny_from_scalar_x_only,
    evaluate_isogeny_x_only,
    random_isogeny_x_only
)
from richelot_isogenies import (
    compute_richelot_chain, 
    compute_Li_from_richelot_chain
)
from compression import (
    compress_curve_and_two_torsion_basis,
    decompress_curve_and_two_torsion_basis
)
from utilities import (
    random_matrix,
    mask_torsion_points,
    canonical_matrix,
    integer_to_bytes,
    bytes_to_integer,
    BiDLP,
    DLP_power_two,
    weil_pairing_pari,
    optimised_strategy
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
    def __init__(self, params):
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

        # l^b = m1^2 T1 + m2^2 T2
        assert (
            self.l_power == m1**2 * self.d1 * self.dA1 + m2**2 * self.d2 * self.dA2
        )
        # T1*T2 = d1*dA*d2
        assert T1 * T2 == self.d1 * self.dA * self.d2

        # gi are scalars which are used for scaling the points
        # used in the (l,l)-chain kernel
        self.g1 = m2 * self.dA2 * self.d2
        self.g2 = m1 * self.d1

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
        self.E0.set_order((self.p+1)**2) 

        # For fast 2^k torsion basis and dlogs
        self.cofactor = (self.p+1)//self.l_power
        self.elligator_tables = precompute_elligator_tables(self.F)
        self.window = params["window"]
        
        # Pre-compute canonical torsion basis for E0[l^b] and E0[d1]
        # self.Pb, self.Qb = torsion_basis(self.E0, self.l_power)
        self.Pb, self.Qb = entangled_torsion_basis(self.E0, self.elligator_tables, self.cofactor)
        self.Pd1, self.Qd1 = torsion_basis(self.E0, self.d1)

        # Compute the byte length of integers for compression
        # and decompression
        self.p_byte_len = (self.p.nbits() + 7) // 8
        self.l_power_byte_len = (self.b + 7) // 8
        self.pk_bytes = 2*self.p_byte_len + 3*self.l_power_byte_len

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

        NOTE: this assumes that the masking matrices are diagonal 
        and that l=2
        """
        if self.l != 2:
            raise NotImplementedError("_extract_Mb is only implemented for l=2")

        # Extract out bytes from shake
        byte_len = (self.b + 7) // 8
        alpha = bytes_to_integer(shake.read(byte_len))

        # Ensure all values are odd
        mask = self.b % 8
        if mask:
            alpha = alpha >> mask
        alpha = (alpha >> 1 << 1) + 1
        beta = inverse_mod(alpha, self.l_power)

        # Construct the matrix in the ring Z/2^bZ
        R = Zmod(self.l_power)
        M = Matrix(R, 2, 2, [alpha, 0, 0, beta])
        return M

    def G(self, x, X):
        """
        Hash function: Z_d2 x M -> Z_d1
        Denoted G() in the paper
        """
        # Convert elements to bytes
        x_bytes = integer_to_bytes(x)

        # TODO: we could only extract the diagonal elements of this
        # matrix, but using all elements is fine
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
            EA, R, S, self.b, self.elligator_tables, self.cofactor, self.window, self.p_byte_len, self.l_power_byte_len
        )
        return ERS_bytes
    
    def decompress_public_key(self, pk_bytes):
        """
        The public key is decompressed by parsing the bytes to get byte
        representations of A_bytes, for the Montgomery coefficient of EA
        and the tuple of integers (aR, bR) and (aS, bS) such that 
        R = [a]P + [b]Q for the canonical basis of EA[l^b] 
        """
        assert len(pk_bytes) == self.pk_bytes, "Length of compressed public key is invalid"
        preimage_data = (self.Pb, self.Qb, self.dA)
        EA, R, S = decompress_curve_and_two_torsion_basis(
            self.F, pk_bytes, preimage_data, self.b, self.elligator_tables, self.cofactor, self.window, self.p_byte_len, self.l_power_byte_len
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
            E1, R1, S1, self.b, self.elligator_tables, self.cofactor, self.window, self.p_byte_len, self.l_power_byte_len
        )

        # Compress (E2, R2, S2)
        ERS2_bytes = compress_curve_and_two_torsion_basis(
            E2, R2, S2, self.b, self.elligator_tables, self.cofactor, self.window, self.p_byte_len, self.l_power_byte_len
        )

        return ERS1_bytes + ERS2_bytes
    
    def decompress_ciphertext(self, ct_bytes):
        """
        Decompression of the ciphertext is essentially two copies of the 
        decompression of the public key, where (Ei, Ri, Si) are decompressed
        from the supplied byte string
        """
        assert len(ct_bytes) == 2*self.pk_bytes, "Length of compressed ciphertext is invalid"

        # Split bytes into (E1, R1, S1) and (E2, R2, S2)
        ERS1_bytes = ct_bytes[:self.pk_bytes]
        ERS2_bytes = ct_bytes[self.pk_bytes:]

        # Grab preimages of points which we will
        # decompress
        _, R_preimage, S_preimage = self.pk
        preimage_data_1 = (self.Pb, self.Qb, self.d1)
        preimage_data_2 = (R_preimage, S_preimage, self.d2)

        # Decompress (E1, R1, S1)
        E1, R1, S1 = decompress_curve_and_two_torsion_basis(
            self.F, ERS1_bytes, preimage_data_1, self.b, self.elligator_tables, self.cofactor, self.window, self.p_byte_len, self.l_power_byte_len
        )

        # Decompress (E2, R2, S2)
        E2, R2, S2 = decompress_curve_and_two_torsion_basis(
            self.F, ERS2_bytes, preimage_data_2, self.b, self.elligator_tables, self.cofactor, self.window, self.p_byte_len, self.l_power_byte_len
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
        imPb, imQb = evaluate_isogeny_x_only(phi_1, self.Pb, self.Qb, self.l_power, self.d1)

        # Compute phi_2, the action of phi_2 on the public key points R, S 
        # E2 = EA / <P + [s1] Q>, EA[d2] = <P, Q>
        phi_2, E2 = isogeny_from_scalar_x_only(
            EA, self.d2, s2
        )
        imR, imS = evaluate_isogeny_x_only(phi_2, R, S, self.l_power, self.d2)

        # Mask torsion points with the matrix B
        R1, S1 = mask_torsion_points(B, imPb, imQb)
        R2, S2 = mask_torsion_points(B, imR, imS)

        return (E1, R1, S1, E2, R2, S2)

    # ========================================== #
    #   Helper functions for trapdoor_inverse    #
    # ========================================== #

    def _recover_si(self, L1, L2, imPd1b, imQd1b, PAd2, QAd2, PA_prime_d2, QA_prime_d2):
        """
        Helper function to recover the secret integers s1, s2 from the
        output of the (l,l)-chain during the inverse trapdoor function
        calculation
        """
        # Scale points to get order d1        
        P_d1 = self.l_power * imPd1b
        Q_d1 = self.l_power * imQd1b
        ePQ_d1 = weil_pairing_pari(P_d1, Q_d1, self.d1)
        s1 = compute_canonical_kernel(
            self.c2 * L1, self.c2 * L2, self.d1, basis=(P_d1, Q_d1), ePQ=ePQ_d1
        )

        # Compute torsion basis
        ePQ_d2 = weil_pairing_pari(PA_prime_d2, QA_prime_d2, self.d2)

        # Compute dlogs and reuse ePQ_d2
        a1, b1 = BiDLP(self.c1 * L1, PA_prime_d2, QA_prime_d2, self.d2, ePQ=ePQ_d2)
        a2, b2 = BiDLP(self.c1 * L2, PA_prime_d2, QA_prime_d2, self.d2, ePQ=ePQ_d2)

        # Compute points to recover kernel from
        imPAd2 = a1 * PAd2 + b1 * QAd2
        imQAd2 = a2 * PAd2 + b2 * QAd2

        # Recompute the kernel generators        
        s2 = compute_canonical_kernel(imPAd2, imQAd2, self.d2)

        return s1, s2

    def _recover_B(self, L1, L2, imPd1b, imQd1b):
        """
        Helper function to recover the masking matrix B from the
        output of the (l,l)-chain during the inverse trapdoor function
        calculation
        """
        # Remove odd order from L1 and L2
        phi_1_dual_R_scaled = self.d1d2 * L1
        phi_1_dual_S_scaled = self.d1d2 * L2

        # Remove d1 order from imPd1b
        phi_A1_Pb = self.clear_d1 * imPd1b
        phi_A1_Qb = self.clear_d1 * imQd1b
    
        # Recover the first coefficient of the matrix by
        # solving the discrete log
        alpha = DLP_power_two(
            phi_1_dual_R_scaled + phi_1_dual_S_scaled, 
            phi_A1_Pb, 
            phi_A1_Qb, 
            self.b,
            self.window
        )

        # Scale out by the constant
        # inv_scalar = (d1**2 * d2 * m1)^(-1) MOD l^b
        alpha = self.inv_scalar * alpha

        # We use unitary matricies where beta = 1/alpha
        # this saves one dlog
        beta  = inverse_mod(alpha, self.l_power)

        # Construct B from alpha and beta
        B = Matrix(Zmod(self.l_power), 2, 2, [alpha, 0, 0, beta])
        B = canonical_matrix(B)

        return B

    def _recover_si_and_B(self, L1, L2, imPd1b, imQd1b, PAd2, QAd2, PA_prime_d2, QA_prime_d2):
        """
        Wrapper function which recovers s1, s2 and the masking matrix B given
        the output of the (l,l)-chain during the inverse trapdoor function
        calculation
        """
        s1, s2 = self._recover_si(L1, L2, imPd1b, imQd1b, PAd2, QAd2, PA_prime_d2, QA_prime_d2)
        B = self._recover_B(L1, L2, imPd1b, imQd1b)
        return s1, s2, B

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
            raise ValueError("Cannot invert without first running keygen. \
                             Question: what was c encrypted with?")

        # Extract data from sk
        A, imPd1b, imQd1b, PAd2, QAd2, EA_prime, PA_prime_d2, QA_prime_d2 = self.sk

        # Extract data from c
        E1, R1, S1, E2, R2, S2 = c

        # Unmask R2, S2
        R2_prime, S2_prime = mask_torsion_points(A.inverse(), R2, S2)

        # Compute the torsion basis needed for the evaluations
        Pd1_1, Qd1_1 = torsion_basis(E1, self.d1)
        Pd2_2, Qd2_2 = torsion_basis(E2, self.d2)

        # ker(Φ) = <(P1, P2), (Q1, Q2)>
        #        = <[m2 d2 dA]R1, [m1 d1]R2), [m2 d2 dA]S1, [m1 d1]S2>
        glue_P1 = self.g1 * R1
        glue_Q1 = self.g1 * S1
        glue_P2 = self.g2 * R2_prime
        glue_Q2 = self.g2 * S2_prime

        # Prepare points to map through Phi
        L1_1, L1_2 = Pd1_1 + R1, Pd2_2
        L2_1, L2_2 = Qd1_1 + S1, Qd2_2

        # Gluing map assumes Weierstrass model, so map everything with 
        # an isomorphism from Montgomery to the short Weierstrass model
        # TODO: Can we modify the gluing map for Montgomery curves?
        E1, (glue_P1, glue_Q1, L1_1, L2_1) = montgomery_to_weierstrass_model(
            E1, (glue_P1, glue_Q1, L1_1, L2_1))
        E2, (glue_P2, glue_Q2, L1_2, L2_2) = montgomery_to_weierstrass_model(
            E2, (glue_P2, glue_Q2, L1_2, L2_2)
        )

        # Package the kernel generators of (l^b, l^b)-chain
        ker_Phi = (glue_P1, glue_Q1, glue_P2, glue_Q2)

        # Compute the (l^b, l^b)-chain
        Phi = compute_richelot_chain(ker_Phi, self.b, self.strategy)

        # Compute L1 and L2 by evaluating Phi on the following points
        # L1 = Phi(L1_1, L1_2)
        # L1 = Phi(L2_1, L2_2)
        L1_preimage = (L1_1, L1_2)
        L2_preimage = (L2_1, L2_2)
        L1, L2 = compute_Li_from_richelot_chain(
            EA_prime, Phi, L1_preimage, L2_preimage
        )

        # Recover the secrets s1, s2 and B from the data of L1 and L2
        s1, s2, B = self._recover_si_and_B(L1, L2, imPd1b, imQd1b, PAd2, QAd2, PA_prime_d2, QA_prime_d2)
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

        sk = (A, imPd1b, imQd1b, PAd2, QAd2, EA_prime)
        Where:
            A is the secret masking matrix
            <imPd1b, imQd1b> = EA_prime[d1 * l^b]
            EA_prime is the intermediate supersingular curve
            <PAd2, QAd2> = EA[d2]

            Precomputed torsion basis:
            <PA_prime_d2, QA_prime_d2> = EA_prime[d2]
        """
        # Random masking matrix
        A = random_matrix(self.l, self.b)

        # Compute a random isogeny phi_A of degree dA1 * dA2
        # The isogeny phi_A is equal to phi_A2_tilde \circ phi_A1_tilde
        # We do this in two steps

        # Step one: compute the isogeny phi_A1 : E0 -> EA_prime of degree dA1 together with phi_A1(Pb)  and phi_A1(Qb)
        # We also need the dual of this isogeny on E0_prime[d1*l^b]
        phi_A1, EA_prime = random_isogeny_x_only(
            self.E0, self.dA1
        )
        imPd1b, imQd1b = evaluate_isogeny_x_only(phi_A1, self.Pb + self.Pd1, self.Qb + self.Qd1, self.l_power*self.d1, self.dA1)

        # Remove the odd-order from the points
        imPb = self.clear_d1 * imPd1b
        imQb = self.clear_d1 * imQd1b

        # Step two: compute the isogeny phi_A2 : EA_prime -> EA of degree dA2 together with the images of 
        # phi_A2(imPb + Pd2_Aprime) and phi_A2(imQb + Qd2_Aprime)
        Pd2_Aprime, Qd2_Aprime = torsion_basis(EA_prime, self.d2)
        phi_A2, EA = random_isogeny_x_only(
            EA_prime, self.dA2
        )
        imPbd2, imQbd2 = evaluate_isogeny_x_only(phi_A2, imPb + Pd2_Aprime, imQb + Qd2_Aprime, self.l_power*self.d2, self.dA2)

        # Remove the odd-order from the points
        imPb = self.clear_d2 * imPbd2
        imQb = self.clear_d2 * imQbd2
        
        # PAd2 = phi_1_tilde \circ phi_A2_tilde (Pd1_tilde) , similarly QAd2
        # Remove the even order from the points
        PAd2 = self.clear_lb * imPbd2
        QAd2 = self.clear_lb * imQbd2

        # Mask torsion points using the matrix A
        R, S = mask_torsion_points(A, imPb, imQb)

        # Precompute torsion basis for later
        PA_prime_d2, QA_prime_d2 = torsion_basis(EA_prime, self.d2)

        # Set key pair
        self.pk = (EA, R, S)
        self.pk_compressed = self.compress_public_key(EA, R, S)
        self.sk = (A, imPd1b, imQd1b, PAd2, QAd2, EA_prime, PA_prime_d2, QA_prime_d2)

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
        R = random_matrix(self.l, self.b)

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