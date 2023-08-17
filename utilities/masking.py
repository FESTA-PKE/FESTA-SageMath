from sage.all import Matrix, ZZ

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
    b = 4 * beta
    aa = R(b * b + 1)
    a = aa.sqrt()  # TODO: fast sqrt mod 2^k?
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
    R = a00 * P + a01 * Q
    S = a10 * P + a11 * Q

    return R, S
