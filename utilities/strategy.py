# ================================================ #
#     Compute optimised strategy for (2,2)-chain   #
# ================================================ #


def optimised_strategy(n, mul_c=0.175):
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
    but while these costs are so noisy from SageMath benchmarks we
    find there's little use in overoptimising here to save a few ms.

    Once a proper, constant-time implementation has been developed
    we can re-address these costs and find the true optimal strategy
    of the (2,2)-chain for decryption.
    """

    eval_c = 1.000
    mul_c = mul_c

    S = {1: []}
    C = {1: 0}
    for i in range(2, n + 1):
        b, cost = min(
            ((b, C[i - b] + C[b] + b * mul_c + (i - b) * eval_c) for b in range(1, i)),
            key=lambda t: t[1],
        )
        S[i] = [b] + S[i - b] + S[b]
        C[i] = cost

    return S[n]
