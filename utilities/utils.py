# Sage imports
from sage.all import (
    ZZ,
    GF,
    proof,
    cached_method,
)

# ====================================== #
#  Dumb check to see if there are zeros  #
# ====================================== #


def have_zero_value(*args):
    return any([x == 0 for x in args])


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
    return n.to_bytes(byte_len, "big")


def bytes_to_integer(b):
    """
    Represent bytes as an integer in big-endian
    """
    return ZZ(int.from_bytes(b, "big"))


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
    p = 2**127 - 1  # Arbitrary large prime
    to_patch = [GF(3), GF(3**2), GF(p), GF(p**2)]
    for x in to_patch:
        type(x).vector_space = cached_method(type(x).vector_space)


# ====================== #
#     Print Debugging    #
# ====================== #


def verbose_print(msg, verbose=False):
    if verbose:
        print(msg)


def print_info(str, banner="="):
    """
    Print information with a banner to help
    with visibility during debug printing
    """
    print(banner * 80)
    print(f"{str}".center(80))
    print(banner * 80)
