# FESTA 

A proof of concept implementation of the isogeny-based PKE FESTA[^1], accompanying 
the research paper 
[FESTA: Fast Encryption from Supersingular Torsion Attacks](https://eprint.iacr.org/2023/660)
by Andrea Basso, Luciano Maino and Giacomo Pope.

[^1]: *party* in Italian :it: :tada:

## FESTA-SageMath

The class `FESTA` implements the PKE and has three main functions:

- `FESTA.keygen()`: generate a keypair for a user
- `FESTA.encrypt(pk, m)`: given a user's public key `pk` encrypt a message `m`.
- `FESTA.decrypt(c)`: given a ciphertext recover the message `m`

Note a user's keypair is stored within the `FESTA` object when the `keygen()` method is called, so 
`sk` does not need to be passed as a parameter to `decrypt(c)`.

### Example Usage

```python
from festa import FESTA
from parameters import festa_params_128

# Initialisation 
alice = FESTA(festa_params_128)
bob   = FESTA(festa_params_128)

# Keygen
alice.keygen()

# Encryption
pk = alice.export_public_key()
m = randint(0, 2**128 - 1)
c = bob.encrypt(pk, m)

# Decryption
assert alice.decrypt(c) == m
```

**SageMath version**: This code was developed and tested using SageMath version 9.8, the most recent stable version (at the time of writing).

**Requirements**: for the OAEP transform, we use `SHAKE` imported from `pycryptodome` to extract random bytes. This can be installed using:
```
sage -pip install -r requirements.txt
```

## Project Overview

FESTA has been implemented as a class in the file [`festa.py`](festa.py). This file contains the full 
implementation of the FESTA trapdoor evaluation and inversion, as well as the functions for the OAEP 
transform to build the PKE and functions to compress and decompress both the public key and ciphertext. 
An effort has been made to keep the comments and doc-strings verbose. The FESTA parameters for all security 
levels, including small toy parameters for testing, are defined in [`parameters.py`](parameters.py). The file
[`parameter_generator.py`](parameter_generator.py) has the code which we used for searching for optimal parameters
targeting all security levels.

For the remainder of this section, we break down the auxiliary files which contain helper functions for various steps of the FESTA protocol.

The file [`supersingular.py`](supersingular.py) contains functions to compute points of a given order $D$ and the generation of canonical 
torsion bases $E[D] = \langle P, Q \rangle$. We include the optimised `entangled_torsion_basis()` following 
[https://ia.cr/2017/1143](https://ia.cr/2017/1143) for computing the torsion basis of $E[2^b]$.
Additionally, this file contains `compute_canonical_kernel()` as described in algorithm six of the 
FESTA paper. 

The file [`compression.py`](compression.py) contains helper functions for the main compression and decompression algorithms 
of FESTA, such as the compression of a curve to a single field element and the compression of a point $P \in E[2^b]$ to a bytestring 
representing three integers $(a,b,c)$ in $\mathbb{Z}/2^b\mathbb{Z}$. We need only three integers as our scaling matricies are unitary
and compatibility of the Weil pairing with isogenies means we can always efficiently recover the fourth. This is handled by the function
`recover_lost_scalar()`.

Our isogenies between elliptic curves are implemented using x-only formula. In [`kummer_line.py`](kummer_line.py) we define classes for 
the Kummer line of Montgomery curves (`KummerLine`) and points on the Kummer line (`KummerPoint`). Essentially, these are a standard 
implementation of x-only Montgomery curves with projective $(X : Z)$ coordinates. We additionally implement x-only isogenies between 
these Kummer lines in [`kummer_isogeny.py`](kummer_isogeny.py) using the Costello-Hisil-Renes formula for small prime degree isogenies
and the VéluSqrt formula for the large prime degree isogenies. The main class `KummerLineIsogeny` accepts a kernel generator of composite 
order and handles factorisation of the isogeny into prime-degree parts.

In several places in the FESTA implementation, we cannot work with x-only formula as we need the affine coordinates for addition of various 
points and the computation of the isogeny between elliptic products. As such, the use of x-only formula for the isogenies requires that it's 
possible to lift the images from the Kummer line back to the curve. 
The file [`isogenies_x_only.py`](isogenies_x_only.py) contains helper functions which do exactly this. They take as input Montgomery curves 
and torsion bases, maps everything to the Kummer line, computes and evaluates isogenies and finally lifts everything to the image curve.

There are two files which are used to implement the $(2^b,2^b)$ isogeny between elliptic products. In [`divisor_arithmetic.py`](divisor_arithmetic.py), 
we implement efficient addition and doubling laws for divisors of Jacobians of genus two hyperelliptic curves. These have been derived by generalising 
a result of [Costello and Lauter](https://eprint.iacr.org/2011/306) to recover algorithms which only require base field operations, avoiding the need 
for slower arithmetic in the polynomial ring. The $(2^b,2^b)$ isogeny itself is implemented in [`richelot_isogenies.py`](richelot_isogenies.py) and has 
been adapted from the previous work on implementing both the [Castryck-Decru](https://github.com/jack4818/Castryck-Decru-SageMath) and 
[MMPPW](https://github.com/Breaking-SIDH/direct-attack) SIDH attacks.

Finally, the file [`utilities.py`](utilities.py) contains a little of everything else. This includes functions to help with the scaling matrices, 
solving discrete logarithms using Weil/Tate pairings, computing whether an element has order exactly $D$ and a function to compute optimal strategies 
based off the SIDH sparse isogeny chains. We additionally include some small optimised formula for computing the square root in $\mathbb{F}\_{p^2}$
using that $p \equiv 3 \pmod 4$. We also implement a function which computes $f^{-1} \pmod g$ where $f,g$ are univariate polynomials 
in $\mathbb{F}\_{p^2}[X]$ and $g$ is of degree two or four. 
Our implementation solves this using linear algebra rather than `xgcd` and is about 2-3x faster
than the generic implementation of SageMath.

## Command Line FESTA

If you want to play around with FESTA and compare run-times, you can use the
file `example_festa.sage` with the following arguments:

```
sage example_festa.sage [--128, --192, --256, --toy, --circulant]
```

- By default, the 128-bit security parameters are selected. To access other parameters:
  - The flag `--192` selects the parameters aiming for 192-bit security 
  - The flag `--256` selects the parameters aiming for 256-bit security 
  - The flag `--toy` selects small toy parameters suitable for debugging
- By default, the masking matrices are diagonal, unitary, invertible matrices.
  - The flag `--circulant` selected the matrices to be circulant, unitary, invertible matrices instead.

#### Example Output

```
User: % sage example_festa.sage
================================================================================
                               Running FESTA_128                                
================================================================================
================================================================================
                           Keygen took: 4.853 seconds
================================================================================
--------------------------------------------------------------------------------
                        Compressed public key: 561 bytes
--------------------------------------------------------------------------------
================================================================================
                          Encrypt took: 3.513 seconds
================================================================================
--------------------------------------------------------------------------------
                       Compressed ciphertext: 1122 bytes                        
--------------------------------------------------------------------------------
================================================================================
                          Decrypt took: 10.102 seconds
================================================================================
```

**Note**: the above output was generated using a single performance core of an Apple M1 PRO CPU, clocked at 3.2 GH.
