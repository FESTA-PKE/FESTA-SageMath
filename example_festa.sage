#!/usr/bin/env sage
"""
To run an example of FESTA with the 128-bit parameter set

```
sage example_festa.sage
```

If you want to test FESTA with the smaller, faster toy parameters, 
run the file with the --toy flag. For example:

```
sage example_festa.sage --toy
```

Additionally, higher security levels of FESTA can be run
with

```
sage example_festa.sage [--128, --192, --256]
```
"""

# Python imports
import sys
import time
from random import randint

# Local imports
from festa import FESTA
from parameters.params import parameter_sets
from utilities.utils import print_info, speed_up_sagemath

# Sage goes vroom!
# Sets up a monkey-patch to cache vector_space which
# dramatically helps with performance of the genus-two
# arithmetic and isogenies
speed_up_sagemath()

# Default is FESTA_128 with diagonal masking matrices
SECURITY = "128"
DIAG = True

for arg in sys.argv[1:]:
    if arg.lower() in ["--toy", "-t"]:
        SECURITY = "TOY"
    elif arg.lower() in ["--192", "-II"]:
        SECURITY = "192"
    elif arg.lower() in ["--256", "-V"]:
        SECURITY = "256"
    elif arg.lower() in ["-c", "--circulant"]:
        DIAG = False

NAME = "FESTA_" + SECURITY

# Initialise Alice and Bob
params = parameter_sets[NAME]
alice = FESTA(params, diag=DIAG)
bob = FESTA(params, diag=DIAG)

# Generate a random message
m = randint(0, 2**alice.lambda_security - 1)

# Start the test!
print_info(f"Running {NAME}")

# ============================== #
#              Keygen            #
# ============================== #
t0 = time.time()
alice.keygen()
keygen_time = time.time() - t0
print_info(f"Keygen took: {keygen_time:.3f} seconds")

# ============================== #
#           Encryption           #
# ============================== #
pk = alice.export_public_key()
print_info(f"Compressed public key: {len(pk)} bytes", banner="-")
t0 = time.time()
c = bob.encrypt(pk, m)
encrypt_time = time.time() - t0
print_info(f"Encrypt took: {encrypt_time:.3f} seconds")
print_info(f"Compressed ciphertext: {len(c)} bytes", banner="-")

# ============================== #
#           Decryption           #
# ============================== #
t0 = time.time()
new_m = alice.decrypt(c)
decrypt_time = time.time() - t0
print_info(f"Decrypt took: {decrypt_time:.3f} seconds")

# Now we know how long everything took, let's make
# sure it works!
assert m == new_m, "Decryption failed!"
