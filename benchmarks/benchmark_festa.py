"""
This file runs FESTA in the same way as `example_festa.sage` but 
additionally uses cProfile to output profiling statistics of 
keygen, encryption and decryption. This is used primarily to 
study which functions consume the most time in each step to allow
us to identify areas which deserve the most attention for optimisation

The benchmarks can be run with the command:

```
sage -python benchmarks/benchmark_festa.py
```

Running this with python has the additional benefit that we can
skip all debugging asserts using the -O flag:

```
sage -python -O benchmarks/benchmark_festa.py
```

NOTE:

If when you run the file this way you get an error about modules not being
found the most likely fix is that the $PYTHONPATH variable needs to be set.

As a one time fix, you can run:

```
export PYTHONPATH=${PWD}
```

Or, you can add this to the sagerc file, so it happens automatically.
"""

# Python imports
import cProfile
import pstats
import sys
import time
from random import randint

# Sage imports
from sage.all import randint, set_random_seed

# Local imports
from festa import FESTA
from parameters.params import parameter_sets
from utilities.utils import print_info, speed_up_sagemath

# Sage goes vroom!
speed_up_sagemath()

# Default is FESTA_128
SECURITY = "128"

for arg in sys.argv[1:]:
    if arg.lower() in ["--toy", "-t"]:
        SECURITY = "TOY"
    elif arg.lower() in ["--192", "-II"]:
        SECURITY = "192"
    elif arg.lower() in ["--256", "-V"]:
        SECURITY = "256"

PKE = FESTA
NAME = "FESTA_" + SECURITY

print_info(f"Benchmarking {NAME}")

N_KG = 10
N_Enc = 10
N_Dec = 10

# Initialise Alice and Bob
params = parameter_sets[NAME]
alice = PKE(params)
precomputation_time = time.time()
bob = PKE(params)
print_info(f"Precomputation took: {time.time() - precomputation_time:5f}")

# Generate a deterministic random message
set_random_seed(0)
m = randint(0, alice.lambda_security - 1)

# ========================= #
#      Start Profiling     #
# ========================= #
# Start the profiler
pr = cProfile.Profile()
pr.enable()

# ===================== #
#        Keygen         #
# ===================== #
keygen_time = time.time()
pr_keygen = cProfile.Profile()
pr_keygen.enable()
#
for _ in range(N_KG):
    alice.keygen()
#
pr_keygen.disable()
pr_keygen.dump_stats("festa_keygen.cProfile")
print_info(f"Keygen took: {(time.time() - keygen_time)/N_KG:.3f} seconds")
p = pstats.Stats("festa_keygen.cProfile")
p.strip_dirs().sort_stats("cumtime").print_stats(30)


# ===================== #
#        Encrypt        #
# ===================== #
pk = alice.export_public_key()
encrypt_time = time.time()
pr_encrypt = cProfile.Profile()
pr_encrypt.enable()
#
for _ in range(N_Enc):
    c = bob.encrypt(pk, m)
#
pr_encrypt.disable()
pr_encrypt.dump_stats("festa_encrypt.cProfile")
print_info(f"Encrypt took: {(time.time() - encrypt_time)/N_Enc:.3f} seconds")
p = pstats.Stats("festa_encrypt.cProfile")
p.strip_dirs().sort_stats("cumtime").print_stats(30)


# ===================== #
#        Decrypt        #
# ===================== #
decrypt_time = time.time()
pr_decrypt = cProfile.Profile()
pr_decrypt.enable()
#
for _ in range(N_Dec):
    new_m = alice.decrypt(c)
#
pr_decrypt.disable()
pr_decrypt.dump_stats("festa_decrypt.cProfile")
print_info(f"Decrypt took: {(time.time() - decrypt_time)/N_Dec:.3f} seconds")
p = pstats.Stats("festa_decrypt.cProfile")
p.strip_dirs().sort_stats("cumtime").print_stats(60)

print(f"Recovered message is equal to message: {m == new_m}")
