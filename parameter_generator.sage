import sys

# Global params
c_min = 2
c_max = 10 # checks all c's in [c_min, ..., c_max]
b_max = 500

proof.all(False)

# lam = 128
lam = 64
smooth_bound = 2^22

for arg in sys.argv[1:]:
    if arg.lower() in ["--192", "--II"]:
        lam = 192
        smooth_bound = 2^24

    elif  arg.lower() in ["--256", "--V"]:
        lam = 256
        smooth_bound = 2^26
        b_max = 3000

B_smooth = prod(primes(smooth_bound))

def find_d1(T1):
    
    d = 1
    dlog = 0

    for p, e in factor(T1):
        d *= p^e
        dlog += e*math.log(p, 2)
        if dlog > 2*lam:
            break

    if dlog < 2*lam:
        return -1

    # we go over once again to remove
    # smaller factors in case we overshot
    for p, e in factor(d)[::-1]:
        if dlog - e*math.log(p, 2) > 2*lam and e%2==0:
            d = d // p^e
            dlog -= e*math.log(p, 2)

    assert dlog > 2*lam

    return d

def find_dAd2(T2dash, res):

    dA = 1
    dAlog = 0
    d2 = 1
    d2log = 0

    for p, e in factor(T2dash):
        d2 *= p^e
        d2log += e*math.log(p, 2)

        if d2log > 2*lam:
            break
    
    if d2log < 2*lam:
        return -1,-1,-1
    
    for p, e in factor(d2)[::-1]:
        if d2 % p^e != 0:
            continue
            
        if d2log - e*math.log(p, 2) > 2*lam and e%2==0:
            d2 = d2 // p^e
            d2log -= e*math.log(p, 2)
            
    T2dash = (T2dash//d2)*res
    
    for p, e in factor(T2dash):
        dA *= p^e
        dAlog += e*math.log(p, 2)

        if dAlog > 2*lam:
            break
    
    if dAlog < 2*lam:
        return -1,-1,-1
    
    for p, e in factor(d2)[::-1]:

            
        if dAlog - e*math.log(p, 2) > 2*lam and e%2==0:
            dA = d2 // p^e
            dAlog -= e*math.log(p, 2)

    dAsf = gcd(dA,B_smooth)
    return dA, d2, dAsf


def make_prime(p):
    for i in range(50000):
        f = 2*i + 1
        if (p*f - 1).is_prime():
            return p*f - 1, f

def log_factored(p_):
    logp = 0
    for p, e in factor(p_):
        logp += e*math.log(p, 2)
    
    return logp

def B_smooth_part(n):
    """
    Compute the B-smooth part of an
    integer
    """
    g = gcd(n, B_smooth)
    b_smooth = g
    while g != 1:
        n = n // g
        g = gcd(n, g)
        b_smooth = g * b_smooth 

    return b_smooth

def update_globals():
    global b,n,N
    b += 1
    n *= 2
    N *= 4

solutions = {}
lenSol = 0

for c in range(c_min, c_max+1):
    print(f"Checking c = {c}...", end=" ")
    # QF used to solve eqn.
    P = 2**c - 1
    QF = gp.Qfb(1, 0, P)

    # Set bounds from c. 
    b = c // 2
    n = 2**b
    N = 2**(2*b)

    for _ in range(b_max):
        # Find a solution of 
        # x^2 + P*y^2 = M = 2^(2b)
        sol = QF.qfbsolve(N)

        # print(sol)
        # No solution found
        if not sol:
            # Update variables
            update_globals()
            continue


        # Compute the B-smooth parts
        x, y = list(map(ZZ, sol))



        T1 = B_smooth_part(n - x)
        T2 = B_smooth_part(n + x)
        m12 = (n - x) // T1
        m22 = (n + x) // T2

        if not (m12.is_square() and m22.is_square()):
            update_globals()
            break
            
        m1 = m12.sqrt()
        m2 = m22.sqrt()

        for j in range(2):
            if j == 1:
                T1, T2 = T2, T1
                m1, m2 = m2, m1
            
            if T1*T2>6*lam:
                d1 = find_d1(T1)
                res = T1//d1
                dA, d2, dAsf = find_dAd2(T2,res)
                dA1 = gcd(res, dA)
                dA2 = gcd(T2,dA)
                
                if d1 == -1 or dA == -1 or d2 == -1:
                    continue

                logp = (d1 * d2* dAsf).nbits() + b # approximately, we don't compute a prime at this stage

                m1 *= sqrt(T1 // (d1*dA1))
                m2 *= sqrt(T2 // (dA2 *d2) ) 

                T1 = d1*dA1
                T2 = d2*dA2
                
                if gcd(m1,T1)==1 and gcd(m2,T2)==1 and (m1 in ZZ) and (m2 in ZZ):
                    solutions[(b, c, T1, T2, m1, m2, d1, dA, d2, dAsf)] = logp
                
        update_globals()

    print(f"(found {len(solutions) - lenSol} solutions)")
    lenSol = len(solutions)

# In case we find nothing
if len(solutions) == 0:
    exit("No solutions found")

print(f"\n\n\nBest solution")



sorted_solutions = sorted(solutions.items(), key=lambda item: item[1])

for best, _ in sorted_solutions[:5]:
    b, c, T1, T2, m1, m2, d1, dA, d2, dAsf = best


    # p = d1*dAd2sf*2^(b+1)
    p = d1*dAsf*d2*2^(b+1)
    p, f = make_prime(p)
    logp = p.nbits()
    print("="*80)
    print(f"Found a solution for {b = }, {c = }")
    print(f"{factor(T1) = }")
    print(f"{factor(T2) = }")
    print(f"{m1 = }")
    print(f"{m2 = }")
    print(f"{f = }")
    print("Is m1^2*T1 + m2^2*T2 == 2^(b+1)?", m1^2*T1 + m2^2*T2 == 2^(b+1))
    print("Is T1*T2 == d1*dA*d2?", T1*T2 == d1*dA*d2)
    print()
    print("Concrete smoothness bits for dA:", factor(dA)[-1][0].nbits())
    print("Concrete smoothness bits for d1d2:", factor(d1*d2)[-1][0].nbits())
    dA1, dA2 = gcd(T1,dA), gcd(T2,dA) 
    print(f"{factor(d1) = }")
    print(f"{factor(dA) = }")
    print(f"{factor(d2) = }")
    print(f"{factor(dA1) = }")
    print(f"{factor(dA2) = }")   
    print(f"{math.log(d1, 2) = }")
    print(f"{math.log(dA, 2) = }")
    print(f"{math.log(d2, 2) = }")
    print(f"{p = }")
    print(f"{logp = }")
    print("="*80)