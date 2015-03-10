from math import floor

# returns a dictionary of prime factors as indices and the number of their
# appearances in the factorization as the value
def primeFactors(n):
   
    factors = {}
    d = 2
    while n > 1:
        while n % d == 0:
            if d in factors:
                factors[d] += 1
            else:
                factors[d] = 1    
            
            n /= d
        d = d + 1
    return factors

# the extended euclidian algorithm
def euclidian(a, b):
    if a == 0:
        return (b, 0, 1)
    else:
        g, y, x = euclidian(b % a, a)
        return (g, x - floor(b / a) * y, y)


# returns the inverse of a mod m
def modInv(a, m):
    g, x, y = euclidian(a, m)
    if g != 1:
        raise Exception('modular inverse does not exist')
    else:
        return x % m


# Apply the Chinese Remainder theorem to a dictionary of congruences
# Congruences is a dictionary of congruences under various moduli
# of the form {moduli, a} where x is congruent to a mod moduli
def CRT( congruences ):
    modSum = 0
    modulus = 1
    for cong1 in congruences:
        modulus *= cong1
        nextTerm = congruences[cong1]
        for cong2 in congruences:
            if cong1 != cong2:
                 nextTerm *= cong2*modInv(cong2,cong1)
        modSum += nextTerm
    return (modSum % modulus)
        

def pohligHellmanSub( G, n, alpha, beta, q, c) :
    betas = [beta]
    ans = []
    j = 0
    while (j <= c - 1):
        delta = (betas[j]**(n/(q**(j+1))) % (G))
        i = 0
        while (delta != ((alpha**(i*n/q)) % (G)) ):
            i += 1
        ans.append(i)
        betas.append((betas[j] * alpha**(-i*(q**j)) % (G)))
        j+= 1
    return ans


def pohligHellman ( G, n, alpha, beta):
    primeDecomp = primeFactors( n )
    congruences = {}
    for prime in primeDecomp:
        power = primeDecomp[prime]
        congruences[prime**power] = 0
        
        print "Prime power " + str( prime ) + "^" + str( power )
        coeffs = pohligHellmanSub ( G, n, alpha, beta, prime, power )
        print coeffs
        for i in range(len(coeffs)):
            congruences[prime**power] += coeffs[i]*prime**i 
    return CRT(congruences)

print pohligHellman ( 131**2, 131*130, 2, 103 )
print pohligHellman ( 29, 28, 2, 18 )
        


