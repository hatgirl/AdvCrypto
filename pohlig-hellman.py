def primeFactors(n):
    """Returns all the prime factors of a positive integer"""
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

def egcd(a, b):
    if a == 0:
        return (b, 0, 1)
    else:
        g, y, x = egcd(b % a, a)
        return (g, x - (b // a) * y, y)

def modinv(a, m):
    g, x, y = egcd(a, m)
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
                 nextTerm *= cong2*modinv(cong2,cong1)
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
        
        for i in range(len(coeffs)):
            congruences[prime**power] += coeffs[i]*prime**i 
    return CRT(congruences)

print pohligHellman ( 131**2, 131*130, 2, 103 )
print pohligHellman ( 29, 28, 2, 18 )
        


