# This is a set of functions to explore Galois Fields and their uses in
# cryptography.
#
# author - Jessica Sorrell
# date - 3/8/15
# For Advanced Cryptography at RIT with Professor Radziszowski
#
# Polynomials are represented by arrays
# Index i of the polynomial array represents the coefficient of x^i

from math import sqrt


# modPolyMultPrecomp calculates some useful congruences for multiplication mod
# a polynomial. Use if performing several multiplications within the same field
def modPolyMultPrecomp (p, m):
    n = len(p) - 1
    #define n equivalences between polynomial with degree deg(p) <= k <= 2deg(p)
    #so that we may use them easy polynomial reduction mod p
    equivalences = [[(m - p[i]) % m for i in range(n)]]
    for i in range (1,n):
        equivalences.append([0 for x in range(n)])
        for j in range (n-1):
            equivalences[i][j+1] = (equivalences[i-1][j] + \
                                   equivalences[i-1][-1]*equivalences[0][j+1])%m
        equivalences[i][0] = (equivalences[i-1][-1]*equivalences[0][0]) % m
    return equivalences 


# modPolyMult multiplies two polynomials given as arrays, modulo another
# polynomial, primitive in Z_m 
def modPolyMult (f, g, p, m, equivalencies):
    
    n = len(p) - 1
    product = [0 for x in range(2*n - 1)]
    #multiply the two polynomials
    for i in range (n):
        for j in range(n):
            product[i+j] = (product[i+j] + f[i]*g[j]) % m
    #reduce mod p            
    for i in range(n,2*n - 1):
        equiv = equivalencies[i - (n)]       
        for j in range(len(equiv)):            
            product[j] = (product[j] + product[i]*equiv[j]) % m
    return product[:n]        


# modular polynomial exponentiation
def modPolyExp (f, e, p, m, equivalencies):
    result = [0 for i in range(len(p)-1)]
    result[0] = 1
    for i in range (e):
        result = modPolyMult (f, result, p, m, equivalencies)
        #wfile.write("Result " + str(i) + ": " + str(result) + "\n")
    return result


# is this polynomial 1?
def isIdentity (f):
    identityFlag = (f[0] == 1)
    if (identityFlag):
        for i in range(1, len(f)):
            if (f[i] != 0):
                identityFlag = False
    return identityFlag
        

#Calculations for question 3 on hw3

# We're working in GF(131^2) so we need to work mod an irreducible
# polynomial of degree 2 in Z_131
# We've chosen x^2 + 1
modP = [1, 0, 1]
equivs = modPolyMultPrecomp (modP, 131)

Calculate the order of every dang ol element in the field
orderlist = []
for i in range(1,131):
    for j in range(131):
        f = [j, i, 0]
        order = 1
        g = f
        while (isIdentity(g) != True):
            g = modPolyMult(f, g, modP, 131, equivs)
            order += 1
        if (order == 17160):
           print ( "$" + str(f[1]) + "x + " + str(f[0]) + "$ is a generator! \\\\" )
        orderlist.append(order)


#Determine how many elements of each order there are
unique = set()
for order in orderlist:
    unique.add(order)

for order in unique:
    count = orderlist.count(order)
    print "There are " + str(count) + " elements of order " + str(order)  




#Shanks for polynomials
def shanks(n, alpha, beta, p, m):

    equivs = modPolyMultPrecomp(p, m)
    s = int(sqrt(n))
    
    aMJs = {}
    bAIs = {}
    logResult = 0
    for j in range (s):
        aMJs[j] = modPolyExp( alpha, s*j, p, m, equivs)

    for i in range (s):
        alphaInv = modPolyExp( alpha, (n-1), p, m, equivs)
        alphaNegI = modPolyExp( alphaInv, i, p, m, equivs)
        bAIs[i] = modPolyMult( beta, alphaNegI, p, m, equivs)

    for index in aMJs:
       for otherIndex in bAIs:
           if aMJs[index] == bAIs[otherIndex]:
               logResult = (s*index + otherIndex) % n
    return logResult

print (shanks( 17160, [3, 1, 0], [101, 1, 0], [1, 0, 1], 131))


        


