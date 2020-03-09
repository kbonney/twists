

def twistedCoeff(E,p,a,b,c):
    #gotta define bracket and chi (legendre)
    #code up BQF treatment as matrix?
    if (b % p == 0 and a % p**4 == 0):
        P_1 = p**(1 - E.weight) * chi(b)
        P_2 = 0
        for r in mresidues(p, isprime = True):
            P_2 += chi(r) * E.coeff(bracket(S,A))
        
        P = P_1 * P_2
        return P
    
    if (maxdivides(p, b) == True and maxdivides(p**4,a) == True):
        P



#

