#!/usr/bin/env sage -python
import pprint
import json
import sys
import functionPile
from fractions import Fraction



#if(len(sys.argv) != 5):
    #print"Please run this script as 'python3 eisenCoeffMaker.py a b c k' where a b c are the coefficients of the BQF index and k is the weight of the Eisenstein series"
#    exit()
#a = int(sys.argv[1])
#b = int(sys.argv[2])
#c = int(sys.argv[3])
#k = int(sys.argv[4])

#this is where we throw everything together
def makeCoeff(a,b,c,k):

    #calculate the discriminant
    disc = 4 * a * c - b * b
    
    #calculate the zeta functions
    z1 = zeta(1 - k)
    z2 = zeta(3 - (2 * k))

    #make a list of nonzero coefficients of our BQF
    list = []
    for i in [a, b, c]:
        if i != 0:
            list.append(i)
            
    sum = 0
    
    #iterate through divisors of our BQF (this forms the core sum)
    for d in range(1, min(list) + 1):
        if a % d == 0 and b % d == 0 and c % d == 0:
            #add term as given in McCarthy for each divisor d
            POW = (d ** (k-1))
            HFUN = functionPile.H(k-1, disc / (d ** 2))
            sum += POW * HFUN
    #finish off the calculation by multiplying by 2/(Z(1 - k) * Z(3 - 32))        
    value = 2 * sum / (z1 * z2) 
    return value

#tada!


    
#print('The coefficient indexed by (' + str(a) + ', ' + str(b) + ', ' + str(c) +
#') for the Siegel Eisenstein series of weight ' + str(k) +   ' is ' + str(makeCoeff(a,b,c,k)))


##########################################



def giveReps(D):
    A = BinaryQF_reduced_representatives(D)
    return A

def genJson(D, k):
    A = dict()
    for y in giveReps(D):
        a = y[0]
        b = y[1]
        c = y[2]
        A[str(y)] = str(makeCoeff(a,b,c,k))
    return A

pprint.pprint(genJson(-128,4))
#or alternatively
#print(json.dumps(genJson(D, k), indent=1))