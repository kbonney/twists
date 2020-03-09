import sys
import json
import math


def BQFsplitter(a,b,c):
    D = dict()
    for n in range(0,a+1):
        a1 = a - n
        for m in range(0,c+1):
            b1 = a - m
            ran = math.floor(2*math.sqrt(a1 * b1))
            for k in range(-ran,ran+1):
                S1cand = (a1,k,b1)
                b2 = b-k
                S2cand = (n,b2,m)
                if (b2)**2 <= 4*m*n:
                    D[str(S1cand) + " : " + str(k**2 - 4*a1*b1)] = str(S2cand) + " : " + str(b2**2 - 4*m*n)
    return D

a = int(sys.argv[1])
b = int(sys.argv[2])    
c = int(sys.argv[3])

print(json.dumps(BQFsplitter(a,b,c),indent=1))

