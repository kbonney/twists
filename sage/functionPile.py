from sage.functions.transcendental import zeta
from sage.quadratic_forms.special_values import quadratic_L_function__exact as lfun
import numpy as np
import math
from fractions import Fraction


#########################################################

#here we have a simple prime checker
#utilzes theorem that says n composite ==> n has prime divisor less than sqrt(n)
#so the contrapositive says if n has no prime divisor less than sqrt(n) then n is prime

def isPrime(n):
    if (n < 2) :
        return False
    for i in range(2, n + 1) :
        if (i * i <= n and n % i == 0) :
            return False
    return True


#########################################################

#here we sum over the divisors of n raised to the kth power
#this is known as the divisor function

def sigmaDivisor(k, n):
    y = 0
    #could cut this to only go to n+1/2 if n gets big
    #~~it probably will~~
    for x in range(1,int(n+1)):
        #we check if x is a divisor of n
        if n % x == 0:
            #if it is a divisor we add its kth power to the value
            y += x ** k
    return int(y)


#########################################################

#found this function online, but it seems to work and the code makes sense so we use it

def mobius(N) :
    # Base Case
    if (N == 1) :
        return 1
    # For a prime factor i
    # check if i^2 is also
    # a factor.
    p = 0
    for i in range(1, N + 1) :
        if (N % i == 0 and
                isPrime(i)) :
            # Check if N is
            # divisible by i^2
            if (N % (i * i) == 0) :
                return 0
            else :
                # i occurs only once,
                # increase f
                p = p + 1
    # All prime factors are
    # contained only once
    # Return 1 if p is even
    # else -1
    if(p % 2 != 0) :
        return -1
    else :
        return 1


#########################################################

#fundamental discriminant finder that Jake came up with, much easier to comprehend.

def fundDisc(N):
    sign = np.sign(N)
    N = abs(N)
    #strip out squares
    for x in range(2, int(math.floor(N ** 0.5)) + 1):
        if 0 == (N % (x**2)):
            return fundDisc(sign * N/(x**2))
    #re-sign the integer
    temp = sign * N
    if 1 == temp % 4:
        return temp
    return 4 * temp

def isFundDisc(N):
    if N == fundDisc(N):
        return True
    else:
        return False


#########################################################

#quick integer sqrt finder, honestly not sure how it works but
#focused on other things right now, will come back once the important
#code is squared away

def isqrt(n):
    x = n
    y = (x + 1) / 2
    while y < x:
        x = y
        y = (x + n / x) / 2
    return x


#########################################################

# kronecker symbol calculator, algorithm taken from Cohen's computational number theory book

def kronecker(a, b):
    if 0 == b:
        if abs(a) != 1:
            return 0
        else:
            return 1
    #remove 2s
    if 0 == a % 2 and 0 == b % 2:
        return 0
    v = 0
    while 0 == b % 2:
        b /= 2
        v += 1
    if 0 == v % 2:
        k = 1
    else:
        if 1==a%8 or 7==a%8:
            k = 1
        elif 3==a%8 or 5==a%8:
            k=-1
    if b < 0:
        b *= -1
        if a < 0:
            k *= -1
    #finished?
    return kronecker2(a,b,k)

def kronecker2(a,b,k):
    if 0 == a:
        if 1 < b:
            return 0
        if 1 == b:
            return k
    v = 0
    while 0 == a % 2:
        v += 1
        a /= 2
    if 1 == v % 2:
        if 1==b%8 or 7==b%8:
            k *= 1
        elif 3==b%8 or 5==b%8:
            k*=-1
    #apply reciprocity
    if 3 == a % 4 and 3 == b % 4:
        k *= (-1)
    r = int(abs(a))
    a = b % r
    b = r
    return kronecker2(a,b,k)


#########################################################


#this function returns divisors of n that take the form p^n, where p is a prime and n is one less than the maximal power of p dividing n

def findPrimes(n):
    #create our set which will contain the prime powers of n with power greater than 1
    #Note: we cannot create a empty set, otherwise it becomes a dictionary, so we start with 1 here and remove it later
    #itereate through integers between 2 and half of n
    primePowers = {1}
    nAb = abs(n)
    for i in range(2, int((nAb/2)+1)):
        #check if i is prime and divides n
        if isPrime(i) and nAb % i == 0:
            #print(i)
            #find largest power (greater than 1) of i that divides n
            y = 1
            while nAb % i ** (y + 1) == 0:
                y += 1
            #when i ** (y+1) does not divide we add i ** y, the greatest power of i dividing n, to our prime powers list
            if y % 2 == 0:
                primePowers.add(i ** y)
            if y % 2 == 1:
                primePowers.add(i ** (y-1))
    #here we remove 1, which was added earlier due to syntax restraints
    primePowers.remove(1)
    return primePowers


#########################################################


# getFundDisc() finds D0 a fundamental discriminant where D=D0f^2
# Jake made a better function. This is obsolete. 

def getFundDisc(D):
    #the way i have constructed the algorithm makes -4 a weird case, so we deal with it seperately for now
    if D == -4:
        return -4
    #here we define D_0, the fundamental discriminant in n
    D_0 = D
    #we iterate through the prime powers given by findPrimes
    for primep in findPrimes(D):
        #we divide out these prime powers to get a squarefree divisor of n
        D_0 = D_0 / (primep)
    #we check the fundamental discriminant condition concerning modulo n (see any standard definition of a fundamental disc)
    #and adjust accordingly
    if D_0 % 4 == 1:
        return D_0
    elif D % 4 == 0:
        return D_0 * 4


#########################################################


def H(k1,N):
    #here we assign variables in a way that meshes with notation in McCarthy and Raum
    k = k1 + 1
    N = -N
    s = 2 - k
    #when N = 0 we calculate the value with the zeta func
    if N == 0:
        x = zeta((2 * s) - 1)
        return x
    #we check if N is 0 or 1 mod 4
    if N % 4 == 0 or N % 4 == 1:
        #find our kronecker character
        D0 = fundDisc(N)
        #find f
        f = isqrt(int(N / D0))
        #initialize our value
        hSum = 0
        #iterate through the divisors of f
        for d in range(1, f + 1):
            if f % d == 0:
                #when d divides f we add the appropriate term to the sum
                MOB = mobius(d)
                KRON = kronecker(D0, d)
                POW = d ** -s
                SIG = sigmaDivisor(1 - (2 * s), int(f / d))
                #print('    adding ' + str(MOB) + ' * ' + str(KRON) + ' * ' + str(POW) + ' * ' + str(SIG) + ' to H func sum.')
                hSum += MOB * KRON * POW * SIG
        LF = lfun(s, D0)
        fCalc = hSum * LF
        return fCalc
    #otherwise the value is 0
    else:
        return 0


#########################################################