

Discs = (1, 5, 8, 12, 13, 17, 21, 24, 28, 29, 33, −3, −4, −7, −8, −11, −15, −19, −20, −23, −24, −31)

def giveReps(D):
    A = BinaryQF_reduced_representatives(D)
    return A

def genJson(D, k):
    A = dict()
    for x in Discs:
        for y in giveReps(D)
        a = y[0]
        b = y[1]
        c = y[2]
            A[str(y)] = str(makeCoeff(k,a,b,c))

            #edit