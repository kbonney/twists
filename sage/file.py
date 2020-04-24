import functionPile as fp
A=fp.SPMF(4)
#print(A.coeff(1,0,1))
#print(A.coeff(1,0,2))
print(A.coeff(81,0,1))
#print(A.coeff(1,0,81))
B = twistForms(A,3)
print(B.coeff(81,0,1))