load('fpMKIII.sage')
A1=SPMF(4)
A2=SPMF(4)
A3=SPMF(4)
A4=multiplyForms(A1,A2)
A5=multiplyForms(A3,A4)
A6=scaleForms(A5,-1)
B=SPMF(12)
C=addForms(A6,B)
D=twistForms(C,3)
X=BinaryQF(81,1,1)
A = matrix([[81,1/2],[1/2,1]])
S = matrix([[1,3],[0,-1/3]])
