import numpy as np
import json
import sys

## TODO:

#finish the reduced form function using thm 2.8 in Cox

#in the alphaCalc, **store each calculated alpha form with its reduced form**
#don't add to the alpha forms if the reduced form is already in the solution set

#CONVERT ALL DISCRIMINANTS TO NEGATIVE
#they're stored in the JSON as positive
#i'm manually flipping them around in the part of the main function where i call alphaCalc
#it's ugly

#make handling of bqf's universal
#most functions referring to bqfs currently take lists [a, b, c]
#functions referring to S or Sprime take matrices [[a, b/2],[b/2, c]]

#clean up testSValues
#this has gotten messy from trying to find bqfs that actually work for case 1
#have it use alphacalc on all the discriminants that satisfy the conditions for representing p**4
#use case1Inverse on those and save only the sensible outputs
#note: this is not necessarily comprehensive

#create vars
weight = 0
level = 0
degree = 0
chosenPrime = 0

#make sure user inputs correctly
if(len(sys.argv) != 3):
    print("Please run this script as 'python coefficientSolver.py filename prime'")
    exit()
chosenPrime = int(sys.argv[2])
filename = sys.argv[1]

#open json file
def readJSON():
    #reads the JSON
    with open(filename) as f:
        data = json.load(f)
        return data

#
def calcAlphaForm(disc):
    #hardcoded search limit of 500
    solutions = []
    for b in range(disc%2, 500, 2):
        if(b**2 != disc):
            if disc%(chosenPrime**4) == (b**2)%(chosenPrime**4):
                if((disc-b**2)%(4*(chosenPrime**4)) == 0):
                    solutions.append([chosenPrime**4, b, -1*(disc-b**2)/(4*(chosenPrime**4))])
            elif disc%(chosenPrime**4) == (-b**2)%(chosenPrime**4):
                if((disc+b**2)%(4*(chosenPrime**4)) == 0):
                    solutions.append([chosenPrime**4, -1*b, -1*(disc+b**2)/(4*(chosenPrime**4))])
    return solutions

def reduce(bqf):
    #use thm 2.8
    print("bam")

def bracket(Sprime, other):
    # other^T*Sprime*other
    return np.matmul(np.matmul(np.transpose(other), Sprime), other)

def inverseBracket(S, other):
    #invert bracket, return S'
    return np.matmul(np.matmul(np.linalg.inv(np.transpose(other)), S), np.linalg.inv(other))

def legendre(input, p):
    #returns legendre (input/p)
    if (input%p == 0):
        return 0
    for x in range (1, p):
        #check if input is quad residue mod p
        if((x*x)%p == input%p):
            return 1
    return -1

def case1Inverse(S, p):
    for x in range(1, p):
        #if():
            otherVal = -float(x)*(1.0/float(p))
            otherMat = np.matrix([[1.0, otherVal], [0.0, float(p)]])
            #store the integer matrices somehow
            print(inverseBracket(S, otherMat))
    return

def case1(Sprime, p, weight):
    #calculates the coefficient corresponding to some S'
    #NOTE: Need to access the respective S coefficients somehow!!
    print("bam")

def testSValues():
    count = 0
    for disc in data["Fourier_coefficients"]:
        if (legendre(-1*int(disc), chosenPrime**4) == 1 or legendre(-1*int(disc), 6*chosenPrime**4) == 1):
            count = count+1
            #calculate coeffs mx^2 + bxy + cy^2
            alphaForm = calcAlphaForm(-1*int(disc))
            if(alphaForm != []):
                print(str(disc) + ": (" + str(alphaForm[0]) + ", " + str(alphaForm[1]) + ", " + str(alphaForm[2]) + ")")
            #alphaForm = reduce(alphForm)
            #relate to a reduced form
            for bqf in data["Fourier_coefficients"][disc]:
                #print(bqf)
                cleanedBqf = bqf.strip('()').replace(" ", "")
                values = cleanedBqf.split(",")
                #NOTE: DO WE DIVIDE THE B-VALUES BY 2???
                bqfMatrix = [[float(values[0]), float(values[1])/2.0],[float(values[1])/2.0, float(values[2])]]
                #print(bqfMatrix)
                #case1Inverse(bqfMatrix, chosenPrime)
        else:
            print(disc)
    print(count)
data = readJSON()

testSValues()
