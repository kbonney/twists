import numpy as np
import json
import sys

## TODO:

#make sure you import level and degree!!

#finish the reduced form function using thm 2.8 in Cox

#make the alphacalc more efficient! probably using some kind of reduction or equivalence classes

#make handling of bqf's universal
#most functions referring to bqfs currently take lists [a, b, c]
#corresponding to ax^2+bxy+cy^2

#clean up testSValues
#this has gotten messy from trying to find bqfs that actually work for case 1
#have it use alphacalc on all the discriminants that satisfy the conditions for representing p**4
#use case1Inverse on those and save only the sensible outputs
#note: this is not necessarily comprehensive

#finish implementation of case1 coefficient calculation
#this is gonna involve pairing off the S bqfs that match their S'.
#Then, if you have enough matching S bqfs, calculate the coefficient for S'

weight = 0
level = 0
degree = 0
chosenPrime = 0

if(len(sys.argv) != 3):
    print("Please run this script as 'python coefficientSolver.py filename prime'")
    exit()
chosenPrime = int(sys.argv[2])
filename = sys.argv[1]


def readJSON():
    #reads the JSON
    with open(filename) as f:
        data = json.load(f)
        return data

def calcAlphaForms(disc):
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
    newA = bqf[0]
    newB = bqf[1]
    newC = bqf[2]
    #
    if newA > newC:
        newC = bqf[0]
        newA = bqf[2]
        newB = newB * -1
    if(abs(newB) > newA):
        while newB != 0 and (abs(newB) > abs(newB - (newB/abs(newB))*abs(2*newA)) or ((abs(newB) == newA) and newB < 0)):
            newC = newA-(newB/abs(newB))*(newA/abs(newA))*newB+newC
            newB = newB - (newB/abs(newB))*abs(2*newA)
            if newA > newC or (newA == newC and newB < 0):
                newCHolder = newC
                newC = newA
                newA = newCHolder
                newB = newB*-1
    return [newA, newB, newC]

def legendre(input, p):
    #returns legendre (input/p)
    if (input%p == 0):
        return 0
    for x in range (1, p):
        #check if input is quad residue mod p
        if((x*x)%p == input%p):
            return 1
    return -1

def case1Inverse(bqf, p, b):
    #takes a list of the form [a, b, c], a prime p, and a value b and returns
    #another bqf according to the inverse bracket in case 1
    newBQF = []
    newBQF.append(bqf[0])
    newBQF.append((2*b*bqf[0]+p*bqf[1])/(p**2))
    newBQF.append(((b**2)*bqf[0]+b*p*bqf[1]+(p**2)*bqf[2])/(p**4))
    return newBQF

def case1(Sprime, Scoeffs):
    print(Sprime)
    print(Scoeffs)
    denominator = legendre(Sprime[1], chosenPrime)*chosenPrime**(weight-1)
    numerator = 0
    for b in range(1, chosenPrime):
        numerator = numerator + legendre(b, chosenPrime)*Scoeffs[b-1][3]
    #calculates the coefficient corresponding to some S'
    #NOTE: Need to access the respective S coefficients somehow!!
    print(str(numerator) + "/" + str(denominator))

def findCase1SPrime(disc):
    #first check whether alpha can be a multiple of p**4
    sols = []
    if (legendre(disc, chosenPrime**4) != 1 and legendre(disc, 6*chosenPrime**4) != 1):
        return sols
    #calculates all the forms with the appropriate discriminant and an alpha value equal to p^4 (NOT COMPREHENSIVE; alpha can be multiples of p^4)
    alphaForms = calcAlphaForms(disc)
    if(alphaForms != []):
        #for each bqf that might have an appropriate corresponding S'
        for possibleBQF in alphaForms:
            #check whether the S' is going to have integer values and whether it meets the conditions for case 1
            if(possibleBQF[1]%chosenPrime == 0):
                bValue = 0
                for b in range(1, chosenPrime):
                    if (((b**2)*possibleBQF[0])+(b*chosenPrime*possibleBQF[1])+((chosenPrime**2)*possibleBQF[2]))%(chosenPrime**4) == 0:
                        bValue = b
                        break
                if(bValue != 0):
                    correspondingBQF = case1Inverse(possibleBQF, chosenPrime, bValue)
                    inserted = 0
                    for sBFQ in sols:
                        if sBFQ[0] == correspondingBQF:
                            sBFQ.append([possibleBQF, b])
                            inserted = 1
                            break
                    if(inserted == 0):
                        sols.append([correspondingBQF, [possibleBQF, b]])
    return sols

def testSValues():

    for disc in data["Fourier_coefficients"]:

        #load up our bqfs from this discriminant
        coefficientSet = []
        for BQF in data["Fourier_coefficients"][disc]:
            #print(bqf)
            cleanedBQF = BQF.strip('()').replace(" ", "")
            coefficientSet.append([list(map(int, cleanedBQF.split(","))), int(data["Fourier_coefficients"][disc][BQF])])

        discVal = int(disc)*-1

        #case1 check to find bqs
        possibleSolutions = findCase1SPrime(discVal)
        if(possibleSolutions != []):
            print(possibleSolutions)
        for solution in possibleSolutions:
            reducedSol = [[0,0,0]]*(chosenPrime-1)
            if(len(solution) == chosenPrime):
                print("Discriminant: " + str(discVal))
                print(solution)
                for bqf in solution[1:]:
                    reducedSol[bqf[1]-1] = reduce(bqf[0])
                print(reducedSol)
            for x in range(len(reducedSol)):
                for y in range(len(coefficientSet)):
                    if reducedSol[x] == coefficientSet[y][0]:
                        reducedSol[x].append(coefficientSet[y][1])
                        break
            if([0,0,0] not in reducedSol):
                case1(solution[0], reducedSol)
            #here we can retrieve coefficients from the old form by comparing values
            #to our reduced bqfs which we should have calculated above

data = readJSON()
weight = int(data["weight"])
testSValues()
