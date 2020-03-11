#!/usr/local/src/SageMath/sage
#usage : sage makeEisenJson.sage k N
import json
import sys
import functionPile as fP



if(len(sys.argv) != 3):
    print("Please run this script as 'python3 eisenCoeffMaker.py k N'")
    exit()
k = int(sys.argv[1])
N = int(sys.argv[2])

print(json.dumps(fP.genJson(k,N),indent=1))