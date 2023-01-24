import os
import numpy as np
import subprocess 
from sys import platform
import random
import time 

##functions
def genEdge(n,w):
    e0 = random.randint(0,n-1)
    e1 = random.randint(0,n-1)
    while (e0==e1):
        e1 =  random.randint(0,n-1)
    return [e0,e1,w]


def modifySimDat(simDat: dict, rate: float):
    #mutate the sim data
    nNodes = len(simDat["nodeProperties"])
    nEdges = len(simDat["edgeProperties"])
    nChangeNode = min(np.random.poisson(rate*nNodes),nNodes)
    nChangeEdge = min(np.random.poisson(rate*nEdges),nEdges)
    nAdd = min(np.random.poisson(rate*nEdges),nEdges)
    nRemove = min(np.random.poisson(rate*nEdges),nEdges)
    
    for ind in np.random.choice(np.arange(0,nNodes),replace=False,size=nChangeNode):
        simDat["nodeProperties"][ind][0] = random.randint(10,100) 
        simDat["nodeProperties"][ind][1] = 1

    for ind in np.random.choice(np.arange(0,nEdges),replace=False,size=nChangeEdge):
        simDat["edgeProperties"][ind] = genEdge(nNodes,random.randint(-4,4))

    for ind in range(nAdd):
        simDat["edgeProperties"].append(genEdge(nNodes,random.randint(-4,4)))
        
    for ind in sorted(np.random.randint(0,nEdges-1,size=nRemove),reverse=True):
        del simDat["edgeProperties"][ind]



def runSims(nSims):
    if platform == "linux" or platform == "linux2":
        runString = f"./sim {nSims}"
    else: ##windows
        runString = f"sim.exe {nSims}"
    p = subprocess.Popen(runString,stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    out, err = p.communicate()
    results = out.decode('UTF-8').split(',')
    timeTaken = results[0]
    data = np.array([[i,float(val)] for i,val in enumerate(results[1:-1])])
    data = data[data[:, 1].argsort()]

    print(f"CPP time taken: {timeTaken }")
    return np.array(data)
    
def getSimDictFromBin(loc,simNo):
    dat = np.fromfile(loc+f"simData{simNo}.bin",dtype='i4')
    nNodes = dat[0]
    nEdges = dat[1]
    nodeOffset = 2;
    edgeOffset = nodeOffset + 2 * nNodes;

    simDesc = {"nodeProperties": [],
                "edgeProperties": []}

    for i in range(0,nNodes):
        simDesc["nodeProperties"].append([dat[nodeOffset+ 2*i],dat[nodeOffset + 1+ 2*i]] )

    for i in range(0,nEdges):
        simDesc["edgeProperties"].append([dat[edgeOffset+ 3*i], dat[edgeOffset+ 1+ 3*i], dat[edgeOffset+2+ 3*i]])
    return simDesc
    

def saveSimDescriptBin(simDescript,simNo,loc):
    n = len(simDescript["nodeProperties"])
    m = len(simDescript["edgeProperties"])

    nodes = np.array([nd for nd in simDescript["nodeProperties"]],dtype=int).flatten()
    edges = np.array([ed for ed in simDescript["edgeProperties"]],dtype=int).flatten()
    saveArr = np.zeros(2+len(nodes) + len(edges), dtype = 'i4') 
    saveArr[0] = n
    saveArr[1] = m
    for i,val in enumerate(nodes):
        saveArr[2+i] = val

    for i,val in enumerate(edges):
        saveArr[2 +len(nodes) +i] = val

    saveArr.tofile(loc+f"simData{simNo}.bin",format='i4')

def genFromBest(top,decimationFactor,mutationRates):
        newNets = []
        mutRates = []
        for j in top:
            for i in range(decimationFactor):
                newNets.append(getSimDictFromBin("",j))
                mutRates.append(mutationRates[i])
        return newNets,mutRates




##########################################################################
###Running Paramters
nIterations = 100000
nNodes = 50
nEdges = 300
popSize = 500
decimationFactor = 20
assert popSize%decimationFactor == 0 
numberToSave = int(popSize/decimationFactor )
nGenerations = 10
mutationRates = np.linspace(0.1,0.001,decimationFactor)





##########################################################################
##Hardcode some parameters (most important one is filePath)
fin = open("main.cpp", "rt")
new_script_name = f"thisSim.cpp"
fout = open(new_script_name, "wt", newline='')

for line in fin:		
    # Find and change the job name
    if line.strip()==f"constexpr int nIterations;":
        op = f"constexpr int nIterations = {nIterations};"
        fout.write(op+"\n")

    elif line.strip()==f"const std::string filePath;":
        if platform == "linux" or platform == "linux2":
            fp = os.getcwd()
            filePath = f"const std::string filePath = \"{fp}/\"; " 

        else: ##windows
            fp = os.getcwd().replace("\\","\\\\")
            filePath = f"const std::string filePath = \"{fp}\\\\\"; " 


        fout.write(filePath+"\n")

    else:
        fout.write(line)

fin.close()
fout.close()


##########################################################################
##Compile c file
if platform == "linux" or platform == "linux2":
    runString = f"g++ {new_script_name} -o sim -std=c++17 -O3"
else: ##windows
    runString = f"g++ .\{new_script_name} -o sim.exe -std=c++17 -static -O3"

p = subprocess.Popen(runString,stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
out, err = p.communicate()
print("Compiling")
print(out)
print(err)
time.sleep(1)

#generate initial population
for simNo in range(popSize):
    simDescript = { 
    "nodeProperties": [[random.randint(10,100),1] for _ in range(nNodes)], 
    "edgeProperties": [genEdge(nNodes,random.randint(-5,5)) for _ in range(nEdges)]}
    saveSimDescriptBin(simDescript,simNo,"")



#run simulations
for gen in range(nGenerations):
    costs = runSims(popSize);

    top = [int(i) for i in costs[0:numberToSave,0]] #determine top x nets
    lowCost = costs[0][1]

    #save best net as an example
    bestInGen = np.fromfile(f"simData{top[0]}.bin",dtype='i4')
    bestInGen.tofile(f"best_gen_{gen}.bin",format='i4')

    newNets,mutRats = genFromBest(top,decimationFactor,mutationRates)
    
    for nN,mr in zip(newNets,mutRats):
        modifySimDat(nN,mr)
    
    for i,nN in enumerate(newNets):
        saveSimDescriptBin(nN,i,"")

    if gen==0:
        np.save(f"lowCosts.npy",np.array([lowCost]))
    else:
        lcs = [ i for i in np.load("lowCosts.npy")]
        lcs.append(lowCost)
        np.save(f"lowCosts.npy",np.array(lcs))

    print(gen,lowCost);
